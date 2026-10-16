# Ahringer multiome pipeline

This pipeline re-implements the single-nucleus multiome (10x GEX + ATAC) analysis from:

> *A lineage-resolved multimodal single-cell atlas reveals the genomic dynamics of early C. elegans development*
> (Ahringer lab), bioRxiv 10.1101/2024.12.02.626321

It is built from the authors' supplementary code (`paper_original_scripts/`, including the master command log
`sc_project_pipeline.20241202.txt`). The paper's scripts hard-code 7 experiments, 3 batches and their Cambridge cluster paths. Here they
are driven by a sample sheet and submitted to Slurm stage by stage. It starts from raw FASTQs.

## 1. Installation (once)

```bash
# from your workstation
rsync -avz multiome_pipeline izblisbon.izb.unibe.ch:/mnt/meister.data/pmeister/Claude/
ssh izblisbon.izb.unibe.ch
cd /mnt/meister.data/pmeister/Claude/multiome_pipeline
```

1. **Cell Ranger ARC**: download it from 10x Genomics (it needs a licence click-through, so the pipeline can't fetch it). The paper used
   2.0.1; 2.0.2 is fine. Set `CELLRANGER_ARC_DIR` in `config/config.sh`.
2. Edit `config/config.sh`: set `PROJECT_DIR` (output location), the Slurm account/partition and the resources.
3. Fill in `config/samples.tsv` (see below).
4. On the login node (needs internet), run `bash setup.sh`. It:
   * creates 4 mamba environments (`envs/*.yaml`): `mo_tools`, `mo_macs`, `mo_r` (Seurat **4**), `mo_meme`
   * creates the `PROJECT_DIR` layout, links `external_data/` and `scripts/`, and extracts the 10x GEX whitelist
   * downloads the WormBase WS285 genome, canonical GTF and GFF3 (EBI mirror)

   Parts can be re-run on their own: `bash setup.sh envs|project|reference`.

### Sample sheet (`config/samples.tsv`, tab-separated)

One row per sample (one 10x multiome library). Several FASTQ runs of the same library go in the same row,
comma-separated; never split one sample over several rows (the array jobs would write to the same folder).

| column | meaning |
|---|---|
| `sample_id` | short id, preferably **letters/digits only** (it prefixes barcodes, e.g. `ctrl1_ACGT...`); never use `-` |
| `gex_fastq_dir`, `gex_fastq_prefix` | GEX FASTQ folder and bcl2fastq sample name(s); comma-separate several runs (`exp031_rna,exp031_more_rna`) |
| `atac_fastq_dir`, `atac_fastq_prefix` | same for ATAC |
| `batch` | processing batch: nuclei that were prepared / captured / sequenced together. Give every independently processed replicate its own batch. Used for batch-enriched cluster detection (stage 7) and SoupX soup ranges |
| `stage_group` | the biological condition of the sample (developmental stage, strain, treatment). **Batches are only compared with other batches of the same `stage_group`**, so differences between conditions are never treated as batch artifacts |
| `idr_rep` | `1` or `2`: pseudo-replicate for IDR (you need samples in both) |
| `utr_ext` | `yes` if the sample is used for 3' UTR extension and bulk scATAC peaks (paper: early-embryo samples only) |
| `barhop`, `ambient` | UMI thresholds for EmptyDrops/SoupX. Leave `NA` until stage 4, then fill them in (see below) |

FASTQs must follow the bcl2fastq naming `<prefix>_S<n>[_L00<n>]_R{1,2}_001.fastq.gz`.

#### Choosing `batch` and `stage_group` (batch-enriched cluster removal, stage 7)

Stage 7 removes WNN clusters whose nuclei come mostly from one batch, on the assumption that a real cell type
is found in every replicate of the same condition. Within each `stage_group` it compares, for every cluster,
the share of each batch's own nuclei that fall in that cluster (so batches of different size are comparable).
A cluster is flagged if one batch's share is more than `BATCH_FOLD` (4) times another's. A stage group is only
compared for a cluster when it contributes meaningfully to it: at least `BATCH_MIN_GROUP_NUCLEI` (30) nuclei
**and** at least `BATCH_MIN_GROUP_FRACTION` (20 %) of the cluster. A cluster that is specific to one condition
is therefore not removed because a few stray nuclei of the other condition happen to sit in one replicate.

Rules of thumb:
* one `batch` per independently processed replicate; replicates of the same condition share a `stage_group`
* conditions you want to compare (strains, treatments, stages) get **different** `stage_group`s. If two
  conditions share a stage group, a cell type that is lost or gained in one of them looks like a batch effect
  and is removed
* a stage group with only one batch is never checked, so with a single replicate per condition nothing is flagged
* leave `LATE_STAGE_GROUPS` (`config/config.sh`) at `late` unless you want the paper's extra rule: clusters
  made up of more than `LATE_MAX_FRACTION` (80 %) nuclei of those stage groups are flagged
* for a first analysis set `BATCH_REMOVAL_MODE=report` in `config/config.sh`: flagged clusters are listed but kept.
  Check `plots/batch_removal/merged_seurat.<annotation>.batch_enrichment_table.tsv` (nuclei per batch, reason for
  flagging, top marker genes) and `batch_enrichment.pdf`, then switch back to `remove` if they are artifacts.
  `PROTECT_CLUSTER_MARKERS` lists marker genes whose clusters are never removed

Example: two replicates of a control strain and two of a test strain, in which some cell types may be lost

| sample_id | batch | stage_group | idr_rep |
|---|---|---|---|
| ctrl1 | ctrl_b1 | control | 1 |
| ctrl2 | ctrl_b2 | control | 2 |
| test1 | test_b1 | test | 1 |
| test2 | test_b2 | test | 2 |

`ctrl1` is compared only with `ctrl2` and `test1` only with `test2`. A cluster present in both control replicates
but missing from the test strain is kept; a cluster found in only one replicate of a strain is flagged.

#### `barhop` and `ambient`

Set them after stage 4 from `plots/barcode_rank_qc/*.pdf` (histogram of UMIs per barcode):
* `barhop`: a value in the gap between the large spike of near-empty barcodes (barcode hopping, sequencing
  errors; a few to ~40 UMIs) and the next peak. Barcodes below it are discarded.
* `ambient`: the upper end of the second peak, which is made of empty droplets containing only ambient RNA
  (the soup). Barcodes between `barhop` and `ambient` define the soup for EmptyDrops and SoupX; only barcodes
  above `ambient` can be called nuclei. Do not set it so high that it cuts into the tail of real nuclei.
* PMW941_1 (L3): `barhop` = 100, `ambient` = 300. Paper (embryos): 20–100 / 150–400.

Check `plots/emptydrops/*ambient_P.pdf` after stage 5: the EmptyDrops P-values of the soup barcodes should be
roughly uniform; a large spike near 0 means `ambient` includes barcodes that are not ambient RNA.

## 2. Running

```bash
./run_pipeline.sh <stage> [--dry-run] [--after JOBID]
```
Each stage submits its Slurm jobs, with dependencies inside the stage. Start the next stage once the previous one has finished
(`squeue -u $USER`; logs are in `PROJECT_DIR/logs`). Use `--after <jobid>` to chain stages that need no manual step.
`--dry-run` prints the `sbatch` commands without submitting anything.

| stage | what it does | original scripts |
|---|---|---|
| `1` reference | annotation tables, filtered gene set, promoter/5'UTR BEDs, operons; `cellranger-arc mkref`; STAR index WS285 | pipeline "annotation preprocessing", `STARsolo_indexing.sh` |
| `2` align | `cellranger-arc count`; fragment files without `-1` and with sample prefix; STARsolo (WS285); per-cell UMI dedup of the BAM | `cellranger_arc_launcher.sh`, `atac_fragments_tsv_parser*`, `atac_fragment_sample_barcode*`, `STAR_launcher.sh`, `bam_barcode_split*`, `STAR_dedup*`, `merge_dedup_bam*` |
| `3` utr_extension | MACS2 on pooled scATAC; 3' extension of transcripts (100 bp bins up to 500 bp, never across an accessible site); extended GTF; STAR index | `macs2_cellranger_fragments_all*`, `gtf_3_utr_extension*.py`, `coveragebed_launcher*` |
| `4` starsolo | STARsolo on the extended annotation (unique + EM multimappers); barcode UMI plots | `STAR_launcher.sh`, `STARsolo_output_parser*` |
| **manual** | open `plots/barcode_rank_qc/*.pdf` and set `barhop` (top of the barcode-hopping peak) and `ambient` (EmptyDrops `lower`) per sample. Paper values: 20–100 / 150–400 | |
| `5` sample_qc | EmptyDrops (FDR<0.001); per sample: mito <5 %, ATAC counts >200 and <5 MAD, scDblFinder, <10,000 UMIs, SoupX (P-cell clusters excluded from the estimate, soupQuantile 0.8), SCT, ATAC assay on bulk peaks | `emptydrops*`, `seurat_soupx*` |
| `6` merge | merge samples, SCT/PCA clustering, ATAC counts, fragments per cluster | `merge_multiple_seurat_objects*` |
| `7` wnn_annotation | MACS2 per cluster; non-exonic peak set; TF-IDF/LSI and WNN clustering; removal of batch-enriched clusters and nuclei with ≤200 ATAC features; re-clustering; marker-based annotation of AB and E/MS/C/D sub-clusters | `macs2_clusters.from_signac*`, `*cluster_peaks_step_2*`, `batch_enriched_clusters_removal*`, `lineage_specific_annotation*` |
| **manual (optional)** | inspect `plots/lineage_annotation`, `plots/batch_removal`; fill `config/late_cluster_reassignment.tsv` for late clusters; adjust `UNASSIGN_CELL_TYPES` | `manual_clusters_reassignment.20240430.txt` (not published) |
| `8` cell_types | 2-cell, Z2/Z3, P2/P3/P4 (m1m2 element accessibility), late clusters, lineages; fragments per cell type/lineage | `final_cell_type_assignment*` |
| `9` peaks_idr | MACS2 per cell type and lineage; pre-IDR peak set; pseudo-replicate IDR per cell type; final peak set (IDR<0.05); embryo stage metadata and FRiP; random-forest early/late cell-cycle calls | `ATAC_frags_for_IDR.R`, `idr_atac.sh`, `idr_vs_all_macs*`, `seurat_integrate_cluster_peaks*`, `cell_cycle_random_forest.R` |
| `10` tracks (optional) | strand-specific RNA bigWigs per cell type (UMI-deduplicated); ATAC CPM bigWigs averaged per embryo stage | `bam_barcode_split.nosuffix.sh`, `merge_replicate_bam*`, `dedup_launcher.clusters*`, `combine_bw_tracks_stage*` |
| `11` transient_genes (optional) | pseudo-bulk DESeq2, pre- vs post-gastrulation lineages (padj<0.001, LFC>2, lfcSE<0.7) | `transient_genes_DESeq*` |
| `motifs <bed> <name>` (optional) | MEME-ChIP (zoops, 6–15 bp, 2nd-order background from all accessible sites), logos, FIMO | `memechip*` |

Main outputs (in `PROJECT_DIR`):
* `seurat_objects/wt_all.all_samples.all_cells.WS285_extended_no_ovlp.final_peakset_IDR_cell_cycle.rds`: final object
  (assays RNA, SCT, WNN, final_peakset, m1m2_elements; metadata cell_type, cell_lineage, stage, FRiP, cell_cycle_rf)
* `cluster_specific_peaks/round_2_all_samples/WS285_extended_no_ovlp/all_peaks.bed`: final accessible sites;
  `cell_type/final_peakset/*.peaks.bed` per cell type; `cell_type/bw/*.bw` ATAC tracks
* `species/elegans/gene_annotation/*extended_no_ovlp.3prime_extended.rnaseq_corrected.gtf`: 3'-extended annotation
* `plots/<step>/`: all QC and annotation plots

## 3. Differences from the original code

Changes that affect results only where the original would have failed or behaved inconsistently:
* **Software versions** (conda availability): Seurat 4.3.0.1, Signac 1.10.0 (paper 1.6.0), scDblFinder 1.12 (1.8.0),
  DropletUtils 1.18 (1.14.2), R 4.2. STAR 2.7.10a, MACS2 2.2.7.1, IDR 2.0.4.2 and MEME 5.5.2 are identical.
* Any number of samples/batches. The batch rule is generalised (see the header of `scripts/batch_enriched_clusters_removal.R`).
  With the paper's design it reproduces the original rule.
* IDR pseudo-replicates come from `idr_rep`. The paper split by `orig.ident`, which after cleanup reads `expXXX_clean`, so `exp` is used instead.
* Total fragment counts (used for FRiP) are matched to cells by barcode. The paper assigned them by position.
* Final annotation: unannotated nuclei are detected by absence from the marker tables, not by the `exp` prefix. P-cell sub-clusters are
  ranked for any number of sub-clusters. The late-cluster table is optional.
* Marker binarisation: if the second-derivative method fails for a marker, 25 %/60 % of the maximum is used as a fallback, with a warning.
  The same applies to the SoupX P-cell step: if it fails, no clusters are excluded.
* IDR failures for one cell type produce an empty set for that type instead of stopping the pipeline.
* Debug prints that would crash on other data were removed from `gtf_3_utr_extension.bin_selection`. The gene filter uses a set (same output, much faster).
* The per-transcript coverage of the *annotated* transcripts (input of the extension bin selection) is not in the published command log.
  It is recomputed here the same way as for the extended bins (whole-locus BED6, strand-specific, deduplicated reads).
* No date stamps in output file names.

## 4. Not included

These parts of the paper need data other than the multiome FASTQs, or files that were not published:
* bulk RNA-seq of mutants (kallisto/DESeq2), ChIP-seq/CUT&Tag (bwa, YAPC)
* promoter annotation from nuclear RNA-seq (lcap "jump" test). The Serizay accessible-site file in the supplement has no gene column,
  so the "rescued promoter" step can't be reproduced either
* inherited vs zygotic promoters (needs germline ATAC from Serizay et al.) and WormCAT enrichment (annotation CSV not provided)

## 5. Caveats

* The marker tables (`external_data/markers_UMAP_annotation.*.tsv`), P-granule genes and m1m2 elements are specific to the
  **early embryo** (up to about 100 cells). For other stages, the annotation stages (7–8) need new marker tables.
* As in the paper, a transcript is only 3'-extended if its extended locus does not overlap any other coding/lincRNA/pseudogene
  locus, including genes nested in its introns. Many long genes are therefore never extended.
* Thresholds that the authors set by eye (`barhop`/`ambient`, late-cluster table, `UNASSIGN_CELL_TYPES`) must be re-evaluated
  for new data.
* The pipeline has been syntax-checked, but it has not been run end-to-end on real data. Expect to tune resources (`RES_*`)
  and check the plots at each manual step.

