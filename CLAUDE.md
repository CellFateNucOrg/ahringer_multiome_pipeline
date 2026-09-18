# Pipeline context for Claude

## What this pipeline does

This pipeline re-implements the single-nucleus multiome (10x GEX + ATAC) analysis from the Ahringer lab
(bioRxiv 10.1101/2024.12.02.626321). It processes raw FASTQs through alignment, QC, clustering, peak calling,
and cell-type annotation for C. elegans single-nucleus data.

**Current use**: running on L3 larvae data (not embryos as in the paper). One sample (`PMW941_1`) with
one GEX prefix (`TEV0_GEX`) and four ATAC prefixes (`TEV0_ATAC_1`–`TEV0_ATAC_4`).

**Cluster**: SLURM HPC at University of Bern. NFS mount: local `/Volumes/meister.data/` = cluster
`/mnt/meister.data/`. Always delete files from the cluster side, not the local NFS mount.

**Key files**:
- `config/config.sh` — sourced by every job; edit resources and paths here
- `config/samples.tsv` — one row per sample, comma-separated prefixes for multi-FASTQ fields
- `run_pipeline.sh <stage>` — submits SLURM jobs for one stage
- `check_logs.sh [N_hours]` — scans logs/ for PASS/FAIL/RUNNING status

## Stage 2 pipeline structure (align)

Jobs submitted and their dependency chain:

```
cellranger-arc count ──┬──► atac_fragments (skip if outputs exist)
                       │
STARsolo ──────────────┴──► star_sort ──► star_dedup
```

All jobs have skip logic — resubmitting stage 2 is safe; only jobs with missing outputs will run.

## Known issues and fixes applied

### samples.tsv format
One row per sample with **comma-separated** ATAC prefixes in `atac_fastq_prefix`. Multiple rows for one
sample causes a race condition (N array tasks all writing to the same output directory).

### cellranger-arc count
- `--create-bam true` is required in cellranger-arc 2.0+ — added to `stages/jobs/cellranger_count.sh`.
- If a SLURM job is cancelled, Martian may leave a `_lock` file. Delete it on the cluster and resubmit;
  cellranger resumes from its checkpoint automatically.

### STARsolo
- 10x Multiome GEX R1 reads are 29 bp, not 28 bp. Fix in `config.sh`:
  `STAR_EXTRA_OPTS="--soloBarcodeReadLength 0"`
- `--outSAMtype BAM Unsorted` (used to split alignment from sorting) is **incompatible** with:
  - `--outWigType wiggle` — removed; not used downstream
  - `CB` and `UB` in `--outSAMattributes` — removed; dedup uses `CR`/`UR` instead
- BAM sorting OOM-killed at 96 GB. Fix: separate `star_sort.sh` job at 128 GB using `samtools sort`
  with local `/scratch/$USER/` via `SCRATCH_DIR` in `config.sh`.

### star_sort.sh (new job)
Runs `samtools sort` + `samtools index`, then generates forward/reverse strand CPM bigWigs via
`bamCoverage`. Outputs: `star_WS285_Aligned.sortedByCoord.out.bam` and `bw/<id>.star_WS285_{fwd,rev}.bw`.
Deletes the unsorted BAM after sorting to save space. Cleans up scratch dir on exit.

### sinto / pkg_resources error
bioconda's `sinto` package has a broken `pkg_resources` dependency. Fix:
```bash
pip install --force-reinstall --no-deps sinto
```
**Do not** use `pip install --force-reinstall sinto` (without `--no-deps`) — this tries to rebuild
`umi_tools` from source and fails. `umi_tools` must always be installed via conda.
The fix is also applied in `setup.sh` after env creation.

### Resource settings (tuned for this dataset)
| Variable | CPUs | RAM | Time |
|---|---|---|---|
| `RES_starsolo` | 12 | 64 GB | 36 h |
| `RES_star_sort` | 8 | 128 GB | 24 h |
| `RES_dedup` | 8 | 96 GB | 24 h |
| `RES_merge` | 4 | 192 GB | 24 h |
| `RES_r_large` | 4 | 192 GB | 48 h |

### Monitoring
- Enable SLURM failure emails: `SLURM_EXTRA="--mail-type=FAIL --mail-user=<email>"` in `config.sh`
- `bash check_logs.sh [N_hours]` — prints PASS/FAIL/RUNNING per log file with first error line shown

## Differences from the paper's scripts
The paper used early embryo data. For L3 larvae:
- Marker tables (stages 7–8) need to be replaced with L3-appropriate markers
- `utr_ext` and `barhop`/`ambient` thresholds must be re-evaluated
- `stage_group` in `samples.tsv` is set to `l3`
