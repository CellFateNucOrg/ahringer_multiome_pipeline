# C. elegans L3 Cell-Type Marker Files

Generated for annotation of L3 snRNA-seq clusters from epi-multiome data.
Compatible with the `scripts/lineage_specific_annotation.R` pipeline format.

## File overview

| File | Genes | Cell types | Description |
|------|-------|-----------|-------------|
| `markers_cao2017_L2.tsv` | 6,147 | 72 | Cao et al. 2017 L2 data (Science) |
| `markers_cengen_L4_herm.tsv` | 3,189 | 168 | CeNGEN L4 hermaphrodite (Taylor et al. 2021) |
| `markers_cengen_L1_herm.tsv` | 3,193 | 179 | CeNGEN L1 hermaphrodite (Taylor et al. 2025 bioRxiv) |
| `markers_cengen_adult_herm.tsv` | 3,584 | 187 | CeNGEN adult hermaphrodite |
| `markers_ghaddar2023_adult.tsv` | 159 | 156 | Ghaddar et al. 2023 adult (Science Advances) |
| `markers_shared_min2studies.tsv` | 2,243 | 160 | **Shared: genes confirmed in ≥2 studies** |
| `markers_expansive_all_studies.tsv` | 8,963 | 381 | **Expansive: all markers from ≥1 study** |
| `cell_type_taxonomy.tsv` | — | 305 | Taxonomy: raw_name → tissue_type.subtype.cell |

## Format

All marker TSV files use the pipeline-compatible format:
- **Rows**: gene names (WormBase gene symbols)
- **Columns**: cell type names
- **Values**: `0` = not expressed; `1` = low expression expected; `2` = high expression expected

Per-study files use raw cell type names from each study.
Combined files (`shared`, `expansive`) use canonical hierarchical names: `tissue_type.tissue_subtype.cell_type`

## Hierarchical taxonomy (combined files)

The three-level hierarchy used in column names of the combined files:

| Tissue type | Subtypes | Example columns |
|-------------|----------|-----------------|
| Neurons | Cholinergic, GABAergic, Dopaminergic, Sensory, Interneuron, Head_motor, Pharyngeal_neuron, Motor | `Neurons.Dopaminergic.CEP`, `Neurons.GABAergic.VD`, `Neurons.Sensory.AFD` |
| Muscle | Body_wall_muscle, Vulval_muscle, Anal_muscle, Intestinal_rectal_muscle, Uterine_muscle, Sex_myoblast, Head_mesodermal_cell | `Muscle.Body_wall_muscle.Body_wall_muscle`, `Muscle.Head_mesodermal_cell.hmc` |
| Intestine | Intestine | `Intestine.Intestine.Intestine`, `Intestine.Intestine.Intestine_anterior` |
| Hypodermis | Hypodermis, Seam_cell | `Hypodermis.Hypodermis.Epidermis_hyp7`, `Hypodermis.Seam_cell.Seam_cell` |
| Pharynx | Pharyngeal_muscle, Pharyngeal_epithelia, Pharyngeal_gland | `Pharynx.Pharyngeal_gland.Pharyngeal_gland`, `Pharynx.Pharyngeal_epithelia.Marginal_cell` |
| Germline | Germline, Sperm, Oocyte | `Germline.Germline.Germline`, `Germline.Sperm.Sperm` |
| Glia | Amphid_sheath, Cephalic_sheath, Phasmid_sheath, GLR, Glia, Socket_cell | `Glia.Amphid_sheath.AMsh`, `Glia.GLR.GLR` |
| Gonad | Gonadal_sheath, Distal_tip, Uterine, Spermatheca, Vulval | `Gonad.Distal_tip.Distal_tip`, `Gonad.Spermatheca.Spermatheca` |
| Coelomocytes | Coelomocytes | `Coelomocytes.Coelomocytes.Coelomocytes` |
| Excretory | Excretory, Excretory_gland, Rectal_gland | `Excretory.Excretory.Excretory_cell` |

## Scoring thresholds

### CeNGEN (L1, L4, adult hermaphrodite)
- Source data: tidy long table with `avg_expr` (mean TPM) and `pct_expr` (% cells expressing)
- Specificity score = avg_expr_in_cell_type / mean_avg_expr_across_all_cell_types
- **Score 2**: gene ranks in top 15 by specificity AND pct_expr ≥ 50% in that cell type
- **Score 1**: gene ranks in top 30 by specificity AND pct_expr ≥ 25% (not already score 2)
- Cell types excluded: Unannotated, Unknown_non-neuronal variants

### Cao 2017 (L2, Tables S6/S7/S8)
- Source: enrichment ratios and p-values from supplementary tables
- S6: tissue-enriched; S7: cell-type enriched; S8: neuron cluster-enriched
- **Score 2**: ratio ≥ 10 AND p < 0.05
- **Score 1**: ratio ≥ 3 AND p < 0.05

### Ghaddar 2023 (adult)
- Source: 10 pre-selected markers per cluster, with marker_score (0–1)
- **Score 2**: top 5 genes per cell group by marker_score
- **Score 1**: remaining 5 genes

## Combined file logic

- `markers_shared_min2studies.tsv`: only genes × cell_type combinations supported by ≥2 independent studies; takes the maximum score across studies
- `markers_expansive_all_studies.tsv`: all genes × cell_type combinations from any study; takes the maximum score

## Data sources

1. **Cao et al. 2017** — *Science* 357:661. L2 stage scRNA-seq (4,000 cells). Supplementary Tables S6 (tissue-enriched), S7 (cell-type enriched), S8 (neuron clusters).
2. **Taylor et al. 2021** — *Cell* 185:2952. CeNGEN L4 hermaphrodite scRNA-seq. Data via CeNGEN pins board (cengenproject/cengen-reference-data, GitHub).
3. **Taylor et al. 2025** — bioRxiv. CeNGEN L1 hermaphrodite. Same source.
4. **CeNGEN adult hermaphrodite** — same pins board.
5. **Ghaddar et al. 2023** — *Science Advances* 9:eadg0506. Adult whole-organism scRNA-seq. Supplementary Data S1.

## Validation spot-checks

| Gene | Expected cell type | Score in shared file |
|------|--------------------|----------------------|
| `dat-1` | Dopaminergic (CEP, ADE, PDE) | 2 at all three |
| `unc-47` | GABAergic (VD, DD, AVL, DVB) | 2 at all four |
| `gcy-8` | AFD sensory neuron | 2 |
| `ttx-1` | AFD sensory neuron | 2 |
| `myo-3` | Body wall muscle | 2 |
| `cat-2` | Dopaminergic (PDE) | 2 |
| `act-5` | Intestine | 2 |

## Usage in the pipeline

Place these files in the `species/elegans/` or `data/external_data/` directory as configured
in your pipeline run script. In `00_run_ahringer_pipeline.sh`, pass the marker file path to
stage 5 (lineage_specific_annotation). The `shared` file is recommended as the primary
reference; use `expansive` for finer-grained annotation.
