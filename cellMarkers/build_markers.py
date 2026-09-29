#!/usr/bin/env python3
"""
Generate C. elegans L3 cell-type marker files from 5 datasets:
  1. Cao 2017 L2 (S6 tissue, S7 cell type, S8 neuron clusters)
  2. CeNGEN L4 hermaphrodite
  3. CeNGEN L1 hermaphrodite
  4. CeNGEN adult hermaphrodite
  5. Ghaddar 2023 adult

Output format: gene × cell_type TSV, values 0/1/2 (not expressed / low / high)
"""

import pandas as pd
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

SCRATCH = Path("/tmp/claude-0/-home-claude/aab2b907-87b8-583f-86b7-bade61b61f50/scratchpad")
UPLOADS = Path("/root/.claude/uploads/aab2b907-87b8-583f-86b7-bade61b61f50")
OUTDIR  = Path("/mnt/user-data/outputs")
OUTDIR.mkdir(parents=True, exist_ok=True)

# =========================================================================
# TAXONOMY  (tissue_type, tissue_subtype, canonical_cell)
# Used to build column headers in the combined/shared/expansive files
# and to cross-reference markers across studies.
# =========================================================================

# Each entry:  raw_name → (tissue_type, subtype, canonical)
# raw_names come from all three levels of all studies
TAXONOMY = {}

def T(raw, tissue, subtype, canonical):
    TAXONOMY[raw] = (tissue, subtype, canonical)

# ---- NEURONS ----
# Dopaminergic
for n in ["CEP","ADE","PDE","Dopaminergic","CEP_ADE_PDE"]:
    T(n,"Neurons","Dopaminergic","Dopaminergic")
# canonical individual cells
T("CEP","Neurons","Dopaminergic","CEP")
T("ADE","Neurons","Dopaminergic","ADE")
T("PDE","Neurons","Dopaminergic","PDE")

# GABAergic motor neurons
for n in ["VD","DD","VD_DD","DVB","AVL","AVL and DVB","GABAergic","GABAergic_neurons"]:
    T(n,"Neurons","GABAergic","GABAergic")
T("VD","Neurons","GABAergic","VD")
T("DD","Neurons","GABAergic","DD")
T("VD_DD","Neurons","GABAergic","VD_DD")
T("DVB","Neurons","GABAergic","DVB")
T("AVL","Neurons","GABAergic","AVL")
T("RIS","Neurons","GABAergic","RIS")

# Cholinergic motor neurons
for n in ["Cholinergic_neurons","Cholinergic"]:
    T(n,"Neurons","Cholinergic","Cholinergic")
for n in ["DA","DA1","DA9","DB","DB01","DB02","VA","VA1","VA12","VB","VB01","VB02",
          "AS","SAA","SAB","AS; DA_VA; DB_VB; and SAB"]:
    T(n,"Neurons","Cholinergic",n)
T("AS; DA_VA; DB_VB; and SAB","Neurons","Cholinergic","AS_DA_VA_DB_VB_SAB")

# Head/neck motor neurons
for n in ["SIA","SIB","SMB","SMD","SMD and subpopulation of RMD_LR and SAA"]:
    T(n,"Neurons","Head_motor",n.split()[0])
for n in ["RMD_DV","RMD_LR","RME_DV","RME_LR","RMF","RMG","RMH","RMD","RME"]:
    T(n,"Neurons","Head_motor",n)
T("SMD and subpopulation of RMD_LR and SAA","Neurons","Head_motor","SMD")

# VC neurons
for n in ["VC","VC_4_5"]:
    T(n,"Neurons","Cholinergic",n)

# Sensory neurons
for n in ["AFD","ASEL","ASER","ASE","ASG","ASH","ASI","ASJ","ASK","ASF","ADF","ADL",
          "AWA","AWB","AWC","AWC_ON","AWC_OFF","BAG","URX","AQR","PQR","URX_AQR_PQR",
          "FLP","PVD","PVD_FLP","PLM","ALM","AVM","PVM","PVM; AVM","OLQ","IL1","IL2",
          "IL2_DV","IL2_LR","URA","URB","URY","PHC","PHA","PHB","OLL",
          "Ciliated_sensory_neurons","Oxygen_sensory_neurons","Touch_receptor_neurons",
          "Touch_receptor","SDQ","SDQ and subpopulation of NSM","SDQ+NSM",
          "IL1 and subpopulation of OLL","PHA_PHB; LUA","URY and subpopulation of PHC"]:
    T(n,"Neurons","Sensory",n.split()[0].split(";")[0])

# Canonical fixes
T("ASE","Neurons","Sensory","ASE"); T("ASEL","Neurons","Sensory","ASEL"); T("ASER","Neurons","Sensory","ASER")
T("AWC_ON","Neurons","Sensory","AWC"); T("AWC_OFF","Neurons","Sensory","AWC")
T("URX_AQR_PQR","Neurons","Sensory","URX_AQR_PQR")
T("PVD_FLP","Neurons","Sensory","PVD_FLP")
T("PLM_ALM","Neurons","Sensory","PLM_ALM")
T("PVM; AVM","Neurons","Sensory","PVM_AVM")
T("IL1 and subpopulation of OLL","Neurons","Sensory","IL1")
T("SDQ and subpopulation of NSM","Neurons","Sensory","SDQ")
T("PHA_PHB; LUA","Neurons","Sensory","PHA_PHB")
T("URY and subpopulation of PHC","Neurons","Sensory","URY")
T("Canal_associated_neurons","Neurons","Sensory","CAN")

# Interneurons
for n in ["AIA","AIB","AIM","AIN","AIY","AIZ","RIA","RIB","RIC","RID","RIF","RIG",
          "RIH","RIR","RIM","RIP","RIS","RIV","AVA","AVB","AVD","AVE","AVF","AVG",
          "AVH","AVJ","AVK","DVA","DVC","PVN","PVP","PVQ","PVR","PVT","PVW",
          "BDU","ALA","AUA","CAN","ADA","ALN","LUA","XXX","PDA","PDB",
          "Other_interneurons","flp-1(+)_interneurons","flp-1(+)"]:
    T(n,"Neurons","Interneuron",n)
T("RIB and RIG","Neurons","Interneuron","RIB_RIG")
T("RIH_RIR","Neurons","Interneuron","RIH_RIR")
T("AVL and DVB","Neurons","Interneuron","AVL_DVB")

# Pharyngeal neurons
for n in ["I1","I2","I3","I4","I5","I6","M1","M2","M3","M4","M5","MC","MI","NSM",
          "M2; M3; M4; I2","I5; I4","M2+M3+M4+I2","Pharyngeal_neurons"]:
    T(n,"Neurons","Pharyngeal_neuron",n.replace("; ","_").replace("+","_"))
T("M2; M3; M4; I2","Neurons","Pharyngeal_neuron","M2_M3_M4_I2")
T("I5; I4","Neurons","Pharyngeal_neuron","I5_I4")

# HSN
T("HSN","Neurons","Motor","HSN")
# SDQ/ALN/PLN neuron cluster
T("SDQ/ALN/PLN","Neurons","Sensory","SDQ_ALN_PLN")
T("ALN","Neurons","Sensory","ALN")
T("PLN","Neurons","Sensory","PLN")
T("PLM","Neurons","Sensory","PLM"); T("ALM","Neurons","Sensory","ALM")

# ---- MUSCLE ----
for n in ["Body_wall_muscle","Body_wall_muscle_anterior","Body_wall_muscle anterior",
          "Body_wall_muscle posterior","Body_wall_muscle middle",
          "Body wall muscle anterior","Body wall muscle middle","Body wall muscle posterior"]:
    T(n,"Muscle","Body_wall_muscle","Body_wall_muscle")
for n in ["Vulval_muscle","vm1","vm2","vm1 (vulval muscle)","vm2 (vulval muscle)"]:
    T(n,"Muscle","Vulval_muscle",n.split()[0])
T("vm1 (vulval muscle)","Muscle","Vulval_muscle","vm1")
T("vm2 (vulval muscle)","Muscle","Vulval_muscle","vm2")
for n in ["Anal_muscle","Anal muscle"]:
    T(n,"Muscle","Anal_muscle","Anal_muscle")
for n in ["Intestinal/rectal_muscle","Intestinal_rectal_muscle"]:
    T(n,"Muscle","Intestinal_rectal_muscle","Intestinal_rectal_muscle")
for n in ["Uterine_muscle","Uterine muscle","Uterine muscle ","Unassigned sex-specific muscle"]:
    T(n,"Muscle","Uterine_muscle","Uterine_muscle")
T("Sex_myoblasts","Muscle","Sex_myoblast","Sex_myoblast")
T("hmc","Muscle","Head_mesodermal_cell","hmc")

# ---- INTESTINE ----
for n in ["Intestine","Intestine anterior","Intestine middle","Intestine posterior",
          "Intestinal/rectal_muscle",  # covered above
          "cat-4(+)ptps-1(+) intestine anterior"]:
    T(n,"Intestine","Intestine",n.replace(" ","_").replace("-","_").replace("(","").replace(")","").replace("+","plus"))
T("Intestine","Intestine","Intestine","Intestine")
T("Intestine anterior","Intestine","Intestine","Intestine_anterior")
T("Intestine middle","Intestine","Intestine","Intestine_middle")
T("Intestine posterior","Intestine","Intestine","Intestine_posterior")
T("cat-4(+)ptps-1(+) intestine anterior","Intestine","Intestine","Intestine_subpopulation")
T("Intestinal-rectal valve","Intestine","Intestine","Intestinal_rectal_valve")

# ---- HYPODERMIS ----
for n in ["Epidermis","Seam_cell","Seam cells","Seam cells (bus+) and excretory duct subpopulation",
          "Seam cells (grd+)","Uterine seam cells","hyp7","hyp4_to_hyp6",
          "Hypodermis head","hyp4_to_hyp6 and Tail hypodermis subpopulation",
          "Non-seam_hypodermis"]:
    T(n,"Hypodermis","Hypodermis","Hypodermis")
T("Epidermis","Hypodermis","Hypodermis","Epidermis_hyp7")  # CeNGEN "Epidermis" = hyp7
T("Seam_cell","Hypodermis","Seam_cell","Seam_cell")
T("Seam cells (bus+) and excretory duct subpopulation","Hypodermis","Seam_cell","Seam_cell")
T("Seam cells (grd+)","Hypodermis","Seam_cell","Seam_cell")
T("Non-seam_hypodermis","Hypodermis","Hypodermis","Hypodermis")
T("hyp7 (hypodermis)","Hypodermis","Hypodermis","hyp7")
T("hyp7","Hypodermis","Hypodermis","hyp7")
T("hyp4_to_hyp6","Hypodermis","Hypodermis","hyp4_to_hyp6")
T("hyp4_to_hyp6 and Tail hypodermis subpopulation","Hypodermis","Hypodermis","hyp4_to_hyp6")
T("Hypodermis head","Hypodermis","Hypodermis","Hypodermis_head")
T("Uterine seam cells","Hypodermis","Seam_cell","Uterine_seam_cell")

# ---- PHARYNX ----
for n in ["Pharyngeal_muscle","pm1_pm2","pm1 pm2","pm3_pm4_pm5","pm6_pm7",
          " pm3_pm4_pm5 and pm6_pm7 (pharyngeal muscle)",
          "pm1_pm2 (pharyngeal muscle)"]:
    T(n,"Pharynx","Pharyngeal_muscle","Pharyngeal_muscle")
T(" pm3_pm4_pm5 and pm6_pm7 (pharyngeal muscle)","Pharynx","Pharyngeal_muscle","pm3_to_pm7")
T("pm1_pm2 (pharyngeal muscle)","Pharynx","Pharyngeal_muscle","pm1_pm2")
T("Pharyngeal_muscle","Pharynx","Pharyngeal_muscle","Pharyngeal_muscle")
T("Pharyngeal_epithelia","Pharynx","Pharyngeal_epithelia","Pharyngeal_epithelia")
T("Pharyngeal_epithelia","Pharynx","Pharyngeal_epithelia","Pharyngeal_epithelia")
T("Pharyngeal epithelia","Pharynx","Pharyngeal_epithelia","Pharyngeal_epithelia")
T("e1_e3 (pharyngeal epithelium)","Pharynx","Pharyngeal_epithelia","e1_e3")
T("e2 (pharyngeal epithelium)","Pharynx","Pharyngeal_epithelia","e2")
T("Marginal_cell","Pharynx","Pharyngeal_epithelia","Marginal_cell")
T("Marginal cells","Pharynx","Pharyngeal_epithelia","Marginal_cell")
for n in ["Pharyngeal_gland_cell","g1A","g1P","g2","Pharyngeal_gland",
          "g1A (pharyngeal gland) ","g1P (pharyngeal gland)","g2 (pharyngeal gland) "]:
    T(n,"Pharynx","Pharyngeal_gland","Pharyngeal_gland")
T("g1A (pharyngeal gland) ","Pharynx","Pharyngeal_gland","g1A")
T("g1P (pharyngeal gland)","Pharynx","Pharyngeal_gland","g1P")
T("g2 (pharyngeal gland) ","Pharynx","Pharyngeal_gland","g2")
T("Pharyngeal_gland_cell","Pharynx","Pharyngeal_gland","Pharyngeal_gland")

# ---- GERMLINE ----
for n in ["Germline","Sperm","Oocytes","Differentiated germ","Meiotic germ cells",
          "Apoptotic germ cells","Mitotic germ cells","Mitotic germ cells ",
          "Spermatocytes","Mature sperm","Syncitial pachytene spermatocytes ",
          "Syncytial pachytene spermatocytes"]:
    T(n,"Germline","Germline","Germline")
T("Germline","Germline","Germline","Germline")
T("Sperm","Germline","Sperm","Sperm")
T("Oocytes","Germline","Oocyte","Oocyte")
T("Differentiated germ","Germline","Germline","Differentiated_germ")
T("Meiotic germ cells","Germline","Germline","Meiotic_germ")
T("Apoptotic germ cells","Germline","Germline","Apoptotic_germ")
T("Mitotic germ cells","Germline","Germline","Mitotic_germ")
T("Mitotic germ cells ","Germline","Germline","Mitotic_germ")
T("Spermatocytes","Germline","Sperm","Spermatocytes")
T("Mature sperm","Germline","Sperm","Mature_sperm")
T("Syncitial pachytene spermatocytes ","Germline","Sperm","Pachytene_spermatocytes")
T("Seminal vesicle (male)","Germline","Sperm","Seminal_vesicle")

# ---- GLIA ----
for n in ["Glia_1","Glia_2","Glia_3","Glia_4","Glia_5","AMsh","CEPsh","PHsh","GLR",
          "Amphid and phasmid sheath","Am/PH_sheath_cells","Cephalic sheath"]:
    T(n,"Glia","Glia",n)
T("AMsh","Glia","Amphid_sheath","AMsh")
T("CEPsh","Glia","Cephalic_sheath","CEPsh")
T("PHsh","Glia","Phasmid_sheath","PHsh")
T("GLR","Glia","GLR","GLR")
T("Amphid and phasmid sheath","Glia","Sheath_cell","Amphid_phasmid_sheath")
T("Am/PH_sheath_cells","Glia","Sheath_cell","Amphid_phasmid_sheath")
T("Cephalic sheath","Glia","Cephalic_sheath","Cephalic_sheath")

# Socket cells
for n in ["AMso","PHso","Arcade_cell","Amphid socket","Phasmid socket",
          "Cephalic and inner labial socket","Socket_cells"]:
    T(n,"Glia","Socket_cell",n)
T("AMso","Glia","Socket_cell","AMso")
T("PHso","Glia","Socket_cell","PHso")
T("Arcade_cell","Glia","Socket_cell","Arcade_cell")
T("Arcade cells and rectal gland subpopulation","Glia","Socket_cell","Arcade_cell")
T("Amphid socket","Glia","Socket_cell","Amphid_socket")
T("Phasmid socket","Glia","Socket_cell","Phasmid_socket")
T("Cephalic and inner labial socket","Glia","Socket_cell","Cephalic_socket")
T("Socket_cells","Glia","Socket_cell","Socket_cell")

# ---- GONAD / REPRODUCTIVE ----
T("Gonadal_sheath_cell","Gonad","Gonadal_sheath","Gonadal_sheath")
T("Gonad","Gonad","Gonad","Gonad")
T("sh1 (gonadal sheath distal)","Gonad","Gonadal_sheath","sh1")
T("sh2 (gonadal sheath distal)","Gonad","Gonadal_sheath","sh2")
T("sh3_sh4 (gonadal sheath proximal)","Gonad","Gonadal_sheath","sh3_sh4")
T("sh5 (gonadal sheath proximal)","Gonad","Gonadal_sheath","sh5")
T("Unassigned sheath cells","Gonad","Gonadal_sheath","Gonadal_sheath")
T("Distal_tip_cell","Gonad","Distal_tip","Distal_tip")
T("Distal tip","Gonad","Distal_tip","Distal_tip")
T("Distal_tip_cells","Gonad","Distal_tip","Distal_tip")
for n in ["Uterine_cell","Dorsal uterine cell","Unassigned uterine cells",
          "Uterine toroid","Uterine-vulval cells"]:
    T(n,"Gonad","Uterine","Uterine_cell")
T("Uterine_cell","Gonad","Uterine","Uterine_cell")
T("Dorsal uterine cell","Gonad","Uterine","Uterine_cell")
T("Uterine toroid","Gonad","Uterine","Uterine_toroid")
T("Uterine-vulval cells","Gonad","Uterine","Uterine_vulval")
T("Spermatheca","Gonad","Spermatheca","Spermatheca")
T("Spermathecal-uterine_junction_or_uterine_toroid","Gonad","Spermatheca","Spermathecal_uterine_junction")
T("Spermatheca bag distal","Gonad","Spermatheca","Spermatheca_bag_distal")
T("Spermatheca bag proximal","Gonad","Spermatheca","Spermatheca_bag_proximal")
T("Spermatheca neck distal","Gonad","Spermatheca","Spermatheca_neck_distal")
T("Spermatheca neck distal-most","Gonad","Spermatheca","Spermatheca_neck_distalmost")
T("Spermatheca-Uterine junction","Gonad","Spermatheca","Spermathecal_uterine_junction")
T("Somatic_gonad_precursors","Gonad","Gonad","Somatic_gonad_precursor")
T("Vulval_cells","Gonad","Vulval","Vulval_cell")
T("Vulval cells","Gonad","Vulval","Vulval_cell")
T("Vulval_precursors","Gonad","Vulval","Vulval_precursor")
T("VC_4_5","Neurons","Cholinergic","VC_4_5")
T("Unassigned sex-specific muscle","Gonad","Vulval","Unassigned_sex_muscle")

# ---- COELOMOCYTES ----
T("Coelomocyte","Coelomocytes","Coelomocytes","Coelomocytes")
T("Coelomocytes","Coelomocytes","Coelomocytes","Coelomocytes")

# ---- EXCRETORY ----
T("Excretory_cell","Excretory","Excretory","Excretory_cell")
T("Excretory cells","Excretory","Excretory","Excretory_cell")
T("Excretory gland","Excretory","Excretory_gland","Excretory_gland")
T("Excretory_gland_cell","Excretory","Excretory_gland","Excretory_gland")
T("Rectal_gland","Excretory","Rectal_gland","Rectal_gland")
T("Excretory_cells","Excretory","Excretory","Excretory_cell")

# ---- PHARYNX (extra) ----
T("Pharynx","Pharynx","Pharynx","Pharynx")
T("Pharyngeal_epithelia","Pharynx","Pharyngeal_epithelia","Pharyngeal_epithelia")

# Broad tissue names from Cao S6
for n in ["Body_wall_muscle","Glia","Gonad","Hypodermis","Intestine","Neurons","Pharynx"]:
    if n not in TAXONOMY:
        T(n, n if n != "Body_wall_muscle" else "Muscle",
          n if n != "Body_wall_muscle" else "Body_wall_muscle",
          n if n != "Body_wall_muscle" else "Body_wall_muscle")

def get_taxonomy(raw):
    """Return (tissue_type, subtype, canonical) for a raw cell type name."""
    if raw in TAXONOMY:
        return TAXONOMY[raw]
    # Fallback: clean up the name
    clean = raw.replace(" ","_").replace("/","_").replace("-","_")
    return ("Other", "Other", clean)

def canonical_col(raw):
    t, s, c = get_taxonomy(raw)
    return f"{t}.{s}.{c}"

# =========================================================================
# STEP 1: Process CeNGEN datasets (L4, L1, adult)
# =========================================================================
CENGEN_EXCLUDE = {'Unannotated','Unknown_non-neuronal_3','Unknown_non_neuronal',
                  'Unknown_non_neuronal_2','Unknown_non_neuronal_3','Embryonic cells'}

def process_cengen(csv_path, study_name):
    print(f"Processing {study_name}...")
    df = pd.read_csv(csv_path)
    df = df[~df['cell_type'].isin(CENGEN_EXCLUDE)].copy()

    # Compute gene-level mean across all cell types for specificity
    gene_mean = df.groupby('gene')['avg_expr'].mean().rename('gene_mean')
    df = df.merge(gene_mean, on='gene', how='left')
    df['specificity'] = df['avg_expr'] / (df['gene_mean'].clip(lower=1))

    cell_types = sorted(df['cell_type'].unique())

    # For each cell type: top markers by specificity, filtered by pct_expr
    # Score 2: pct_expr >= 50 AND specificity top 15
    # Score 1: pct_expr >= 25 AND specificity top 30 (not already score 2)
    rows = {}  # gene → {cell_type: score}

    for ct, grp in df.groupby('cell_type'):
        grp_sorted = grp.sort_values('specificity', ascending=False)
        high = grp_sorted[grp_sorted['pct_expr'] >= 50].head(15)['gene'].tolist()
        low  = grp_sorted[grp_sorted['pct_expr'] >= 25].head(30)['gene'].tolist()
        low  = [g for g in low if g not in high]
        for g in high:
            rows.setdefault(g, {})[ct] = 2
        for g in low:
            rows.setdefault(g, {})[ct] = rows.get(g, {}).get(ct, 1)

    # Build matrix
    mat = pd.DataFrame(0, index=sorted(rows.keys()), columns=cell_types, dtype=np.int8)
    for gene, ct_scores in rows.items():
        for ct, score in ct_scores.items():
            if ct in mat.columns:
                mat.loc[gene, ct] = score

    # Drop all-zero columns
    mat = mat.loc[:, mat.max() > 0]
    # Drop all-zero rows
    mat = mat.loc[mat.max(axis=1) > 0]

    print(f"  {study_name}: {mat.shape[0]} genes × {mat.shape[1]} cell types")
    return mat, cell_types

mat_L4, ct_L4 = process_cengen(SCRATCH/"L4_herm.csv",    "CeNGEN_L4")
mat_L1, ct_L1 = process_cengen(SCRATCH/"L1_herm.csv",    "CeNGEN_L1")
mat_ad, ct_ad = process_cengen(SCRATCH/"adult_herm.csv", "CeNGEN_adult")

# =========================================================================
# STEP 2: Process Cao 2017 (S6 tissue + S7 cell type + S8 neuron clusters)
# =========================================================================
print("Processing Cao 2017...")
cao_file = UPLOADS / "f20d7413-aam8940_cao_sm_tables_s1_to_s14.xlsx"

def score_cao(ratio, pval):
    if pval > 0.05: return 0
    if ratio >= 10: return 2
    if ratio >= 3:  return 1
    return 0

# S6: tissue-enriched markers
s6 = pd.read_excel(cao_file, sheet_name="Table S6", header=1)
s6['score'] = s6.apply(lambda r: score_cao(r['ratio'], r['pval']), axis=1)
s6 = s6[s6['score'] > 0]
cao_tissue_cts = sorted(s6['max.tissue'].unique())

# S7: cell-type enriched
s7 = pd.read_excel(cao_file, sheet_name="Table S7", header=1)
s7['score'] = s7.apply(lambda r: score_cao(r['ratio'], r['pval']), axis=1)
s7 = s7[s7['score'] > 0]
cao_celltype_cts = sorted(s7['max.cell.type'].unique())

# S8: neuron clusters
s8 = pd.read_excel(cao_file, sheet_name="Table S8", header=1)
s8['score'] = s8.apply(lambda r: score_cao(r['ratio'], r['pval']), axis=1)
s8 = s8[s8['score'] > 0]
cao_cluster_cts = sorted(s8['max.cluster'].unique())

# Combine all Cao cell types
cao_cts = list(dict.fromkeys(cao_tissue_cts + cao_celltype_cts + cao_cluster_cts))

# Build Cao matrix
all_cao_genes = sorted(set(s6['symbol'].tolist() + s7['symbol'].tolist() + s8['symbol'].tolist()))
mat_cao = pd.DataFrame(0, index=all_cao_genes, columns=cao_cts, dtype=np.int8)

for _, row in s6.iterrows():
    g = row['symbol']; ct = row['max.tissue']; sc = int(row['score'])
    if g in mat_cao.index: mat_cao.loc[g, ct] = max(mat_cao.loc[g, ct], sc)

for _, row in s7.iterrows():
    g = row['symbol']; ct = row['max.cell.type']; sc = int(row['score'])
    if g in mat_cao.index and ct in mat_cao.columns:
        mat_cao.loc[g, ct] = max(mat_cao.loc[g, ct], sc)

for _, row in s8.iterrows():
    g = row['symbol']; ct = row['max.cluster']; sc = int(row['score'])
    if g in mat_cao.index and ct in mat_cao.columns:
        mat_cao.loc[g, ct] = max(mat_cao.loc[g, ct], sc)

mat_cao = mat_cao.loc[mat_cao.max(axis=1) > 0, mat_cao.max() > 0]
print(f"  Cao2017: {mat_cao.shape[0]} genes × {mat_cao.shape[1]} cell types")

# =========================================================================
# STEP 3: Process Ghaddar 2023
# =========================================================================
print("Processing Ghaddar 2023...")
ghad_file = UPLOADS / "cd540a3a-adg0506_Data_S1.xlsx"
g = pd.read_excel(ghad_file)
g = g[g['Final annotation'].notna()].copy()

ghad_cts = sorted(g['Final annotation'].unique())
all_ghad_genes = sorted(g['gene_short_name'].unique())
mat_ghad = pd.DataFrame(0, index=all_ghad_genes, columns=ghad_cts, dtype=np.int8)

for ct, grp in g.groupby('Final annotation'):
    grp_sorted = grp.sort_values('marker_score', ascending=False).reset_index(drop=True)
    for i, row in grp_sorted.iterrows():
        gene = row['gene_short_name']
        score = 2 if i < 5 else 1   # top 5 → score 2; rest → score 1
        if gene in mat_ghad.index:
            mat_ghad.loc[gene, ct] = max(mat_ghad.loc[gene, ct], score)

mat_ghad = mat_ghad.loc[mat_ghad.max(axis=1) > 0, mat_ghad.max() > 0]
print(f"  Ghaddar2023: {mat_ghad.shape[0]} genes × {mat_ghad.shape[1]} cell types")

# =========================================================================
# STEP 4: Save per-study TSV files (raw column names, values 0/1/2)
# =========================================================================
print("\nSaving per-study TSV files...")

studies = {
    "cao2017_L2":          mat_cao,
    "cengen_L4_herm":      mat_L4,
    "cengen_L1_herm":      mat_L1,
    "cengen_adult_herm":   mat_ad,
    "ghaddar2023_adult":   mat_ghad,
}

for name, mat in studies.items():
    out = OUTDIR / f"markers_{name}.tsv"
    mat.index.name = "gene"
    mat.to_csv(out, sep="\t")
    print(f"  Saved {out.name}: {mat.shape}")

# =========================================================================
# STEP 5: Build combined matrices using canonical taxonomy
# =========================================================================
print("\nBuilding combined (canonical) matrices...")

def add_canonical_cols(mat, study_name):
    """Return a long DataFrame: gene, canonical_col, score, study"""
    rows = []
    for ct in mat.columns:
        canon = canonical_col(ct)
        for gene, score in mat[ct].items():
            if score > 0:
                rows.append({"gene": gene, "cell_type": canon, "score": score, "study": study_name})
    return pd.DataFrame(rows)

long_dfs = []
for study, mat in [("cao2017_L2",mat_cao), ("cengen_L4",mat_L4),
                   ("cengen_L1",mat_L1), ("cengen_adult",mat_ad),
                   ("ghaddar2023",mat_ghad)]:
    long_dfs.append(add_canonical_cols(mat, study))

long_all = pd.concat(long_dfs, ignore_index=True)

# For each (gene, cell_type): count studies that support it, take max score
agg = (long_all.groupby(['gene','cell_type'])
       .agg(n_studies=('study','nunique'), max_score=('score','max'))
       .reset_index())

all_canon_cts = sorted(agg['cell_type'].unique())
all_canon_genes = sorted(agg['gene'].unique())

# --- Shared markers: supported by ≥ 2 studies ---
shared_entries = agg[agg['n_studies'] >= 2]
shared_cts = sorted(shared_entries['cell_type'].unique())
shared_genes = sorted(shared_entries['gene'].unique())
mat_shared = pd.DataFrame(0, index=shared_genes, columns=sorted(all_canon_cts), dtype=np.int8)
for _, row in shared_entries.iterrows():
    mat_shared.loc[row['gene'], row['cell_type']] = int(row['max_score'])
# Drop all-zero cols
shared_cts_present = mat_shared.columns[mat_shared.max() > 0].tolist()
mat_shared = mat_shared[shared_cts_present]
mat_shared = mat_shared.loc[mat_shared.max(axis=1) > 0]
print(f"  Shared (≥2 studies): {mat_shared.shape[0]} genes × {mat_shared.shape[1]} cell types")

# --- Expansive markers: all (≥1 study), but cell types with ≥2 studies get priority ---
# Use max_score for all entries; keep cell types covered by ≥1 study
mat_exp = pd.DataFrame(0, index=sorted(agg['gene'].unique()),
                       columns=sorted(agg['cell_type'].unique()), dtype=np.int8)
for _, row in agg.iterrows():
    mat_exp.loc[row['gene'], row['cell_type']] = int(row['max_score'])
mat_exp = mat_exp.loc[:, mat_exp.max() > 0]
mat_exp = mat_exp.loc[mat_exp.max(axis=1) > 0]
print(f"  Expansive (≥1 study): {mat_exp.shape[0]} genes × {mat_exp.shape[1]} cell types")

# Save combined files
mat_shared.index.name = "gene"
mat_exp.index.name = "gene"
(mat_shared.to_csv(OUTDIR / "markers_shared_min2studies.tsv", sep="\t"))
(mat_exp.to_csv(OUTDIR / "markers_expansive_all_studies.tsv", sep="\t"))
print(f"  Saved markers_shared_min2studies.tsv")
print(f"  Saved markers_expansive_all_studies.tsv")

# =========================================================================
# STEP 6: Save taxonomy reference table
# =========================================================================
print("\nSaving taxonomy reference...")
tax_rows = []
for raw, (tissue, subtype, canon) in TAXONOMY.items():
    tax_rows.append({"raw_name": raw, "tissue_type": tissue,
                     "tissue_subtype": subtype, "canonical_cell_type": canon,
                     "canonical_col": f"{tissue}.{subtype}.{canon}"})
tax_df = pd.DataFrame(tax_rows).drop_duplicates().sort_values(['tissue_type','tissue_subtype','canonical_cell_type'])
tax_df.to_csv(OUTDIR / "cell_type_taxonomy.tsv", sep="\t", index=False)
print(f"  Saved cell_type_taxonomy.tsv ({len(tax_df)} entries)")

# =========================================================================
# STEP 7: Summary statistics
# =========================================================================
print("\n=== SUMMARY ===")
for study, mat in studies.items():
    n2 = (mat == 2).sum().sum()
    n1 = (mat == 1).sum().sum()
    print(f"  {study}: {mat.shape[0]} genes, {mat.shape[1]} cell types, {n2} high markers, {n1} low markers")
print(f"  Shared:   {mat_shared.shape[0]} genes × {mat_shared.shape[1]} cell types")
print(f"  Expansive:{mat_exp.shape[0]} genes × {mat_exp.shape[1]} cell types")

# Show which cell types are in shared file
print("\nCell types in shared file:")
for ct in sorted(mat_shared.columns):
    n = (mat_shared[ct] > 0).sum()
    print(f"  {ct}: {n} marker genes")

print("\nDone!")
