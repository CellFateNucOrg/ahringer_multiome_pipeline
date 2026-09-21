"""Add common gene names to STARsolo features.tsv and store the unique+multimapper (EM) matrix in <raw>/um_reads.
Same output as the paper's STARsolo_output_parser.20231108.py, but re-runnable (no interactive gzip prompts)
and robust to genes without a common name."""
import gzip
import os
import shutil
import sys

sample_name = sys.argv[1]
gene_annotation = sys.argv[2]
wbid_gene_name = sys.argv[3]

raw_star_path = f"data/sc/rnaseq/{sample_name}/star_{gene_annotation}/star_{gene_annotation}_Solo.out/GeneFull/raw"
out_dir = os.path.join(raw_star_path, "um_reads")
os.makedirs(out_dir, exist_ok=True)

names = {}
for line in open(wbid_gene_name):
    G = line.rstrip("\n").split("\t")
    if len(G) >= 2:
        names[G[0]] = G[1]

with gzip.open(os.path.join(out_dir, "features.tsv.gz"), "wt") as new_feats:
    for line in open(os.path.join(raw_star_path, "features.tsv")):
        G = line.rstrip("\n").split("\t")
        G[1] = names.get(G[1], G[1])
        new_feats.write("\t".join(G) + "\n")

for src, dst in (("barcodes.tsv", "barcodes.tsv.gz"), ("UniqueAndMult-EM.mtx", "matrix.mtx.gz")):
    with open(os.path.join(raw_star_path, src), "rb") as fin, gzip.open(os.path.join(out_dir, dst), "wb") as fout:
        shutil.copyfileobj(fin, fout)
