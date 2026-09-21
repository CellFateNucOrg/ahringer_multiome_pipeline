import sys

# identical output to the paper's script; genes kept in a set instead of a list (the list lookup is extremely slow)
input_file=sys.argv[1]
gene_types=sys.argv[2]
output_file=sys.argv[3]
genes=sys.argv[4]

allowed_types = set(gene_types.split(","))

genes_to_retain=set()
for line in open(genes):
    if line.rstrip().split("\t")[1] in allowed_types:
        genes_to_retain.add(line.rstrip().split("\t")[0])

new=open(output_file, "w")
for line in open(input_file):
    if line[0] == "#":
        continue
    G=line.rstrip().split("\t")
    K = G[8].split('"')
    gene_name=K[1]
    if gene_name in genes_to_retain:
        new.write(line)

new.close()
