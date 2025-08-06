"""
    split ligand file into smaller files with max 15000 ligands each
"""

import os

# max. input molecule per output file
n = 15000

with open(snakemake.input[0], "r", encoding="utf-8") as infile:
    intext = infile.read()  # read the entire file into memory

nligands = intext.split("ENDMDL")
# number of output files
file_number = (len(nligands) // n) + 1

for j in range(file_number):
    if not os.path.exists(snakemake.output[0]):
        os.makedirs(snakemake.output[0], exist_ok=True)
    with open(
        os.path.join(snakemake.output[0], "".join([str(j), ".pdbqt"])),
        "w",
        encoding="utf-8",
    ) as outfile:
        count = 0
        for i in nligands:
            if n * j <= count < n * (j + 1):
                outfile.write(i)
                outfile.write("ENDMDL")
            count += 1
