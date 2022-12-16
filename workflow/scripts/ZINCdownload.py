"""download ligands in pdbqt format from ZINC database"""
import os
import requests

out = str(snakemake.output)
weight_log = out.split("/")[-2]
attr = out.rsplit("/", 1)[-1]

file_name = "".join([attr[:6], ".xaa.pdbqt.gz"])
url = os.path.join("http://files.docking.org/3D", weight_log, attr[2:6], file_name)

folder = os.path.join(
    snakemake.config["INPUT_DIR"], "ZINC", weight_log, "".join([attr[:6], ".pdbqt.gz"])
)

r = requests.get(url, timeout=60)
with open(folder, "wb") as outfile:
    outfile.write(r.content)

hashed = folder + ".sh2"
shell: "sha256sum {folder} > {hashed}"
