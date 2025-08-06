"""download ligand subsets in mol2 format from ZINC database"""

import requests

out = snakemake.output[0]
subset = snakemake.params.sub

url = "".join(
    ["https://zinc15.docking.org/substances/subsets/", subset, ".mol2?count=all"]
)

r = requests.get(url, timeout=60)
with open(out, "wb") as outfile:
    outfile.write(r.content)
