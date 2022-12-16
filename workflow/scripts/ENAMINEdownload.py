"""download from ENAMINE database"""

import os
import requests

for collection in snakemake.config["ENAMINE_INPUT"]:
    url = os.path.join(snakemake.config["ENAMINE_URL"], collection, ".zip")
    folder = os.path.join(
        snakemake.config["INPUT_DIR"], "ENAMINE", os.path.dirname(collection)
    )
    r = requests.get(url, allow_redirects=True, timout=60)
    with open(
        os.path.join(snakemake.config["INPUT_DIR"], "ENAMINE", collection), "wb"
    ) as outfile:
        outfile.write(r.content)
    # hashed = (folder+".sh2")
    # shell('sha256sum {folder} > {hashed}')
