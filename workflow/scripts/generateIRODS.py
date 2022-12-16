"""Generating iRODS template file"""

import json
import pandas as pd

target1 =  snakemake.config["TARGETS"]
targets = snakemake.config["RESCREENING_TARGETS"]

name = snakemake.config["EXPERIMENT_NAME"]

results = snakemake.input[0]
outfile = snakemake.output[0]

target_dict = {}

target_dict["primary_targets"] = target1
target_dict["rescreening_targets"] = targets

out_dict = {
    "@context": "http://schema.org",
    "@id": "https://something",
    "@type": "ScreeningResults",
    "dct:conformsTo": "not defined yet",
    "name": "",
    "url": "",
    "description": "<enter your discription here>",
}

out_dict["name"] = name

out_dict.update(target_dict)

union = pd.read_csv(results)
union = union[0:10].to_dict()

results = {"results": union}

out_dict.update(results)

with open(outfile, "w", encoding = "utf-8") as f:
    json.dump(out_dict, f, sort_keys = False, indent = 4)
