"""merge all docking output files together"""

import shutil
import subprocess

outFile = snakemake.output[0]
inFiles = snakemake.input

# TODO: change to ZIPFILE to support better result portability
with open(outFile, "wb") as outF:
    for file in inFiles:
        try:
            out = subprocess.check_output(
                ["gunzip", "-t", file]
            )  # test if gzip file is valid before merging

            # TODO: avoid overwriting the file descriptor!
            with open(file, "rb") as sourceFile:
                shutil.copyfileobj(sourceFile, outF)
        except subprocess.CalledProcessError:
            print(f"Erros occured during the docking of {file}")
