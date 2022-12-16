"""
preparation of ligand input file for VinaLC containing the path
 to every ligand with same weigth + logP
 """
import os

input_directory = snakemake.input.in_dir
output_directory = snakemake.config["OUTPUT_DIR"]

if "ZINC" in snakemake.config["DATABASE"]:
    weightLog = input_directory[-2:]
    for weight in snakemake.config["ZINC_INPUT"]["WEIGHT"]:
        for log in snakemake.config["ZINC_INPUT"]["LOGP"]:
            for react in snakemake.config["ZINC_INPUT"]["REACT"]:
                for purchase in snakemake.config["ZINC_INPUT"]["PURCHASE"]:
                    for ph in snakemake.config["ZINC_INPUT"]["PH"]:
                        for charge in snakemake.config["ZINC_INPUT"]["CHARGE"]:
                            fname = weightlog + react + purchase + ph + charge + ".pdbqt"
                            ligand_file = os.path.join(input_directory, fname)
                            with open(os.path.join(snakemake.output.library ),"r", encoding='utf-8') as file_object:
                                file_object.write(ligand_file + "\n")
