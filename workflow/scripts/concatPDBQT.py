"""read all files in input folder and concat them together"""
import os
from snakemake.shell import shell

input_directory = snakemake.params.in_dir
out_file = snakemake.output[0]

if not os.path.exists(os.path.join(input_directory, "/convertedPDBQT")):
    os.makedirs(os.path.join(input_directory, "/convertedPDBQT"))

with open(out_file, "a", encoding="utf-8") as out:
    for filename in os.listdir(input_directory):
        if os.path.isfile(os.path.join(input_directory, filename)):
            if filename.endswith(".pdbqt"):
                with open(
                    os.path.join(input_directory, filename), "r", encoding="utf-8"
                ) as inFile:
                    (name, ext) = os.path.splitext(filename)
                    out.write("MODEL\n")
                    out.write("REMARK Name = " + name + "\n")
                    out.write(inFile.read())
                    out.write("ENDMDL\n")
            else:
                newFilename = os.path.join(
                    input_directory, "/convertedPDBQT/", filename, ".pdbqt"
                )
                in_babel = os.path.join(input_directory, "/", filename)
                shell("obabel -isdf {in_babel} -opdbqt -O {newFilename}")
                with open(newFilename, encoding="utf-8") as newFile:
                    out.write(newFile.read())
with open(out_file, "r", encoding="utf-8") as file:
    counter = 0
    for line in file:
        if "ATOM" in line:
            xyz = line.split()[5:8]
            if xyz == ["0.000", "0.000", "0.000"]:
                counter += 1
        if "ENDMDL" in line:
            counter = 0
        if counter >= 2:
            print("Please provide 3d coordinates in input files")
            sys.exit("Please provide 3d coordinates in input files")
            break
