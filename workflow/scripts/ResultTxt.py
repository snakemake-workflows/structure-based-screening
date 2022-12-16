"""safes table with best screening results"""

from itertools import tee
from pathlib import Path
import os
import pandas as pd

infile  = snakemake.input[0]
outfile = snakemake.output[0]

linesep = os.linesep

def pairwise(iterable):
    """
    :param iterable:
    :return: a, b
    """
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

def extract_pdbqt(infile):
    """
    :param infile: path to infile
    :param outfile: path to outfile
    :return: receptor name, IDs and enthalpies
    """
    target_IDs  = []
    enthalpies  = []
    references  = set()
    with open(infile, encoding='utf-8') as to_parse:
        # this line reads: 'REMARK RECEPTOR path/to/target.pdbqt'
        headerline = to_parse.readline()
        # we want to extract the target name WITHOUT the suffix
        receptor   = Path(headerline.split(' ')[-1]).stem
        # we then proceed extracting the enthalpy per ligand
        # we iterate two line per iteration, because the file is structured, like:
        #
        # REMARK VINA RESULT:  Enthalpy 1 Enthalpy 2 Enthalpy 3
        # REMARK  Name = <name>
        for i, j in pairwise(to_parse):
            if 'Name' in j:
                ID = j.split(' = ')[-1].strip(linesep)
                if ID not in references:
                    references.add(ID)
                    target_IDs.append(ID)
                    enthalpies.append(i.split()[3])

    return receptor, target_IDs, enthalpies

if __name__ == '__main__':
    receptor, target_IDs, enthalpies = extract_pdbqt(infile)
    # putting into data frame:
    out_df = pd.DataFrame(enthalpies, index = target_IDs, columns = [receptor])
    # finally, save the output
    out_df.to_csv(outfile)
