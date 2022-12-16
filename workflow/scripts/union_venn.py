''''safes table with union ligands of all rescreening results'''

from pathlib import Path
import itertools
import os
import venn
import pandas as pd

linesep = os.linesep

rescreen_results = snakemake.input.re_results

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b= itertools.tee(iterable)
    next(b, None)
    return zip(a, b)

def group3(iterable):
    '''iterates over 3 lines at once'''
    a, b, c = itertools.tee(iterable, 3)
    next(b, None)
    next(c, None)
    next(c, None)
    return zip(a, b, c)

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

def extract_pdbqt_rescreening(infile):
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
        for a, _, j in group3(to_parse):
            if 'Name' in j:
                ID = j.split(' = ')[-1].strip(linesep)
                if ID not in references:
                    references.add(ID)
                    target_IDs.append(ID)
                    enthalpies.append(a.split()[3])

    return receptor, target_IDs, enthalpies

def unionVenn(bestLigandFile, outFile, reBest):
    '''outputs union table of '''
    ID = []
    value = []
    receptor, ID, value = extract_pdbqt(bestLigandFile)

    out_df = pd.DataFrame(value, index = ID, columns = [receptor])

    inputArray = []

    for inFile in reBest:
        inputArray.append(inFile)

    for i in inputArray: #join all rescreening results to first screening
        receptor, ID, value = extract_pdbqt_rescreening(i)
        df = pd.DataFrame(value, index = ID, columns = [receptor])
        out_df = out_df.join(df, how = 'inner')
    out_df.to_csv(outFile)

def getLigands(result_file):
    """ take results file and outputs ligand names and receptor name """
    ligandNames = []
    with open(result_file, encoding = 'utf-8') as result:
        for line in result:
            if 'REMARK RECEPTOR' in line:
                receptor = line.split('/')[-1]
            if 'Name' in line:
                ligandNames.append(line.split()[-1])
    return ligandNames, receptor
def makeVenn(result_list):
    """takes rescreening results and outputs a venn diagram"""
    ligand_lists = []
    receptor_list = []
    for result in result_list:
        ligands,receptor = getLigands(result)
        ligand_lists.append(ligands)
        receptor_list.append(receptor)
    venn.venn(ligand_lists, snakemake.output[1], receptor_list, figsize=(16,16))

#rescreen_results = str(rescreen_results).split()


if 1 < len(rescreen_results) <= 4:
    makeVenn(rescreen_results)
else:
    open(snakemake.output[1], 'a', encoding = 'utf-8').close()
unionVenn( snakemake.input[0], snakemake.output[0], rescreen_results)
