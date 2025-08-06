"""preparation of proteins specified in configfile"""

import os
from Bio.PDB import PDBParser, PDBIO


def removeChains(model, chainlist):
    """
    removes chains not specified in chainlist from model

    Parameters
    ----------
    model :
        Description of parameter `model`.
    chainlist : string
        Description of parameter `chainlist`.

    """

    residue_to_remove = []
    chain_to_remove = []
    for chain in model:
        for residue in chain:
            if residue.id[0] != " ":
                residue_to_remove.append((chain.id, residue.id))
        if not chain:
            chain_to_remove.append(chain.id)
    for chain in model:
        if chain.get_id() not in chainlist:
            chain_to_remove.append(chain.get_id())

    for residue in residue_to_remove:
        model[residue[0]].detach_child(residue[1])
    for chain in chain_to_remove:
        model.detach_child(chain)


def prepareRec(inputfile, outputfile, target):
    """
    select chains to delete depending on config definition
    """
    print(target)
    ID = target.split(",")
    chains = ID[1].split(" ")
    parser = PDBParser()  # MMCIFParser()
    structure = parser.get_structure(ID[0], inputfile)
    model = structure[0]
    removeChains(model, chains)

    io = PDBIO()
    io.set_structure(structure)
    out = outputfile
    print("printing outfile")
    io.save(out)


head, tail = os.path.split(snakemake.input[0])
filename = tail.split(".")[0]

if any(filename in target for target in snakemake.config["TARGETS"]):
    print("filename in targets")
    prepareRec(snakemake.input[0], snakemake.output[0], snakemake.config["TARGETS"][0])

elif filename in str(snakemake.config["RESCREENING_TARGETS"]):
    for target in snakemake.config["RESCREENING_TARGETS"]:
        if filename in target:
            prepareRec(snakemake.input[0], snakemake.output[0], target)
