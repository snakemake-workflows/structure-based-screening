# An HPC conformant Structure Based Screening Workflow

This document descibes a Workflow for screening massive amounts of small molecules (ligands) against a target (PDB-File)


## Optaining the Audodock Tool Suite

In order to prepare the workflow run, users need to define a so-called 'grid file' to describe the region of a protein target of interest (e.g. a receptor binding domain or an allosteric binding domain). 

To prepare the grid file we recommend using the Autodock Tool Suite (ADT), which can be obtained via the [MGLTools download site](https://ccsb.scripps.edu/mgltools/downloads/). Please download a recent version of the tool suite for your operating system, unpack and install according to the instructions.

Then start ADT. 


Note that, the workflow is designed to perform direct docking respectively screening. Here, a putative binding domain should not exceed 9.000 Å³ as the workflow employs VinaLC (xxx link xxx), which in turn is based on Vina (xxx link xxx). These programms document this recommandation in their documentation (xxx link xxx). To run 'blind docking' experiments, consult the documentation, too.
