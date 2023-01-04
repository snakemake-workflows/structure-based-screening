import glob
import re
import os
import requests
from snakemake.logging import logger


def getTranches():
    '''return traches from parsing last log file'''
    log_files = glob.glob('.snakemake/log/*')
    sorted_files = sorted(log_files, key=os.path.getmtime)
    tranch_list = []
    pattern =  '[A-K][A-K][ABCEGI][ABCDEF][RMLH][NMLOP]'
    with open(sorted_files[-2]) as log:
        matches = re.findall(pattern, log.read())
        return list(set(matches))

def library(wildcards):
    if DATABASE[0]=="ZINC": # ZINC database selected

        if SUBSET=="TRANCHES": # Tranches selected
            out = []
            r_zinc = requests.get("http://files.docking.org/3D/")
            if True: #r_zinc.status_code != 200:  #test if ZINC database is available
                logger.info("The ZINC database is not accessible right now. Perhaps it is temporarily down?")
                user_input = input("Have you already run this workflow in the current folder with the same input data?(y/n) \n")
                if user_input == "y":
                    logger.info("Trying to proceed without accessing the ZINC database (http://files.docking.org)")
                    tranch_list = getTranches() #get list with 6 letter code specifing ZINC tranches

                    receptorID = config["TARGETS"][0].split(',')[0]
                    database = config["DATABASE"]

                    for i in tranch_list:
                        wl = i[0:2]
                        entry = f"{OUTPUT_DIR}/output/{receptorID}/{receptorID}_{database}_{wl}_{i}.pdbqt.gz"
                        out.append(entry)
                    if not out:
                        logger.error("No tranche parameter found in last log file.")
                        sys.exit(1)
                    else:
                        return out
                else:
                    logger.error("Aborting the snakemake run for now, as the data from the ZINC database (http://files.docking.org) are temporarily not available. Please try at a later time.")
                    sys.exit(1)
            rawOut= expand(path.join(OUTPUT_DIR, "output", "{receptorID}",
                "{receptorID}_{database}_{w}{l}_{w}{l}{r}{p}{ph}{c}.pdbqt.gz"),
                receptorID = config["TARGETS"][0].split(',')[0],
                database = config["DATABASE"],
                w  = config["ZINC_INPUT"]["WEIGHT"],
                l  = config["ZINC_INPUT"]["LOGP"],
                r  = config["ZINC_INPUT"]["REACT"],
                p  = config["ZINC_INPUT"]["PURCHASE"],
                ph = config["ZINC_INPUT"]["PH"],
                c  = config["ZINC_INPUT"]["CHARGE"],
                dataset = (config["ZINC_INPUT"]["WEIGHT"]) + (config["ZINC_INPUT"]["LOGP"]))
            for i in rawOut:
                weighLog = i.split("_")[-2]
                restAttr = (i.split("_")[-1])[2:6]
                url = "".join(['http://files.docking.org/3D/', weighLog,'/', restAttr, '/', weighLog,restAttr,'.xaa.pdbqt.gz'])
                r = requests.get(url)
                if r.status_code == 200:
                    out.append(i)
            if not out:
                logger.error("All selected tranches are empty; select other parameters!")
                sys.exit(1)
            else:
                return out
        else: # if not tranches select subset

            out = expand(path.join(OUTPUT_DIR, "output", "{receptorID}",
                "{receptorID}_{database}_subsets_{subset}.pdbqt.gz"),
                  database = config["DATABASE"],
                  subset = config["SUBSET"],
                  receptorID = config["TARGETS"][0].split(',')[0])

            url = "https://zinc15.docking.org/substances/subsets/"+SUBSET+".mol2?count=all"
            r = requests.get(url)
            if r.status_code == 200: #test if subset is valid
                return out

            r_zinc = requests.get("https://zinc15.docking.org/")
            if r_zinc.status_code != 200: #test if ZINC database is available
                logger.info("The ZINC database is not accessible right now. Perhaps it is temporarily down?")
                #if ZINC not available, but dataset is already downloaded --> continue
                subset_dir = path.join(INPUT_DIR,config["SUBSET"]+'.mol2')
                if os.path.isfile(subset_dir):
                    return out
                else:
                    logger.error("Subset is not availiable in the specified data folder. \n Abort snakemake run, try again later")
                    sys.exit(1)

            else:
                logger.error( "Invalid subset name!")
                sys.exit(1)

    else: # not ZINC database --> local input data
        best = expand(path.join(OUTPUT_DIR, "output", "{receptorID}",
              "{receptorID}_{database}_{dataset}_local.pdbqt.gz"),
              database = config["DATABASE"],
              dataset  = config["LOC_DATA"],
              receptorID = config["TARGETS"][0].split(',')[0])
        return best

rule makeHistogram:
    input:
        path.join(OUTPUT_DIR, "results","{receptorID}.pdbqt.gz")
    output:
        report(path.join(OUTPUT_DIR,"results","{receptorID}_hist.png"), category="Histogram")
    envmodules:
        config["PYPLOT"]
    script:
        "../scripts/makeHistogram.py"

rule bestLigands:
    input:
        library
    output:
        path.join(OUTPUT_DIR, "results","{receptorID}.pdbqt.gz")
    script:
        "../scripts/mergeOutput.py"

rule dockingResults:
    input:
        path.join(OUTPUT_DIR, "results","{receptorID}.pdbqt.gz")
    output:
        path.join(OUTPUT_DIR, "results", "{receptorID}_{percentage}.pdbqt")
    envmodules:
        config["PYTHON"]
    threads: config["DOCKING_RESULTS"]["threads"]
    resources:
        account = config["ACCOUNT"],
        resources("DOCKING_RESULTS")
    script:
        "../scripts/sortResult.py"

rule dockingResultsTxt:
    input:
        path.join(OUTPUT_DIR, "results", "{receptorID}_{percentage}.pdbqt")
    output:
        path.join(OUTPUT_DIR, "results", "{receptorID}_{percentage}.csv")
    wildcard_constraints:
        receptorID="[^/]+",
        percentage="[^/]+"
    script:
        "../scripts/ResultTxt.py"

rule removeDuplicateLigands:
    input:
        path.join(OUTPUT_DIR, "results", "{receptorID}_{percentage}.pdbqt")
    output:
        path.join(OUTPUT_DIR, "rescreening", "unique",  "{receptorID}_{percentage}.pdbqt")
    shell:
        "sed '/MODEL [2-9]/,/ENDMDL/d' {input} > {output}"

checkpoint split2:
    input:
        path.join(OUTPUT_DIR, "rescreening", "unique",  "{receptorID}_{percentage}.pdbqt")
    output:
        directory(os.path.join(TMP_DIR, "rescreening_ligands_{percentage}", "{receptorID}"))
    script:
        "../scripts/splitFile.py"

rule prepareLigands2:
    input:
        ligands = path.join(TMP_DIR, "rescreening_ligands_{percentage}", "{receptorID}","{i}.pdbqt")
    output:
        ligands = path.join(OUTPUT_DIR, "rescreening_{percentage}", "{name}_{receptorID}", "{i}.txt")
    shell:
        "echo {input.ligands} > {output.ligands}"

rule prepareSecondDocking:
    input:
        grid = path.join(OUTPUT_DIR,"grid","{name}_grid.txt"),
        receptor = path.join(PREPARED_DIR,"receptor", "{name}.pdbqt")
    output:
        grid = path.join(OUTPUT_DIR, "rescreening_{percentage}", "{name}_{receptorID}", "{name}.grd"),
        receptor = path.join(OUTPUT_DIR, "rescreening_{percentage}", "{name}_{receptorID}", "{name}.rec")
    shell:
        """
        cp {input.grid} {output.grid}
        echo {input.receptor} > {output.receptor}
        """

rule docking2:
    input:
        ligands = path.join(OUTPUT_DIR, "rescreening_{percentage}", "{name}_{receptorID}", "{i}.txt"),
        grid = path.join(OUTPUT_DIR, "rescreening_{percentage}", "{name}_{receptorID}", "{name}.grd"),
        receptor = path.join(OUTPUT_DIR, "rescreening_{percentage}", "{name}_{receptorID}", "{name}.rec")
    output:
        path.join(OUTPUT_DIR, "rescreening_{percentage}", "{name}_{receptorID}","{name}.rec_{i}.txt.pdbqt.gz")
    params:
        dir = path.join(OUTPUT_DIR, "rescreening_{percentage}"),
        gridfile = path.join(config["GRID_DIR"], "{name}.gpf"),
        cutOff = config["CUTOFF_VALUE"]
    resources:
        account = config["ACCOUNT"],
        mpi = True,
        resources("DOCKING_RESULTS")

    envmodules:
        config["VINALC"]
    shell:
        """(
        cd {params.dir}/{wildcards.name}_{wildcards.receptorID}
        space=$(grep 'spacing' {params.gridfile} | cut -f2 -d' ')
        #TODO: remove the srun statement once snakemake is supporting SLURM
        srun --cpu-bind=rank vinalc --recList {wildcards.name}.rec --ligList {wildcards.i}.txt --geoList {wildcards.name}.grd --granularity $space --useScoreCF --scoreCF {params.cutOff}
        cd -
        )"""

def aggregate_in2(wildcards):
    checkpoint_output = checkpoints.split2.get(**wildcards).output[0]
    files_names = expand(path.join(OUTPUT_DIR, "rescreening_{{percentage}}", "{{name}}_{{receptorID}}","{{name}}.rec_{i}.txt.pdbqt.gz"),
        i = glob_wildcards(os.path.join(checkpoint_output, '{i}.pdbqt')).i)
    return files_names

rule mergeDocking2:
    input:
        unpack(aggregate_in2)
    output:
        path.join(OUTPUT_DIR, "output", "rescreening_{percentage}", "{name}_{receptorID}","{name}.pdbqt.gz")
    shell:
        '''
            cat {input} >> {output}
        '''

rule dockingResults2:
    input:
        path.join(OUTPUT_DIR, "output", "rescreening_{percentage}","{name}_{receptorID}","{name}.pdbqt.gz")
    output:
        path.join(OUTPUT_DIR,"output","rescreening_{percentage}","{name}_{receptorID}","{name}_best.pdbqt")
    envmodules:
        config["PYTHON"]
    script:
        "../scripts/sortResult.py"

rule makeVenn:
    input:
        best = path.join(OUTPUT_DIR, "results", "{receptorID}_{percentage}.pdbqt"),
        re_results = expand(path.join(OUTPUT_DIR,"output","rescreening_{percentage}","{name}_{receptorID}","{name}_best.pdbqt"),
            percentage = config["RESULT_NUMBER"],
            receptorID = config["TARGETS"][0].split(',')[0],
            name = targetList)
    output:
        report(path.join(OUTPUT_DIR,"results","rescreening_{percentage}","{receptorID}","union.csv"), category="Rescreening"),
        report(path.join(OUTPUT_DIR,"results","rescreening_{percentage}","{receptorID}","venn.png"), category="Rescreening")
    script:
        "../scripts/union_venn.py"
