import glob
import re
import os
import requests
from snakemake.logging import logger
import builtins
import importlib
from urllib.parse import urlparse

# all rules in this file are local rules
localrules: dockingResults, dockingResultsTxt, bestLigands, makeHistogram, mergeDocking

def url_reachable(url):
    """
    test for reachable URL
    """
    try:
        r = requests.head(url, allow_redirects=True, timeout=5)
        # Accept common codes: 200 OK, 301/302 redirects, 403 Forbidden (some servers block HEAD)
        return r.status_code in (200, 301, 302, 403)
    except Exception:
        return False


def check_zinc_url(url):
    """Robust check for ZINC availability.

    Behaviour:
    - If config contains `ZINC_IGNORE_CHECK` and it's truthy, return True.
    - If a user override function `zinc_available(url)` exists in
      `workflow.scripts.user_checks` or `user_checks`, call it and use its result.
    - Otherwise try http and https variants, allow redirects and accept 200/301/302/403.
    """
    # allow config to skip checks (useful for clusters/non-interactive runs)
    try:
        if config.get("ZINC_IGNORE_CHECK"):
            return True
    except Exception:
        pass

    # user override hook: look for zinc_available(url) in user_checks
    for modname in ("workflow.scripts.user_checks", "user_checks"):
        try:
            mod = importlib.import_module(modname)
            if hasattr(mod, "zinc_available"):
                try:
                    return bool(mod.zinc_available(url))
                except Exception:
                    # if the user function errors, fall back to default checks
                    logger.warning(f"user zinc_available in {modname} raised an exception; falling back to default checks")
        except Exception:
            continue

    # Try provided url, and fallback to https/http variants
    variants = [url]
    parsed = urlparse(url)
    if parsed.scheme == "http":
        variants.insert(0, url.replace("http://", "https://", 1))
    elif parsed.scheme == "https":
        variants.insert(0, url.replace("https://", "http://", 1))

    for u in variants:
        if url_reachable(u):
            return True
    return False


def getTranches():
    """return traches from parsing last log file"""
    log_files = glob.glob(".snakemake/log/*")
    sorted_files = sorted(log_files, key=os.path.getmtime)
    tranch_list = []
    pattern = "[A-K][A-K][ABCEGI][ABCDEF][RMLH][NMLOP]"
    with open(sorted_files[-2]) as log:
        matches = re.findall(pattern, log.read())
        return list(set(matches))


def library_files(wildcards):
    if DATABASE[0] == "ZINC":  # ZINC database selected

        if SUBSET == "TRANCHES":  # Tranches selected
            out = []
            # Use ZINC_MIRROR from config if available
            zinc_mirror = config.get("ZINC_MIRROR", "files.docking.org")
            if not zinc_mirror.startswith("http://") and not zinc_mirror.startswith("https://"):
                zinc_mirror = "http://" + zinc_mirror
            zinc_mirror = zinc_mirror.rstrip("/")
            zinc_test_url = f"{zinc_mirror}/3D/"
            
            # test for ZINC reachability (robust):
            if not check_zinc_url(zinc_test_url):
                logger.info(
                    f"The ZINC database mirror ({zinc_mirror}) is not accessible right now. Perhaps it is temporarily down?"
                )
                user_input = __import__("builtins").input(
                    "Have you already run this workflow in the current folder with the same input data?(y/n) \n"
                )
                if user_input == "y":
                    logger.info(
                        "Trying to proceed without accessing the ZINC database (http://files.docking.org)"
                    )
                    tranch_list = (
                        getTranches()
                    )  # get list with 6 letter code specifing ZINC tranches

                    receptorID = config["TARGETS"][0].split(",")[0]
                    database = config["DATABASE"]

                    for i in tranch_list:
                        wl = i[0:2]  # first 2 chars are weight+logP = dataset
                        full_name = i  # all 6 chars = name
                        entry = f"docking/{receptorID}/{receptorID}_{database}_{wl}_{full_name}.pdbqt.gz"
                        out.append(entry)
                    if not out:
                        logger.error("No tranche parameter found in last log file.")
                        sys.exit(1)
                    else:
                        return out
                else:
                    logger.error(
                        f"Aborting the snakemake run for now, as the data from the ZINC database ({zinc_mirror}) are temporarily not available. Please try at a later time."
                    )
                    sys.exit(1)
            rawOut = expand(
                path.join(                    
                    "docking",
                    "{receptorID}",
                    "{receptorID}_{database}_{dataset}_{name}.pdbqt.gz",
                ),
                receptorID=config["TARGETS"][0].split(",")[0],
                database=config["DATABASE"],
                dataset=[w+l for w in config["ZINC_INPUT"]["WEIGHT"] for l in config["ZINC_INPUT"]["LOGP"]],
                name=[w+l+r+p+ph+c for w in config["ZINC_INPUT"]["WEIGHT"] 
                      for l in config["ZINC_INPUT"]["LOGP"]
                      for r in config["ZINC_INPUT"]["REACT"]
                      for p in config["ZINC_INPUT"]["PURCHASE"]
                      for ph in config["ZINC_INPUT"]["PH"]
                      for c in config["ZINC_INPUT"]["CHARGE"]],
            )
            for i in rawOut:
                weighLog = i.split("_")[-2]
                restAttr = (i.split("_")[-1])[2:6]
                url = f"{zinc_mirror}/3D/{weighLog}/{restAttr}/{weighLog}{restAttr}.xaa.pdbqt.gz"
                try:
                    r = requests.get(url, allow_redirects=True, timeout=10)
                    # Accept common positive/redirect/forbidden codes and treat them as present
                    if r.status_code in (200, 301, 302, 403):
                        out.append(i)
                except Exception:
                    # treat as not present and continue
                    continue
            if not out:
                logger.error(
                    "All selected tranches are empty; select other parameters!"
                )
                sys.exit(1)
            else:
                return out
        else:  # if not tranches select subset

            out = expand(
                path.join(
                    "docking",
                    "{receptorID}",
                    "{receptorID}_{database}_subsets_{subset}.pdbqt.gz",
                ),
                database=config["DATABASE"],
                subset=config["SUBSET"],
                receptorID=config["TARGETS"][0].split(",")[0],
            )

            url = (
                "https://zinc15.docking.org/substances/subsets/"
                + SUBSET
                + ".mol2?count=all"
            )
            try:
                r = requests.get(url, allow_redirects=True, timeout=10)
                if r.status_code == 200:  # test if subset is valid
                    return out
            except Exception as e:
                logger.warning(f"Could not connect to ZINC to validate subset: {e}")

            try:
                r_zinc = requests.get("https://zinc15.docking.org/", allow_redirects=True, timeout=10)
                if r_zinc.status_code != 200:  # test if ZINC database is available
                    logger.info(
                        "The ZINC database is not accessible right now. Perhaps it is temporarily down?"
                    )
                    # if ZINC not available, but dataset is already downloaded --> continue
                    subset_dir = path.join(INPUT_DIR, config["SUBSET"] + ".mol2")
                    if os.path.isfile(subset_dir):
                        return out
                    else:
                        logger.error(
                            "Subset is not availiable in the specified data folder. \n Abort snakemake run, try again later"
                        )
                        sys.exit(1)
                else:
                    logger.error("Invalid subset name!")
                    sys.exit(1)
            except Exception as e:
                logger.info(
                    f"The ZINC database is not accessible right now (error: {e}). Perhaps it is temporarily down?"
                )
                # if ZINC not available, but dataset is already downloaded --> continue
                subset_dir = path.join(DATABASE, config["SUBSET"] + ".mol2")
                if os.path.isfile(subset_dir):
                    return out
                else:
                    logger.error(
                        "Subset is not availiable in the specified data folder. \n Abort snakemake run, try again later"
                    )
                    sys.exit(1)

    else:  # not ZINC database --> local input data
        best = expand(
            path.join(
                "docking",
                "{receptorID}",
                "{receptorID}_{database}_{dataset}_local.pdbqt.gz",
            ),
            database=config["DATABASE"],
            dataset=config["LOC_DATA"],
            receptorID=config["TARGETS"][0].split(",")[0],
        )
        return best


rule makeHistogram:
    input:
        path.join("results", "{receptorID}.pdbqt.gz"),
    output:
        report(
            path.join("results", "{receptorID}_hist.png"),
            category="Histogram",
        ),
    log:
        "logs/makeHistogram_{receptorID}.log",
    conda:
        "../envs/plotting.yml"
    envmodules:
        config["PYPLOT"],
    script:
        "../scripts/makeHistogram.py"


rule bestLigands:
    input:
        library_files,
    output:
        path.join("results", "{receptorID}.pdbqt.gz"),
    log:
        "logs/bestLigands_{receptorID}.log",
    script:
        "../scripts/mergeOutput.py"


rule dockingResults:
    input:
        path.join("results", "{receptorID}.pdbqt.gz"),
    output:
        path.join("results", "{receptorID}_{percentage}.pdbqt"),
    envmodules:
        config["PYTHON"],
    log:
        "logs/dockingResults_{receptorID}_{percentage}.log",
    threads: 1
    script:
        "../scripts/sortResult.py"


rule dockingResultsTxt:
    input:
        path.join("results", "{receptorID}_{percentage}.pdbqt"),
    output:
        path.join("results", "{receptorID}_{percentage}.csv"),
    log:
        "logs/dockingResultsTxt_{receptorID}_{percentage}.log",
    conda:
        "../envs/simple_pandas.yml"
    wildcard_constraints:
        receptorID="[^/]+",
        percentage="[^/]+",
    script:
        "../scripts/ResultTxt.py"


rule removeDuplicateLigands:
    input:
        path.join("results", "{receptorID}_{percentage}.pdbqt"),
    output:
        path.join(
            "rescreening", "unique", "{receptorID}_{percentage}.pdbqt"
        ),
    log:
        "logs/removeDuplicateLigands_{receptorID}_{percentage}.log",
    shell:
        "sed '/MODEL [2-9]/,/ENDMDL/d' {input} > {output}"


checkpoint split2:
    input:
        path.join(
            "rescreening", "unique", "{receptorID}_{percentage}.pdbqt"
        ),
    output:
        temp(directory(
            os.path.join("scratch", "rescreening_ligands_{percentage}", "{receptorID}")
        )),
    log:
        "logs/split2_{receptorID}_{percentage}.log",
    script:
        "../scripts/splitFile.py"


rule prepareLigands2:
    input:
        ligands=path.join(
            "scratch", "rescreening_ligands_{percentage}", "{receptorID}", "{i}.pdbqt"
        ),
    output:
        ligands=path.join(
            "rescreening_{percentage}", "{name}_{receptorID}", "{i}.txt"
        ),
    log:
        "logs/prepareLigands2_{receptorID}_{percentage}_{name}_{i}.log",
    shell:
        "echo {input.ligands} > {output.ligands}"


rule prepareSecondDocking:
    input:
        grid=path.join("grid", "{name}_grid.txt"),
        receptor=path.join("prepared", "receptor", "{name}.pdbqt"),
    output:
        grid=path.join(
            "rescreening_{percentage}", "{name}_{receptorID}", "{name}.grd"
        ),
        receptor=path.join(
           "rescreening_{percentage}", "{name}_{receptorID}", "{name}.rec"
        ),
    log:
        "logs/prepareSecondDocking_{name}_{receptorID}_{percentage}.log",
    run:
        import shutil

        shutil.copy(input.grid, output.grid)
        shutil.copy(input.receptor, output.receptor)


rule docking2:
    input:
        ligands=path.join(
            "rescreening_{percentage}", "{name}_{receptorID}", "{i}.txt"
        ),
        grid=path.join(
            "rescreening_{percentage}", "{name}_{receptorID}", "{name}.grd"
        ),
        receptor=path.join(
            "rescreening_{percentage}", "{name}_{receptorID}", "{name}.rec"
        ),
    output:
        path.join(
            "rescreening_{percentage}",
            "{name}_{receptorID}",
            "{name}.rec_{i}.txt.pdbqt.gz",
        ),
    log:
        "logs/docking2_{name}_{receptorID}_{percentage}_{i}.log",
    params:
        dir=path.join("rescreening_{percentage}"),
        gridfile=path.join(config["GRID_DIR"], "{name}.gpf"),
        cutOff=config["CUTOFF_VALUE"],
    conda:
        "../envs/vinalc.yml"
    envmodules:
        config["VINALC"],
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
    files_names = expand(
        path.join(
            "rescreening",
            "rescreening_{{percentage}}",
            "{{name}}_{{receptorID}}",
            "{{name}}.rec_{i}.txt.pdbqt.gz",
        ),
        i=glob_wildcards(os.path.join(checkpoint_output, "{i}.pdbqt")).i,
    )
    return files_names


rule mergeDocking2:
    input:
        unpack(aggregate_in2),
    output:
        path.join(
            "rescreening",
            "rescreening_{percentage}",
            "{name}_{receptorID}",
            "{name}.pdbqt.gz",
        ),
    log:
        "logs/mergeDocking2_{name}_{receptorID}_{percentage}.log",
    shell:
        """
            cat {input} >> {output}
        """


rule dockingResults2:
    input:
        path.join(
            "rescreening",
            "rescreening_{percentage}",
            "{name}_{receptorID}",
            "{name}.pdbqt.gz",
        ),
    output:
        path.join(
            "resreening",
            "rescreening_{percentage}",
            "{name}_{receptorID}",
            "{name}_best.pdbqt",
        ),
    envmodules:
        config["PYTHON"],
    log:
        "logs/dockingResults2_{name}_{receptorID}_{percentage}.log",
    script:
        "../scripts/sortResult.py"


rule makeVenn:
    input:
        best=path.join("results", "{receptorID}_{percentage}.pdbqt"),
        re_results=expand(
            path.join(
                "rescreening",
                "rescreening_{percentage}",
                "{name}_{receptorID}",
                "{name}_best.pdbqt",
            ),
            percentage=config["RESULT_NUMBER"],
            receptorID=config["TARGETS"][0].split(",")[0],
            name=targetList,
        ),
    output:
        report(
            path.join(
                "results",
                "rescreening_{percentage}",
                "{receptorID}",
                "union.csv",
            ),
            category="Rescreening",
        ),
        report(
            path.join(
                "rescreening",
                "results",
                "rescreening_{percentage}",
                "{receptorID}",
                "venn.png",
            ),
            category="Rescreening",
        ),
    conda:
        "../envs/plotting.yml"
    log:
        "logs/makeVenn_{receptorID}_{percentage}.log",
    script:
        "../scripts/union_venn.py"
