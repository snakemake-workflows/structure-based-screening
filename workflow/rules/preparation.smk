ruleorder: convertMol2 > gunzip


# most rules in this file are local rules:
localrules:
    targetProtein,
    getZINCdata,
    getZINCSubsets,
    convertMol2,
    mergeLocalInput,
    ENAMINEdownload,
    SDFToPDBQT,
    prepareReceptor,
    makeReceptorPDBQT,
    gunzip,
    cleanLigands,
    split,
    prepareGeometry,
    prepareLibrary,
    prepareDocking,


rule targetProtein:
    output:
        path.join("PDB", "receptor", "{receptorID}.pdb.gz"),
    log:
        "logs/targetProtein/{receptorID}.log",
    shell:
        "curl -o {output} {config[TARGET_URL]}/{wildcards.receptorID}.pdb.gz &> {log} 2>&1"


rule getZINCdata:
    output:
        temp(path.join(DATABASE, "{dataset}", "{name}.pdbqt.gz")),
    log: "logs/downloadZINC/{dataset}_{name}.log",
    message:
        "Downloading ZINC data for {wildcards.name} from ZINC database {wildcards.dataset}...",
    script:
        "../scripts/ZINCdownload.py"


rule getZINCSubsets:
    params:
        sub=config["SUBSET"],
    output:
        expand(
            path.join(DATABASE, "subsets", "{subset}.mol2"),
            subset=config["SUBSET"],
        ),
    script:
        "../scripts/downloadZINCsubsets.py"


rule convertMol2:
    input:
        path.join("ZINC", "subsets", "{subset}.mol2"),
    output:
        temp(path.join("scratch", "unzipped", "ZINC", "subsets", "{subset}.pdbqt")),
    conda:
        "../envs/openbabel.yml"
    envmodules:
        config["OPENBABEL"],
    shell:
        "obabel {input} -opdbqt -O {output}"


rule mergeLocalInput:
    params:
        in_dir=LOCAL_INPUT,
    output:
        temp(path.join("scratch", "unzipped", "{database}", "{dataset}", "total.pdbqt")),
    conda:
        "../envs/openbabel.yml"
    envmodules:
        config["OPENBABEL"],
    script:
        "../scripts/concatPDBQT.py"


rule ENAMINEdownload:
    output:
        expand(
            path.join("ENAMINE", "{ENAMINE_collection}"),
            ENAMINE_collection=config["ENAMINE_INPUT"],
        ),
    script:
        "../scripts/ENAMINEdownload.py"


rule SDFToPDBQT:
    input:
        path.join("scratch", "unzipped", "{database}", "{dataset}", "{name}.sdf"),
    output:
        temp(
            path.join("scratch", "unzipped", "{database}", "{dataset}", "{name}.pdbqt")
        ),
    conda:
        "../envs/openbabel.yml"
    envmodules:
        config["OPENBABEL"],
    shell:
        "obabel {input} -opdbqt -O {output}"


rule prepareReceptor:
    input:
        rules.targetProtein.output,
    output:
        temp(path.join("scratch", "PDB", "receptor", "{receptorID}.pdb")),
    conda:
        "../envs/biopython.yml"
    envmodules:
        config["BIOPYTHON"],
    log:
        "logs/prepareReceptor/{receptorID}.log",
    script:
        "../scripts/prepareReceptor.py"


rule makeReceptorPDBQT:
    input:
        path.join("scratch", "PDB", "receptor", "{receptorID}.pdb"),
    output:
        path.join("prepared", "receptor", "{receptorID}.pdbqt"),
    conda:
        "../envs/openbabel.yml"
    envmodules:
        config["OPENBABEL"],
    log:
        "logs/makeReceptorPDBQT/{receptorID}.log",
    shell:
        "obabel -ipdb {input} -opdbqt -O {output} -xr > {log} 2>&1"


rule gunzip:
    input:
        path.join("{database}", "{dataset}", "{name}.{filetype}.gz"),
    output:
        temp(
            path.join(
                "scratch", "unzipped", "{database}", "{dataset}", "{name}.{filetype}"
            )
        ),
    conda:
        "../envs/basic.yml"
    log:
        "logs/gunzip/{database}_{dataset}_{name}_{filetype}.log",
    run:
        import gzip, shutil, os
        try:
            with gzip.open(input[0], "rb") as src, open(output[0], "wb") as dst:
                shutil.copyfileobj(src, dst)
        except Exception as e:
            # log the error if a log file is defined, then create an empty output (touch)
            try:
                with open(log[0], "a") as lf:
                    lf.write(str(e) + "\n")
            except Exception:
                pass
            open(output[0], "wb").close()


rule cleanLigands:
    input:
        path.join("scratch", "unzipped", "{database}", "{dataset}", "{name}.pdbqt"),
    output:
        temp(path.join("scratch", "cleaned", "{database}", "{dataset}", "{name}.pdbqt")),
    log:
        "logs/cleanLigands/{database}_{dataset}_{name}.log",
    shell:
        "awk '/MODEL/ {{s=1}} s {{a[++c]=$0}} /Si/ {{p=1}} /ENDMDL/ {{if (!p) for (i=1;i<=c;i++) print a[i]; delete a;s=p=c=0;next}} !s' {input} > {output}"


checkpoint split:
    input:
        path.join("scratch", "cleaned", "{database}", "{dataset}", "{name}.pdbqt"),
    output:
        temp(
            directory(
                os.path.join(
                    "scratch", "prepared", "{database}", "{dataset}", "{name}_dir"
                )
            )
        ),
    log:
        "logs/split/{database}_{dataset}_{name}.log",
    script:
        "../scripts/splitFile.py"


rule energyMin:
    input:
        path.join(
            "scratch", "prepared", "{database}", "{dataset}", "{name}_dir", "{i}.pdbqt"
        ),
    output:
        path.join("minimized", "{database}", "{dataset}", "{name}", "{i}.pdbqt"),
    params:
        algorithm=config["ENERGY_MIN_ALGORITHM"],
        steps=config["STEPS"],
        forcefield=config["FORCEFIELD"],
        convergence=config["CONVERGENCE_CRITERIA"],
    threads: 8
    log:
        "logs/energyMin/{database}_{dataset}_{name}_{i}.log",
    conda:
        "../envs/openbabel.yml"
    envmodules:
        config["OPENBABEL"],
    shell:
        "obabel -ipdbqt {input} -opdbqt -O {output} --minimize --{params.algorithm} --steps {params.steps} --ff {params.forcefield} --crit {params.convergence} > {log} 2>&1"


rule prepareGeometry:
    input:
        path.join(config["GRID_DIR"], "{receptorID}.gpf"),
    output:
        path.join("grid", "{receptorID}_grid.txt"),
    log:
        "logs/prepareGeometry/{receptorID}.log",
    run:
        grid_params = []

        with open(input[0], "r") as f:
            for line in f:
                # Match lines starting with 'npts' or 'gridcenter'
                if line.startswith(("npts", "gridcenter")):
                    # Extract fields 2-4 (space-separated values after the first field)
                    parts = line.strip().split()
                    if len(parts) >= 4:
                        grid_params.append(" ".join(parts[1:4]))

                        # Reverse the order (equivalent to 'tac')
        grid_params.reverse()

        # Write to output file with space separation and trailing newline
        with open(output[0], "w") as f:
            f.write(" ".join(grid_params) + " \n")


rule prepareLibrary:
    input:
        path.join("minimized", "{database}", "{dataset}", "{receptorID}", "{i}.pdbqt"),
    output:
        library=path.join("library", "{database}_{dataset}_{receptorID}_{i}.txt"),
    log:
        "logs/prepareLibrary/{database}_{dataset}_{receptorID}_{i}.log",
    shell:
        "echo {input} > {output}"


rule prepareDocking:
    input:
        rules.makeReceptorPDBQT.output,
    output:
        path.join("receptor", "{receptorID}.txt"),
    log:
        "logs/prepareDocking/{receptorID}.log",
    shell:
        "echo {input} > {output}"
