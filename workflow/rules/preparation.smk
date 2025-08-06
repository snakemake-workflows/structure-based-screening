ruleorder: convertMol2 > gunzip


rule targetProtein:
    output:
        path.join(INPUT_DIR, "PDB", "receptor", "{receptorID}.pdb.gz"),
    shell:
        "curl -o {output} {config[TARGET_URL]}/{wildcards.receptorID}.pdb.gz"


rule getZINCdata:
    output:
        path.join(INPUT_DIR, "ZINC", "{dataset}", "{name}.pdbqt.gz"),
    script:
        "../scripts/ZINCdownload.py"


rule getZINCSubsets:
    params:
        sub=config["SUBSET"],
    output:
        expand(
            path.join(INPUT_DIR, "ZINC", "subsets", "{subset}.mol2"),
            subset=config["SUBSET"],
        ),
    script:
        "../scripts/downloadZINCsubsets.py"


rule convertMol2:
    input:
        path.join(INPUT_DIR, "ZINC", "subsets", "{subset}.mol2"),
    output:
        path.join(TMP_DIR, "unzipped", "ZINC", "subsets", "{subset}.pdbqt"),
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
        path.join(TMP_DIR, "unzipped", "{database}", "{dataset}", "local.pdbqt"),
    envmodules:
        config["OPENBABEL"],
    script:
        "../scripts/concatPDBQT.py"


rule ENAMINEdownload:
    output:
        expand(
            path.join(INPUT_DIR, "ENAMINE", "{ENAMINE_collection}"),
            ENAMINE_collection=config["ENAMINE_INPUT"],
        ),
    script:
        "../scripts/ENAMINEdownload.py"


rule SDFToPDBQT:
    input:
        path.join(TMP_DIR, "unzipped", "{database}", "{dataset}", "{name}.sdf"),
    output:
        path.join(TMP_DIR, "unzipped", "{database}", "{dataset}", "{name}.pdbqt"),
    conda:
        "../envs/openbabel.yml"
    envmodules:
        config["OPENBABEL"],
    shell:
        "obabel {input} -opdbqt -O {output}"


rule prepareReceptor:
    input:
        path.join(TMP_DIR, "unzipped", "PDB", "receptor", "{name}.pdb"),
    output:
        path.join(TMP_DIR, "PDB", "receptor", "{name}.pdb"),
    conda:
        "../envs/biopython.yml"
    envmodules:
        config["BIOPYTHON"],
    script:
        "../scripts/prepareReceptor.py"


rule makeReceptorPDBQT:
    input:
        path.join(TMP_DIR, "PDB", "receptor", "{name}.pdb"),
    output:
        path.join(PREPARED_DIR, "receptor", "{name}.pdbqt"),
    conda:
        "../envs/openbabel.yml"
    envmodules:
        config["OPENBABEL"],
    shell:
        "obabel -ipdb {input} -opdbqt -O {output} -xr"


rule gunzip:
    input:
        path.join(INPUT_DIR, "{database}", "{dataset}", "{name}.{filetype}.gz"),
    output:
        path.join(TMP_DIR, "unzipped", "{database}", "{dataset}", "{name}.{filetype}"),
    shell:
        "gunzip < {input} > {output} || touch {output}"


rule cleanLigands:
    input:
        path.join(TMP_DIR, "unzipped", "{database}", "{dataset}", "{name}.pdbqt"),
    output:
        path.join(TMP_DIR, "cleaned", "{database}", "{dataset}", "{name}.pdbqt"),
    shell:
        "awk '/MODEL/ {{s=1}} s {{a[++c]=$0}} /Si/ {{p=1}} /ENDMDL/ {{if (!p) for (i=1;i<=c;i++) print a[i]; delete a;s=p=c=0;next}} !s' {input} > {output}"


checkpoint split:
    input:
        path.join(TMP_DIR, "cleaned", "{database}", "{dataset}", "{name}.pdbqt"),
    output:
        directory(
            os.path.join(TMP_DIR, "prepared", "{database}", "{dataset}", "{name}_dir")
        ),
    script:
        "../scripts/splitFile.py"


rule energyMin:
    input:
        path.join(
            TMP_DIR, "prepared", "{database}", "{dataset}", "{name}_dir", "{i}.pdbqt"
        ),
    output:
        path.join(MIN_DIR, "{database}", "{dataset}", "{name}", "{i}.pdbqt"),
    params:
        algorithm=config["ENERGY_MIN_ALGORITHM"],
        steps=config["STEPS"],
        forcefield=config["FORCEFIELD"],
        convergence=config["CONVERGENCE_CRITERIA"],
    threads: config["ENERGY_MIN"]["threads"]
    resources:
        partition=config["ENERGY_MIN"]["partition"],
        runtime=config["ENERGY_MIN"]["runtime"],
        mem_mb=config["ENERGY_MIN"]["mem_mb"],
    conda:
        "../envs/openbabel.yml"
    envmodules:
        config["OPENBABEL"],
    shell:
        "obabel -ipdbqt {input} -opdbqt -O {output} --minimize --{params.algorithm} --steps {params.steps} --ff {params.forcefield} --crit {params.convergence}"


rule prepareGeometry:
    input:
        path.join(config["GRID_DIR"], "{receptorID}.gpf"),
    output:
        path.join(OUTPUT_DIR, "grid", "{receptorID}_grid.txt"),
    shell:
        "egrep 'npts|gridcenter' {input} |cut -f2-4 -d' '| tac |tr '\n' ' ' > {output} && sed -i -e '$a\ ' {output}"


rule prepareLibrary:
    input:
        path.join(MIN_DIR, "{database}", "{dataset}", "{name}", "{i}.pdbqt"),
    output:
        library=path.join(OUTPUT_DIR, "library", "{database}_{dataset}_{name}_{i}.txt"),
    shell:
        "echo {input} > {output}"


rule prepareDocking:
    input:
        path.join(PREPARED_DIR, "receptor", "{receptorID}.pdbqt"),
    output:
        path.join(OUTPUT_DIR, "receptor", "{receptorID}.txt"),
    shell:
        "echo {input} > {output}"
