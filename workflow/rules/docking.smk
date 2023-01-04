def get_spacing(gridfile):
    with open(gridfile) as infile:
        for line in gridfile:
            match = re.findall("spacing", line)
            if match:
                 spacing = match.split(" ")[1]
                 return spacing
            

rule docking:
    input:
        receptor = path.join(OUTPUT_DIR,"receptor","{receptorID}.txt"),
        geometry = path.join(OUTPUT_DIR,"grid","{receptorID}_grid.txt"),
        ligands = path.join(OUTPUT_DIR,"library","{database}_{dataset}_{name}_{i}.txt")
    output:
        path.join(OUTPUT_DIR, "output","{receptorID}","{dataset}","{receptorID}.txt_{database}_{dataset}_{name}_{i}.txt.pdbqt.gz")
    envmodules:
        config["VINALC"]
    params:
        dir = path.join(OUTPUT_DIR, "output"),
        gridfile = path.join(config["GRID_DIR"], "{receptorID}.gpf"),
        errorDir = path.join(OUTPUT_DIR, "errorLigands.txt"),
        space = get_space(params.gridfile)
    resources:
        partition = config["DOCKING"]["partition"],
        runtime = config["DOCKING"]["time"],
        constraint = config["DOCKING"]["constraint"],
        mpi = "srun",
        mem_mb_per_cpu = config["DOCKING"]["mem_per_cpu"],
        ntasks = config["DOCKING"]["ntasks"]
    shell:
        """(
        mkdir -p {params.dir}/{wildcards.receptorID}/{wildcards.dataset}
        cd {params.dir}/{wildcards.receptorID}/{wildcards.dataset}
        cp {input.receptor} .
        cp {input.geometry} .
        cp {input.ligands} .
        space=$(grep "spacing" {params.gridfile} | cut -f2 -d" ")
        {resources.mpi} vinalc --recList {wildcards.receptorID}.txt --ligList {wildcards.database}_{wildcards.dataset}_{wildcards.name}_{wildcards.i}.txt --geoList {wildcards.receptorID}_grid.txt --granularity $space
        cd -
        )"""

def aggregate_in(wildcards):
    checkpoint_output = checkpoints.split.get(**wildcards).output[0]
    files_names = expand(path.join(OUTPUT_DIR,"output","{{receptorID}}", "{{dataset}}", "{{receptorID}}.txt_{{database}}_{{dataset}}_{{name}}_{i}.txt.pdbqt.gz"),
        i = glob_wildcards(os.path.join(checkpoint_output, "{i}.pdbqt")).i)
    return files_names

rule mergeDocking:
    input:
        unpack(aggregate_in)
    output:
        path.join(OUTPUT_DIR, "output", "{receptorID}", "{receptorID}_{database}_{dataset}_{name}.pdbqt.gz")
    script:
        "../scripts/mergeOutput.py"
