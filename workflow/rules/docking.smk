localrules:
    prepare_docking_local,
    prepare_docking_ligand,


from snakemake.exceptions import WorkflowError


def get_spacing(gridfile):
    """Return spacing as float parsed from a .gpf grid file.

    Raise WorkflowError if the file does not contain a usable spacing value.
    """
    try:
        with open(gridfile) as fh:
            for line in fh:
                if "spacing" in line.lower():
                    parts = line.strip().split()
                    # expect: spacing <number>
                    try:
                        return float(parts[-1])
                    except Exception:
                        raise WorkflowError(
                            f"unable to convert spacing to float from grid file: {gridfile}"
                        )
    except FileNotFoundError:
        raise WorkflowError(f"grid file not found: {gridfile}")

    # nothing found
    raise WorkflowError(f"no spacing found in grid file: {gridfile}")


rule prepare_docking_local:
    """Prepare files for docking: copy receptor, geometry and ligand list into the per-job output directory.

    This runs as a regular (non-MPI) job so we can perform setup steps safely.
    """
    input:
        receptor=rules.makeReceptorPDBQT.output,
        geometry=path.join("grid", "{receptorID}_grid.txt"),
    output:
        temp(path.join("docking", "{receptorID}", "{dataset}", "{receptorID}.txt")),
        temp(path.join("docking", "{receptorID}", "{dataset}", "{receptorID}_grid.txt")),
    message:
        (
            f"  Copying receptor from {str(input.receptor)} to {str(output[0])}; "
            f"Copying geometry from {str(input.geometry)} to {str(output[1])}"
        )
    run:
        import shutil

        shutil.copy(str(input.receptor), str(output[0]))
        shutil.copy(str(input.geometry), str(output[1]))


rule prepare_docking_ligand:
    """Copy a single ligand-list entry into the per-job output directory.

    This is kept as a separate rule because ligand files carry additional
    wildcards (database, name, i) and must not be mixed with the per-receptor
    outputs in `prepare_docking_local`.
    """
    input:
        ligands=path.join("library", "{database}_{dataset}_{name}_{i}.txt"),
    output:
        temp(
            path.join(
                "docking",
                "{receptorID}",
                "{dataset}",
                "{database}_{dataset}_{name}_{i}.txt",
            )
        ),
    params:
        directory=path.join("docking"),
    run:
        import shutil
        import os

        outdir = os.path.join(params.directory, wildcards.receptorID, wildcards.dataset)
        os.makedirs(outdir, exist_ok=True)

        shutil.copy(input.ligands, output[0])


rule docking:
    """MPI docking rule: invoke the MPI program (vinalc) as a single command.

    The preparation step is handled by `prepare_docking_local` which creates the
    necessary files in output/{receptorID}/{dataset}/ so this rule can
    run the MPI program only.
    """
    input:
        rec=rules.prepare_docking_local.output[0],
        geo=rules.prepare_docking_local.output[1],
        lig=rules.prepare_docking_ligand.output.lig,
    output:
        path.join(
            "docking",
            "{receptorID}",
            "{dataset}",
            "{receptorID}_{database}_{dataset}_{name}_{i}.pdbqt.gz",
        ),
    conda:
        "../envs/vinalc.yml"
    envmodules:
        config["VINALC"],
    params:
        # get spacing from the receptor's .gpf at runtime using wildcards
        space=lambda wildcards: get_spacing(
            os.path.join(config["GRID_DIR"], f"{wildcards.receptorID}.gpf")
        ),
    log:
        "logs/docking/{receptorID}_{dataset}_{database}_{name}_{i}.log",
    resources:
        mpi="mpiexec",
    shell:
        (
            "cd docking/{wildcards.receptorID}/{wildcards.dataset} ; "
            "{resources.mpi} vinalc --recList {wildcards.receptorID}.txt "
            "--ligList {wildcards.database}_{wildcards.dataset}_{wildcards.name}_{wildcards.i}.txt "
            "--geoList {wildcards.receptorID}_grid.txt --granularity {params.space} "
        )


def aggregate_in(wildcards):
    checkpoint_output = checkpoints.split.get(**wildcards).output[0]
    files_names = expand(
        path.join(
            "docking",
            "{{receptorID}}",
            "{{dataset}}",
            "{{receptorID}}_{{database}}_{{dataset}}_{{name}}_{i}.pdbqt.gz",
        ),
        i=glob_wildcards(os.path.join(checkpoint_output, "{i}.pdbqt")).i,
    )
    return files_names


rule mergeDocking:
    input:
        unpack(aggregate_in),
    output:
        temp(
            path.join(
                "docking",
                "{receptorID}",
                "{receptorID}_{database}_{dataset}_{name}.pdbqt.gz",
            )
        ),
    script:
        "../scripts/mergeOutput.py"
