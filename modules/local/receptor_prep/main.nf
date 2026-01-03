process RECEPTOR_PREP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openbabel:3.1.1--2' :
        'quay.io/biocontainers/openbabel:3.1.1--2' }"

    input:
    tuple val(meta), path(receptor_pdb)

    output:
    tuple val(meta), path("*.pdbqt"), emit: pdbqt
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ph = params.receptor_ph ?: 7.4
    def remove_water = params.remove_water != false ? true : false
    def remove_heteroatoms = params.remove_heteroatoms ?: false
    """
    # Step 1: Clean PDB file - remove water and optionally other heteroatoms
    if [ "${remove_heteroatoms}" = "true" ]; then
        # Remove all HETATM records (waters, ligands, ions, etc.)
        grep -v "^HETATM" ${receptor_pdb} > cleaned_receptor.pdb
        echo "Removed all heteroatoms including waters, ligands, and ions"
    elif [ "${remove_water}" = "true" ]; then
        # Remove only water molecules (HOH, WAT, H2O, TIP, SPC residues)
        grep -v "^\\(HETATM\\|ATOM\\).*\\(HOH\\|WAT\\|H2O\\|TIP\\|SPC\\|SOL\\)" ${receptor_pdb} > cleaned_receptor.pdb
        echo "Removed water molecules"
    else
        cp ${receptor_pdb} cleaned_receptor.pdb
        echo "No cleaning performed"
    fi

    # Count atoms before and after cleaning
    echo "Original atoms: \$(grep -c "^\\(ATOM\\|HETATM\\)" ${receptor_pdb} || echo 0)"
    echo "Cleaned atoms: \$(grep -c "^\\(ATOM\\|HETATM\\)" cleaned_receptor.pdb || echo 0)"

    # Step 2: Convert PDB to PDBQT using Open Babel
    # Add only polar hydrogens at specified pH, compute Gasteiger charges
    # Note: -p adds polar hydrogens (for H-bonding), -xr marks as rigid receptor
    obabel \\
        cleaned_receptor.pdb \\
        -O ${prefix}_receptor.pdbqt \\
        -xr \\
        -p ${ph} \\
        --partialcharge gasteiger \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        openbabel: \$(obabel --version 2>&1 | head -1 | sed 's/Open Babel //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_receptor.pdbqt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        openbabel: \$(obabel --version 2>&1 | head -1 | sed 's/Open Babel //')
    END_VERSIONS
    """
}
