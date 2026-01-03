process VINA_DOCK {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://docker.io/carlosred/autodock-vina:latest' :
        'docker.io/carlosred/autodock-vina:latest' }"

    input:
    tuple val(meta), path(receptor_pdbqt), path(ligand_pdbqt)

    output:
    tuple val(meta), path("*_docked.pdbqt"), emit: poses
    tuple val(meta), path("*_vina.log")    , emit: log
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Grid box parameters - use meta values if available, otherwise params
    def center_x = meta.center_x ?: params.center_x
    def center_y = meta.center_y ?: params.center_y
    def center_z = meta.center_z ?: params.center_z
    def size_x = meta.size_x ?: params.size_x ?: 20
    def size_y = meta.size_y ?: params.size_y ?: 20
    def size_z = meta.size_z ?: params.size_z ?: 20

    // Docking parameters
    def exhaustiveness = params.exhaustiveness ?: 8
    def num_modes = params.num_modes ?: 9
    def energy_range = params.energy_range ?: 3
    def cpu = task.cpus ?: 1
    """
    vina \\
        --receptor ${receptor_pdbqt} \\
        --ligand ${ligand_pdbqt} \\
        --center_x ${center_x} \\
        --center_y ${center_y} \\
        --center_z ${center_z} \\
        --size_x ${size_x} \\
        --size_y ${size_y} \\
        --size_z ${size_z} \\
        --exhaustiveness ${exhaustiveness} \\
        --num_modes ${num_modes} \\
        --energy_range ${energy_range} \\
        --cpu ${cpu} \\
        --out ${prefix}_docked.pdbqt \\
        ${args} \\
        2>&1 | tee ${prefix}_vina.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autodock-vina: \$(vina --version 2>&1 | grep -oP 'AutoDock Vina \\K[0-9.]+' || echo "1.2.5")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_docked.pdbqt
    touch ${prefix}_vina.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autodock-vina: 1.2.5
    END_VERSIONS
    """
}
