process PDB_DOWNLOAD {
    tag "$meta.id"
    label 'process_single'

    // Use the same container as RECEPTOR_PREP (openbabel) which has curl/wget available
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openbabel:3.1.1--2' :
        'quay.io/biocontainers/openbabel:3.1.1--2' }"

    input:
    tuple val(meta), val(pdb_id)

    output:
    tuple val(meta), path("*.pdb"), emit: pdb
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pdb_id_upper = pdb_id.toUpperCase()
    """
    # Download PDB file from RCSB (try wget first, then curl)
    PDB_URL="https://files.rcsb.org/download/${pdb_id_upper}.pdb"

    if command -v wget &> /dev/null; then
        wget -q -O ${prefix}.pdb "\$PDB_URL"
        DOWNLOAD_TOOL="wget"
    elif command -v curl &> /dev/null; then
        curl -fsSL "\$PDB_URL" -o ${prefix}.pdb
        DOWNLOAD_TOOL="curl"
    else
        echo "ERROR: Neither wget nor curl is available" >&2
        exit 1
    fi

    # Verify the file was downloaded and is valid
    if [ ! -s ${prefix}.pdb ]; then
        echo "ERROR: Failed to download PDB ${pdb_id_upper}" >&2
        exit 1
    fi

    # Check if it's a valid PDB file (should start with HEADER or have ATOM records)
    if ! grep -q "^ATOM\\|^HETATM\\|^HEADER" ${prefix}.pdb; then
        echo "ERROR: Downloaded file does not appear to be a valid PDB file" >&2
        exit 1
    fi

    echo "Successfully downloaded PDB ${pdb_id_upper} using \$DOWNLOAD_TOOL"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        download_tool: \$DOWNLOAD_TOOL
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pdb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        download_tool: stub
    END_VERSIONS
    """
}
