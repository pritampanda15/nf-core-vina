process SPLIT_SDF {
    tag "split_library"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://docker.io/ndonyapour/meeko:latest' :
        'docker.io/ndonyapour/meeko:latest' }"

    input:
    path(sdf_library)

    output:
    path("ligands/*.sdf"), emit: ligands
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3

    import os
    from rdkit import Chem

    os.makedirs('ligands', exist_ok=True)

    # Read all molecules from SDF
    suppl = Chem.SDMolSupplier('${sdf_library}', removeHs=False)

    count = 0
    for i, mol in enumerate(suppl):
        if mol is None:
            print(f"Warning: Could not read molecule {i}")
            continue

        # Get molecule name from properties or generate one
        name = mol.GetProp('_Name') if mol.HasProp('_Name') and mol.GetProp('_Name').strip() else f"ligand_{i:06d}"
        # Clean name for filename
        name = "".join(c if c.isalnum() or c in '-_' else '_' for c in name)

        # Write individual SDF
        writer = Chem.SDWriter(f'ligands/{name}.sdf')
        writer.write(mol)
        writer.close()
        count += 1

    print(f"Split {count} molecules from library")

    # Write versions
    import rdkit
    with open("versions.yml", "w") as f:
        f.write(f'"${task.process}":\\n')
        f.write(f'    rdkit: {rdkit.__version__}\\n')
    """

    stub:
    """
    mkdir -p ligands
    touch ligands/stub_ligand.sdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rdkit: stub
    END_VERSIONS
    """
}
