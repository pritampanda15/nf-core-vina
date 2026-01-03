process LIGAND_PREP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://docker.io/ndonyapour/meeko:latest' :
        'docker.io/ndonyapour/meeko:latest' }"

    input:
    tuple val(meta), path(ligand)

    output:
    tuple val(meta), path("*.pdbqt"), emit: pdbqt
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ph = params.ligand_ph ?: 7.4
    def ligand_ext = ligand.getExtension().toLowerCase()
    """
    #!/usr/bin/env python3

    import sys
    from pathlib import Path
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    # Determine input format and read molecule
    ligand_file = "${ligand}"
    prefix = "${prefix}"

    if ligand_file.endswith('.sdf') or ligand_file.endswith('.mol'):
        suppl = Chem.SDMolSupplier(ligand_file, removeHs=False)
        mol = next(suppl)
    elif ligand_file.endswith('.mol2'):
        mol = Chem.MolFromMol2File(ligand_file, removeHs=False)
    elif ligand_file.endswith('.pdb'):
        mol = Chem.MolFromPDBFile(ligand_file, removeHs=False)
    elif ligand_file.endswith('.smi') or ligand_file.endswith('.smiles'):
        with open(ligand_file, 'r') as f:
            smiles = f.read().strip().split()[0]
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
    else:
        raise ValueError(f"Unsupported ligand format: {ligand_file}")

    if mol is None:
        raise ValueError(f"Failed to parse ligand: {ligand_file}")

    # Add hydrogens if not present
    if mol.GetNumAtoms() == mol.GetNumHeavyAtoms():
        mol = Chem.AddHs(mol, addCoords=True)

    # Generate 3D coordinates if not present
    if mol.GetNumConformers() == 0:
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)

    # Prepare molecule with Meeko
    preparator = MoleculePreparation()
    mol_setup_list = preparator.prepare(mol)

    if not mol_setup_list:
        raise ValueError("Meeko preparation failed")

    # Write PDBQT
    pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(mol_setup_list[0])

    if not is_ok:
        raise ValueError(f"PDBQT writing failed: {error_msg}")

    with open(f"{prefix}_ligand.pdbqt", "w") as f:
        f.write(pdbqt_string)

    print(f"Successfully prepared ligand: {prefix}_ligand.pdbqt")

    # Write versions
    import meeko
    import rdkit
    with open("versions.yml", "w") as f:
        f.write(f'"${task.process}":\\n')
        f.write(f'    meeko: {meeko.__version__}\\n')
        f.write(f'    rdkit: {rdkit.__version__}\\n')
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_ligand.pdbqt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meeko: stub
        rdkit: stub
    END_VERSIONS
    """
}
