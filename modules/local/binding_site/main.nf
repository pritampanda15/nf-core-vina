process BINDING_SITE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'quay.io/biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(pdb_file)

    output:
    tuple val(meta), path("*_binding_site.json"), emit: binding_site
    tuple val(meta), path("*_binding_site.txt") , emit: summary
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ligand_id = meta.ligand_id ?: params.ligand_id ?: ''
    def padding = params.box_padding ?: 5
    """
    #!/usr/bin/env python3
    import json
    import sys

    # Configuration
    ligand_id = "${ligand_id}".strip()
    padding = ${padding}
    prefix = "${prefix}"

    # Common molecules to exclude (waters, ions, buffers)
    EXCLUDE = {
        'HOH', 'WAT', 'H2O', 'TIP', 'SPC', 'SOL',  # Water
        'NA', 'CL', 'MG', 'ZN', 'CA', 'K', 'FE', 'MN', 'CU', 'NI', 'CO',  # Ions
        'SO4', 'PO4', 'NO3', 'CO3',  # Common anions
        'GOL', 'EDO', 'PEG', 'PGE', 'MPD',  # Glycols/PEG
        'ACE', 'NH2', 'NME',  # Caps
        'DMS', 'BME', 'DTT',  # Reducing agents
        'TRS', 'EPE', 'MES', 'HEP',  # Buffers
        'ACT', 'FMT', 'ACY',  # Small organics
        'IMD', 'IPA', 'EOH', 'MOH'  # Solvents
    }

    def parse_pdb(filename, ligand_filter=None):
        \"\"\"Parse PDB file and extract ligand coordinates.\"\"\"
        ligands = {}

        with open(filename, 'r') as f:
            for line in f:
                if line.startswith('HETATM'):
                    res_name = line[17:20].strip()

                    # Skip excluded molecules
                    if res_name in EXCLUDE:
                        continue

                    # If specific ligand requested, filter by it
                    if ligand_filter and res_name != ligand_filter:
                        continue

                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])

                        if res_name not in ligands:
                            ligands[res_name] = []
                        ligands[res_name].append((x, y, z))
                    except (ValueError, IndexError):
                        continue

        return ligands

    def calculate_binding_site(coords, padding=5):
        \"\"\"Calculate center and box size from coordinates.\"\"\"
        if not coords:
            return None

        xs = [c[0] for c in coords]
        ys = [c[1] for c in coords]
        zs = [c[2] for c in coords]

        center_x = sum(xs) / len(xs)
        center_y = sum(ys) / len(ys)
        center_z = sum(zs) / len(zs)

        # Box size = ligand extent + padding on each side
        size_x = (max(xs) - min(xs)) + 2 * padding
        size_y = (max(ys) - min(ys)) + 2 * padding
        size_z = (max(zs) - min(zs)) + 2 * padding

        # Minimum box size of 20 Angstroms
        size_x = max(20, size_x)
        size_y = max(20, size_y)
        size_z = max(20, size_z)

        return {
            'center_x': round(center_x, 2),
            'center_y': round(center_y, 2),
            'center_z': round(center_z, 2),
            'size_x': round(size_x, 0),
            'size_y': round(size_y, 0),
            'size_z': round(size_z, 0)
        }

    # Parse PDB file
    ligand_filter = ligand_id if ligand_id else None
    ligands = parse_pdb("${pdb_file}", ligand_filter)

    if not ligands:
        print("ERROR: No ligands found in PDB file", file=sys.stderr)
        if ligand_filter:
            print(f"  Searched for ligand: {ligand_filter}", file=sys.stderr)
        sys.exit(1)

    # If multiple ligands found and no filter, use the largest one
    if len(ligands) > 1:
        print(f"Found {len(ligands)} ligands: {', '.join(ligands.keys())}")
        largest = max(ligands.keys(), key=lambda k: len(ligands[k]))
        print(f"Using largest ligand: {largest}")
        selected_ligand = largest
        coords = ligands[largest]
    else:
        selected_ligand = list(ligands.keys())[0]
        coords = ligands[selected_ligand]

    print(f"Ligand: {selected_ligand} ({len(coords)} atoms)")

    # Calculate binding site
    binding_site = calculate_binding_site(coords, padding)

    if not binding_site:
        print("ERROR: Could not calculate binding site", file=sys.stderr)
        sys.exit(1)

    binding_site['ligand_id'] = selected_ligand
    binding_site['atom_count'] = len(coords)
    binding_site['padding'] = padding

    # Write JSON output
    with open(f"{prefix}_binding_site.json", 'w') as f:
        json.dump(binding_site, f, indent=2)

    # Write human-readable summary
    with open(f"{prefix}_binding_site.txt", 'w') as f:
        f.write(f"Binding Site Detection Summary\\n")
        f.write(f"{'='*40}\\n")
        f.write(f"Ligand ID: {selected_ligand}\\n")
        f.write(f"Atom count: {len(coords)}\\n")
        f.write(f"Padding: {padding} Angstroms\\n\\n")
        f.write(f"Grid Box Center:\\n")
        f.write(f"  center_x: {binding_site['center_x']}\\n")
        f.write(f"  center_y: {binding_site['center_y']}\\n")
        f.write(f"  center_z: {binding_site['center_z']}\\n\\n")
        f.write(f"Grid Box Size:\\n")
        f.write(f"  size_x: {binding_site['size_x']}\\n")
        f.write(f"  size_y: {binding_site['size_y']}\\n")
        f.write(f"  size_z: {binding_site['size_z']}\\n")

    print(f"Center: ({binding_site['center_x']}, {binding_site['center_y']}, {binding_site['center_z']})")
    print(f"Box size: ({binding_site['size_x']}, {binding_site['size_y']}, {binding_site['size_z']})")

    # Write versions
    with open('versions.yml', 'w') as f:
        f.write('"${task.process}":\\n')
        f.write('    python: ' + sys.version.split()[0] + '\\n')
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '{"center_x": 0, "center_y": 0, "center_z": 0, "size_x": 20, "size_y": 20, "size_z": 20}' > ${prefix}_binding_site.json
    echo "Stub output" > ${prefix}_binding_site.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
