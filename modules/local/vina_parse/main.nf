process VINA_PARSE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'biocontainers/python:3.11' }"

    input:
    tuple val(meta), path(vina_log), path(docked_pdbqt)

    output:
    tuple val(meta), path("*_scores.csv"), emit: scores
    tuple val(meta), path("*_summary.txt"), emit: summary
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3

    import re
    import csv
    import sys

    sample_id = "${meta.id}"
    log_file = "${vina_log}"
    pdbqt_file = "${docked_pdbqt}"
    prefix = "${prefix}"

    # Parse Vina log file for scores
    scores = []
    parsing_results = False

    with open(log_file, 'r') as f:
        for line in f:
            line = line.strip()
            # Look for the results table header
            if 'mode' in line and 'affinity' in line:
                parsing_results = True
                continue
            if parsing_results:
                # Parse result lines: mode affinity(kcal/mol) rmsd_lb rmsd_ub
                match = re.match(r'^\\s*(\\d+)\\s+(-?[\\d.]+)\\s+([\\d.]+)\\s+([\\d.]+)', line)
                if match:
                    mode = int(match.group(1))
                    affinity = float(match.group(2))
                    rmsd_lb = float(match.group(3))
                    rmsd_ub = float(match.group(4))
                    scores.append({
                        'sample_id': sample_id,
                        'mode': mode,
                        'affinity_kcal_mol': affinity,
                        'rmsd_lower_bound': rmsd_lb,
                        'rmsd_upper_bound': rmsd_ub
                    })
                elif line.startswith('Writing'):
                    parsing_results = False

    # Write CSV with all scores
    with open(f"{prefix}_scores.csv", 'w', newline='') as csvfile:
        fieldnames = ['sample_id', 'mode', 'affinity_kcal_mol', 'rmsd_lower_bound', 'rmsd_upper_bound']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for score in scores:
            writer.writerow(score)

    # Write summary
    with open(f"{prefix}_summary.txt", 'w') as f:
        f.write(f"Docking Summary for {sample_id}\\n")
        f.write("=" * 50 + "\\n\\n")
        if scores:
            best = scores[0]
            f.write(f"Best binding affinity: {best['affinity_kcal_mol']:.2f} kcal/mol\\n")
            f.write(f"Number of poses: {len(scores)}\\n\\n")
            f.write("All poses:\\n")
            f.write("-" * 50 + "\\n")
            f.write(f"{'Mode':>5} {'Affinity':>12} {'RMSD_LB':>10} {'RMSD_UB':>10}\\n")
            f.write("-" * 50 + "\\n")
            for s in scores:
                f.write(f"{s['mode']:>5} {s['affinity_kcal_mol']:>12.2f} {s['rmsd_lower_bound']:>10.2f} {s['rmsd_upper_bound']:>10.2f}\\n")
        else:
            f.write("No docking results found.\\n")

    print(f"Parsed {len(scores)} poses for {sample_id}")

    # Write versions
    with open("versions.yml", "w") as f:
        f.write(f'"${task.process}":\\n')
        f.write(f'    python: {sys.version.split()[0]}\\n')
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_scores.csv
    touch ${prefix}_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
