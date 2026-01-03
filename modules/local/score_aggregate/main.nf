process SCORE_AGGREGATE {
    tag "aggregate"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'biocontainers/python:3.11' }"

    input:
    path(score_files)

    output:
    path "docking_results_all.csv"    , emit: all_scores
    path "docking_results_summary.csv", emit: best_scores
    path "docking_results_mqc.csv"    , emit: multiqc
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3

    import csv
    import glob
    import sys
    from collections import defaultdict

    # Collect all score files
    score_files = sorted(glob.glob("*_scores.csv"))

    # Aggregate all scores
    all_scores = []
    best_scores = {}

    for score_file in score_files:
        with open(score_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                all_scores.append(row)
                sample_id = row['sample_id']
                affinity = float(row['affinity_kcal_mol'])
                mode = int(row['mode'])

                # Track best score per sample
                if sample_id not in best_scores or affinity < best_scores[sample_id]['affinity_kcal_mol']:
                    best_scores[sample_id] = {
                        'sample_id': sample_id,
                        'best_mode': mode,
                        'affinity_kcal_mol': affinity,
                        'rmsd_lower_bound': float(row['rmsd_lower_bound']),
                        'rmsd_upper_bound': float(row['rmsd_upper_bound']),
                        'total_poses': 0
                    }
                best_scores[sample_id]['total_poses'] += 1

    # Write all scores
    with open('docking_results_all.csv', 'w', newline='') as f:
        fieldnames = ['sample_id', 'mode', 'affinity_kcal_mol', 'rmsd_lower_bound', 'rmsd_upper_bound']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for score in all_scores:
            writer.writerow(score)

    # Write summary (best scores only)
    with open('docking_results_summary.csv', 'w', newline='') as f:
        fieldnames = ['sample_id', 'best_mode', 'affinity_kcal_mol', 'rmsd_lower_bound', 'rmsd_upper_bound', 'total_poses']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for sample_id in sorted(best_scores.keys()):
            writer.writerow(best_scores[sample_id])

    # Write MultiQC compatible table
    with open('docking_results_mqc.csv', 'w', newline='') as f:
        f.write("# id: 'docking_results'\\n")
        f.write("# section_name: 'Docking Results'\\n")
        f.write("# description: 'AutoDock Vina docking results summary'\\n")
        f.write("# format: 'csv'\\n")
        f.write("# plot_type: 'table'\\n")
        fieldnames = ['Sample', 'Best Affinity (kcal/mol)', 'Total Poses']
        writer = csv.writer(f)
        writer.writerow(fieldnames)
        for sample_id in sorted(best_scores.keys()):
            bs = best_scores[sample_id]
            writer.writerow([sample_id, f"{bs['affinity_kcal_mol']:.2f}", bs['total_poses']])

    print(f"Aggregated {len(all_scores)} poses from {len(best_scores)} samples")

    # Write versions
    with open("versions.yml", "w") as f:
        f.write(f'"${task.process}":\\n')
        f.write(f'    python: {sys.version.split()[0]}\\n')
    """

    stub:
    """
    touch docking_results_all.csv
    touch docking_results_summary.csv
    touch docking_results_mqc.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
