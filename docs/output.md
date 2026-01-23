# nf-core/moleculardocking: Output

## Introduction

This document describes the output produced by the pipeline. Most of the results are molecular docking outputs from AutoDock Vina, with an aggregated summary in the MultiQC report.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Binding Site Detection](#binding-site-detection) - Auto-detect docking box (optional)
- [Receptor Preparation](#receptor-preparation) - Convert receptor PDB to PDBQT
- [Ligand Preparation](#ligand-preparation) - Convert ligand to PDBQT
- [Molecular Docking](#molecular-docking) - AutoDock Vina docking
- [Score Parsing](#score-parsing) - Extract binding affinities
- [Aggregated Results](#aggregated-results) - Combined results from all samples
- [MultiQC](#multiqc) - Aggregate report with docking summary
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Binding Site Detection

<details markdown="1">
<summary>Output files</summary>

- `binding_site/`
  - `*_binding_site.json`: JSON file containing detected binding site coordinates and box dimensions.
  - `*_binding_site.txt`: Human-readable summary of the detected binding site.

</details>

When `--auto_binding_site` is enabled, the pipeline automatically detects the docking box from co-crystallized ligands in the receptor PDB file. This step:

- Parses HETATM records from the PDB file
- Excludes common non-ligand molecules (water, ions, buffers, solvents)
- Calculates the centroid of the largest ligand as the box center
- Determines box dimensions based on ligand extent plus padding (default: 5 Å)

**JSON output format:**
```json
{
  "ligand_id": "MK1",
  "atom_count": 45,
  "center_x": 13.1,
  "center_y": 22.5,
  "center_z": 5.6,
  "size_x": 25,
  "size_y": 22,
  "size_z": 20,
  "padding": 5
}
```

### Receptor Preparation

<details markdown="1">
<summary>Output files</summary>

- `receptor_prep/`
  - `*_receptor.pdbqt`: Prepared receptor files in PDBQT format with Gasteiger charges and hydrogens added at specified pH.

</details>

The receptor preparation step uses Open Babel to convert PDB files to PDBQT format. This includes:

- Adding hydrogen atoms at the specified pH (default: 7.4)
- Computing Gasteiger partial charges
- Converting to AutoDock PDBQT format

### Ligand Preparation

<details markdown="1">
<summary>Output files</summary>

- `ligand_prep/`
  - `*_ligand.pdbqt`: Prepared ligand files in PDBQT format with proper torsion tree setup.

</details>

The ligand preparation step uses Meeko (RDKit-based) to convert various ligand formats to PDBQT:

- Supports SDF, MOL2, PDB, and SMILES input formats
- 3D coordinate generation for SMILES (if needed)
- Proper rotatable bond detection
- Torsion tree setup for flexible docking

### Molecular Docking

<details markdown="1">
<summary>Output files</summary>

- `docking/`
  - `*_docked.pdbqt`: Docked ligand poses in PDBQT format, ranked by binding affinity.
  - `*_docking.log`: AutoDock Vina log file with docking scores and RMSD values.

</details>

The docking step runs AutoDock Vina with the prepared receptor and ligand files. Output includes:

- Multiple docked poses (up to `num_modes` poses)
- Binding affinity predictions (kcal/mol) - more negative = stronger binding
- RMSD values for pose clustering

#### Interpreting Docking Scores

| Affinity (kcal/mol) | Interpretation |
|---------------------|----------------|
| < -10               | Very strong binding |
| -10 to -8           | Strong binding |
| -8 to -6            | Moderate binding |
| -6 to -4            | Weak binding |
| > -4                | Very weak binding |

### Score Parsing

<details markdown="1">
<summary>Output files</summary>

- `scores/`
  - `*_scores.csv`: CSV file with all docking poses and their scores for each sample.
  - `*_summary.txt`: Human-readable summary of docking results.

</details>

The score parsing step extracts binding affinities from Vina output:

**CSV format:**
```csv
sample_id,mode,affinity_kcal_mol,rmsd_lower_bound,rmsd_upper_bound
dock_1,1,-8.5,0.0,0.0
dock_1,2,-8.2,1.5,2.3
dock_1,3,-7.9,2.1,3.5
```

### Aggregated Results

<details markdown="1">
<summary>Output files</summary>

- `results/`
  - `docking_results_all.csv`: All docking poses from all samples combined.
  - `docking_results_summary.csv`: Best pose for each sample with summary statistics.
  - `docking_results_mqc.csv`: MultiQC-compatible summary table.

</details>

The aggregation step combines results from all samples into convenient summary files:

**Summary CSV format:**
```csv
sample_id,best_mode,affinity_kcal_mol,rmsd_lower_bound,rmsd_upper_bound,total_poses
dock_1,1,-8.5,0.0,0.0,9
dock_2,1,-7.2,0.0,0.0,9
dock_3,1,-9.1,0.0,0.0,9
```

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: A standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: Directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: Directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. The report includes:

- **Docking Results Table**: Summary of best binding affinities per sample
- **Software Versions**: All tools and versions used in the pipeline
- **Workflow Summary**: Pipeline parameters and execution details

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

## Output directory structure

```
results/
├── binding_site/               # Only when --auto_binding_site is enabled
│   ├── sample1_binding_site.json
│   └── sample1_binding_site.txt
├── receptor_prep/
│   ├── sample1_receptor.pdbqt
│   └── sample2_receptor.pdbqt
├── ligand_prep/
│   ├── sample1_ligand.pdbqt
│   └── sample2_ligand.pdbqt
├── docking/
│   ├── sample1_docked.pdbqt
│   ├── sample1_docking.log
│   ├── sample2_docked.pdbqt
│   └── sample2_docking.log
├── scores/
│   ├── sample1_scores.csv
│   ├── sample1_summary.txt
│   ├── sample2_scores.csv
│   └── sample2_summary.txt
├── results/
│   ├── docking_results_all.csv
│   ├── docking_results_summary.csv
│   └── docking_results_mqc.csv
├── multiqc/
│   ├── multiqc_report.html
│   └── multiqc_data/
└── pipeline_info/
    ├── execution_report.html
    ├── execution_timeline.html
    ├── execution_trace.txt
    └── software_versions.yml
```
