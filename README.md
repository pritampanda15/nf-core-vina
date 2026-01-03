<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-vina_logo_dark.png">
    <img alt="nf-core/vina" src="docs/images/nf-core-vina_logo_light.png">
  </picture>
</h1>

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/nf-core/vina)
[![GitHub Actions CI Status](https://github.com/nf-core/vina/actions/workflows/nf-test.yml/badge.svg)](https://github.com/nf-core/vina/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/vina/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/vina/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/vina/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.5.1-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.5.1)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/vina)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23vina-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/vina)[![Follow on Bluesky](https://img.shields.io/badge/bluesky-%40nf__core-1185fe?labelColor=000000&logo=bluesky)](https://bsky.app/profile/nf-co.re)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/vina** is a high-throughput, containerized molecular docking pipeline using **AutoDock Vina** with modern ligand and receptor preparation tools. It is designed for structure-based virtual screening campaigns in drug discovery.

![workflow](docs/images/nf-core-vina.png)

### Key Features

- **Modern toolchain only** - No legacy Python 2.7 dependencies or AutoDockTools
- **Receptor preparation** using Open Babel (PDB to PDBQT with polar hydrogens and Gasteiger charges)
- **Ligand preparation** using Meeko/RDKit (SDF/MOL2/SMILES to PDBQT)
- **Virtual screening mode** - Screen thousands of ligands against multiple receptors efficiently
- **Multi-molecule SDF support** - Automatically split and process compound libraries
- **Auto binding site detection** - Automatically detect docking box from co-crystallized ligands in PDB
- **Docking** with AutoDock Vina >= 1.2
- **Chemically correct** - Proper protonation at specified pH, Gasteiger charges, rotatable bond handling
- **Fully containerized** - Docker, Singularity, Podman, Conda support
- **HPC ready** - SLURM, AWS Batch configurations included

### What This Pipeline Does NOT Use

This pipeline explicitly **avoids deprecated tooling**:

- AutoDockTools (ADT) / MGLTools
- Python 2.7
- `prepare_ligand4.py`, `prepare_receptor4.py`
- AutoGrid (grid defined by box coordinates instead)

### Pipeline Summary

1. **Binding Site Detection** (optional) - Auto-detect docking box from co-crystallized ligand in PDB structure
2. **Receptor Preparation** - Convert PDB to PDBQT using Open Babel with polar hydrogen addition and Gasteiger charges
3. **Ligand Preparation** - Convert SDF/MOL2/SMILES to PDBQT using Meeko/RDKit with proper torsion tree setup
4. **Library Splitting** (virtual screening) - Split multi-molecule SDF files into individual ligands
5. **Molecular Docking** - Run AutoDock Vina with user-defined or auto-detected docking box
6. **Score Parsing** - Extract binding affinities and generate per-sample CSV reports
7. **Aggregation** - Combine all results into summary tables ranked by affinity
8. **Reporting** - Generate MultiQC HTML report with docking results

### Pipeline Diagram

```mermaid
flowchart LR
    subgraph Input[" "]
        A[/"Samplesheet"/]
        B[/"Receptor (PDB)"/]
        C[/"Ligand (SDF)"/]
    end

    subgraph Prep["Preparation"]
        D["BINDING_SITE"]
        E["RECEPTOR_PREP"]
        F["LIGAND_PREP"]
    end

    subgraph Dock["Docking"]
        G["VINA_DOCK"]
    end

    subgraph Analysis["Analysis"]
        H["VINA_PARSE"]
        I["SCORE_AGGREGATE"]
        J["MULTIQC"]
    end

    subgraph Output[" "]
        K[/"Poses (PDBQT)"/]
        L[/"Scores (CSV)"/]
        M[/"Report (HTML)"/]
    end

    A --> B & C
    B -.->|optional| D
    D -.-> E
    B --> E
    C --> F
    E & F --> G
    G --> H & K
    H --> I & L
    I --> J --> M

    style D stroke:#ef6c00,stroke-dasharray: 5 5
```

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your receptor-ligand pairs:

`samplesheet.csv`:

```csv
sample,receptor,ligand,center_x,center_y,center_z
dock_1,/path/to/receptor.pdb,/path/to/ligand1.sdf,10.5,20.3,15.0
dock_2,/path/to/receptor.pdb,/path/to/ligand2.mol2,10.5,20.3,15.0
dock_3,/path/to/receptor.pdb,/path/to/ligand3.smi,10.5,20.3,15.0
```

Each row represents a docking job with a receptor-ligand pair and optional docking box coordinates.

Now, you can run the pipeline using:

```bash
nextflow run nf-core/vina \
   -profile docker \
   --input samplesheet.csv \
   --outdir results \
   --center_x 10.5 \
   --center_y 20.3 \
   --center_z 15.0
```

Or with per-sample box coordinates defined in the samplesheet:

```bash
nextflow run nf-core/vina \
   -profile docker \
   --input samplesheet.csv \
   --outdir results
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/vina/usage) and the [parameter documentation](https://nf-co.re/vina/parameters).

## Pipeline output

The pipeline generates the following outputs:

```
results/
├── binding_site/           # Auto-detected binding site coordinates (if enabled)
├── receptor_prep/          # Prepared receptor PDBQT files
├── ligand_prep/            # Prepared ligand PDBQT files
├── docking/                # Docked poses and Vina logs
├── scores/                 # Per-sample score CSVs
├── results/                # Aggregated results
│   ├── docking_results_all.csv       # All poses from all samples
│   └── docking_results_summary.csv   # Best pose per sample
├── multiqc/                # MultiQC HTML report
└── pipeline_info/          # Execution reports
```

For more details about the output files and reports, please refer to the [output documentation](https://nf-co.re/vina/output).

## Credits

nf-core/vina was originally written by Pritam Kumar Panda.

We thank the following projects and their developers:

- [AutoDock Vina](https://vina.scripps.edu/) - Scripps Research
- [Open Babel](https://openbabel.org/) - Open Babel development team
- [Meeko](https://github.com/forlilab/Meeko) - Forli Lab, Scripps Research
- [RDKit](https://www.rdkit.org/) - RDKit community

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#vina` channel](https://nfcore.slack.com/channels/vina) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
