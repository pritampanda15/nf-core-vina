# nf-core/vina: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/vina/usage](https://nf-co.re/vina/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

**nf-core/vina** is a high-throughput molecular docking pipeline using AutoDock Vina with modern ligand and receptor preparation tools. This pipeline is designed for structure-based virtual screening campaigns.

### Important: Modern Toolchain Only

This pipeline uses exclusively modern tools and does **NOT** use:

- AutoDockTools (ADT) / MGLTools
- Python 2.7
- `prepare_ligand4.py`, `prepare_receptor4.py`
- AutoGrid (grid defined by box coordinates instead)

## Samplesheet input

You will need to create a samplesheet with information about the receptor-ligand pairs you would like to dock before running the pipeline. Use this parameter to specify its location.

```bash
--input '[path to samplesheet file]'
```

### Samplesheet format

The samplesheet is a CSV file with the following columns:

| Column     | Required | Description                                      |
|------------|----------|--------------------------------------------------|
| `sample`   | Yes      | Unique sample identifier (no spaces)              |
| `receptor` | Yes      | Path to receptor PDB file                         |
| `ligand`   | Yes      | Path to ligand file (SDF, MOL2, PDB, or SMILES)   |
| `center_x` | No       | X coordinate of docking box center (Angstroms)    |
| `center_y` | No       | Y coordinate of docking box center (Angstroms)    |
| `center_z` | No       | Z coordinate of docking box center (Angstroms)    |
| `size_x`   | No       | Box size in X dimension (default: 20 Angstroms)   |
| `size_y`   | No       | Box size in Y dimension (default: 20 Angstroms)   |
| `size_z`   | No       | Box size in Z dimension (default: 20 Angstroms)   |

### Minimal samplesheet

If you're using global docking box parameters (via CLI), you only need three columns:

```csv title="samplesheet.csv"
sample,receptor,ligand
dock_1,/path/to/receptor.pdb,/path/to/aspirin.sdf
dock_2,/path/to/receptor.pdb,/path/to/ibuprofen.mol2
dock_3,/path/to/receptor.pdb,/path/to/caffeine.smi
```

### Full samplesheet with per-sample box coordinates

For different docking boxes per sample:

```csv title="samplesheet.csv"
sample,receptor,ligand,center_x,center_y,center_z,size_x,size_y,size_z
dock_1,/path/to/receptor1.pdb,/path/to/ligand1.sdf,10.5,20.3,15.0,20,20,20
dock_2,/path/to/receptor2.pdb,/path/to/ligand2.mol2,5.2,10.1,8.5,25,25,25
dock_3,/path/to/receptor1.pdb,/path/to/ligand3.smi,10.5,20.3,15.0,30,30,30
```

### Supported ligand formats

The pipeline supports the following ligand input formats:

- **SDF** (.sdf) - Recommended for prepared ligands
- **MOL2** (.mol2) - Commonly exported from molecular modeling software
- **PDB** (.pdb) - Protein Data Bank format
- **SMILES** (.smi, .smiles) - Will generate 3D coordinates automatically

## Defining the docking box

The docking box (search space) can be defined in two ways:

### Option 1: Global parameters (CLI)

Use command-line parameters to set the same docking box for all samples:

```bash
nextflow run nf-core/vina \
   --input samplesheet.csv \
   --outdir results \
   --center_x 10.5 \
   --center_y 20.3 \
   --center_z 15.0 \
   --size_x 20 \
   --size_y 20 \
   --size_z 20
```

### Option 2: Per-sample in samplesheet

Define docking box coordinates in the samplesheet for each sample. This is useful when docking to different binding sites or receptors.

Per-sample values override global parameters.

### Option 3: Auto-detect from co-crystallized ligand

If your receptor PDB contains a co-crystallized ligand (from X-ray crystallography), the pipeline can automatically detect the binding site:

```bash
nextflow run nf-core/vina \
   --input samplesheet.csv \
   --outdir results \
   --auto_binding_site true \
   -profile docker
```

The pipeline will:
1. Parse the receptor PDB for HETATM records (excluding waters, ions, buffers)
2. Calculate the centroid of the co-crystallized ligand
3. Determine box dimensions based on ligand extent + padding
4. Use these coordinates for docking

#### Specifying a specific ligand

If the PDB contains multiple ligands, specify which one to use:

```bash
nextflow run nf-core/vina \
   --input samplesheet.csv \
   --outdir results \
   --auto_binding_site true \
   --ligand_id MK1 \
   -profile docker
```

#### Adjusting box padding

The default padding around the ligand is 5 Angstroms. Adjust if needed:

```bash
--box_padding 10  # Larger box for more flexibility
```

#### When to use auto binding site detection

- **Use when**: Your receptor PDB comes from a co-crystal structure with a known ligand bound
- **Don't use when**: You want to dock to a different site, or the PDB has no ligand

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/vina \
   --input ./samplesheet.csv \
   --outdir ./results \
   --center_x 10.5 \
   --center_y 20.3 \
   --center_z 15.0 \
   -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

### Example with all docking parameters

```bash
nextflow run nf-core/vina \
   --input ./samplesheet.csv \
   --outdir ./results \
   --center_x 10.5 \
   --center_y 20.3 \
   --center_z 15.0 \
   --size_x 25 \
   --size_y 25 \
   --size_z 25 \
   --exhaustiveness 32 \
   --num_modes 9 \
   --energy_range 3 \
   --receptor_ph 7.4 \
   --ligand_ph 7.4 \
   -profile docker
```

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/vina -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
center_x: 10.5
center_y: 20.3
center_z: 15.0
size_x: 25
size_y: 25
size_z: 25
exhaustiveness: 32
num_modes: 9
energy_range: 3
receptor_ph: 7.4
ligand_ph: 7.4
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

## Docking parameters explained

### Exhaustiveness

The `exhaustiveness` parameter controls how thoroughly the docking algorithm searches the conformational space:

- **8** (default): Quick screening, suitable for initial virtual screening
- **32**: Production quality, good balance of speed and accuracy
- **64+**: High accuracy studies, significantly slower but more reliable

Higher values increase the probability of finding the global minimum but also increase computation time proportionally.

### Number of modes

The `num_modes` parameter sets the maximum number of binding poses to output:

- Each pose is ranked by predicted binding affinity
- Output is limited to poses within `energy_range` of the best pose
- Default is 9 poses

### Energy range

The `energy_range` parameter (in kcal/mol) filters poses:

- Only poses within this energy difference from the best pose are reported
- Default is 3 kcal/mol
- Increase for more pose diversity

## Virtual Screening Mode

For high-throughput virtual screening with many ligands against one or more receptors, use the **screening mode**.

### Concept

In screening mode, the pipeline:
1. Prepares each unique receptor once
2. Splits your ligand library (multi-molecule SDF file) into individual molecules
3. Automatically creates all receptor × ligand docking combinations
4. Efficiently parallelizes thousands of docking jobs

This is ideal for screening compound libraries (1000s of ligands) against protein targets.

### Samplesheet for Screening Mode

Create a simpler samplesheet with just the receptor(s) and docking box coordinates:

```csv title="screening_samplesheet.csv"
sample,receptor,ligand,center_x,center_y,center_z,size_x,size_y,size_z
1HSG_pocket,receptors/1HSG.pdb,dummy.sdf,16.0,25.0,4.0,20,20,20
2RH1_pocket,receptors/2RH1.pdb,dummy.sdf,10.5,15.2,20.3,25,25,25
```

Note: The `ligand` column can contain any placeholder file (e.g., `dummy.sdf`) in screening mode as ligands come from the library file.

### Running Virtual Screening

```bash
nextflow run nf-core/vina \
   --input screening_samplesheet.csv \
   --outdir results \
   --screening_mode true \
   --ligand_library /path/to/compound_library.sdf \
   -profile docker
```

### Example: 1000 Ligands × 2 Receptors

Given a ligand library with 1000 compounds and 2 receptors in the samplesheet, the pipeline will:
- Prepare 2 receptors (once each)
- Prepare 1000 ligands
- Run 2000 docking jobs (1000 × 2)
- Aggregate all results into ranked score tables

### Preparing Ligand Libraries

Your compound library should be a multi-molecule SDF file. Tools like RDKit or Open Babel can combine individual molecules:

```python
from rdkit import Chem

writer = Chem.SDWriter('library.sdf')
for smiles in smiles_list:
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    writer.write(mol)
writer.close()
```

Or using Open Babel:
```bash
obabel *.sdf -O combined_library.sdf
```

## Receptor Preparation Options

### Water Removal (Default: ON)

By default, the pipeline removes water molecules (HOH, WAT, H2O, etc.) from receptor PDB files:

```bash
--remove_water true   # Default - removes water molecules
--remove_water false  # Keep water molecules
```

### Heteroatom Removal

To remove ALL heteroatoms including waters, ligands, ions, and cofactors:

```bash
--remove_heteroatoms true  # Remove all HETATM records
```

Use this when you want a clean protein-only structure for blind docking.

## HPC and Cloud execution

### SLURM clusters

Use the provided SLURM profile:

```bash
nextflow run nf-core/vina \
   --input samplesheet.csv \
   --outdir results \
   --center_x 10.5 \
   --center_y 20.3 \
   --center_z 15.0 \
   -profile singularity,slurm
```

### AWS Batch

```bash
nextflow run nf-core/vina \
   --input s3://bucket/samplesheet.csv \
   --outdir s3://bucket/results \
   --center_x 10.5 \
   --center_y 20.3 \
   --center_z 15.0 \
   -profile awsbatch \
   -work-dir s3://bucket/work
```

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

Available profiles:

- `docker` - Use Docker containers
- `singularity` - Use Singularity containers
- `podman` - Use Podman containers
- `conda` - Use Conda environments
- `slurm` - Use SLURM scheduler
- `awsbatch` - Use AWS Batch
- `test` - Run with test data

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

```bash
nextflow run nf-core/vina -resume
```

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
