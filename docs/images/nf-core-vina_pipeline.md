# nf-core/vina Pipeline Diagram

```mermaid
flowchart LR
    subgraph Input
        A[/"Samplesheet\n(CSV)"/]
        B[/"Receptor\n(PDB)"/]
        C[/"Ligand\n(SDF/MOL2/SMI)"/]
        D[/"Ligand Library\n(Multi-SDF)"/]
    end

    subgraph Optional
        E{"BINDING_SITE\n(Auto-detect)"}
    end

    subgraph Preparation
        F["RECEPTOR_PREP\nOpen Babel"]
        G["LIGAND_PREP\nMeeko/RDKit"]
        H["SPLIT_SDF\n(Screening mode)"]
    end

    subgraph Docking
        I["VINA_DOCK\nAutoDock Vina"]
    end

    subgraph Analysis
        J["VINA_PARSE"]
        K["SCORE_AGGREGATE"]
        L["MULTIQC"]
    end

    subgraph Output
        M[/"Docked Poses\n(PDBQT)"/]
        N[/"Score Tables\n(CSV)"/]
        O[/"MultiQC Report\n(HTML)"/]
    end

    A --> B & C
    B -.-> E
    E -.-> F
    B --> F
    C --> G
    D --> H --> G
    F --> I
    G --> I
    I --> J & M
    J --> K & N
    K --> L
    L --> O

    style E stroke-dasharray: 5 5
```

## Pipeline Steps

| Step | Module | Tool | Description |
|------|--------|------|-------------|
| 1 | BINDING_SITE | Python | Auto-detect docking box from co-crystallized ligand (optional) |
| 2 | RECEPTOR_PREP | Open Babel | Convert PDB to PDBQT, add polar hydrogens, compute charges |
| 3 | LIGAND_PREP | Meeko/RDKit | Convert SDF/MOL2/SMILES to PDBQT with torsion tree |
| 4 | SPLIT_SDF | Python | Split multi-molecule SDF for virtual screening (optional) |
| 5 | VINA_DOCK | AutoDock Vina | Perform molecular docking |
| 6 | VINA_PARSE | Python | Extract binding affinities and poses |
| 7 | SCORE_AGGREGATE | Python | Combine and rank all docking results |
| 8 | MULTIQC | MultiQC | Generate HTML summary report |
