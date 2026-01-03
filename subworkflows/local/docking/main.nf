//
// DOCKING: Subworkflow for molecular docking with AutoDock Vina
//
// Supports two modes:
// 1. Standard mode: Each samplesheet row is one docking job
// 2. Screening mode: Cross-product of receptors × ligands for virtual screening
//

include { RECEPTOR_PREP   } from '../../../modules/local/receptor_prep/main'
include { LIGAND_PREP     } from '../../../modules/local/ligand_prep/main'
include { VINA_DOCK       } from '../../../modules/local/vina_dock/main'
include { VINA_PARSE      } from '../../../modules/local/vina_parse/main'
include { SCORE_AGGREGATE } from '../../../modules/local/score_aggregate/main'
include { SPLIT_SDF       } from '../../../modules/local/split_sdf/main'

workflow DOCKING {

    take:
    ch_samplesheet  // channel: [ val(meta), path(receptor), path(ligand) ]

    main:

    ch_versions = Channel.empty()

    //
    // Check if screening mode is enabled (ligand_library provided)
    //
    if (params.screening_mode && params.ligand_library) {
        //
        // VIRTUAL SCREENING MODE
        // Split ligand library and create receptor × ligand combinations
        //

        // Split SDF library into individual ligands
        ch_ligand_library = Channel.fromPath(params.ligand_library)
        SPLIT_SDF ( ch_ligand_library )
        ch_versions = ch_versions.mix(SPLIT_SDF.out.versions)

        // Get unique receptors from samplesheet (avoid re-preparing same receptor)
        ch_receptors = ch_samplesheet
            .map { meta, receptor, ligand ->
                def receptor_id = meta.id
                [ receptor_id, meta, receptor ]
            }
            .unique { it[0] }  // Unique by receptor_id
            .map { receptor_id, meta, receptor -> [ meta, receptor ] }

        // Prepare receptors
        RECEPTOR_PREP ( ch_receptors )
        ch_versions = ch_versions.mix(RECEPTOR_PREP.out.versions.first())

        // Prepare each ligand from the split library
        ch_ligands = SPLIT_SDF.out.ligands
            .flatten()
            .map { ligand_file ->
                def ligand_name = ligand_file.baseName
                def meta = [ id: ligand_name ]
                [ meta, ligand_file ]
            }

        LIGAND_PREP ( ch_ligands )
        ch_versions = ch_versions.mix(LIGAND_PREP.out.versions.first())

        // Create cross-product of receptors × ligands
        ch_docking_input = RECEPTOR_PREP.out.pdbqt
            .combine(LIGAND_PREP.out.pdbqt)
            .map { receptor_meta, receptor_pdbqt, ligand_meta, ligand_pdbqt ->
                // Create combined meta with receptor + ligand info
                def combined_meta = [
                    id: "${receptor_meta.id}_${ligand_meta.id}",
                    receptor_id: receptor_meta.id,
                    ligand_id: ligand_meta.id,
                    center_x: receptor_meta.center_x,
                    center_y: receptor_meta.center_y,
                    center_z: receptor_meta.center_z,
                    size_x: receptor_meta.size_x,
                    size_y: receptor_meta.size_y,
                    size_z: receptor_meta.size_z
                ]
                [ combined_meta, receptor_pdbqt, ligand_pdbqt ]
            }

    } else {
        //
        // STANDARD MODE
        // Each samplesheet row is one docking job
        //

        // Prepare receptors
        ch_receptor = ch_samplesheet.map { meta, receptor, ligand -> [ meta, receptor ] }
        RECEPTOR_PREP ( ch_receptor )
        ch_versions = ch_versions.mix(RECEPTOR_PREP.out.versions.first())

        // Prepare ligands
        ch_ligand = ch_samplesheet.map { meta, receptor, ligand -> [ meta, ligand ] }
        LIGAND_PREP ( ch_ligand )
        ch_versions = ch_versions.mix(LIGAND_PREP.out.versions.first())

        // Combine receptor and ligand for docking (joined by meta.id)
        ch_docking_input = RECEPTOR_PREP.out.pdbqt
            .join(LIGAND_PREP.out.pdbqt)
            .map { meta, receptor_pdbqt, ligand_pdbqt ->
                [ meta, receptor_pdbqt, ligand_pdbqt ]
            }
    }

    //
    // MODULE: Run AutoDock Vina docking
    //
    VINA_DOCK ( ch_docking_input )
    ch_versions = ch_versions.mix(VINA_DOCK.out.versions.first())

    //
    // MODULE: Parse Vina results
    //
    ch_parse_input = VINA_DOCK.out.log
        .join(VINA_DOCK.out.poses)
        .map { meta, log, poses ->
            [ meta, log, poses ]
        }

    VINA_PARSE ( ch_parse_input )
    ch_versions = ch_versions.mix(VINA_PARSE.out.versions.first())

    //
    // MODULE: Aggregate all scores
    //
    ch_all_scores = VINA_PARSE.out.scores
        .map { meta, scores -> scores }
        .collect()

    SCORE_AGGREGATE ( ch_all_scores )
    ch_versions = ch_versions.mix(SCORE_AGGREGATE.out.versions)

    emit:
    receptor_pdbqt = RECEPTOR_PREP.out.pdbqt      // channel: [ val(meta), path(pdbqt) ]
    ligand_pdbqt   = LIGAND_PREP.out.pdbqt        // channel: [ val(meta), path(pdbqt) ]
    poses          = VINA_DOCK.out.poses          // channel: [ val(meta), path(pdbqt) ]
    vina_logs      = VINA_DOCK.out.log            // channel: [ val(meta), path(log) ]
    scores         = VINA_PARSE.out.scores        // channel: [ val(meta), path(csv) ]
    summaries      = VINA_PARSE.out.summary       // channel: [ val(meta), path(txt) ]
    all_scores     = SCORE_AGGREGATE.out.all_scores    // channel: path(csv)
    best_scores    = SCORE_AGGREGATE.out.best_scores   // channel: path(csv)
    multiqc_scores = SCORE_AGGREGATE.out.multiqc       // channel: path(csv)
    versions       = ch_versions                       // channel: [ path(versions.yml) ]
}
