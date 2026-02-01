//
// DOCKING: Subworkflow for molecular docking with AutoDock Vina
//
// Supports two modes:
// 1. Standard mode: Each samplesheet row is one docking job
// 2. Screening mode: Cross-product of receptors × ligands for virtual screening
//
// Receptor input supports:
// - Local PDB file paths (e.g., /path/to/receptor.pdb)
// - PDB IDs for automatic download from RCSB (e.g., 5KIR)
//

include { RECEPTOR_PREP          } from '../../../modules/local/receptor_prep/main'
include { LIGAND_PREP            } from '../../../modules/local/ligand_prep/main'
include { MOLECULARDOCKING_DOCK  } from '../../../modules/local/moleculardocking_dock/main'
include { MOLECULARDOCKING_PARSE } from '../../../modules/local/moleculardocking_parse/main'
include { SCORE_AGGREGATE        } from '../../../modules/local/score_aggregate/main'
include { SPLIT_SDF              } from '../../../modules/local/split_sdf/main'
include { BINDING_SITE           } from '../../../modules/local/binding_site/main'
include { PDB_DOWNLOAD           } from '../../../modules/local/pdb_download/main'

//
// Function to check if a string is a PDB ID (4-character alphanumeric)
//
def isPdbId(receptor) {
    return receptor ==~ /^[a-zA-Z0-9]{4}$/
}

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

        //
        // Handle PDB ID vs file path for receptors
        // PDB IDs (4-character alphanumeric) are downloaded from RCSB
        //

        // Branch the samplesheet to handle PDB IDs vs local files
        ch_samplesheet
            .branch {
                meta, receptor, ligand ->
                    pdb_id: isPdbId(receptor.toString())
                    local_file: true
            }
            .set { ch_receptor_type }

        // Multicast PDB ID channel for download and join
        ch_receptor_type.pdb_id
            .map { meta, receptor, ligand -> [ meta, receptor.toString(), ligand ] }
            .multiMap { meta, pdb_id, ligand ->
                for_download: [ meta, pdb_id ]
                for_join: [ meta, ligand ]
            }
            .set { ch_pdb_ids_multi }

        // Download PDB structures for PDB IDs
        PDB_DOWNLOAD ( ch_pdb_ids_multi.for_download )
        ch_versions = ch_versions.mix(PDB_DOWNLOAD.out.versions.first().ifEmpty([]))

        // Merge downloaded PDBs with ligand info
        ch_downloaded = PDB_DOWNLOAD.out.pdb
            .join(ch_pdb_ids_multi.for_join)
            .map { meta, pdb_file, ligand -> [ meta, pdb_file, ligand ] }

        // Combine local files with downloaded PDBs
        ch_samplesheet_with_receptors = ch_receptor_type.local_file.mix(ch_downloaded)

        // Check if auto binding site detection is needed
        if (params.auto_binding_site) {
            //
            // AUTO BINDING SITE DETECTION
            // Detect binding site from co-crystallized ligand in PDB
            //

            // Branch samples into those with and without coordinates
            // Note: nf-schema sets missing optional fields to empty lists [], not null
            // So we need to check if the value is a number (not null and not empty list)
            ch_samplesheet_with_receptors
                .branch {
                    meta, receptor, ligand ->
                        def hasCoords = (meta.center_x instanceof Number) &&
                                        (meta.center_y instanceof Number) &&
                                        (meta.center_z instanceof Number)
                        has_coords: hasCoords
                        needs_detection: !hasCoords
                }
                .set { ch_branched }

            // Samples that already have coordinates - pass through directly
            ch_with_coords = ch_branched.has_coords

            // For samples needing detection, use multiMap to create both channels at once
            // This avoids double consumption of the needs_detection channel
            ch_branched.needs_detection
                .multiMap { meta, receptor, ligand ->
                    for_binding_site: [ meta, receptor ]
                    for_join: [ meta.id, meta, receptor, ligand ]
                }
                .set { ch_needs_detection_multi }

            // Run binding site detection
            BINDING_SITE ( ch_needs_detection_multi.for_binding_site )
            ch_versions = ch_versions.mix(BINDING_SITE.out.versions.first().ifEmpty([]))

            // Join binding site results with original sample data using meta.id as key
            ch_detected = BINDING_SITE.out.binding_site
                .map { meta, json_file -> [ meta.id, json_file ] }
                .join(ch_needs_detection_multi.for_join)
                .map { id, json_file, meta, receptor, ligand ->
                    def json = new groovy.json.JsonSlurper().parse(json_file)
                    def updated_meta = meta + [
                        center_x: json.center_x,
                        center_y: json.center_y,
                        center_z: json.center_z,
                        size_x: json.size_x,
                        size_y: json.size_y,
                        size_z: json.size_z,
                        detected_ligand: json.ligand_id
                    ]
                    [ updated_meta, receptor, ligand ]
                }

            // Merge samples with and without detection
            ch_samplesheet_updated = ch_with_coords.mix(ch_detected)

        } else {
            ch_samplesheet_updated = ch_samplesheet_with_receptors
        }

        // Prepare receptors
        ch_receptor = ch_samplesheet_updated.map { meta, receptor, ligand -> [ meta, receptor ] }
        RECEPTOR_PREP ( ch_receptor )
        ch_versions = ch_versions.mix(RECEPTOR_PREP.out.versions.first())

        // Prepare ligands
        ch_ligand = ch_samplesheet_updated.map { meta, receptor, ligand -> [ meta, ligand ] }
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
    MOLECULARDOCKING_DOCK ( ch_docking_input )
    ch_versions = ch_versions.mix(MOLECULARDOCKING_DOCK.out.versions.first())

    //
    // MODULE: Parse docking results
    //
    ch_parse_input = MOLECULARDOCKING_DOCK.out.log
        .join(MOLECULARDOCKING_DOCK.out.poses)
        .map { meta, log, poses ->
            [ meta, log, poses ]
        }

    MOLECULARDOCKING_PARSE ( ch_parse_input )
    ch_versions = ch_versions.mix(MOLECULARDOCKING_PARSE.out.versions.first())

    //
    // MODULE: Aggregate all scores
    //
    ch_all_scores = MOLECULARDOCKING_PARSE.out.scores
        .map { meta, scores -> scores }
        .collect()

    SCORE_AGGREGATE ( ch_all_scores )
    ch_versions = ch_versions.mix(SCORE_AGGREGATE.out.versions)

    emit:
    receptor_pdbqt = RECEPTOR_PREP.out.pdbqt               // channel: [ val(meta), path(pdbqt) ]
    ligand_pdbqt   = LIGAND_PREP.out.pdbqt                 // channel: [ val(meta), path(pdbqt) ]
    poses          = MOLECULARDOCKING_DOCK.out.poses       // channel: [ val(meta), path(pdbqt) ]
    docking_logs   = MOLECULARDOCKING_DOCK.out.log         // channel: [ val(meta), path(log) ]
    scores         = MOLECULARDOCKING_PARSE.out.scores     // channel: [ val(meta), path(csv) ]
    summaries      = MOLECULARDOCKING_PARSE.out.summary    // channel: [ val(meta), path(txt) ]
    all_scores     = SCORE_AGGREGATE.out.all_scores        // channel: path(csv)
    best_scores    = SCORE_AGGREGATE.out.best_scores       // channel: path(csv)
    multiqc_scores = SCORE_AGGREGATE.out.multiqc           // channel: path(csv)
    versions       = ch_versions                           // channel: [ path(versions.yml) ]
}
