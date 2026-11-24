//
// spoQC (>= 0.1) quality control.
//
// Every spoQC step runs as one SPOQC_STEP task. Steps exchange data through spoQC's tmp and report
// folders; each task gets the output folders of the steps it depends on and merges them (as links).
// The dependencies below follow the spoQC metric/prior definitions (needs_metrics / needs_priors).

include { CELL_QC_ANNDATA } from '../../../modules/local/spoQC/cellanndata/main'

include { SPOQC_STEP as SPOQC_ANNOTATION        } from '../../../modules/local/spoQC/step/main'
include { SPOQC_STEP as SPOQC_GENERALQC         } from '../../../modules/local/spoQC/step/main'
include { SPOQC_STEP as SPOQC_DOUBLETQC         } from '../../../modules/local/spoQC/step/main'
include { SPOQC_STEP as SPOQC_VOIDQC            } from '../../../modules/local/spoQC/step/main'
include { SPOQC_STEP as SPOQC_CELLQC            } from '../../../modules/local/spoQC/step/main'
include { SPOQC_STEP as SPOQC_HQCR_IDENT        } from '../../../modules/local/spoQC/step/main'
include { SPOQC_STEP as SPOQC_HQPR_METRICES     } from '../../../modules/local/spoQC/step/main'
include { SPOQC_STEP as SPOQC_HQPR_CLUSTERING   } from '../../../modules/local/spoQC/step/main'
include { SPOQC_STEP as SPOQC_HQPR_REFINEMENT   } from '../../../modules/local/spoQC/step/main'
include { SPOQC_STEP as SPOQC_AMBIENTQC         } from '../../../modules/local/spoQC/step/main'
include { SPOQC_STEP as SPOQC_HQTR_METRICES     } from '../../../modules/local/spoQC/step/main'
include { SPOQC_STEP as SPOQC_HQTR_QV           } from '../../../modules/local/spoQC/step/main'
include { SPOQC_STEP as SPOQC_HQTR_AC           } from '../../../modules/local/spoQC/step/main'
include { SPOQC_STEP as SPOQC_HQTR_CLUSTERING   } from '../../../modules/local/spoQC/step/main'
include { SPOQC_STEP as SPOQC_HQTR_REFINEMENT   } from '../../../modules/local/spoQC/step/main'
include { SPOQC_STEP as SPOQC_HQCR_CELLTYPE     } from '../../../modules/local/spoQC/step/main'
include { SPOQC_STEP as SPOQC_HQPR_BOUNDING_BOX } from '../../../modules/local/spoQC/step/main'
include { SPOQC_STEP as SPOQC_HQPR_CELLTYPE     } from '../../../modules/local/spoQC/step/main'
include { SPOQC_STEP as SPOQC_HQTR_BOUNDING_BOX } from '../../../modules/local/spoQC/step/main'
include { SPOQC_STEP as SPOQC_HQTR_CELLTYPE     } from '../../../modules/local/spoQC/step/main'
include { SPOQC_STEP as SPOQC_COMBINE_MASKS     } from '../../../modules/local/spoQC/step/main'
include { SPOQC_STEP as SPOQC_TRANSCRIPTQC      } from '../../../modules/local/spoQC/step/main'
include { SPOQC_STEP as SPOQC_CELLCYCLEQC       } from '../../../modules/local/spoQC/step/main'
include { SPOQC_STEP as SPOQC_MODELQC           } from '../../../modules/local/spoQC/step/main'
include { SPOQC_STEP as SPOQC_WHOLE_SLIDE_QC    } from '../../../modules/local/spoQC/step/main'
include { SPOQC_STEP as SPOQC_ANALYSIS_OVERVIEW } from '../../../modules/local/spoQC/step/main'
include { SPOQC_STEP as SPOQC_ANALYSIS_CATEGORY } from '../../../modules/local/spoQC/step/main'
include { SPOQC_STEP as SPOQC_ANALYSIS_CLUSTER  } from '../../../modules/local/spoQC/step/main'
include { SPOQC_STEP as SPOQC_FINAL_REPORT      } from '../../../modules/local/spoQC/step/main'

// -------------------------------------------------------------------------------------------
// Helpers: build the [meta, spatialdata, annotation, [upstream out folders]] input tuples
// -------------------------------------------------------------------------------------------

// ch_base: [meta, spatialdata, annotation]; ups: list of channels [meta, out]
def withUpstream(ch_base, List ups) {
    def ch = ch_base.map { meta, sd, ann -> [meta, sd, ann, []] }
    ups.each { up ->
        ch = ch.join(up, by: 0).map { meta, sd, ann, dirs, d -> [meta, sd, ann, dirs + [d]] }
    }
    return ch
}

// Same per [meta, staining]: ups_st are [meta, staining, out], ups_sample ([meta, out]) are
// broadcast to every staining of the sample. Returns [meta, spatialdata, annotation, dirs, staining].
def withUpstreamSt(ch_base_st, List ups_st, List ups_sample) {
    def ch = ch_base_st.map { meta, st, sd, ann -> [[meta, st], meta, st, sd, ann, []] }
    ups_st.each { up ->
        ch = ch.join(up.map { meta, st, d -> [[meta, st], d] }, by: 0)
            .map { key, meta, st, sd, ann, dirs, d -> [key, meta, st, sd, ann, dirs + [d]] }
    }
    ups_sample.each { up ->
        ch = ch.map { key, meta, st, sd, ann, dirs -> [meta, key, st, sd, ann, dirs] }
            .combine(up, by: 0)
            .map { meta, key, st, sd, ann, dirs, d -> [key, meta, st, sd, ann, dirs + [d]] }
    }
    return ch.map { _key, meta, st, sd, ann, dirs -> [meta, sd, ann, dirs, st] }
}

// Split [meta, spatialdata, annotation, dirs, staining] into the process inputs.
def stIn(ch) {
    return ch.multiMap { meta, sd, ann, dirs, st ->
        data: [meta, sd, ann, dirs]
        st: st
    }
}

workflow SPOQC {

    take:
    ch_sd                   // channel: [ val(meta), path(spatialdata) ]
    ch_annotation_src       // channel: [ val(meta), path(annotation) or [] ]
    ch_stainings            // channel: [ val(meta), val(staining) ] - one item per sample x staining

    main:

    // ---------------------------------------------------------------------------------------
    // Annotation: use the samplesheet annotation if given, else spoQC's unsupervised annotation
    // ---------------------------------------------------------------------------------------

    ch_ann_in = ch_sd.join(ch_annotation_src, by: 0, remainder: true)
        .map { meta, sd, ann -> [meta, sd, ann ?: []] }
        .branch { _meta, _sd, ann ->
            given: ann
            missing: true
        }

    SPOQC_ANNOTATION(
        ch_ann_in.missing.map { meta, sd, _ann -> [meta, sd, [], []] },
        "annotation",
        "",
    )

    ch_base = ch_ann_in.given
        .mix(ch_ann_in.missing.join(SPOQC_ANNOTATION.out.annotation, by: 0).map { meta, sd, _a, ann -> [meta, sd, ann] })
    // Steps that do not use the annotation start from the SpatialData directly, so they run in
    // parallel with the annotation instead of waiting for it.
    ch_base_noann = ch_sd.map { meta, sd -> [meta, sd, []] }
    ch_base_st = ch_base.combine(ch_stainings, by: 0).map { meta, sd, ann, st -> [meta, st, sd, ann] }
    ch_base_st_noann = ch_sd.combine(ch_stainings, by: 0).map { meta, sd, st -> [meta, st, sd, []] }

    // ---------------------------------------------------------------------------------------
    // HQCR: cell metrics -> hqcr_ident  (final HQCR probability map)
    // ---------------------------------------------------------------------------------------

    SPOQC_GENERALQC(withUpstream(ch_base, []), "generalqc", "")
    SPOQC_DOUBLETQC(withUpstream(ch_base, []), "doubletqc", "")
    SPOQC_CELLQC(withUpstream(ch_base_noann, []), "cellqc", "")
    SPOQC_VOIDQC(withUpstream(ch_base_noann, [SPOQC_DOUBLETQC.out.out]), "voidqc", "")   // void needs doublet_score

    SPOQC_HQCR_IDENT(
        withUpstream(ch_base, [SPOQC_GENERALQC.out.out, SPOQC_DOUBLETQC.out.out, SPOQC_CELLQC.out.out, SPOQC_VOIDQC.out.out]),
        "hqcr_ident",
        "",
    )

    // ---------------------------------------------------------------------------------------
    // HQPR (per staining): metrics -> clustering -> refinement  (final HQPR probability map)
    // ---------------------------------------------------------------------------------------

    ch_hqpr_met_in = stIn(withUpstreamSt(ch_base_st_noann, [], []))
    SPOQC_HQPR_METRICES(ch_hqpr_met_in.data, "hqpr_metrices", ch_hqpr_met_in.st)
    ch_hqpr_met = SPOQC_HQPR_METRICES.out.out_st.map { meta, st, d -> [meta, st, d] }

    ch_hqpr_clu_in = stIn(withUpstreamSt(ch_base_st_noann, [ch_hqpr_met], []))
    SPOQC_HQPR_CLUSTERING(ch_hqpr_clu_in.data, "hqpr_clustering", ch_hqpr_clu_in.st)
    ch_hqpr_clu = SPOQC_HQPR_CLUSTERING.out.out_st.map { meta, st, d -> [meta, st, d] }

    ch_hqpr_ref_in = stIn(withUpstreamSt(ch_base_st_noann, [ch_hqpr_clu], []))
    SPOQC_HQPR_REFINEMENT(ch_hqpr_ref_in.data, "hqpr_refinement", ch_hqpr_ref_in.st)
    ch_hqpr_ref = SPOQC_HQPR_REFINEMENT.out.out_st.map { meta, st, d -> [meta, st, d] }

    // ---------------------------------------------------------------------------------------
    // HQTR: metrics, QV, ambient -> AC -> clustering -> refinement  (final HQTR probability map)
    // ---------------------------------------------------------------------------------------

    SPOQC_AMBIENTQC(withUpstream(ch_base_noann, []), "ambientqc", "")
    SPOQC_HQTR_METRICES(withUpstream(ch_base_noann, []), "hqtr_metrices", "")
    SPOQC_HQTR_QV(withUpstream(ch_base_noann, []), "hqtr_qv", "")
    SPOQC_HQTR_AC(withUpstream(ch_base_noann, [SPOQC_AMBIENTQC.out.out]), "hqtr_ac", "")   // ac needs global_moran_I
    SPOQC_HQTR_CLUSTERING(
        withUpstream(ch_base_noann, [SPOQC_HQTR_METRICES.out.out, SPOQC_HQTR_QV.out.out, SPOQC_HQTR_AC.out.out]),
        "hqtr_clustering",
        "",
    )
    SPOQC_HQTR_REFINEMENT(withUpstream(ch_base_noann, [SPOQC_HQTR_CLUSTERING.out.out]), "hqtr_refinement", "")

    if (params.spoqc_cell_anndata) {
        ch_export_stainings = ch_stainings.groupTuple(by: 0)
        ch_export_masks = SPOQC_HQCR_IDENT.out.out
            .mix(SPOQC_HQTR_REFINEMENT.out.out, ch_hqpr_ref.map { meta, _st, directory -> [meta, directory] })
            .groupTuple(by: 0)
        ch_export_input = ch_sd.join(ch_export_stainings, by: 0, remainder: true)
            .join(ch_export_masks, by: 0, remainder: true)
            .map { meta, spatialdata, stainings, masks -> [meta, spatialdata, masks ?: [], stainings ?: []] }
        CELL_QC_ANNDATA(ch_export_input, params.spoqc_cell_h5ad)
    }

    ch_masks = SPOQC_HQCR_IDENT.out.out.map { meta, d -> [meta, 'hqcr', '', d] }
        .mix(ch_hqpr_ref.map { meta, st, d -> [meta, 'hqpr', st, d] })
        .mix(SPOQC_HQTR_REFINEMENT.out.out.map { meta, d -> [meta, 'hqtr', '', d] })


    // -----------------------------------------------------------------------------------
    // Cell type refinements, bounding boxes, combined masks
    // -----------------------------------------------------------------------------------

    ch_cell_metrics = [SPOQC_GENERALQC.out.out, SPOQC_DOUBLETQC.out.out, SPOQC_CELLQC.out.out, SPOQC_VOIDQC.out.out]

    SPOQC_HQCR_CELLTYPE(withUpstream(ch_base, ch_cell_metrics + [SPOQC_HQCR_IDENT.out.out]), "hqcr_celltype", "")
    SPOQC_HQTR_BOUNDING_BOX(withUpstream(ch_base_noann, [SPOQC_HQTR_REFINEMENT.out.out]), "hqtr_bounding_box", "")
    SPOQC_HQTR_CELLTYPE(withUpstream(ch_base, ch_cell_metrics + [SPOQC_HQTR_REFINEMENT.out.out]), "hqtr_celltype", "")

    ch_bbox_in = stIn(withUpstreamSt(ch_base_st_noann, [ch_hqpr_ref], []))
    SPOQC_HQPR_BOUNDING_BOX(ch_bbox_in.data, "hqpr_bounding_box", ch_bbox_in.st)
    ch_hqpr_bbox = SPOQC_HQPR_BOUNDING_BOX.out.out_st.map { meta, st, d -> [meta, st, d] }

    ch_ct_in = stIn(withUpstreamSt(ch_base_st, [ch_hqpr_ref], ch_cell_metrics))
    SPOQC_HQPR_CELLTYPE(ch_ct_in.data, "hqpr_celltype", ch_ct_in.st)
    ch_hqpr_ct = SPOQC_HQPR_CELLTYPE.out.out_st.map { meta, st, d -> [meta, st, d] }

    ch_cm_in = stIn(withUpstreamSt(ch_base_st_noann, [ch_hqpr_ref],
        [SPOQC_HQCR_IDENT.out.out, SPOQC_HQTR_REFINEMENT.out.out]))
    SPOQC_COMBINE_MASKS(ch_cm_in.data, "combine_masks", ch_cm_in.st)
    ch_cm = SPOQC_COMBINE_MASKS.out.out_st.map { meta, st, d -> [meta, st, d] }

    // -----------------------------------------------------------------------------------
    // Extra QC
    // -----------------------------------------------------------------------------------

    SPOQC_TRANSCRIPTQC(withUpstream(ch_base, []), "transcriptqc", "")
    SPOQC_CELLCYCLEQC(withUpstream(ch_base, []), "cellcycleqc", "")
    SPOQC_MODELQC(withUpstream(ch_base_noann, []), "modelqc", "")
    SPOQC_WHOLE_SLIDE_QC(withUpstream(ch_base_noann, []), "whole_slide_qc", "")

    // -----------------------------------------------------------------------------------
    // Analysis (all stainings of a sample) and final report
    // -----------------------------------------------------------------------------------

    ch_all_st = ch_hqpr_met.mix(ch_hqpr_clu, ch_hqpr_ref, ch_hqpr_bbox, ch_hqpr_ct, ch_cm)
        .map { meta, _st, d -> [meta, d] }
        .groupTuple(by: 0)

    ch_analysis_ups = [
        SPOQC_GENERALQC.out.out, SPOQC_DOUBLETQC.out.out, SPOQC_CELLQC.out.out, SPOQC_VOIDQC.out.out,
        SPOQC_HQCR_IDENT.out.out, SPOQC_HQCR_CELLTYPE.out.out,
        SPOQC_AMBIENTQC.out.out, SPOQC_HQTR_METRICES.out.out, SPOQC_HQTR_QV.out.out, SPOQC_HQTR_AC.out.out,
        SPOQC_HQTR_CLUSTERING.out.out, SPOQC_HQTR_REFINEMENT.out.out, SPOQC_HQTR_BOUNDING_BOX.out.out,
        SPOQC_HQTR_CELLTYPE.out.out,
    ]
    ch_analysis_in = withUpstream(ch_base, ch_analysis_ups)
        .join(ch_all_st, by: 0, remainder: true)     // samples without stainings have no HQPR outputs
        .map { meta, sd, ann, dirs, st_dirs -> [meta, sd, ann, dirs + (st_dirs ?: [])] }

    SPOQC_ANALYSIS_OVERVIEW(ch_analysis_in, "analysis_overview", "")
    SPOQC_ANALYSIS_CATEGORY(ch_analysis_in, "analysis_category", "")
    SPOQC_ANALYSIS_CLUSTER(ch_analysis_in, "analysis_cluster", "")

    ch_report_in = withUpstream(ch_base, ch_analysis_ups + [
            SPOQC_TRANSCRIPTQC.out.out, SPOQC_CELLCYCLEQC.out.out, SPOQC_MODELQC.out.out, SPOQC_WHOLE_SLIDE_QC.out.out,
            SPOQC_ANALYSIS_OVERVIEW.out.out, SPOQC_ANALYSIS_CATEGORY.out.out, SPOQC_ANALYSIS_CLUSTER.out.out,
        ])
        .join(ch_all_st, by: 0, remainder: true)     // samples without stainings have no HQPR outputs
        .map { meta, sd, ann, dirs, st_dirs -> [meta, sd, ann, dirs + (st_dirs ?: [])] }

    SPOQC_FINAL_REPORT(ch_report_in, "final_report", "")
    ch_report = SPOQC_FINAL_REPORT.out.out

    emit:
    masks  = ch_masks       // channel: [ val(meta), val(modality), val(staining), path(out) ] - final probability maps
    report = ch_report      // channel: [ val(meta), path(out) ]
}
