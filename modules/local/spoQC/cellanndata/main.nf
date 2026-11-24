process CELL_QC_ANNDATA {
    tag "${meta.id}"
    label 'process_low'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cb8fc03fa657c164c5d83f075578bbb5d9c10f1178165f94e94f33c67efca1a1/data' :
        'community.wave.seqera.io/library/spatialdata-io_spatialdata:b264928c30680e87' }"

    input:
    tuple val(meta), path(spatialdata, stageAs: "input/*"), path(masks, stageAs: "mask??/*"), val(stainings)

    output:
    tuple val(meta), path("*_qc.zarr"), emit: spatialdata
    tuple val("${task.process}"), val('anndata'), eval("python3 -c 'import anndata; print(anndata.__version__)'"), topic: versions, emit: versions_anndata

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    def staining_args = stainings.collect { "--staining ${it as Integer}" }.join(' ')
    """
    export NUMBA_CACHE_DIR="\$PWD/numba-cache"
    mkdir -p "\$NUMBA_CACHE_DIR"
    spoqc_cell_anndata.py \\
        --spatialdata "${spatialdata}" \\
        --spoqc-output mask??/* \\
        ${staining_args} \\
        --output "${prefix}_qc.zarr"
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    mkdir -p "${prefix}_qc.zarr"
    """
}
