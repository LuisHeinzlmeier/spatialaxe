process SPOQC_STEP {
    tag "${meta.id}:${step}${staining != null && staining != '' ? ':' + staining : ''}"
    label 'spoqc'

    container "${params.spoqc_container}"

    input:
    tuple val(meta), path(spatialdata, stageAs: "input/*"), path(annotation, stageAs: "annotation/*"), path(upstream, stageAs: "upstream_??/*")
    val(step)
    val(staining)

    output:
    tuple val(meta), path("out")                                                     , emit: out
    tuple val(meta), val(staining), path("out")                                      , emit: out_st
    tuple val(meta), path("new/**")                                                  , emit: new_files, optional: true
    tuple val(meta), path("out/report/annotation/unsupervised_cell_annotation.tsv")  , emit: annotation, optional: true
    tuple val("${task.process}"), val('spoqc'), eval("python3 -m spoqc --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1 || echo unknown"), topic: versions, emit: versions_spoqc

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("SPOQC_STEP does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args           = task.ext.args ?: ''
    def plots_arg      = params.spoqc_plots ? '' : '--no-plots'
    def annotation_arg = annotation ? "-a ${annotation}" : ''
    def staining_arg   = (staining != null && staining != '') ? "--staining ${staining}" : ''  // 0 is a valid channel
    """
    export MPLCONFIGDIR=\$PWD/tmp NUMBA_CACHE_DIR=\$PWD/tmp
    mkdir -p tmp out/spoQC_tmp out/report

    # spoQC steps communicate through one tmp folder and one report folder: every step adds files
    # and reads what earlier steps wrote. Merge the outputs of all upstream steps as symlinks
    # (nothing is copied) so this step sees the same layout as in a single spoQC run.
    for up in upstream_*/out; do
        [ -d "\$up" ] || continue
        src=\$(readlink -f "\$up")
        cp -rsn "\$src/spoQC_tmp/." out/spoQC_tmp/ 2>/dev/null || true
        cp -rsn "\$src/report/." out/report/ 2>/dev/null || true
    done
    # spoQC rewrites this file in every step; never write through a link into an upstream task.
    rm -f out/report/staining_log.txt

    python3 -m spoqc \\
        -i ${spatialdata} \\
        -o out \\
        -t out/spoQC_tmp \\
        -s ${step} \\
        -n ${task.cpus} \\
        --dataset ${meta.id} \\
        ${annotation_arg} \\
        ${staining_arg} \\
        ${plots_arg} \\
        ${args}

    # spoQC writes its figures to <output_dir>/report/. Fail loudly if the annotation is missing,
    # otherwise every downstream step would silently be skipped.
    if [ "${step}" = "annotation" ] && [ ! -s out/report/annotation/unsupervised_cell_annotation.tsv ]; then
        echo "ERROR: spoQC did not write out/report/annotation/unsupervised_cell_annotation.tsv" >&2
        exit 1
    fi

    # Files this step created are the real files in out/ (inherited ones are symlinks).
    # They are linked into new/ for publishing.
    find out -type f | while read -r f; do
        mkdir -p "new/\$(dirname "\${f#out/}")"
        ln -s "\$(readlink -f "\$f")" "new/\${f#out/}"
    done
    """

    stub:
    """
    mkdir -p out/spoQC_tmp out/report/${step} new/report/${step}
    touch out/report/${step}/done.txt
    ln -s \$PWD/out/report/${step}/done.txt new/report/${step}/done.txt
    if [ "${step}" = "annotation" ]; then
        mkdir -p out/report/annotation
        printf 'Barcode\\tCluster\\n' > out/report/annotation/unsupervised_cell_annotation.tsv
    fi
    """
}
