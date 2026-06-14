process MIXMHC2_PREPARE_INPUT {

    tag "$sample_id"

    input:
    tuple val(sample_id), path(source_file), val(hla_dir), val(input_source)

    output:
    tuple val(sample_id), path("${sample_id}_mixmhc2_input.csv"), path("${sample_id}_mixmhc2_pairs.tsv"), val(input_source)

    script:
    """
    PYTHONPATH=${projectDir}/scr python3 -m pipeline.mixmhc2_prepare \
        --input_file ${source_file} \
        --output_file ${sample_id}_mixmhc2_input.csv \
        --pairs_file ${sample_id}_mixmhc2_pairs.tsv \
        --hla_dir ${hla_dir} \
        --default_sample ${sample_id} \
        --input_source ${input_source}
    """
}

process MIXMHC2_RUN {

    tag "$sample_id"

    input:
    tuple val(sample_id), path(pairs_file)

    output:
    tuple val(sample_id), path("${sample_id}_mixmhc2_predictions.csv")

    script:
    """
    MIXMHC2_COMMAND="${params.mixmhc2_command}"

    PYTHONPATH=${projectDir}/scr python3 -m pipeline.mixmhc2_run \
        --input_file ${pairs_file} \
        --output_file ${sample_id}_mixmhc2_predictions.csv \
        --command "\$MIXMHC2_COMMAND" \
        --sample_id ${sample_id}
    """
}

process MIXMHC2_MERGE {

    tag "$sample_id"

    publishDir "output", mode: 'copy'

    input:
    tuple val(sample_id), path(prepared_csv), path(prediction_csv)

    output:
    tuple val(sample_id), path("${sample_id}_mixmhc2_scored.csv")

    script:
    """
    PYTHONPATH=${projectDir}/scr python3 -m pipeline.mixmhc2_merge \
        --input_file ${prepared_csv} \
        --prediction_file ${prediction_csv} \
        --output_file ${sample_id}_mixmhc2_scored.csv
    """
}
