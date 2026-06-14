process GENERATE_MIXMHC2_PEPTIDES {

    tag "$sample_id"

    input:
    tuple val(sample_id), path(transcript_csv), path(hla_file)

    output:
    tuple val(sample_id), path("${sample_id}_peptides.csv")

    script:
    """
    PYTHONPATH=${projectDir}/scr python3 ${projectDir}/scr/pipeline/peptide_cli.py \
        --input_file ${transcript_csv} \
        --output_file ${sample_id}_peptides.csv \
        --hla_file ${hla_file} \
        --min_len ${params.mixmhc2_min_len} \
        --max_len ${params.mixmhc2_max_len}
    """
}
