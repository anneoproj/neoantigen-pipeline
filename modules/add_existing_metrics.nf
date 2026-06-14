process ADD_EXISTING_METRICS {

    tag "$sample_id"

    input:
    tuple val(sample_id), path(candidate_csv), path(mhcflurry_file), path(mhcnuggets_file), path(pssmhcpan_file)

    output:
    tuple val(sample_id), path("${sample_id}_existing_metrics.csv")

    script:
    """
    PYTHONPATH=${projectDir}/scr python3 -m pipeline.add_existing_metrics \
        --input_file ${candidate_csv} \
        --mhcflurry_file ${mhcflurry_file} \
        --mhcnuggets_file ${mhcnuggets_file} \
        --pssmhcpan_file ${pssmhcpan_file} \
        --output_file ${sample_id}_existing_metrics.csv
    """
}
