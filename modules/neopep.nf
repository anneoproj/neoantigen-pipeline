process PREPARE_NEOPEP {

    tag "${input_file.baseName}"

    input:
    path input_file

    output:
    tuple val("neopep"), path("Neopep_data_org_expanded.csv"), path("Neopep_data_org_mhcflurry_input.csv"), path("Neopep_data_org_netmhcpan_batches")

    script:
    """
    PYTHONPATH=${projectDir}/scr python3 -m pipeline.neopep_prepare \
        --input_file ${input_file} \
        --output_dir . \
        --chunk_size 1000 \
        --chunksize 5
    """
}

process ADD_NEOPEP_SCORES {

    tag "$sample_id"

    input:
    tuple val(sample_id), path(expanded_csv), path(netmhcpan_files), path(mhcflurry_file), path(mhcnuggets_file), path(pssmhcpan_file)

    output:
    tuple val(sample_id), path("${sample_id}_neopep_scored.csv")

    script:
    """
    PYTHONPATH=${projectDir}/scr python3 -m pipeline.add_neopep_scores \
        --input_file ${expanded_csv} \
        --netmhcpan_files ${netmhcpan_files} \
        --mhcflurry_file ${mhcflurry_file} \
        --mhcnuggets_file ${mhcnuggets_file} \
        --pssmhcpan_file ${pssmhcpan_file} \
        --output_file ${sample_id}_neopep_scored.csv
    """
}
