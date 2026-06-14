#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { BUILD_TRANSCRIPTS } from './modules/build_transcripts.nf'
include { GENERATE_PEPTIDES } from './modules/build_peptides.nf'
include { BUILD_PEPTIDE_TXT } from './modules/build_peptides_txt.nf'
include { MAKE_FASTA; NETCHOP } from './modules/netchop.nf'
include { SPLIT_NETMHCPAN_BATCHES; NETMHCPAN } from './modules/netmhcpan.nf'
include { NETMHCPAN as NEOPEP_NETMHCPAN } from './modules/netmhcpan.nf'
include { BUILD_MHCFLURRY_INPUT; MHCFLURRY } from './modules/mhcflurry.nf'
include { BUILD_MHCFLURRY_INPUT as BUILD_EXISTING_MHCFLURRY_INPUT; MHCFLURRY as EXISTING_MHCFLURRY } from './modules/mhcflurry.nf'
include { MHCFLURRY as NEOPEP_MHCFLURRY } from './modules/mhcflurry.nf'
include { GENERATE_MIXMHC2_PEPTIDES } from './modules/build_mixmhc2_peptides.nf'
include { MIXMHC2_PREPARE_INPUT; MIXMHC2_RUN; MIXMHC2_MERGE } from './modules/mixmhc2.nf'
include { MHCNUGGETS; PSSMHCPAN } from './modules/additional_predictors.nf'
include { MHCNUGGETS as EXISTING_MHCNUGGETS; PSSMHCPAN as EXISTING_PSSMHCPAN } from './modules/additional_predictors.nf'
include { MHCNUGGETS as NEOPEP_MHCNUGGETS; PSSMHCPAN as NEOPEP_PSSMHCPAN } from './modules/additional_predictors.nf'
include { ENSURE_DOCKER_IMAGES; ENSURE_METRIC_DOCKER_IMAGES; ENSURE_NEOPEP_DOCKER_IMAGES } from './modules/docker_images.nf'
include { ADD_FINAL_SCORE } from './modules/add_scores.nf'
include { ADD_EXISTING_METRICS } from './modules/add_existing_metrics.nf'
include { PREPARE_NEOPEP; ADD_NEOPEP_SCORES } from './modules/neopep.nf'


workflow {
    if (params.mixmhc2_only) {
        if (params.mixmhc2_input_glob) {
            mixmhc2_source_files = Channel.fromPath(params.mixmhc2_input_glob, checkIfExists: true)
                .map { path -> tuple(path.baseName, path, "custom_peptides") }
        } else if (params.mixmhc2_peptides_input) {
            mixmhc2_source_files = Channel.value(file(params.mixmhc2_peptides_input))
                .map { path -> tuple(path.baseName, path, "custom_peptides") }
        } else if (params.maf_glob) {
            mixmhc2_source_files = Channel.fromPath(params.maf_glob)
                .map { file -> tuple(file.baseName, file) }
                | BUILD_TRANSCRIPTS
                .map { sample_id, transcript_csv ->
                    def hla = file("${params.hla_dir}/${sample_id}_hla.txt")
                    tuple(sample_id, transcript_csv, hla)
                }
                | GENERATE_MIXMHC2_PEPTIDES
                .map { sample_id, peptides_csv -> tuple(sample_id, peptides_csv, "maf_generated") }
        } else {
            throw new IllegalArgumentException("mixmhc2_only mode requires --mixmhc2_input_glob, --mixmhc2_peptides_input, or a valid --maf_glob")
        }

        mixmhc2_prepared = mixmhc2_source_files
            .map { sample_id, source_file, source_type ->
                tuple(sample_id, source_file, params.hla_dir, source_type)
            }
            | MIXMHC2_PREPARE_INPUT

        mixmhc2_predictions = mixmhc2_prepared
            .map { sample_id, prepared_csv, pairs_file, source_type -> tuple(sample_id, pairs_file) }
            | MIXMHC2_RUN

        mixmhc2_scored = mixmhc2_prepared
            .join(mixmhc2_predictions)
            .map { sample_id, prepared_csv, pairs_file, source_type, prediction_file ->
                tuple(sample_id, prepared_csv, prediction_file)
            }
            | MIXMHC2_MERGE

    } else if (params.neopep_input) {
        neopep_images_ready = ENSURE_NEOPEP_DOCKER_IMAGES()
        neopep_prepared = PREPARE_NEOPEP(Channel.value(file(params.neopep_input)))

        neopep_netmhcpan_batch_files = neopep_prepared.netmhcpan_batches
            .flatMap { sample_id, batch_dir ->
                batch_dir.toFile().listFiles()
                    .findAll { batch_file -> batch_file.name.endsWith(".txt") }
                    .collect { batch_file ->
                        def match = batch_file.name =~ /^(HLA-[^_]+)__chunk_(\d+)\.txt$/
                        if (!match.find()) {
                            throw new IllegalArgumentException("Unexpected Neopep NetMHCpan batch file: ${batch_file.name}")
                        }
                        tuple(sample_id, match.group(1), match.group(2), file(batch_file.toString()))
                    }
            }

        neopep_netmhcpan_out = neopep_netmhcpan_batch_files
            .combine(neopep_images_ready)
            .map { sample_id, hla, chunk_id, pep, ready -> tuple(sample_id, hla, chunk_id, pep) }
            | NEOPEP_NETMHCPAN

        neopep_mhcflurry_out = neopep_prepared.mhcflurry_input
            .combine(neopep_images_ready)
            .map { sample_id, mhcflurry_input_csv, ready -> tuple(sample_id, mhcflurry_input_csv) }
            | NEOPEP_MHCFLURRY

        neopep_mhcnuggets_out = neopep_prepared.expanded
            .combine(neopep_images_ready)
            .map { sample_id, expanded_csv, ready -> tuple(sample_id, expanded_csv) }
            | NEOPEP_MHCNUGGETS

        neopep_pssmhcpan_out = neopep_prepared.expanded
            .combine(neopep_images_ready)
            .map { sample_id, expanded_csv, ready -> tuple(sample_id, expanded_csv) }
            | NEOPEP_PSSMHCPAN

        neopep_scores = neopep_prepared.expanded
            .join(neopep_netmhcpan_out.groupTuple())
            .join(neopep_mhcflurry_out)
            .join(neopep_mhcnuggets_out)
            .join(neopep_pssmhcpan_out)
            .map { sample_id, expanded_csv, netmhcpan_files, mhcflurry_file, mhcnuggets_file, pssmhcpan_file ->
                tuple(sample_id, expanded_csv, netmhcpan_files, mhcflurry_file, mhcnuggets_file, pssmhcpan_file)
            }
            | ADD_NEOPEP_SCORES

    } else if (params.existing_candidates_glob) {
        metric_images_ready = ENSURE_METRIC_DOCKER_IMAGES()

        existing_candidates = Channel.fromPath(params.existing_candidates_glob, checkIfExists: true)
            .map { candidate ->
                def match = candidate.baseName =~ /^(PT\d+)/
                if (!match.find()) {
                    throw new IllegalArgumentException("Could not infer sample id from existing candidate file: ${candidate.name}")
                }
                tuple(match.group(1), candidate)
            }

        existing_candidates.view { "Existing candidate input: $it" }

        existing_mhcflurry_inputs = existing_candidates | BUILD_EXISTING_MHCFLURRY_INPUT

        existing_mhcflurry_out = existing_mhcflurry_inputs
            .combine(metric_images_ready)
            .map { sample_id, mhcflurry_input_csv, ready -> tuple(sample_id, mhcflurry_input_csv) }
            | EXISTING_MHCFLURRY

        existing_mhcnuggets_out = existing_candidates
            .combine(metric_images_ready)
            .map { sample_id, candidate_csv, ready -> tuple(sample_id, candidate_csv) }
            | EXISTING_MHCNUGGETS

        existing_pssmhcpan_out = existing_candidates
            .combine(metric_images_ready)
            .map { sample_id, candidate_csv, ready -> tuple(sample_id, candidate_csv) }
            | EXISTING_PSSMHCPAN

        existing_metrics = existing_candidates
            .join(existing_mhcflurry_out)
            .join(existing_mhcnuggets_out)
            .join(existing_pssmhcpan_out)
            .map { sample_id, candidate_csv, mhcflurry_file, mhcnuggets_file, pssmhcpan_file ->
                tuple(sample_id, candidate_csv, mhcflurry_file, mhcnuggets_file, pssmhcpan_file)
            }
            | ADD_EXISTING_METRICS

    } else {
        docker_images_ready = ENSURE_DOCKER_IMAGES()

        maf_files = Channel.fromPath(params.maf_glob)
            .map { file -> tuple(file.baseName, file) }

        maf_files.view { "MAF input: $it" }

        transcripts = maf_files | BUILD_TRANSCRIPTS

        peptide_inputs = transcripts
            .map { sample_id, transcript_csv ->
                def hla = file("${params.hla_dir}/${sample_id}_hla.txt")
                tuple(sample_id, transcript_csv, hla)
            }

        peptides = peptide_inputs | GENERATE_PEPTIDES
        peptides_txt = peptides | BUILD_PEPTIDE_TXT

        mhcflurry_inputs = peptides | BUILD_MHCFLURRY_INPUT

        mhcflurry_out = mhcflurry_inputs
            .combine(docker_images_ready)
            .map { sample_id, mhcflurry_input_csv, ready -> tuple(sample_id, mhcflurry_input_csv) }
            | MHCFLURRY

        mhcnuggets_out = peptides
            .combine(docker_images_ready)
            .map { sample_id, peptides_csv, ready -> tuple(sample_id, peptides_csv) }
            | MHCNUGGETS

        pssmhcpan_out = peptides
            .combine(docker_images_ready)
            .map { sample_id, peptides_csv, ready -> tuple(sample_id, peptides_csv) }
            | PSSMHCPAN

        fasta_ch = (transcripts | MAKE_FASTA)
            .flatMap { sample_id, fasta_files ->
                def files = fasta_files instanceof Collection ? fasta_files : [fasta_files]
                files.collect { fasta -> tuple(sample_id, fasta) }
            }

        netchop_out = fasta_ch
            .combine(docker_images_ready)
            .map { sample_id, fasta, ready -> tuple(sample_id, fasta) }
            | NETCHOP

        netmhcpan_batches = peptides_txt
            .map { sample_id, pep ->
                def hla = file("${params.hla_dir}/${sample_id}_hla.txt")
                tuple(sample_id, pep, hla)
            }
            | SPLIT_NETMHCPAN_BATCHES

        netmhcpan_batch_files = netmhcpan_batches
            .flatMap { sample_id, batch_files ->
                def files = batch_files instanceof Collection ? batch_files : [batch_files]
                files.collect { batch_file ->
                    def match = batch_file.name =~ /^(HLA-[^_]+)__chunk_(\d+)\.txt$/
                    if (!match.find()) {
                        throw new IllegalArgumentException("Unexpected NetMHCpan batch file: ${batch_file.name}")
                    }
                    tuple(sample_id, match.group(1), match.group(2), batch_file)
                }
            }

        netmhcpan_out = netmhcpan_batch_files
            .combine(docker_images_ready)
            .map { sample_id, hla, chunk_id, pep, ready -> tuple(sample_id, hla, chunk_id, pep) }
            | NETMHCPAN

        final_scores = peptides
            .join(netchop_out.groupTuple())
            .join(netmhcpan_out.groupTuple())
            .join(mhcflurry_out)
            .join(mhcnuggets_out)
            .join(pssmhcpan_out)
            .map { sample_id, peptides_csv, netchop_dirs, netmhcpan_files, mhcflurry_file, mhcnuggets_file, pssmhcpan_file ->
                def expression_file = file("${params.expression_dir}/${sample_id}_kallisto_expressions.csv")
                tuple(sample_id, peptides_csv, netchop_dirs, netmhcpan_files, mhcflurry_file, mhcnuggets_file, pssmhcpan_file, expression_file)
            }
            | ADD_FINAL_SCORE
    }
}
