# Neoantigen Discovery Pipeline

This Nextflow pipeline processes somatic mutations and ranks neoantigen candidates using peptide-MHC binding, NetChop cleavage scores, and expression-derived weighting.

## Workflow

1. Build mutated protein sequences from each sample MAF.
2. Generate mutation-overlapping 8-14mer peptides.
3. Run NetChop per transcript FASTA.
4. Split NetMHCpan work by HLA allele and peptide chunk.
5. Combine NetChop, NetMHCpan, and expression data into one final scored candidate table.

The final published result is written to `output/`.

## Inputs

By default, place one sample MAF per file in:

```text
data/mafs/<sample_id>.csv
```

Place the matching HLA alleles in:

```text
data/netmhcpan_input/<sample_id>_hla.txt
```

Place the matching expression table in:

```text
data/expressions/<sample_id>_kallisto_expressions.csv
```

For example, sample `PT1` uses:

```text
data/mafs/PT1.csv
data/netmhcpan_input/PT1_hla.txt
data/expressions/PT1_kallisto_expressions.csv
```

The sample id is taken from the MAF filename, so the HLA and expression filenames must use the same sample id.

The input locations can also be provided as run parameters:

```bash
nextflow run main.nf \
    --maf_glob "data/mafs/*.csv" \
    --hla_dir data/netmhcpan_input \
    --expression_dir data/expressions
```

## Requirements

Create the Python environment and install dependencies with uv:

```bash
uv sync
```

The pipeline also requires Nextflow and Docker. Download NetChop and NetMHCpan manually from (https://services.healthtech.dtu.dk/services/NetChop-3.1/) and (https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/) axxordingly, then place the downloaded archives in the Docker build directories:

```text
docker/netchop/netChop.tar.gz
docker/netmhcpan/netMHCpan.tar.gz
```

The workflow builds the local Docker images automatically from `docker/netchop` and `docker/netmhcpan`.

## Run

```bash
nextflow run main.nf
```

NetMHCpan chunking can be tuned at runtime:

```bash
nextflow run main.nf \
    --netmhcpan_chunk_size 2000 \
    --netmhcpan_max_forks 4
```

Increase `netmhcpan_max_forks` if the machine has enough CPU and memory for more concurrent NetMHCpan containers. Decrease `netmhcpan_chunk_size` if individual NetMHCpan tasks are still too large.

## Output

Each sample produces:

```text
output/<sample_id>_neoantigen_candidates.csv
```

The final table contains peptide metadata, HLA allele, NetChop score, IC50, expression score, component scores, and final `PeptideScore`.

## MixMHC2-only mode

The workflow has an isolated MHC-II branch triggered with `--mixmhc2_only true`.

In MHC-II mode:

- If `--mixmhc2_peptides_input` or `--mixmhc2_input_glob` is set, MixMHC2 runs directly on those files.
- If neither is set, peptides are generated from MAF files using `mixmhc2_min_len`/`mixmhc2_max_len` (default 12-21).
- If `HLA` is present in the peptide input, it is used directly.
- If `HLA` is missing, MixMHC2 uses `data/netmhcpan_input/<sample>_hla.txt` by sample and builds the cartesian product.

Supported peptide input formats:

- `TXT`: one peptide per line (minimum valid format)
- `CSV/TSV`: required `Peptide`; optional `HLA`, `Sample`, `Source`/`Origin`, `Protein`, `Gene`

Output files keep input columns and add:

- `MixMHC2pred_Score`
- `MixMHC2pred_PercentileRank`
- `InputSource` (`custom_peptides` or `maf_generated`)

Run example:

```bash
nextflow run main.nf \
    --mixmhc2_only true \
    --mixmhc2_input_glob "data/mixmhc2_inputs/*.tsv" \
    --mixmhc2_command "mixmhc2predictor --input {input_file} --output {output_file}"
```

### Smoke-test fixtures

Use these local fixtures to validate MixMHC2-only mode with custom inputs:

```text
tests/fixtures/mixmhc2/
├── hla/PT1_hla.txt
├── hla/PT2_hla.txt
├── input/PT1_custom_hla.csv
├── input/PT2_custom_no_hla.tsv
├── input/PT1_custom_txt.txt
└── mock_mixmhc2predictor.py
```

- `PT1_custom_hla.csv` includes explicit HLA and should work even without fallback HLA lookup.
- `PT2_custom_no_hla.tsv` has no HLA column and is expanded using `PT2_hla.txt`.
- `PT1_custom_txt.txt` is plain TXT format.

Example checks:

```bash
nextflow run main.nf \
  --mixmhc2_only true \
  --mixmhc2_input_glob "tests/fixtures/mixmhc2/input/PT1_custom_hla.csv" \
  --hla_dir ${PWD}/tests/fixtures/mixmhc2/hla \
  --mixmhc2_command "python3 ${PWD}/tests/fixtures/mixmhc2/mock_mixmhc2predictor.py --input {input_file} --output {output_file}"

nextflow run main.nf \
  --mixmhc2_only true \
  --mixmhc2_input_glob "tests/fixtures/mixmhc2/input/*.tsv" \
  --hla_dir ${PWD}/tests/fixtures/mixmhc2/hla \
  --mixmhc2_command "python3 ${PWD}/tests/fixtures/mixmhc2/mock_mixmhc2predictor.py --input {input_file} --output {output_file}"

nextflow run main.nf \
  --mixmhc2_only true \
  --mixmhc2_peptides_input tests/fixtures/mixmhc2/input/PT1_custom_txt.txt \
  --hla_dir ${PWD}/tests/fixtures/mixmhc2/hla \
  --mixmhc2_command "python3 ${PWD}/tests/fixtures/mixmhc2/mock_mixmhc2predictor.py --input {input_file} --output {output_file}"
```

The output should contain all source columns and:

- `MixMHC2pred_Score`
- `MixMHC2pred_PercentileRank`
- `InputSource` (`custom_peptides`)
```
