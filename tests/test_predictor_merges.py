import os
import sys
import tempfile
import unittest

import pandas as pd

sys.path.insert(0, os.path.abspath("scr"))

from pipeline.add_existing_metrics import process_existing_file
from pipeline.add_final_scores import (
    add_mhcnuggets_scores,
    add_netchop_scores,
    add_pssmhcpan_scores,
)
from pipeline.add_neopep_scores import process_neopep_file
from pipeline.neopep_prepare import prepare_neopep
from scoring.hla import (
    hla_to_colon,
    hla_to_compact_hla,
    hla_to_mhcflurry,
    hla_to_star,
    normalize_hla,
    normalize_peptide,
)
from pipeline.mixmhc2_merge import add_mixmhc2_scores
from pipeline.mixmhc2_prepare import prepare_mixmhc2_inputs
from scoring.mixmhc2 import parse_mixmhc2_output
from scoring.peptide_score import add_peptide_score


class PredictorMergeTests(unittest.TestCase):
    def test_hla_and_peptide_normalization(self):
        self.assertEqual(normalize_hla("HLA-A*01:01"), "A0101")
        self.assertEqual(normalize_hla("HLA-A01:01"), "A0101")
        self.assertEqual(normalize_hla("A0101"), "A0101")
        self.assertEqual(hla_to_mhcflurry("HLA-B08:01"), "HLA-B*08:01")
        self.assertEqual(hla_to_colon("HLA-B08:01"), "HLA-B08:01")
        self.assertEqual(hla_to_star("HLA-B08:01"), "HLA-B*08:01")
        self.assertEqual(hla_to_compact_hla("HLA-B08:01"), "HLA-B0801")
        self.assertEqual(normalize_peptide("AA-C "), "AAC")

    def _write_prediction_csv(self, columns, row):
        handle = tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False)
        handle.write(",".join(columns) + "\n")
        handle.write(",".join(map(str, row)) + "\n")
        handle.close()
        return handle.name

    def test_new_predictor_joins(self):
        df = pd.DataFrame({
            "HLA": ["HLA-A01:01", "HLA-B08:01"],
            "Peptide": ["PEPTIDE", "PEPTIDE"],
        })

        mhcnuggets_path = self._write_prediction_csv(
            ["HLA", "Peptide", "MHCnuggets_IC50"],
            ["HLA-A01:01", "PEPTIDE", 123.4],
        )
        pssmhcpan_path = self._write_prediction_csv(
            ["HLA", "Peptide", "PSSMHCpan_IC50"],
            ["HLA-A0101", "PEPTIDE", 345.6],
        )

        try:
            merged = add_mhcnuggets_scores(df, mhcnuggets_path)
            merged = add_pssmhcpan_scores(merged, pssmhcpan_path)

            self.assertEqual(merged.loc[0, "MHCnuggets_IC50"], 123.4)
            self.assertEqual(merged.loc[0, "PSSMHCpan_IC50"], 345.6)
            self.assertTrue(pd.isna(merged.loc[1, "MHCnuggets_IC50"]))
        finally:
            for path in [mhcnuggets_path, pssmhcpan_path]:
                os.unlink(path)

    def test_netchop_uses_mutation_position_coordinate(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            netchop_dir = os.path.join(temp_dir, "NCOR1_ENST00000268712_tumor")
            os.makedirs(netchop_dir)

            with open(os.path.join(netchop_dir, "netchop.txt"), "w") as handle:
                handle.write(
                    "ENST00000268712\t1275   R  S   0.974995 PT25|NCOR1\n"
                    "ENST00000268712\t1284   S  S   0.701637 PT25|NCOR1\n"
                    "ENST00000268712\t1286   L  S   0.975806 PT25|NCOR1\n"
                    "ENST00000268712\t1295   S  .   0.025248 PT25|NCOR1\n"
                )

            df = pd.DataFrame({
                "Transcript_ID": ["ENST00000268712"],
                "Protein_position": [1285],
                "Peptide_Length": [11],
                "Tumor_start_position": [1276],
                "Peptide": ["ALPRGSPHSEL"],
            })

            scored = add_netchop_scores(df, temp_dir)
            self.assertAlmostEqual(
                scored.loc[0, "NetChop_Tumor"],
                0.701637 * 0.025248,
            )

    def test_existing_metrics_merge_preserves_old_scores(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            input_path = os.path.join(temp_dir, "PT25_candidates.csv")
            output_path = os.path.join(temp_dir, "PT25_with_metrics.csv")

            pd.DataFrame({
                "Sample": ["PT25"],
                "Peptide": ["ALPRGSPHSEL"],
                "HLA": ["HLA-C0102"],
                "IC50": [29.03],
                "ExpressionScore": [1.0],
                "PeptideScore": [0.7963356527286876],
                "peptide_is_validated": [1],
            }).to_csv(input_path, index=False)

            mhcflurry_path = self._write_prediction_csv(
                [
                    "allele",
                    "peptide",
                    "mhcflurry_affinity",
                    "mhcflurry_affinity_percentile",
                    "mhcflurry_processing_score",
                    "mhcflurry_presentation_score",
                    "mhcflurry_presentation_percentile",
                ],
                ["HLA-C*01:02", "ALPRGSPHSEL", 82.27, 0.35, 0.18, 0.66, 0.49],
            )
            mhcnuggets_path = self._write_prediction_csv(
                ["HLA", "Peptide", "MHCnuggets_IC50"],
                ["HLA-C01:02", "ALPRGSPHSEL", 101.0],
            )
            pssmhcpan_path = self._write_prediction_csv(
                ["HLA", "Peptide", "PSSMHCpan_IC50"],
                ["HLA-C0102", "ALPRGSPHSEL", 303.0],
            )

            try:
                process_existing_file(
                    input_file=input_path,
                    mhcflurry_file=mhcflurry_path,
                    mhcnuggets_file=mhcnuggets_path,
                    pssmhcpan_file=pssmhcpan_path,
                    output_file=output_path,
                )

                merged = pd.read_csv(output_path)
                self.assertEqual(merged.loc[0, "PeptideScore"], 0.7963356527286876)
                self.assertEqual(merged.loc[0, "peptide_is_validated"], 1)
                self.assertEqual(merged.loc[0, "MHCflurry_Affinity"], 82.27)
                self.assertEqual(merged.loc[0, "MHCnuggets_IC50"], 101.0)
                self.assertEqual(merged.loc[0, "PSSMHCpan_IC50"], 303.0)
            finally:
                for path in [
                    mhcflurry_path,
                    mhcnuggets_path,
                    pssmhcpan_path,
                ]:
                    os.unlink(path)

    def test_neopep_prepare_expands_mutant_alleles_and_patient_expression(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            input_path = os.path.join(temp_dir, "neopep.tsv")
            output_dir = os.path.join(temp_dir, "prepared")
            with open(input_path, "w") as handle:
                handle.write(
                    "patient\tmutant_seq\tmut_netchop_score_ct\trnaseq_TPM\t"
                    "mutant_best_alleles\tmutant_best_alleles_netMHCpan\t"
                    "mutant_other_significant_alleles_netMHCpan\twt_best_alleles\n"
                )
                handle.write("P1\tPEPA\t0.1\t0\tA0201\tA0201\t\tA0301\n")
                handle.write("P1\tPEPB\t0.2\t10\tA0201,B0702\tB0702\tC0102\tA0301\n")
                handle.write("P1\tPEPC\t0.3\t30\tC0102\t\t\tA0301\n")
                handle.write("P1\tPEPD\t0.4\t100\tB0801\tB0801\t\tA0301\n")

            prepare_neopep(input_path, output_dir, chunk_size=2, chunksize=3)

            expanded = pd.read_csv(
                os.path.join(output_dir, "Neopep_data_org_expanded.csv")
            )
            pepb_hlas = set(expanded.loc[expanded["mutant_seq"] == "PEPB", "HLA"])
            self.assertEqual(pepb_hlas, {"A0201", "B0702", "C0102"})
            self.assertFalse((expanded["HLA"] == "A0301").any())
            self.assertTrue(os.path.exists(os.path.join(output_dir, "netmhcpan_batches")))

    def test_neopep_final_merge_scores_without_ic50_filter(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            input_path = os.path.join(temp_dir, "expanded.csv")
            netmhcpan_path = os.path.join(temp_dir, "netmhcpan.out")
            output_path = os.path.join(temp_dir, "scored.csv")

            pd.DataFrame({
                "patient": ["P1", "P1"],
                "mutant_seq": ["PEPTIDEA", "PEPTIDEB"],
                "rnaseq_TPM": ["10", "20"],
                "HLA": ["A0201", "A0201"],
                "Peptide": ["PEPTIDEA", "PEPTIDEB"],
                "NetChop_Tumor": ["0.5", "0.5"],
                "gene_abundance": ["10", "20"],
                "ExpressionScore": ["1.0", "1.0"],
            }).to_csv(input_path, index=False)

            with open(netmhcpan_path, "w") as handle:
                handle.write("------\n")
                handle.write(
                    "1 HLA-A02:01 PEPTIDEA PEPTIDEA 0 0 0 0 0 PEPTIDEA ID "
                    "0 0 0 0 50 SB\n"
                )
                handle.write(
                    "1 HLA-A02:01 PEPTIDEB PEPTIDEB 0 0 0 0 0 PEPTIDEB ID "
                    "0 0 0 0 700 WB\n"
                )

            mhcflurry_path = self._write_prediction_csv(
                [
                    "allele",
                    "peptide",
                    "mhcflurry_affinity",
                    "mhcflurry_affinity_percentile",
                    "mhcflurry_processing_score",
                    "mhcflurry_presentation_score",
                    "mhcflurry_presentation_percentile",
                ],
                ["HLA-A*02:01", "PEPTIDEA", 12.3, 0.4, 0.1, 0.8, 0.2],
            )
            mhcnuggets_path = self._write_prediction_csv(
                ["HLA", "Peptide", "MHCnuggets_IC50"],
                ["HLA-A02:01", "PEPTIDEA", 11.0],
            )
            pssmhcpan_path = self._write_prediction_csv(
                ["HLA", "Peptide", "PSSMHCpan_IC50"],
                ["HLA-A0201", "PEPTIDEA", 33.0],
            )

            try:
                process_neopep_file(
                    input_file=input_path,
                    netmhcpan_files=[netmhcpan_path],
                    mhcflurry_file=mhcflurry_path,
                    mhcnuggets_file=mhcnuggets_path,
                    pssmhcpan_file=pssmhcpan_path,
                    output_file=output_path,
                )

                scored = pd.read_csv(output_path)
                self.assertEqual(len(scored), 2)
                self.assertEqual(scored.loc[0, "mutant_seq"], "PEPTIDEA")
                self.assertEqual(scored.loc[0, "MHCnuggets_IC50"], 11.0)
                self.assertEqual(scored.loc[0, "PSSMHCpan_IC50"], 33.0)
                self.assertIn("PEPTIDEB", set(scored["mutant_seq"]))
            finally:
                for path in [
                    mhcflurry_path,
                    mhcnuggets_path,
                    pssmhcpan_path,
                ]:
                    os.unlink(path)

    def test_peptide_score_keeps_high_ic50_with_lower_binding_score(self):
        df = pd.DataFrame({
            "IC50": [50.0, 5000.0],
            "NetChop_Tumor": [0.5, 0.5],
            "ExpressionScore": [1.0, 1.0],
        })

        scored = add_peptide_score(df)
        self.assertEqual(len(scored), 2)
        self.assertEqual(scored.loc[0, "IC50"], 50.0)
        high_ic50 = scored.loc[scored["IC50"] == 5000.0].iloc[0]
        self.assertLess(high_ic50["MHCBindingScore"], 0.0)

    def test_mixmhc2_txt_input_with_fallback_hla(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            input_path = os.path.join(temp_dir, "custom.txt")
            with open(input_path, "w") as handle:
                handle.write("PEPTIDE-1\nAA-BC\n")

            hla_dir = temp_dir
            with open(os.path.join(hla_dir, "PT1_hla.txt"), "w") as handle:
                handle.write("HLA-DRB1*15:01,HLA-DPA1*01:03\n")

            prepared_df, pairs_df = prepare_mixmhc2_inputs(
                input_file=input_path,
                output_file=os.path.join(temp_dir, "prepared.csv"),
                pairs_file=os.path.join(temp_dir, "pairs.tsv"),
                hla_dir=hla_dir,
                default_sample="PT1",
                input_source="custom_peptides",
            )

            self.assertEqual(len(prepared_df), 4)
            self.assertEqual(len(pairs_df), 4)
            self.assertTrue(all(prepared_df["InputSource"] == "custom_peptides"))
            self.assertEqual(
                set(prepared_df["Peptide"]),
                {"PEPTIDE1", "AABC"},
            )

    def test_mixmhc2_csv_input_uses_hla_column(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            input_path = os.path.join(temp_dir, "custom.csv")
            with open(input_path, "w") as handle:
                handle.write(
                    "Sample,Peptide,HLA,Protein,Gene,Source\n"
                    "PT1,AA-1,DRB1*15:01,TP53,TP53,viral\n"
                )

            prepared_df, pairs_df = prepare_mixmhc2_inputs(
                input_file=input_path,
                output_file=os.path.join(temp_dir, "prepared.csv"),
                pairs_file=os.path.join(temp_dir, "pairs.tsv"),
                hla_dir=temp_dir,
                default_sample="PT1",
                input_source="custom_peptides",
            )

            self.assertEqual(len(prepared_df), 1)
            self.assertEqual(prepared_df.loc[0, "HLA"], "DRB11501")
            self.assertEqual(prepared_df.loc[0, "Source"], "viral")
            self.assertEqual(pairs_df.loc[0, "Peptide"], "AA1")

    def test_mixmhc2_tsv_input_keeps_optional_origin_and_fallback_hla(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            input_path = os.path.join(temp_dir, "custom.tsv")
            with open(input_path, "w") as handle:
                handle.write("Sample\tPeptide\tOrigin\tSource\n")
                handle.write("PT1\tAA-1\tviral\tcustom\n")
            with open(os.path.join(temp_dir, "PT1_hla.txt"), "w") as handle:
                handle.write("HLA-DPA1*01:03\n")

            prepared_df, pairs_df = prepare_mixmhc2_inputs(
                input_file=input_path,
                output_file=os.path.join(temp_dir, "prepared.csv"),
                pairs_file=os.path.join(temp_dir, "pairs.tsv"),
                hla_dir=temp_dir,
                default_sample="PT1",
                input_source="custom_peptides",
            )

            self.assertEqual(len(prepared_df), 1)
            self.assertEqual(prepared_df.loc[0, "Peptide"], "AA1")
            self.assertEqual(prepared_df.loc[0, "HLA"], "DPA10103")
            self.assertEqual(prepared_df.loc[0, "Source"], "custom")
            self.assertEqual(prepared_df.loc[0, "Origin"], "viral")
            self.assertEqual(pairs_df.loc[0, "HLA"], "DPA1_01_03")

    def test_mixmhc2_merge_returns_original_columns_plus_scores(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            input_path = os.path.join(temp_dir, "input.csv")
            output_path = os.path.join(temp_dir, "output.csv")
            pred_path = os.path.join(temp_dir, "pred.csv")

            pd.DataFrame({
                "Sample": ["PT1"],
                "Peptide": ["AA1"],
                "HLA": ["DRB1_15_01"],
                "Source": ["custom"],
            }).to_csv(input_path, index=False)

            with open(pred_path, "w") as handle:
                handle.write("HLA,Peptide,Score,Percentile\nDRB1_15_01,AA1,0.77,1.2\n")

            add_mixmhc2_scores(
                df=pd.read_csv(input_path),
                prediction_file=pred_path,
            ).to_csv(output_path, index=False)

            merged = pd.read_csv(output_path)
            self.assertTrue({"Sample", "Peptide", "HLA", "Source"}.issubset(merged.columns))
            self.assertEqual(merged.loc[0, "MixMHC2pred_Score"], 0.77)
            self.assertEqual(merged.loc[0, "MixMHC2pred_PercentileRank"], 1.2)

    def test_parse_mixmhc2_output(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as handle:
            handle.write("HLA,Peptide,%Rank,PredScore\nDRB1_15_01,AA1,1.2,0.77\n")
            path = handle.name

        try:
            parsed = parse_mixmhc2_output(path)
            self.assertEqual(parsed.loc[0, "MixMHC2pred_Score"], 0.77)
            self.assertEqual(parsed.loc[0, "MixMHC2pred_PercentileRank"], 1.2)
            self.assertEqual(parsed.loc[0, "Peptide"], "AA1")
        finally:
            os.unlink(path)


if __name__ == "__main__":
    unittest.main()
