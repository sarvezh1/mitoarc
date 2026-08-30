import csv
import gzip
import importlib.machinery
import importlib.util
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def load_python(name, path):
    loader = importlib.machinery.SourceFileLoader(name, str(path))
    spec = importlib.util.spec_from_loader(name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


b12b = load_python("block12b_analysis", ROOT / "benchmark/scripts/block12b_analysis.py")
numt = load_python("numt_evidence", ROOT / "bin/numt_evidence")


def rows(name):
    with (ROOT / "benchmark/summaries/block12b" / name).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class StrategyCTests(unittest.TestCase):
    def test_snv_retained_case(self):
        self.assertEqual(numt.c_status(True, (5, 5, 0, 0)), "RETAINED")

    def test_snv_flagged_case(self):
        self.assertEqual(numt.c_status(True, (5, 0, 2, 3)), "FLAGGED_NUMT_EVIDENCE")

    def test_indel_unresolved_scope_case(self):
        self.assertEqual(numt.c_status(True, ("NA", "NA", "NA", "NA")), "UNRESOLVED_SUPPORT")

    def test_indel_cigar_support_extraction(self):
        insertion = ["r1", "0", "chrM", "95", "60", "6M1I9M", "*", "0", "0", "AAAAAATAAAAAAAAA", "I" * 16]
        deletion = ["r2", "0", "chrM", "95", "60", "6M1D9M", "*", "0", "0", "A" * 15, "I" * 15]
        self.assertTrue(b12b.record_supports_indel(insertion, 100, "A", "AT"))
        self.assertTrue(b12b.record_supports_indel(deletion, 100, "AC", "A"))

    def test_normalized_key_matching(self):
        refs = {"chrM": "ACGTAAAAAACGT"}
        self.assertEqual(
            b12b.bh.normalize_allele("chrM", 6, "AA", "A", refs),
            b12b.bh.normalize_allele("chrM", 5, "AA", "A", refs),
        )

    def test_observed_indel_trace_is_unresolved_despite_support(self):
        trace = rows("block12b_strategy_c_indel_trace.tsv")
        self.assertTrue(all(int(row["bam_support_alt_fragments"]) > 0 for row in trace))
        self.assertTrue(all(row["pipeline_supporting_fragments"] == "NA" for row in trace))
        self.assertTrue(all(row["final_strategy_c_status"] == "UNRESOLVED_SUPPORT" for row in trace))


class CircularMatchedTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = rows("block12b_circular_matched_summary.tsv")

    def test_matched_mode_cardinality(self):
        detail = [row for row in self.data if row["summary_level"] == "TRUTH_ALLELE"]
        pooled = [row for row in self.data if row["summary_level"] == "POOLED"]
        self.assertEqual(len(detail), 216)
        self.assertEqual(len(pooled), 8)
        self.assertEqual({row["mode"] for row in detail}, {"STANDARD_ONLY_MATCHED", "STANDARD_PLUS_SHIFTED_MATCHED"})

    def test_non_boundary_equivalence(self):
        pooled = [row for row in self.data if row["summary_level"] == "POOLED" and row["boundary_class"] == "NON_BOUNDARY"]
        self.assertTrue(all(float(row["recall_delta"]) == 0 for row in pooled))

    def test_boundary_shifted_evidence_inclusion(self):
        pooled = [row for row in self.data if row["summary_level"] == "POOLED" and row["boundary_class"] == "BOUNDARY"]
        self.assertEqual({row["variant_class"] for row in pooled}, {"SNV", "indel"})
        self.assertEqual({row["mode"] for row in pooled}, {"STANDARD_ONLY_MATCHED", "STANDARD_PLUS_SHIFTED_MATCHED"})
        self.assertTrue(all(float(row["recall_delta"]) == 0 for row in pooled))


class MutserveIndelTests(unittest.TestCase):
    def test_raw_indel_parsing(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "calls.vcf.gz"
            with gzip.open(path, "wt") as out:
                out.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                out.write("chrM\t10\t.\tA\tAT,ATT\t.\tPASS\t.\n")
            parsed = b12b.open_vcf(path)
        self.assertEqual([(x["pos"], x["ref"], x["alt"]) for x in parsed], [(10, "A", "AT"), (10, "A", "ATT")])

    def test_normalized_indel_representation(self):
        refs = {"chrM": "ACGTAAAAAACGT"}
        left = b12b.bh.normalize_allele("chrM", 6, "AA", "A", refs)
        self.assertEqual(left, (4, "TA", "T"))

    def test_no_silent_allele_loss(self):
        trace = rows("block12b_mutserve_indel_trace.tsv")
        self.assertEqual(len(trace), 216)
        self.assertTrue(all(row["raw_native_emitted"] == "false" for row in trace))
        self.assertTrue(all(row["canonical_emitted"] == "false" for row in trace))
        self.assertTrue(all(row["loss_stage"] == "CALLER_NOT_EMITTED" for row in trace))

    def test_positive_control_classification(self):
        controls = rows("block12b_mutserve_indel_positive_control.tsv")
        self.assertEqual({row["variant_class"] for row in controls}, {"insertion", "deletion"})
        self.assertTrue(all(row["gatk_emitted"] == "true" for row in controls))
        self.assertTrue(all(row["mutserve_native_emitted"] == "false" for row in controls))
        self.assertTrue(all(row["classification"] == "MUTSERVE2_INDEL_NOT_EMITTED_BY_FROZEN_MODE" for row in controls))


if __name__ == "__main__":
    unittest.main(verbosity=2)
