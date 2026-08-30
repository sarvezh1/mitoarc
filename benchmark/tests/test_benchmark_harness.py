#!/usr/bin/env python3
import argparse
import csv
import gzip
import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "benchmark_harness.py"
SPEC = importlib.util.spec_from_file_location("benchmark_harness", SCRIPT)
bh = importlib.util.module_from_spec(SPEC)
assert SPEC.loader
SPEC.loader.exec_module(bh)


class HarnessTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.reference = self.root / "reference.fa"
        self.reference.write_text(">chrM\nACGTAAAAAACGTACGTACGT\n", encoding="utf-8")
        self.truth = self.root / "truth.tsv"
        self.write_truth([
            {
                "truth_id": "t1", "chrom": "chrM", "pos": "2", "ref": "C", "alt": "T",
                "variant_class": "SNV", "expected_vaf": "0.2", "boundary_class": "NON_BOUNDARY",
                "numt_context": "NO_NUMT_CHALLENGE", "source": "test", "notes": "NA",
            }
        ])
        self.input_bam = self.root / "input.bam"
        self.input_bam.write_bytes(b"fixture")

    def tearDown(self):
        self.temp.cleanup()

    def write_truth(self, rows, path=None):
        path = path or self.truth
        bh.write_tsv(path, bh.TRUTH_FIELDS, rows)
        return path

    def valid_manifest_row(self):
        return {
            "benchmark_id": "d100_vaf020_A_gatk_circ_r1", "sample_id": "S1",
            "dataset_type": "real_validation", "reference": str(self.reference),
            "input_bam_or_fastq": str(self.input_bam), "truth_vcf": str(self.truth),
            "depth_target": "100", "vaf_target": "0.20", "vaf_bin": "VAF_0_10_0_50",
            "variant_class": "SNV", "variant_position": "2", "boundary_class": "NON_BOUNDARY",
            "numt_context": "NO_NUMT_CHALLENGE", "caller": "gatk", "numt_strategy": "A",
            "circular_mode": "standard_plus_shifted", "replicate": "1", "seed": "42",
            "numt_truth": "NA", "numt_design": "NA",
        }

    def write_manifest(self, rows, fields=None):
        path = self.root / "manifest.tsv"
        bh.write_tsv(path, fields or bh.MANIFEST_FIELDS, rows)
        return path

    def assert_manifest_failure(self, field, value):
        row = self.valid_manifest_row()
        row[field] = value
        with self.assertRaises(bh.HarnessError):
            bh.validate_manifest(self.write_manifest([row]))

    def test_valid_manifest(self):
        self.assertEqual(len(bh.validate_manifest(self.write_manifest([self.valid_manifest_row()]))), 1)

    def test_malformed_manifest(self):
        with self.assertRaises(bh.HarnessError):
            bh.validate_manifest(self.write_manifest([self.valid_manifest_row()], bh.MANIFEST_FIELDS[:-1]))

    def test_duplicate_benchmark_id(self):
        row = self.valid_manifest_row()
        with self.assertRaisesRegex(bh.HarnessError, "duplicate benchmark_id"):
            bh.validate_manifest(self.write_manifest([row, row]))

    def test_missing_truth(self):
        self.assert_manifest_failure("truth_vcf", str(self.root / "missing.tsv"))

    def test_generated_input_requires_provenance(self):
        row = self.valid_manifest_row()
        row["dataset_type"] = "downsampled"
        with self.assertRaisesRegex(bh.HarnessError, "missing its required provenance"):
            bh.validate_manifest(self.write_manifest([row]))

    def test_unsupported_caller(self):
        self.assert_manifest_failure("caller", "unknown")

    def test_unsupported_numt_strategy(self):
        self.assert_manifest_failure("numt_strategy", "D")

    def test_invalid_circular_mode(self):
        self.assert_manifest_failure("circular_mode", "circular_magic")

    def test_invalid_vaf(self):
        row = self.valid_manifest_row()
        row["vaf_target"] = "1.2"
        with self.assertRaises(bh.HarnessError):
            bh.validate_manifest(self.write_manifest([row]))

    def test_duplicate_truth_id(self):
        fields, rows = bh.read_tsv(self.truth)
        self.write_truth([rows[0], rows[0]])
        with self.assertRaisesRegex(bh.HarnessError, "duplicate or empty truth_id"):
            bh.validate_truth(self.truth, self.reference)

    def test_truth_ref_mismatch(self):
        fields, rows = bh.read_tsv(self.truth)
        rows[0]["ref"] = "A"
        self.write_truth(rows)
        with self.assertRaisesRegex(bh.HarnessError, "REF mismatch"):
            bh.validate_truth(self.truth, self.reference)

    def test_invalid_coordinate(self):
        fields, rows = bh.read_tsv(self.truth)
        rows[0]["pos"] = "999"
        self.write_truth(rows)
        with self.assertRaisesRegex(bh.HarnessError, "invalid mitochondrial coordinate"):
            bh.validate_truth(self.truth, self.reference)

    def test_invalid_truth_vaf(self):
        fields, rows = bh.read_tsv(self.truth)
        rows[0]["expected_vaf"] = "-0.1"
        self.write_truth(rows)
        with self.assertRaisesRegex(bh.HarnessError, "impossible expected_vaf"):
            bh.validate_truth(self.truth, self.reference)

    def test_indel_left_normalization(self):
        truth = self.write_truth([{
            "truth_id": "indel", "chrom": "chrM", "pos": "5", "ref": "AA", "alt": "A",
            "variant_class": "indel", "expected_vaf": "0.1", "boundary_class": "NON_BOUNDARY",
            "numt_context": "NO_NUMT_CHALLENGE", "source": "test", "notes": "NA",
        }], self.root / "indel.tsv")
        row = bh.validate_truth(truth, self.reference)[0]
        self.assertEqual((row["pos"], row["ref"], row["alt"]), (4, "TA", "T"))

    def test_exact_matching_and_fp_fn(self):
        calls = self.root / "calls.tsv"
        bh.write_tsv(calls, bh.CALL_FIELDS, [
            {"benchmark_id": "b", "caller": "gatk", "chrom": "chrM", "pos": 2, "ref": "C", "alt": "T", "variant_key": "chrM:2:C:T", "filter_status": "PASS", "called_vaf": "0.25", "called_depth": "100", "truth_status": "UNMATCHED"},
            {"benchmark_id": "b", "caller": "gatk", "chrom": "chrM", "pos": 3, "ref": "G", "alt": "A", "variant_key": "chrM:3:G:A", "filter_status": "PASS", "called_vaf": "0.1", "called_depth": "100", "truth_status": "UNMATCHED"},
        ])
        row = self.valid_manifest_row()
        row["benchmark_id"] = "b"
        summary, details = bh.match_one(row, self.truth, self.reference, calls)
        self.assertEqual((summary["tp"], summary["fp"], summary["fn"]), (1, 1, 0))
        self.assertEqual(summary["precision"], "0.5")
        tp = next(item for item in details if item["truth_status"] == "TP")
        self.assertAlmostEqual(float(tp["absolute_vaf_error"]), 0.05)
        self.assertAlmostEqual(float(tp["relative_vaf_error"]), 0.25)

    def test_fn_accounting(self):
        calls = self.root / "empty_calls.tsv"
        bh.write_tsv(calls, bh.CALL_FIELDS, [])
        summary, _ = bh.match_one(self.valid_manifest_row(), self.truth, self.reference, calls)
        self.assertEqual((summary["tp"], summary["fp"], summary["fn"]), (0, 0, 1))
        self.assertEqual(summary["precision"], "NA")
        self.assertEqual(summary["recall"], "0")

    def test_division_by_zero_metrics(self):
        self.assertEqual(bh.safe_metric(0, 0), "NA")
        self.assertEqual(bh.safe_metric(0, 2), "0")

    def test_deterministic_seed(self):
        truth = self.write_truth([{
            "truth_id": "t", "chrom": "chrM", "pos": "10", "ref": "A", "alt": "G",
            "variant_class": "SNV", "expected_vaf": "0.2", "boundary_class": "NON_BOUNDARY",
            "numt_context": "NO_NUMT_CHALLENGE", "source": "test", "notes": "NA",
        }], self.root / "seed_truth.tsv")
        hashes = []
        for name in ["one", "two"]:
            args = argparse.Namespace(
                benchmark_id="seed", sample="S", reference=str(self.reference), truth=str(truth),
                output_prefix=str(self.root / name), depth=100, seed=99, mt_contig="chrM",
                read_length=20, fragment_length=20,
            )
            bh.deterministic_fastq(args)
            hashes.append((bh.sha256(self.root / f"{name}_R1.fastq.gz"), bh.sha256(self.root / f"{name}_R2.fastq.gz")))
        self.assertEqual(hashes[0], hashes[1])

    def test_resource_schema(self):
        self.assertEqual(bh.RESOURCE_FIELDS, [
            "benchmark_id", "attempt", "wall_clock_seconds", "cpu_seconds", "cpu_hours",
            "peak_rss_mb", "max_container_memory_mb", "disk_input_bytes", "disk_work_bytes_peak",
            "disk_output_bytes", "process_count", "cached_process_count", "submitted_process_count",
            "nextflow_status",
        ])

    def test_resume_attempt_numbering(self):
        run = self.root / "run"
        (run / "nxf" / "attempt-001").mkdir(parents=True)
        (run / "nxf" / "attempt-002").mkdir()
        self.assertEqual(bh.latest_attempt(run), 2)

    def test_resume_uses_benchmark_specific_run_name(self):
        run = self.root / "run"
        attempt = run / "nxf" / "attempt-001"
        attempt.mkdir(parents=True)
        (attempt / "command.stdout.log").write_text(
            "Launching `/pipeline/main.nf` [specific_run]\n", encoding="utf-8"
        )
        (attempt / "run_metadata.json").write_text('{"status":"PASS"}\n', encoding="utf-8")
        self.assertEqual(bh.prior_successful_run_name(run), "specific_run")


if __name__ == "__main__":
    unittest.main(verbosity=2)
