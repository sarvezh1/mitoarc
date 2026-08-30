import csv
import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


panel = load("block12_panel", ROOT / "benchmark/scripts/block12_panel.py")


class ControlledBenchmarkTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.sequence = panel.read_fasta(panel.SOURCE_REF)["chrM"]
        cls.truth = panel.panel_template(cls.sequence)

    def test_panel_balance(self):
        self.assertEqual(len(self.truth), 36)
        for target in panel.TARGETS:
            self.assertEqual(sum(r["target_vaf"] == target for r in self.truth), 6)

    def test_panel_no_read_overlap(self):
        positions = [int(r["pos"]) for r in self.truth]
        for i, a in enumerate(positions):
            for b in positions[i + 1:]:
                self.assertGreater(min((a-b) % len(self.sequence), (b-a) % len(self.sequence)), panel.READ_LENGTH)

    def test_reference_and_seed_contract(self):
        refs = panel.read_fasta(ROOT / "benchmark/datasets/block12/reference/block12_reference.fa")
        self.assertEqual(len(refs["chrM"]), 16569)
        self.assertTrue(all(f"chr{i}" in refs for i in range(1, 23)))
        seeds = [panel.dataset_seed(rep, depth) for rep in panel.REPLICATE_SEEDS for depth in panel.DEPTHS]
        self.assertEqual(len(seeds), len(set(seeds)))

    def test_completed_summary_cardinality(self):
        def rows(name):
            with (ROOT / "benchmark/summaries/block12" / name).open() as handle:
                return list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows("block12_truth_panel.tsv")), 648)
        self.assertEqual(len(rows("block12_variant_metrics.tsv")), 4536)
        self.assertEqual(len(rows("block12_numt_observability_audit.tsv")), 432)
        self.assertEqual(len(rows("block12_run_metrics.tsv")), 18)

    def test_observability_denominator(self):
        with (ROOT / "benchmark/summaries/block12/block12_numt_observability_audit.tsv").open() as handle:
            audit = list(csv.DictReader(handle, delimiter="\t"))
        unique = {(r["benchmark_id"], r["truth_id"]): r for r in audit}.values()
        observed = {}
        for depth in panel.DEPTHS:
            q = [r for r in unique if int(r["depth"]) == depth]
            observed[depth] = sum(r["observable_numt_challenge"] == "true" for r in q)
        self.assertEqual(observed, {100: 24, 250: 33, 500: 33, 1000: 36, 2000: 36, 5000: 36})


if __name__ == "__main__":
    unittest.main()
