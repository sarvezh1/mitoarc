#!/usr/bin/env python3
"""Extract frozen controlled-benchmark callsets and compute exact truth metrics."""

from __future__ import annotations

import argparse
import csv
import gzip
import importlib.util
import json
import math
import statistics
import hashlib
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BENCH = ROOT / "benchmark"
MANIFEST = BENCH / "manifests/block12_primary_matrix.tsv"
RUNS = BENCH / "runs"
METRICS = BENCH / "metrics/block12"
SUMMARY = BENCH / "summaries/block12"

spec = importlib.util.spec_from_file_location("benchmark_harness", BENCH / "scripts/benchmark_harness.py")
bh = importlib.util.module_from_spec(spec)
assert spec.loader
spec.loader.exec_module(bh)

CALLSETS = [
    ("gatk", "A", "standard_plus_shifted", "gatk_canonical_A", "variants/gatk/{s}.gatk.filtered.vcf.gz", "vcf"),
    ("gatk", "B", "standard_plus_shifted", "gatk_canonical_B", "variants/gatk/{s}_strategy_b.gatk.filtered.vcf.gz", "vcf"),
    ("gatk", "C", "standard_plus_shifted", "gatk_evidence_C", "variants/gatk/{s}.gatk.filtered.vcf.gz", "vcf"),
    ("gatk", "A", "standard_only", "gatk_standard_raw_intermediate", "variants/gatk/intermediate/{s}.standard.gatk.raw.vcf.gz", "vcf"),
    ("mutserve2", "A", "standard_only", "mutserve2_canonical_A", "variants/mutserve2/{s}.mutserve2.calls.tsv", "mutserve"),
    ("mutserve2", "B", "standard_only", "mutserve2_canonical_B", "variants/mutserve2/{s}_strategy_b.mutserve2.calls.tsv", "mutserve"),
    ("mutserve2", "C", "standard_only", "mutserve2_evidence_C", "variants/mutserve2/{s}.mutserve2.calls.tsv", "mutserve"),
]
ALL_CALL_FIELDS = ["benchmark_id", "sample_id", "replicate", "depth_target", "caller", "numt_strategy",
                   "circular_mode", "callset_source", "chrom", "pos", "ref", "alt", "variant_key",
                   "filter_status", "callable", "called_vaf", "called_depth", "truth_status", "truth_id",
                   "fp_context", "evidence_low", "evidence_ambiguous", "evidence_high"]
VARIANT_FIELDS = ["benchmark_id", "sample_id", "replicate", "seed", "depth_target", "caller", "numt_strategy",
                  "circular_mode", "callset_source", "truth_id", "variant_key", "chrom", "pos", "ref", "alt",
                  "variant_class", "boundary_class", "distance_to_boundary", "numt_context", "vaf_bin",
                  "target_vaf", "achieved_vaf", "truth_status", "called_vaf", "absolute_vaf_error",
                  "relative_vaf_error", "called_depth", "filter_status"]
CALL_METRIC_FIELDS = ["benchmark_id", "replicate", "depth_target", "caller", "numt_strategy", "circular_mode",
                      "callset_source", "vaf_bin", "variant_class", "boundary_class", "numt_context", "tp", "fp",
                      "fn", "precision", "recall", "f1", "numt_associated_fp", "ambiguous_signal_fp"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows({field: row.get(field, "NA") for field in fields} for row in rows)


def manifest_rows(selected: set[str]) -> list[dict]:
    rows = bh.validate_manifest(MANIFEST)
    unknown = selected - {r["benchmark_id"] for r in rows}
    if unknown: raise RuntimeError(f"unknown benchmark IDs: {sorted(unknown)}")
    return [r for r in rows if not selected or r["benchmark_id"] in selected]


def active(filter_status: str, strategy: str) -> bool:
    if strategy == "C":
        return filter_status == "RETAINED"
    return filter_status in {"PASS", ".", "SITE", "CALLED", "RETAINED", "RAW_UNFILTERED"}


def evidence_map(path: Path, caller: str) -> dict[str, dict]:
    result = {}
    for row in read_tsv(path):
        if row["caller"] != caller: continue
        key = f"{row['contig']}:{row['position']}:{row['ref']}:{row['alt']}"
        result[key] = row
    return result


def nominal_target(notes: str) -> str:
    for item in notes.split(";"):
        if item.startswith("nominal_target="): return item.split("=", 1)[1]
    return "NA"


def distance(notes: str, pos: int) -> int:
    for item in notes.split(";"):
        if item.startswith("distance_to_boundary="): return int(item.split("=", 1)[1])
    return min(pos - 1, 16569 - pos)


def metric(n: int, d: int) -> str:
    return "NA" if d == 0 else f"{n / d:.10g}"


def f1_value(precision: str, recall: str) -> str:
    if "NA" in {precision, recall}: return "NA"
    p, r = float(precision), float(recall)
    return "0" if p + r == 0 else f"{2*p*r/(p+r):.10g}"


def extract_one(row: dict) -> tuple[list[dict], list[dict], list[dict]]:
    bid, sample = row["benchmark_id"], row["sample_id"]
    results = RUNS / bid / "results"
    reference = Path(row["reference"])
    refs = bh.parse_fasta(reference)
    truth = bh.validate_truth(Path(row["truth_vcf"]), reference)
    truth_map = {str(t["variant_key"]): t for t in truth}
    challenge_rows = read_tsv(BENCH / "truth/block12/block12_numt_challenges.tsv")
    challenge_map = {f"chrM:{r['pos']}:{r['ref']}:{r['alt']}": r for r in challenge_rows}
    evidence_path = results / "variants/numt" / f"{sample}.numt_variant_evidence.tsv"
    all_calls, variant_rows, call_metrics = [], [], []
    for caller, strategy, circular, source_label, pattern, source_type in CALLSETS:
        source = results / pattern.format(s=sample)
        if not source.is_file(): raise RuntimeError(f"{bid}: missing {source_label}: {source}")
        calls = bh.normalize_vcf(source, bid, caller) if source_type == "vcf" else bh.normalize_mutserve(source, bid)
        ev = evidence_map(evidence_path, caller) if strategy == "C" else {}
        call_map = {}
        for call in calls:
            chrom, pos, ref, alt = str(call["chrom"]), int(call["pos"]), str(call["ref"]), str(call["alt"])
            pos, ref, alt = bh.normalize_allele(chrom, pos, ref, alt, refs)
            key = f"{chrom}:{pos}:{ref}:{alt}"
            call.update(pos=pos, ref=ref, alt=alt, variant_key=key)
            evidence = ev.get(key, {})
            if strategy == "C": call["filter_status"] = evidence.get("strategy_c_status", "NOT_EVALUATED")
            is_active = active(str(call["filter_status"]), strategy)
            if is_active:
                if key in call_map: raise RuntimeError(f"{bid}/{source_label}: duplicate allele {key}")
                call_map[key] = call
            truth_status = ("TP" if is_active else "FILTERED_TRUE") if key in truth_map else (
                "FP" if is_active else "FILTERED_NONTRUTH")
            challenge = challenge_map.get(key)
            fp_context = (challenge["truth_class"] if challenge else "OTHER_NONTRUTH") if key not in truth_map else "NA"
            all_calls.append({"benchmark_id": bid, "sample_id": sample, "replicate": row["replicate"],
                              "depth_target": row["depth_target"], "caller": caller, "numt_strategy": strategy,
                              "circular_mode": circular, "callset_source": source_label, "chrom": chrom, "pos": pos,
                              "ref": ref, "alt": alt, "variant_key": key, "filter_status": call["filter_status"],
                              "callable": str(is_active).lower(), "called_vaf": call["called_vaf"],
                              "called_depth": call["called_depth"], "truth_status": truth_status,
                              "truth_id": truth_map.get(key, {}).get("truth_id", "NA"), "fp_context": fp_context,
                              "evidence_low": evidence.get("low_numt_evidence_fragments", "NA"),
                              "evidence_ambiguous": evidence.get("ambiguous_numt_evidence_fragments", "NA"),
                              "evidence_high": evidence.get("high_numt_evidence_fragments", "NA")})
        common = {"benchmark_id": bid, "sample_id": sample, "replicate": row["replicate"], "seed": row["seed"],
                  "depth_target": row["depth_target"], "caller": caller, "numt_strategy": strategy,
                  "circular_mode": circular, "callset_source": source_label}
        for key, truth_row in truth_map.items():
            call = call_map.get(key)
            achieved = float(truth_row["expected_vaf"])
            called = None if not call or call["called_vaf"] == "NA" else float(call["called_vaf"])
            absolute = None if called is None else abs(called - achieved)
            relative = None if absolute is None or achieved == 0 else absolute / achieved
            variant_rows.append({**common, "truth_id": truth_row["truth_id"], "variant_key": key,
                "chrom": truth_row["chrom"], "pos": truth_row["pos"], "ref": truth_row["ref"], "alt": truth_row["alt"],
                "variant_class": truth_row["variant_class"], "boundary_class": truth_row["boundary_class"],
                "distance_to_boundary": distance(truth_row["notes"], int(truth_row["pos"])),
                "numt_context": truth_row["numt_context"], "vaf_bin": bh.vaf_bin(achieved),
                "target_vaf": nominal_target(truth_row["notes"]), "achieved_vaf": truth_row["expected_vaf"],
                "truth_status": "TP" if call else "FN", "called_vaf": call["called_vaf"] if call else "NA",
                "absolute_vaf_error": "NA" if absolute is None else f"{absolute:.10g}",
                "relative_vaf_error": "NA" if relative is None else f"{relative:.10g}",
                "called_depth": call["called_depth"] if call else "NA", "filter_status": call["filter_status"] if call else "NA"})
        fps = [(key, c) for key, c in call_map.items() if key not in truth_map]
        # Allocate FPs to their challenge stratum when known; OTHER is retained
        # explicitly and never silently assigned to a truth VAF bin.
        strata = defaultdict(lambda: {"tp": 0, "fn": 0, "fp": 0, "numt_fp": 0, "amb_fp": 0})
        for vr in variant_rows[-len(truth):]:
            skey = (vr["vaf_bin"], vr["variant_class"], vr["boundary_class"], vr["numt_context"])
            strata[skey]["tp" if vr["truth_status"] == "TP" else "fn"] += 1
        for key, call in fps:
            challenge = challenge_map.get(key)
            if challenge:
                vb = bh.vaf_bin(float(challenge["target_vaf"]))
                context = "NUMT_LIKE" if challenge["truth_class"] == "NUMT_FALSE" else "AMBIGUOUS_CONTEXT"
                skey = (vb, "SNV", "NON_BOUNDARY", context)
                strata[skey]["fp"] += 1
                strata[skey]["numt_fp" if challenge["truth_class"] == "NUMT_FALSE" else "amb_fp"] += 1
            else:
                strata[("NA", "MIXED", "MIXED", "OTHER_NONTRUTH")]["fp"] += 1
        for skey, counts in sorted(strata.items()):
            p = metric(counts["tp"], counts["tp"] + counts["fp"])
            r = metric(counts["tp"], counts["tp"] + counts["fn"])
            call_metrics.append({**common, "vaf_bin": skey[0], "variant_class": skey[1],
                                 "boundary_class": skey[2], "numt_context": skey[3], **counts,
                                 "precision": p, "recall": r, "f1": f1_value(p, r),
                                 "numt_associated_fp": counts["numt_fp"], "ambiguous_signal_fp": counts["amb_fp"]})
    return all_calls, variant_rows, call_metrics


def extract(args) -> None:
    selected = set(args.benchmark_id or [])
    for row in manifest_rows(selected):
        calls, variants, metrics = extract_one(row)
        out = METRICS / row["benchmark_id"]
        write_tsv(out / "block12_all_calls.tsv", ALL_CALL_FIELDS, calls)
        write_tsv(out / "block12_variant_metrics.tsv", VARIANT_FIELDS, variants)
        write_tsv(out / "block12_call_metrics.tsv", CALL_METRIC_FIELDS, metrics)
        print(f"PASS\t{row['benchmark_id']}\tcall_records={len(calls)}\ttruth_comparisons={len(variants)}\tstrata={len(metrics)}")


def pilot_qc(args) -> None:
    selected = set(args.benchmark_id or [])
    failures = []
    for row in manifest_rows(selected):
        bid, sample = row["benchmark_id"], row["sample_id"]
        base = METRICS / bid
        calls, variants, metrics = read_tsv(base / "block12_all_calls.tsv"), read_tsv(base / "block12_variant_metrics.tsv"), read_tsv(base / "block12_call_metrics.tsv")
        sources = {r["callset_source"] for r in variants}
        if sources != {x[3] for x in CALLSETS}: failures.append(f"{bid}: callset sources {sorted(sources)}")
        if len(variants) != 36 * len(CALLSETS): failures.append(f"{bid}: truth comparisons={len(variants)}")
        if any(r["truth_status"] not in {"TP", "FN"} for r in variants): failures.append(f"{bid}: invalid truth status")
        preflight = RUNS / bid / "results/variants/numt" / f"{sample}.numt_fixture_validation.tsv"
        evidence = RUNS / bid / "results/variants/numt" / f"{sample}.numt_variant_evidence.tsv"
        if not preflight.is_file() or len(read_tsv(preflight)) != 12: failures.append(f"{bid}: NUMT preflight incomplete")
        if not evidence.is_file(): failures.append(f"{bid}: NUMT evidence missing")
        if not metrics: failures.append(f"{bid}: metrics empty")
        print(f"QC\t{bid}\tcallsets={len(sources)}\ttruth_comparisons={len(variants)}\tcall_records={len(calls)}")
    if failures: raise RuntimeError("; ".join(failures))
    print("PILOT_QC_PASS")


def quantiles(values: list[float]) -> dict[str, object]:
    if not values:
        return {"n": 0, "mean": "NA", "median": "NA", "sd": "NA", "min": "NA", "max": "NA", "q25": "NA", "q75": "NA"}
    ordered = sorted(values)
    def q(p):
        x = (len(ordered) - 1) * p
        lo, hi = math.floor(x), math.ceil(x)
        return ordered[lo] if lo == hi else ordered[lo] * (hi - x) + ordered[hi] * (x - lo)
    return {"n": len(values), "mean": f"{statistics.mean(values):.10g}", "median": f"{statistics.median(values):.10g}",
            "sd": "NA" if len(values) < 2 else f"{statistics.stdev(values):.10g}",
            "min": f"{min(values):.10g}", "max": f"{max(values):.10g}",
            "q25": f"{q(.25):.10g}", "q75": f"{q(.75):.10g}"}


def aggregate_counts(records: list[dict]) -> dict[str, object]:
    tp = sum(r.get("truth_status") == "TP" for r in records)
    fn = sum(r.get("truth_status") == "FN" for r in records)
    fp = sum(r.get("truth_status") == "FP" for r in records)
    p, rec = metric(tp, tp + fp), metric(tp, tp + fn)
    errors = [float(r["absolute_vaf_error"]) for r in records if r.get("absolute_vaf_error", "NA") != "NA"]
    return {"tp": tp, "fp": fp, "fn": fn, "precision": p, "recall": rec, "f1": f1_value(p, rec),
            "absolute_vaf_error_mean": "NA" if not errors else f"{statistics.mean(errors):.10g}",
            "absolute_vaf_error_median": "NA" if not errors else f"{statistics.median(errors):.10g}"}


def consolidate(args) -> None:
    rows = manifest_rows(set())
    all_calls, variants, call_metrics, run_metrics, truth_panel, provenance = [], [], [], [], [], []
    for row in rows:
        bid, sample = row["benchmark_id"], row["sample_id"]
        base = METRICS / bid
        all_calls += read_tsv(base / "block12_all_calls.tsv")
        variants += read_tsv(base / "block12_variant_metrics.tsv")
        call_metrics += read_tsv(base / "block12_call_metrics.tsv")
        resources = sorted(base.glob("benchmark_resource_metrics.attempt-*.tsv"))
        normal = [p for p in resources if int(p.stem.rsplit("-", 1)[-1]) == (2 if bid == "block12_R1_d100" else 1)]
        if len(normal) != 1: raise RuntimeError(f"{bid}: expected one primary resource record; got {normal}")
        run_metrics += read_tsv(normal[0])
        truth = bh.validate_truth(Path(row["truth_vcf"]), Path(row["reference"]))
        for item in truth:
            truth_panel.append({"benchmark_id": bid, "sample_id": sample, "replicate": row["replicate"],
                                "depth_target": row["depth_target"], "seed": row["seed"], **item,
                                "target_vaf": nominal_target(str(item["notes"])),
                                "achieved_vaf": item["expected_vaf"],
                                "distance_to_boundary": distance(str(item["notes"]), int(item["pos"]))})
        prov_path = Path(row["input_bam_or_fastq"].split(";")[0]).with_name(sample + ".provenance.tsv")
        prov = read_tsv(prov_path)[0]
        attempt = RUNS / bid / "nxf" / ("attempt-002" if bid == "block12_R1_d100" else "attempt-001") / "run_metadata.json"
        meta = json.loads(attempt.read_text())
        provenance.append({**prov, "replicate": row["replicate"], "depth_target": row["depth_target"],
                           "truth_checksum": sha256(Path(row["truth_vcf"])), "numt_truth_checksum": sha256(Path(row["numt_truth"])),
                           "numt_design_checksum": sha256(Path(row["numt_design"])), "manifest_checksum": meta["manifest_sha256"],
                           "pipeline_checksum": meta["pipeline_checksum"], "nextflow_run_name": meta["nextflow_run_name"],
                           "run_attempt": meta["attempt"], "run_metadata": str(attempt), "run_status": meta["status"]})
    write_tsv(SUMMARY / "block12_all_calls.tsv", ALL_CALL_FIELDS, all_calls)
    write_tsv(SUMMARY / "block12_variant_metrics.tsv", VARIANT_FIELDS, variants)
    write_tsv(SUMMARY / "block12_run_metrics.tsv", list(run_metrics[0]), run_metrics)
    truth_fields = ["benchmark_id", "sample_id", "replicate", "depth_target", "seed"] + bh.TRUTH_FIELDS + ["variant_key", "target_vaf", "achieved_vaf", "distance_to_boundary"]
    write_tsv(SUMMARY / "block12_truth_panel.tsv", truth_fields, truth_panel)
    write_tsv(SUMMARY / "block12_experiment_provenance.tsv", list(provenance[0]), provenance)

    # Explicit designed -> mapped -> extracted -> outcome observability chain.
    audit = []
    challenge_template = {r["case_id"]: r for r in read_tsv(BENCH / "truth/block12/block12_numt_challenges.tsv")}
    for row in rows:
        bid, sample = row["benchmark_id"], row["sample_id"]
        preflight = {r["read_group"].split("X", 1)[0]: r for r in read_tsv(RUNS / bid / "results/variants/numt" / f"{sample}.numt_fixture_validation.tsv")}
        calls = [r for r in all_calls if r["benchmark_id"] == bid]
        for case_id, challenge in challenge_template.items():
            observed = preflight.get(case_id)
            if observed is None: raise RuntimeError(f"{bid}: missing preflight case {case_id}")
            key = f"chrM:{challenge['pos']}:{challenge['ref']}:{challenge['alt']}"
            for caller in ("gatk", "mutserve2"):
                outcomes = {}
                for strategy in ("A", "B", "C"):
                    candidates = [r for r in calls if r["caller"] == caller and r["numt_strategy"] == strategy
                                  and r["variant_key"] == key and r["callset_source"] != "gatk_standard_raw_intermediate"]
                    if not candidates: outcomes[strategy] = "NOT_CALLED"
                    else:
                        c = candidates[0]
                        outcomes[strategy] = "CALLED_FALSE_SIGNAL" if c["callable"] == "true" else f"NOT_RETAINED:{c['filter_status']}"
                extracted = int(observed["extracted_fragments"]) > 0
                achieved_challenge_vaf = int(observed["pairs"]) / int(row["depth_target"])
                audit.append({"benchmark_id": bid, "replicate": row["replicate"], "depth": row["depth_target"],
                    "truth_id": case_id, "requested_vaf": challenge["target_vaf"],
                    "achieved_input_vaf": f"{achieved_challenge_vaf:.10g}", "vaf_bin": bh.vaf_bin(achieved_challenge_vaf),
                    "numt_context": "NUMT_LIKE" if challenge["truth_class"] == "NUMT_FALSE" else "AMBIGUOUS_CONTEXT",
                    "caller": caller, "designed_challenge": "true", "realized_challenge": str(int(observed["fullref_fragments_observed"]) > 0).lower(),
                    "nuclear_alt_generated": "true", "nuclear_alt_mapped": str(observed["has_nuclear_alternative"] == "true").lower(),
                    "nuclear_alt_extracted": str(extracted).lower(), "observable_numt_challenge": str(extracted).lower(),
                    "fullref_fragments_observed": observed["fullref_fragments_observed"], "extracted_fragments": observed["extracted_fragments"],
                    "fullref_primary_contig": observed["fullref_primary_contig"], "fullref_primary_mapq": observed["fullref_primary_mapq"],
                    "strategy_a_outcome": outcomes["A"], "strategy_b_outcome": outcomes["B"], "strategy_c_outcome": outcomes["C"]})
    audit_fields = ["benchmark_id", "replicate", "depth", "truth_id", "requested_vaf", "achieved_input_vaf", "vaf_bin", "numt_context", "caller",
                    "designed_challenge", "realized_challenge", "nuclear_alt_generated", "nuclear_alt_mapped", "nuclear_alt_extracted",
                    "observable_numt_challenge", "fullref_fragments_observed", "extracted_fragments", "fullref_primary_contig",
                    "fullref_primary_mapq", "strategy_a_outcome", "strategy_b_outcome", "strategy_c_outcome"]
    write_tsv(SUMMARY / "block12_numt_observability_audit.tsv", audit_fields, audit)

    # H1: exact true-MT retention plus observable false signal denominator.
    numt_rows = []
    for row in rows:
        bid = row["benchmark_id"]
        for caller in ("gatk", "mutserve2"):
            for strategy in ("A", "B", "C"):
                for context in ("NUMT_LIKE", "AMBIGUOUS_CONTEXT"):
                    for target in map(str, (0.005, 0.015, 0.035, 0.075, 0.25, 0.75)):
                        vr = [r for r in variants if r["benchmark_id"] == bid and r["caller"] == caller and r["numt_strategy"] == strategy
                              and r["callset_source"] != "gatk_standard_raw_intermediate" and r["numt_context"] == context and r["target_vaf"] == target]
                        ar = [r for r in audit if r["benchmark_id"] == bid and r["caller"] == caller and r["numt_context"] == context
                              and abs(float(challenge_template[r["truth_id"]]["target_vaf"]) - float(target)) < 1e-12]
                        tp = sum(x["truth_status"] == "TP" for x in vr); fn = sum(x["truth_status"] == "FN" for x in vr)
                        observable = sum(x["observable_numt_challenge"] == "true" for x in ar)
                        outcome_field = f"strategy_{strategy.lower()}_outcome"
                        fp = sum(x["observable_numt_challenge"] == "true" and x[outcome_field] == "CALLED_FALSE_SIGNAL" for x in ar)
                        unresolved = sum(x["filter_status"] == "UNRESOLVED_SUPPORT" for x in all_calls if x["benchmark_id"] == bid
                                         and x["caller"] == caller and x["numt_strategy"] == strategy and x["fp_context"] == ("NUMT_FALSE" if context == "NUMT_LIKE" else "AMBIGUOUS"))
                        p, rec = metric(tp, tp + fp), metric(tp, tp + fn)
                        achieved_true = vr[0]["achieved_vaf"] if vr else "NA"
                        achieved_false = ar[0]["achieved_input_vaf"] if ar else "NA"
                        numt_rows.append({"summary_level": "REPLICATE", "benchmark_id": bid, "replicate": row["replicate"],
                            "depth": row["depth_target"], "target_vaf": target, "achieved_true_vaf": achieved_true,
                            "achieved_challenge_vaf": achieved_false,
                            "vaf_bin": vr[0]["vaf_bin"] if vr else bh.vaf_bin(float(achieved_false)),
                            "numt_context": context, "caller": caller, "numt_strategy": strategy,
                            "designed_challenges": len(ar), "observable_challenges": observable,
                            "unobservable_challenges": len(ar)-observable, "true_mt_tp": tp, "true_mt_fn": fn,
                            "observable_numt_fp": fp, "precision": p, "recall": rec, "f1": f1_value(p, rec),
                            "unresolved_evidence": unresolved})
    numt_fields = list(numt_rows[0])
    write_tsv(SUMMARY / "block12_numt_summary.tsv", numt_fields, numt_rows)

    # H2 exact boundary/non-boundary comparison from the same reads.
    circular = []
    for keys, group in group_rows([r for r in variants if r["caller"] == "gatk" and r["numt_strategy"] == "A" and
                                   r["callset_source"] in {"gatk_canonical_A", "gatk_standard_raw_intermediate"}],
                                  ["replicate", "depth_target", "target_vaf", "vaf_bin", "variant_class", "boundary_class", "circular_mode", "callset_source"]):
        circular.append(dict(zip(["replicate", "depth", "target_vaf", "vaf_bin", "variant_class", "boundary_class", "circular_mode", "callset_source"], keys),
                             **aggregate_counts(group)))
    write_tsv(SUMMARY / "block12_circular_summary.tsv", list(circular[0]), circular)

    # H4 baseline caller comparison using shared controlled truth and exact FPs.
    caller_summary = []
    baseline_sources = {"gatk": "gatk_canonical_A", "mutserve2": "mutserve2_canonical_A"}
    for caller, source in baseline_sources.items():
        relevant = [r for r in variants if r["caller"] == caller and r["callset_source"] == source]
        for keys, group in group_rows(relevant, ["replicate", "depth_target", "target_vaf", "vaf_bin", "variant_class", "boundary_class", "numt_context"]):
            base = dict(zip(["replicate", "depth", "target_vaf", "vaf_bin", "variant_class", "boundary_class", "numt_context"], keys))
            # Challenge-key FPs are assigned through the call metric layer.
            matching_fp = sum(int(m["fp"]) for m in call_metrics if m["caller"] == caller and m["callset_source"] == source
                              and m["replicate"] == base["replicate"] and m["depth_target"] == base["depth"]
                              and m["vaf_bin"] == base["vaf_bin"] and m["variant_class"] == base["variant_class"]
                              and m["boundary_class"] == base["boundary_class"] and m["numt_context"] == base["numt_context"])
            counts = aggregate_counts(group); counts["fp"] = matching_fp
            counts["precision"] = metric(counts["tp"], counts["tp"] + matching_fp)
            counts["f1"] = f1_value(counts["precision"], counts["recall"])
            caller_summary.append({"summary_level": "REPLICATE", "caller": caller, **base, **counts})
    write_tsv(SUMMARY / "block12_caller_summary.tsv", list(caller_summary[0]), caller_summary)

    # H3 run resources joined to each canonical A/B/C scientific callset.
    resources = {r["benchmark_id"]: r for r in run_metrics}
    depth_rows = []
    for row in rows:
        bid = row["benchmark_id"]
        for caller, strategy, circular_mode, source, _, _ in CALLSETS:
            if source == "gatk_standard_raw_intermediate": continue
            vr = [r for r in variants if r["benchmark_id"] == bid and r["callset_source"] == source]
            calls_active = [r for r in all_calls if r["benchmark_id"] == bid and r["callset_source"] == source and r["callable"] == "true"]
            counts = aggregate_counts(vr + [{"truth_status": "FP"} for r in calls_active if r["truth_status"] == "FP"])
            res = resources[bid]
            depth_rows.append({"summary_level": "REPLICATE", "benchmark_id": bid, "replicate": row["replicate"],
                "depth": row["depth_target"], "caller": caller, "numt_strategy": strategy, "circular_mode": circular_mode,
                **counts, "wall_clock_seconds": res["wall_clock_seconds"], "cpu_hours": res["cpu_hours"],
                "peak_rss_mb": res["peak_rss_mb"], "work_bytes_peak": res["disk_work_bytes_peak"],
                "output_bytes": res["disk_output_bytes"], "runtime_ratio_vs_100x": "NA", "cpu_ratio_vs_100x": "NA",
                "work_storage_ratio_vs_100x": "NA"})
    for record in depth_rows:
        baseline = next(r for r in depth_rows if r["replicate"] == record["replicate"] and r["caller"] == record["caller"]
                        and r["numt_strategy"] == record["numt_strategy"] and r["depth"] == "100")
        for field, out in (("wall_clock_seconds", "runtime_ratio_vs_100x"), ("cpu_hours", "cpu_ratio_vs_100x"),
                           ("work_bytes_peak", "work_storage_ratio_vs_100x")):
            record[out] = "NA" if float(baseline[field]) == 0 else f"{float(record[field])/float(baseline[field]):.10g}"
    write_tsv(SUMMARY / "block12_depth_summary.tsv", list(depth_rows[0]), depth_rows)

    # Long-format descriptive statistics across R1/R2/R3.
    descriptive = []
    def add_descriptive(analysis, records, group_fields, metrics):
        for keys, group in group_rows(records, group_fields):
            base = {"analysis": analysis, **dict(zip(group_fields, keys))}
            for name in metrics:
                values = []
                for item in group:
                    value = item.get(name, "NA")
                    if value not in {"NA", "", None}:
                        try: values.append(float(value))
                        except ValueError: pass
                descriptive.append({**base, "metric": name, **quantiles(values)})
    add_descriptive("H1_NUMT", numt_rows,
                    ["depth", "target_vaf", "vaf_bin", "numt_context", "caller", "numt_strategy"],
                    ["observable_challenges", "unobservable_challenges", "true_mt_tp", "true_mt_fn", "observable_numt_fp", "precision", "recall", "f1", "unresolved_evidence"])
    add_descriptive("H2_CIRCULAR", circular,
                    ["depth", "target_vaf", "vaf_bin", "variant_class", "boundary_class", "circular_mode", "callset_source"],
                    ["tp", "fn", "recall", "absolute_vaf_error_mean"])
    add_descriptive("H3_DEPTH", depth_rows,
                    ["depth", "caller", "numt_strategy", "circular_mode"],
                    ["precision", "recall", "f1", "absolute_vaf_error_mean", "wall_clock_seconds", "cpu_hours", "peak_rss_mb", "work_bytes_peak", "output_bytes", "runtime_ratio_vs_100x", "cpu_ratio_vs_100x", "work_storage_ratio_vs_100x"])
    add_descriptive("H4_CALLER", caller_summary,
                    ["depth", "target_vaf", "vaf_bin", "variant_class", "boundary_class", "numt_context", "caller"],
                    ["tp", "fp", "fn", "precision", "recall", "f1", "absolute_vaf_error_mean"])
    descriptive_fields = ["analysis", "depth", "target_vaf", "vaf_bin", "variant_class", "boundary_class", "numt_context",
                          "caller", "numt_strategy", "circular_mode", "callset_source", "metric", "n", "mean", "median", "sd", "min", "max", "q25", "q75"]
    write_tsv(SUMMARY / "block12_descriptive_statistics.tsv", descriptive_fields, descriptive)
    print(f"PASS\tconsolidated calls={len(all_calls)} variants={len(variants)} audit={len(audit)} provenance={len(provenance)}")


def group_rows(rows: list[dict], fields: list[str]):
    groups = defaultdict(list)
    for row in rows: groups[tuple(row[field] for field in fields)].append(row)
    return sorted(groups.items())


def final_qc(args) -> None:
    errors = []
    expected = {"block12_primary_matrix.tsv", "block12_truth_panel.tsv", "block12_all_calls.tsv", "block12_variant_metrics.tsv",
                "block12_run_metrics.tsv", "block12_numt_summary.tsv", "block12_numt_observability_audit.tsv",
                "block12_circular_summary.tsv", "block12_depth_summary.tsv", "block12_caller_summary.tsv",
                "block12_failures.tsv", "block12_experiment_provenance.tsv"}
    for name in expected:
        if not (SUMMARY / name).is_file(): errors.append(f"missing {name}")
    truth = read_tsv(SUMMARY / "block12_truth_panel.tsv")
    variants = read_tsv(SUMMARY / "block12_variant_metrics.tsv")
    audit = read_tsv(SUMMARY / "block12_numt_observability_audit.tsv")
    runs = read_tsv(SUMMARY / "block12_run_metrics.tsv")
    if len(truth) != 648: errors.append(f"truth rows {len(truth)} != 648")
    if len(variants) != 4536: errors.append(f"variant rows {len(variants)} != 4536")
    if len(audit) != 432: errors.append(f"audit rows {len(audit)} != 432")
    if len(runs) != 18 or any(r["nextflow_status"] != "PASS" for r in runs): errors.append("run metrics incomplete")
    if len({(r["benchmark_id"], r["truth_id"]) for r in truth}) != 648: errors.append("truth IDs not unique per run")
    if len({r["seed"] for r in truth}) != 18: errors.append("dataset seeds not unique")
    if errors: raise RuntimeError("; ".join(errors))
    print("FINAL_MATRIX_QC_PASS\t18_runs\t648_truth_rows\t4536_truth_callset_comparisons\t432_observability_rows")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for command in ("extract", "pilot-qc", "consolidate", "final-qc"):
        item = sub.add_parser(command)
        item.add_argument("--benchmark-id", action="append")
    args = parser.parse_args()
    {"extract": extract, "pilot-qc": pilot_qc, "consolidate": consolidate, "final-qc": final_qc}[args.command](args)


if __name__ == "__main__":
    main()
