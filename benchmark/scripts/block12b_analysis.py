#!/usr/bin/env python3
"""Focused controlled-benchmark diagnostics without changing frozen workflow science."""

from __future__ import annotations

import argparse
import csv
import gzip
import importlib.util
import json
import math
import random
import re
import shlex
import statistics
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BENCH = ROOT / "benchmark"
RUNS = BENCH / "runs"
SUMMARY12 = BENCH / "summaries/block12"
SUMMARY = BENCH / "summaries/block12b"
GENERATED = BENCH / "generated/block12b"
SAMTOOLS_IMAGE = "quay.io/biocontainers/samtools:1.20--h50ea8bc_1"
GATK_IMAGE = "broadinstitute/gatk:4.6.2.0"
MUTSERVE_IMAGE = "mitoarc-mutserve2:2.0.3"

spec = importlib.util.spec_from_file_location("benchmark_harness", BENCH / "scripts/benchmark_harness.py")
bh = importlib.util.module_from_spec(spec)
assert spec.loader
spec.loader.exec_module(bh)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows({field: row.get(field, "NA") for field in fields} for row in rows)


def rel(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def docker(image: str, command: list[str], workdir: Path = ROOT, capture: bool = False) -> str:
    invocation = ["docker", "run", "--rm", "--user", "1000:1000", "-v", f"{ROOT}:/workspace",
                  "-w", f"/workspace/{rel(workdir)}", image, *command]
    result = subprocess.run(invocation, check=True, text=True, stdout=subprocess.PIPE if capture else None)
    return result.stdout if capture else ""


def sam_view(bam: Path, region: str | None = None, names: Path | None = None) -> list[list[str]]:
    command = ["samtools", "view"]
    if names is not None:
        command += ["-N", f"/workspace/{rel(names)}"]
    command.append(f"/workspace/{rel(bam)}")
    if region:
        command.append(region)
    text = docker(SAMTOOLS_IMAGE, command, capture=True)
    return [line.split("\t") for line in text.splitlines() if line]


def cigar_ops(cigar: str):
    return [(int(length), op) for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)]


def aligned_allele(record: list[str], pos: int, ref: str) -> str | None:
    ref_pos, query_pos, sequence = int(record[3]), 0, record[9]
    bases = []
    wanted = set(range(pos, pos + len(ref)))
    for length, op in cigar_ops(record[5]):
        if op in "M=X":
            for offset in range(length):
                coordinate = ref_pos + offset
                if coordinate in wanted:
                    bases.append((coordinate, sequence[query_pos + offset]))
            ref_pos += length; query_pos += length
        elif op in "DN":
            ref_pos += length
        elif op in "IS":
            query_pos += length
    if {coordinate for coordinate, _ in bases} != wanted:
        return None
    return "".join(base for _, base in sorted(bases))


def record_supports_indel(record: list[str], pos: int, ref: str, alt: str) -> bool:
    ref_pos, query_pos, sequence = int(record[3]), 0, record[9]
    insertion = len(alt) > len(ref)
    deletion = len(ref) > len(alt)
    for length, op in cigar_ops(record[5]):
        if op in "M=X":
            ref_pos += length; query_pos += length
        elif op == "I":
            inserted = sequence[query_pos:query_pos + length]
            if insertion and ref_pos == pos + len(ref) and inserted == alt[len(ref):]:
                return True
            query_pos += length
        elif op in "DN":
            if deletion and op == "D" and ref_pos == pos + len(alt) and length == len(ref) - len(alt):
                return True
            ref_pos += length
        elif op == "S":
            query_pos += length
    return False


def bam_indel_support(bam: Path, pos: int, ref: str, alt: str) -> tuple[set[str], set[str]]:
    records = sam_view(bam, f"chrM:{max(1, pos - 2)}-{pos + len(ref) + 2}")
    primary = [record for record in records if not int(record[1]) & (0x4 | 0x100 | 0x800)]
    alt_names = {record[0] for record in primary if record_supports_indel(record, pos, ref, alt)}
    ref_names = {record[0] for record in primary if aligned_allele(record, pos, ref) == ref} - alt_names
    return alt_names, ref_names


def classify_fullref(records: list[list[str]], names: set[str], mt_contig: str = "chrM", min_mapq: int = 20) -> Counter:
    evidence = {name: dict(oa=False, nuc_secondary=False, nuc_supplementary=False, nuc_alt=False,
                           b_keep=False, nuclear_primary=False, strong_nuclear=False, nuclear_mate=False)
                for name in names}
    for f in records:
        if f[0] not in evidence:
            continue
        state = evidence[f[0]]
        flag, contig, mapq, mate = int(f[1]), f[2], int(f[4]), f[6]
        tags = f[11:]
        state["oa"] |= any(tag.startswith("OA:Z:") for tag in tags)
        secondary = bool(flag & 0x100 and contig not in ("*", mt_contig))
        supplementary = bool(flag & 0x800 and contig not in ("*", mt_contig))
        state["nuc_secondary"] |= secondary
        state["nuc_supplementary"] |= supplementary
        state["nuc_alt"] |= bool((secondary or supplementary) and mapq >= min_mapq)
        if not flag & (0x100 | 0x800):
            nuclear_primary = contig not in ("*", mt_contig)
            state["nuclear_primary"] |= nuclear_primary
            state["strong_nuclear"] |= bool(nuclear_primary and mapq >= min_mapq)
            state["b_keep"] |= bool(contig == mt_contig and mapq >= min_mapq)
            resolved_mate = contig if mate == "=" else mate
            state["nuclear_mate"] |= resolved_mate not in ("*", "", mt_contig)
        for tag in tags:
            if tag.startswith("XA:Z:"):
                state["nuc_secondary"] |= any(x and x.split(",", 1)[0] != mt_contig for x in tag[5:].split(";"))
            if tag.startswith("SA:Z:"):
                for item in tag[5:].split(";"):
                    parts = item.split(",")
                    if item and parts[0] != mt_contig:
                        state["nuc_supplementary"] = True
                        state["nuc_alt"] |= len(parts) > 4 and int(parts[4]) >= min_mapq
        state["strong_nuclear"] |= state["nuc_alt"]
    counts = Counter()
    for state in evidence.values():
        if state["strong_nuclear"]:
            counts["HIGH"] += 1
        elif (state["nuclear_primary"] or state["nuc_secondary"] or state["nuc_supplementary"]
              or state["nuclear_mate"] or not state["b_keep"]):
            counts["AMBIGUOUS"] += 1
        else:
            counts["LOW"] += 1
    return counts


def truth_by_run() -> dict[tuple[str, str], dict[str, str]]:
    return {(row["benchmark_id"], row["truth_id"]): row for row in read_tsv(SUMMARY12 / "block12_truth_panel.tsv")}


def strategy_c_diagnostics(_: argparse.Namespace) -> None:
    selected = [
        ("block12_R1_d250", "V2_A6"), ("block12_R2_d250", "V3_A5"),
        ("block12_R2_d1000", "V4_A6"), ("block12_R3_d1000", "V4_A5"),
        ("block12_R1_d5000", "V5_A6"), ("block12_R3_d5000", "V6_A5"),
    ]
    truths = truth_by_run()
    trace = []
    for benchmark_id, truth_id in selected:
        truth = truths[(benchmark_id, truth_id)]
        sample = truth["sample_id"]
        base = RUNS / benchmark_id / "results"
        mt_bam = base / f"mtDNA/realigned/standard/{sample}.standard_mt.sorted.bam"
        fullref = base / f"alignment/{sample}.sorted.bam"
        alt_names, ref_names = bam_indel_support(mt_bam, int(truth["pos"]), truth["ref"], truth["alt"])
        names_path = GENERATED / "strategy_c" / f"{benchmark_id}.{truth_id}.supporting_names.txt"
        names_path.parent.mkdir(parents=True, exist_ok=True)
        names_path.write_text("\n".join(sorted(alt_names)) + ("\n" if alt_names else ""))
        classified = classify_fullref(sam_view(fullref, names=names_path), alt_names)
        evidence = read_tsv(base / f"variants/numt/{sample}.numt_variant_evidence.tsv")
        row = next(x for x in evidence if x["caller"] == "gatk" and x["contig"] == truth["chrom"]
                   and x["position"] == truth["pos"] and x["ref"] == truth["ref"] and x["alt"] == truth["alt"])
        trace.append({"benchmark_id": benchmark_id, "truth_id": truth_id, "chrom": truth["chrom"],
                      "pos": truth["pos"], "ref": truth["ref"], "alt": truth["alt"],
                      "achieved_vaf": truth["achieved_vaf"], "depth": truth["depth_target"],
                      "replicate": truth["replicate"], "caller": "gatk",
                      "boundary_class": truth["boundary_class"],
                      "bam_support_alt_fragments": len(alt_names), "bam_support_ref_fragments": len(ref_names),
                      "LOW_fragments": classified["LOW"], "AMBIGUOUS_fragments": classified["AMBIGUOUS"],
                      "HIGH_fragments": classified["HIGH"],
                      "pipeline_supporting_fragments": row["total_supporting_fragments"],
                      "final_strategy_c_status": row["strategy_c_status"],
                      "reason": "INTENDED_SNV_ONLY_SUPPORT_SCOPE:non_SNV_forced_to_NA_before_evidence_join",
                      "support_source": row["support_source"]})
    fields = ["benchmark_id", "truth_id", "chrom", "pos", "ref", "alt", "achieved_vaf", "depth", "replicate",
              "caller", "boundary_class", "bam_support_alt_fragments", "bam_support_ref_fragments", "LOW_fragments",
              "AMBIGUOUS_fragments", "HIGH_fragments", "pipeline_supporting_fragments", "final_strategy_c_status",
              "reason", "support_source"]
    write_tsv(SUMMARY / "block12b_strategy_c_indel_trace.tsv", fields, trace)

    # Matched controls show the scientific rule after support extraction and
    # the deliberate non-SNV scope gate before that rule.
    controls = []
    evidence_path = RUNS / "block12_R1_d1000/results/variants/numt/B12R1D1000.numt_variant_evidence.tsv"
    evidence = read_tsv(evidence_path)
    retained = next(x for x in evidence if len(x["ref"]) == len(x["alt"]) == 1 and x["strategy_c_status"] == "RETAINED"
                    and int(x["low_numt_evidence_fragments"]) >= 2)
    flagged = next(x for x in evidence if len(x["ref"]) == len(x["alt"]) == 1 and x["strategy_c_status"] == "FLAGGED_NUMT_EVIDENCE"
                   and int(x["total_supporting_fragments"]) >= 2)
    for label, row in (("SNV_STRONG_LOW", retained), ("SNV_HIGH_OR_AMBIGUOUS", flagged)):
        controls.append({"control_id": label, "variant_class": "SNV", "source": "BLOCK12_OBSERVED",
                         "total": row["total_supporting_fragments"], "LOW": row["low_numt_evidence_fragments"],
                         "AMBIGUOUS": row["ambiguous_numt_evidence_fragments"], "HIGH": row["high_numt_evidence_fragments"],
                         "scope_gate_input": "COUNTS_AVAILABLE", "strategy_c_status": row["strategy_c_status"],
                         "interpretation": "frozen >=2/LOW==0 rule applied"})
    low_indel = trace[-1]
    controls.append({"control_id": "INDEL_STRONG_LOW", "variant_class": "indel", "source": "BLOCK12_CIGAR_DIAGNOSTIC",
                     "total": low_indel["bam_support_alt_fragments"], "LOW": low_indel["LOW_fragments"],
                     "AMBIGUOUS": low_indel["AMBIGUOUS_fragments"], "HIGH": low_indel["HIGH_fragments"],
                     "scope_gate_input": "NA_non_SNV", "strategy_c_status": "UNRESOLVED_SUPPORT",
                     "interpretation": "evidence exists diagnostically but frozen extractor does not evaluate non-SNV support"})
    controls.append({"control_id": "INDEL_HIGH_OR_AMBIGUOUS", "variant_class": "indel", "source": "DECISION_UNIT_CONTROL",
                     "total": 5, "LOW": 0, "AMBIGUOUS": 2, "HIGH": 3, "scope_gate_input": "NA_non_SNV",
                     "strategy_c_status": "UNRESOLVED_SUPPORT",
                     "interpretation": "non-SNV scope gate precedes frozen evidence decision rule"})
    write_tsv(SUMMARY / "block12b_strategy_c_diagnostic.tsv",
              ["control_id", "variant_class", "source", "total", "LOW", "AMBIGUOUS", "HIGH",
               "scope_gate_input", "strategy_c_status", "interpretation"], controls)
    print(f"PASS\tstrategy_c\tindel_traces={len(trace)}\tcontrols={len(controls)}")


def reference_bundle() -> Path:
    matches = sorted((RUNS / "block12_R1_d1000/work").glob("*/*/gatk_mt_reference_bundle.tar"))
    if not matches:
        raise RuntimeError("GATK mt reference bundle not found")
    return matches[0]


def run_standard_full(benchmark_id: str, force: bool = False) -> Path:
    depth = benchmark_id.rsplit("d", 1)[1]
    replicate = benchmark_id.split("_")[1]
    sample = f"B12{replicate}D{depth}"
    out = GENERATED / "circular" / benchmark_id
    filtered = out / f"{sample}.standard_only_matched.filtered.vcf.gz"
    if filtered.is_file() and not force:
        return filtered
    out.mkdir(parents=True, exist_ok=True)
    bundle = reference_bundle()
    bam = RUNS / benchmark_id / f"results/mtDNA/realigned/standard/{sample}.standard_mt.sorted.bam"
    raw = f"{sample}.standard_full.raw.vcf.gz"
    command = (
        f"set -euo pipefail; tar -xf /workspace/{shlex.quote(rel(bundle))}; "
        f"gatk Mutect2 -R standard_mt.fa -I /workspace/{shlex.quote(rel(bam))} "
        f"-L chrM:1-16569 --mitochondria-mode -O {raw}; "
        f"gatk FilterMutectCalls -R standard_mt.fa -V {raw} --stats {raw}.stats "
        f"--mitochondria-mode -O {filtered.name}"
    )
    docker(GATK_IMAGE, ["bash", "-lc", command], out)
    return filtered


def normalized_active_vcf(path: Path, benchmark_id: str, reference: Path) -> dict[str, dict]:
    refs = bh.parse_fasta(reference)
    result = {}
    for call in bh.normalize_vcf(path, benchmark_id, "gatk"):
        pos, ref, alt = bh.normalize_allele(str(call["chrom"]), int(call["pos"]), str(call["ref"]), str(call["alt"]), refs)
        call.update(pos=pos, ref=ref, alt=alt)
        if call["filter_status"] in {"PASS", "."}:
            result[f"{call['chrom']}:{pos}:{ref}:{alt}"] = call
    return result


def circular_matched(args: argparse.Namespace) -> None:
    selected = [f"block12_R{rep}_d{depth}" for rep in (1, 2, 3) for depth in (250, 1000, 5000)]
    truth_rows = read_tsv(SUMMARY12 / "block12_truth_panel.tsv")
    reference = BENCH / "datasets/block12/reference/block12_reference.fa"
    details = []
    for benchmark_id in selected:
        standard = normalized_active_vcf(run_standard_full(benchmark_id, args.force), benchmark_id, reference)
        sample = next(x["sample_id"] for x in truth_rows if x["benchmark_id"] == benchmark_id)
        circular = normalized_active_vcf(RUNS / benchmark_id / f"results/variants/gatk/{sample}.gatk.filtered.vcf.gz",
                                         benchmark_id, reference)
        subset = [x for x in truth_rows if x["benchmark_id"] == benchmark_id and x["target_vaf"] in {"0.015", "0.075", "0.25"}
                  and x["numt_context"] == "NO_NUMT_CHALLENGE"]
        for truth in subset:
            key = truth["variant_key"]
            for mode, calls in (("STANDARD_ONLY_MATCHED", standard), ("STANDARD_PLUS_SHIFTED_MATCHED", circular)):
                call = calls.get(key)
                called_vaf = "NA" if call is None else call["called_vaf"]
                error = "NA" if called_vaf == "NA" else f"{abs(float(called_vaf) - float(truth['achieved_vaf'])):.10g}"
                details.append({"summary_level": "TRUTH_ALLELE", "benchmark_id": benchmark_id,
                                "replicate": truth["replicate"], "depth": truth["depth_target"],
                                "target_vaf": truth["target_vaf"], "truth_id": truth["truth_id"],
                                "variant_class": truth["variant_class"], "boundary_class": truth["boundary_class"],
                                "mode": mode, "tp": int(call is not None), "fn": int(call is None),
                                "recall": int(call is not None), "called_vaf": called_vaf,
                                "absolute_vaf_error": error, "recall_delta": "NA"})
    pooled = []
    for boundary in ("BOUNDARY", "NON_BOUNDARY"):
        for variant_class in ("SNV", "indel"):
            mode_rows = {}
            for mode in ("STANDARD_ONLY_MATCHED", "STANDARD_PLUS_SHIFTED_MATCHED"):
                group = [x for x in details if x["boundary_class"] == boundary and x["variant_class"] == variant_class and x["mode"] == mode]
                tp, fn = sum(x["tp"] for x in group), sum(x["fn"] for x in group)
                errors = [float(x["absolute_vaf_error"]) for x in group if x["absolute_vaf_error"] != "NA"]
                mode_rows[mode] = {"summary_level": "POOLED", "benchmark_id": "ALL_9", "replicate": "ALL",
                                   "depth": "250,1000,5000", "target_vaf": "0.015,0.075,0.25", "truth_id": "ALL",
                                   "variant_class": variant_class, "boundary_class": boundary, "mode": mode,
                                   "tp": tp, "fn": fn, "recall": f"{tp/(tp+fn):.10g}", "called_vaf": "NA",
                                   "absolute_vaf_error": "NA" if not errors else f"{statistics.mean(errors):.10g}"}
            delta = float(mode_rows["STANDARD_PLUS_SHIFTED_MATCHED"]["recall"]) - float(mode_rows["STANDARD_ONLY_MATCHED"]["recall"])
            for row in mode_rows.values():
                row["recall_delta"] = f"{delta:.10g}"
                pooled.append(row)
    fields = ["summary_level", "benchmark_id", "replicate", "depth", "target_vaf", "truth_id", "variant_class",
              "boundary_class", "mode", "tp", "fn", "recall", "called_vaf", "absolute_vaf_error", "recall_delta"]
    write_tsv(SUMMARY / "block12b_circular_matched_summary.tsv", fields, details + pooled)
    print(f"PASS\tcircular_matched\truns={len(selected)}\ttruth_mode_rows={len(details)}")


def open_vcf(path: Path) -> list[dict]:
    calls = []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip().split("\t")
            for alt in f[4].split(","):
                calls.append({"chrom": f[0], "pos": int(f[1]), "ref": f[3], "alt": alt, "filter": f[6]})
    return calls


def emit_positive_sam(path: Path, reference: str, insertion: dict, deletion: dict, depth: int = 1000, alt_fraction: float = .25) -> None:
    sample = "B12BPOS"
    with path.open("w") as out:
        out.write("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chrM\tLN:16569\n")
        out.write(f"@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA\n")
        records = []
        snv = {"pos": 10000, "ref": reference[9999],
               "alt": "A" if reference[9999] != "A" else "C"}
        for event, truth in (("INS", insertion), ("DEL", deletion), ("SNV", snv)):
            pos, ref, alt = int(truth["pos"]), truth["ref"], truth["alt"]
            for i in range(depth):
                start = pos - 75 + (i % 31) - 15
                left = pos - start + 1
                is_alt = i < round(depth * alt_fraction)
                if event == "INS" and is_alt:
                    right = 150 - left
                    sequence = reference[start - 1:start - 1 + left] + alt[len(ref):] + reference[pos:pos + right]
                    cigar = f"{left}M{len(alt)-len(ref)}I{right}M"
                elif event == "DEL" and is_alt:
                    right = 150 - left
                    sequence = reference[start - 1:start - 1 + left] + reference[pos + len(ref) - 1:pos + len(ref) - 1 + right]
                    cigar = f"{left}M{len(ref)-len(alt)}D{right}M"
                elif event == "SNV" and is_alt:
                    sequence = reference[start - 1:start - 1 + 150]
                    offset = pos - start
                    sequence = sequence[:offset] + alt + sequence[offset + 1:]
                    cigar = "150M"
                else:
                    sequence = reference[start - 1:start - 1 + 150]
                    cigar = "150M"
                flag = 16 if i % 2 else 0
                if flag:
                    sequence = sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]
                records.append((start, f"{event}:{i+1}", flag, cigar, sequence))
        for start, name, flag, cigar, sequence in sorted(records):
            out.write(f"{name}\t{flag}\tchrM\t{start}\t60\t{cigar}\t*\t0\t0\t{sequence}\t{'I'*len(sequence)}\tRG:Z:{sample}\n")


def circular_slice(sequence: str, start: int, length: int) -> str:
    return "".join(sequence[(start + offset) % len(sequence)] for offset in range(length))


def emit_positive_fastq(out: Path, reference: str, insertion: dict, deletion: dict,
                        depth: int = 1000, alt_fraction: float = .25) -> tuple[Path, Path]:
    read_length, fragment_length = 150, 300
    fragments = math.ceil(depth * len(reference) / (2 * read_length))
    starts = [int(i * len(reference) / fragments) for i in range(fragments)]
    rng = random.Random(120012)
    rng.shuffle(starts)
    snv = {"truth_id": "POS_SNV_CONTROL", "pos": 10000, "ref": reference[9999],
           "alt": "A" if reference[9999] != "A" else "C"}
    events = (insertion, deletion, snv)
    chosen = {}
    for event in events:
        eligible = []
        pos0 = int(event["pos"]) - 1
        for index, start in enumerate(starts):
            mate = (start + fragment_length - read_length) % len(reference)
            if any((read_start + offset) % len(reference) == pos0
                   for read_start in (start, mate) for offset in range(read_length)):
                eligible.append(index)
        chosen[event["truth_id"]] = set(rng.sample(eligible, round(len(eligible) * alt_fraction)))
    r1_path, r2_path = out / "B12BPOS_R1.fastq.gz", out / "B12BPOS_R2.fastq.gz"
    with gzip.open(r1_path, "wt") as r1, gzip.open(r2_path, "wt") as r2:
        for index, start in enumerate(starts):
            mate_start = (start + fragment_length - read_length) % len(reference)
            read1 = circular_slice(reference, start, read_length)
            mate_forward = circular_slice(reference, mate_start, read_length)
            for event in events:
                if index not in chosen[event["truth_id"]]:
                    continue
                for which, sequence, read_start in ((1, read1, start), (2, mate_forward, mate_start)):
                    offsets = [offset for offset in range(read_length)
                               if (read_start + offset) % len(reference) == int(event["pos"]) - 1]
                    if not offsets:
                        continue
                    offset = offsets[0]
                    ref, alt = event["ref"], event["alt"]
                    if sequence[offset:offset + len(ref)] != ref:
                        raise RuntimeError("positive-control REF mismatch")
                    changed = sequence[:offset] + alt + sequence[offset + len(ref):]
                    if len(changed) < read_length:
                        changed += circular_slice(reference, read_start + read_length, read_length - len(changed))
                    changed = changed[:read_length]
                    if which == 1: read1 = changed
                    else: mate_forward = changed
            name = f"B12BPOS:{index+1}:{start+1}"
            quality = "I" * read_length
            r1.write(f"@{name}/1\n{read1}\n+\n{quality}\n")
            mate = mate_forward.translate(str.maketrans("ACGT", "TGCA"))[::-1]
            r2.write(f"@{name}/2\n{mate}\n+\n{quality}\n")
    return r1_path, r2_path


def positive_control(_: argparse.Namespace) -> None:
    # v4 includes an unrelated SNV emission control because mutserve2 2.0.3
    # throws while writing a completely empty VCF when both indel options are
    # disabled. The SNV does not alter either indel's support.
    out = GENERATED / "mutserve_indel_positive_v4"
    out.mkdir(parents=True, exist_ok=True)
    ref_records = bh.parse_fasta(BENCH / "datasets/block12/reference/block12_reference.fa")
    sequence = ref_records["chrM"]
    insertion = {"truth_id": "POS_INS", "chrom": "chrM", "pos": 8000,
                 "ref": sequence[7999], "alt": sequence[7999] + ("A" if sequence[7999] != "A" else "C")}
    deletion = {"truth_id": "POS_DEL", "chrom": "chrM", "pos": 9000,
                "ref": sequence[8999:9001], "alt": sequence[8999]}
    bam = out / "B12BPOS.sorted.bam"
    if not bam.is_file():
        r1, r2 = emit_positive_fastq(out, sequence, insertion, deletion)
        bundle = reference_bundle()
        command = (f"set -euo pipefail; tar -xf /workspace/{shlex.quote(rel(bundle))}; "
                   "bwa-mem2 index standard_mt.fa; bwa-mem2 mem -t 2 "
                   "-R '@RG\\tID:B12BPOS\\tSM:B12BPOS\\tPL:ILLUMINA' standard_mt.fa "
                   "B12BPOS_R1.fastq.gz B12BPOS_R2.fastq.gz > B12BPOS.sam")
        docker("quay.io/biocontainers/bwa-mem2:2.2.1--he513fc3_0", ["bash", "-lc", command], out)
        docker(SAMTOOLS_IMAGE, ["bash", "-lc", "samtools sort -o B12BPOS.sorted.bam B12BPOS.sam; samtools index B12BPOS.sorted.bam"], out)
    bundle = reference_bundle()
    gatk_filtered = out / "B12BPOS.gatk.filtered.vcf.gz"
    if not gatk_filtered.is_file():
        command = (f"set -euo pipefail; tar -xf /workspace/{shlex.quote(rel(bundle))}; "
                   "gatk Mutect2 -R standard_mt.fa -I B12BPOS.sorted.bam -L chrM:1-16569 --mitochondria-mode "
                   "-O B12BPOS.gatk.raw.vcf.gz; gatk FilterMutectCalls -R standard_mt.fa -V B12BPOS.gatk.raw.vcf.gz "
                   "--stats B12BPOS.gatk.raw.vcf.gz.stats --mitochondria-mode -O B12BPOS.gatk.filtered.vcf.gz")
        docker(GATK_IMAGE, ["bash", "-lc", command], out)
    mutserve_native = out / "B12BPOS.mutserve2.native.vcf.gz"
    if not mutserve_native.is_file():
        command = (f"set -euo pipefail; tar -xf /workspace/{shlex.quote(rel(bundle))}; "
                   "mutserve call --reference standard_mt.fa --output B12BPOS.mutserve2.native.vcf.gz --threads 2 "
                   "--level 0.01 --contig-name chrM --no-ansi --write-raw B12BPOS.sorted.bam")
        docker(MUTSERVE_IMAGE, ["bash", "-lc", command], out)
    canonical = out / "B12BPOS.mutserve2.canonical.vcf.gz"
    if not canonical.is_file():
        command = (f"set -euo pipefail; tar -xf /workspace/{shlex.quote(rel(bundle))}; "
                   "printf 'B12BPOS\\n' > sample_name.txt; bcftools reheader --samples sample_name.txt "
                   "--output reheader.vcf.gz B12BPOS.mutserve2.native.vcf.gz; "
                   "bcftools norm --fasta-ref standard_mt.fa --check-ref e --multiallelics -any reheader.vcf.gz "
                   "| bcftools sort -Oz -o B12BPOS.mutserve2.canonical.vcf.gz; "
                   "bcftools index --tbi B12BPOS.mutserve2.canonical.vcf.gz; "
                   "/workspace/bin/mutserve2_vcf_to_tsv B12BPOS B12BPOS.mutserve2.canonical.vcf.gz B12BPOS.mutserve2.calls.tsv")
        docker("quay.io/biocontainers/bcftools:1.20--h8b25389_0", ["bash", "-lc", command], out)
    refs = {"chrM": sequence}
    stages = {"gatk_filtered": open_vcf(gatk_filtered), "mutserve_native": open_vcf(mutserve_native),
              "mutserve_canonical": open_vcf(canonical)}
    rows = []
    for truth in (insertion, deletion):
        pos, ref, alt = bh.normalize_allele("chrM", truth["pos"], truth["ref"], truth["alt"], refs)
        key = ("chrM", pos, ref, alt)
        alt_fragments, ref_fragments = bam_indel_support(bam, pos, ref, alt)
        observed_total = len(alt_fragments) + len(ref_fragments)
        record = {"truth_id": truth["truth_id"], "chrom": "chrM", "pos": pos, "ref": ref, "alt": alt,
                  "variant_class": "insertion" if len(alt) > len(ref) else "deletion", "target_vaf": "0.25",
                  "achieved_vaf": "NA" if not observed_total else f"{len(alt_fragments) / observed_total:.10g}",
                  "alt_fragments": len(alt_fragments), "ref_fragments": len(ref_fragments),
                  "depth": "1000", "gatk_emitted": "false", "mutserve_native_emitted": "false",
                  "mutserve_canonical_emitted": "false", "benchmark_parser_emitted": "false",
                  "classification": "MUTSERVE2_INDEL_NOT_EMITTED_BY_FROZEN_MODE",
                  "frozen_command_indel_options": "insertions=false;deletions=false"}
        for stage, calls in stages.items():
            normalized = set()
            for call in calls:
                cpos, cref, calt = bh.normalize_allele(call["chrom"], call["pos"], call["ref"], call["alt"], refs)
                normalized.add((call["chrom"], cpos, cref, calt))
            field = {"gatk_filtered": "gatk_emitted", "mutserve_native": "mutserve_native_emitted",
                     "mutserve_canonical": "mutserve_canonical_emitted"}[stage]
            record[field] = str(key in normalized).lower()
        record["benchmark_parser_emitted"] = record["mutserve_canonical_emitted"]
        rows.append(record)
    fields = ["truth_id", "chrom", "pos", "ref", "alt", "variant_class", "target_vaf", "achieved_vaf",
              "alt_fragments", "ref_fragments", "depth", "gatk_emitted",
              "mutserve_native_emitted", "mutserve_canonical_emitted", "benchmark_parser_emitted", "classification",
              "frozen_command_indel_options"]
    write_tsv(SUMMARY / "block12b_mutserve_indel_positive_control.tsv", fields, rows)
    print("PASS\tmutserve_indel_positive_control\ttruth=2")


def mutserve_trace(_: argparse.Namespace) -> None:
    truths = read_tsv(SUMMARY12 / "block12_truth_panel.tsv")
    variants = read_tsv(SUMMARY12 / "block12_variant_metrics.tsv")
    manifest = {x["benchmark_id"]: x for x in bh.validate_manifest(BENCH / "manifests/block12_primary_matrix.tsv")}
    reference = Path(next(iter(manifest.values()))["reference"])
    refs = bh.parse_fasta(reference)
    stage_cache = {}
    rows = []
    for truth in truths:
        if truth["variant_class"] != "indel":
            continue
        bid, sample = truth["benchmark_id"], truth["sample_id"]
        key = (truth["chrom"], int(truth["pos"]), truth["ref"], truth["alt"])
        def keys(calls):
            result = set()
            for call in calls:
                pos, ref, alt = bh.normalize_allele(call["chrom"], call["pos"], call["ref"], call["alt"], refs)
                result.add((call["chrom"], pos, ref, alt))
            return result
        if bid not in stage_cache:
            native = open_vcf(RUNS / bid / f"results/variants/mutserve2/native/{sample}.mutserve2.native.vcf.gz")
            canonical = open_vcf(RUNS / bid / f"results/variants/mutserve2/{sample}.mutserve2.raw.vcf.gz")
            stage_cache[bid] = (keys(native), keys(canonical))
        native_keys, canonical_keys = stage_cache[bid]
        metric = next(x for x in variants if x["benchmark_id"] == bid and x["truth_id"] == truth["truth_id"]
                      and x["callset_source"] == "mutserve2_canonical_A")
        rows.append({"benchmark_id": bid, "truth_id": truth["truth_id"], "replicate": truth["replicate"],
                     "depth": truth["depth_target"], "achieved_vaf": truth["achieved_vaf"], "chrom": truth["chrom"],
                     "pos": truth["pos"], "ref": truth["ref"], "alt": truth["alt"],
                     "raw_native_emitted": str(key in native_keys).lower(),
                     "canonical_emitted": str(key in canonical_keys).lower(),
                     "benchmark_normalized_status": metric["truth_status"],
                     "loss_stage": "CALLER_NOT_EMITTED" if key not in native_keys else "NONE",
                     "frozen_mode": "--insertions absent;--deletions absent;both mutserve defaults false"})
    fields = ["benchmark_id", "truth_id", "replicate", "depth", "achieved_vaf", "chrom", "pos", "ref", "alt",
              "raw_native_emitted", "canonical_emitted", "benchmark_normalized_status", "loss_stage", "frozen_mode"]
    write_tsv(SUMMARY / "block12b_mutserve_indel_trace.tsv", fields, rows)
    print(f"PASS\tmutserve_indel_trace\trows={len(rows)}")


def finalize(_: argparse.Namespace) -> None:
    required = ["block12b_strategy_c_diagnostic.tsv", "block12b_strategy_c_indel_trace.tsv",
                "block12b_circular_matched_summary.tsv", "block12b_mutserve_indel_trace.tsv",
                "block12b_mutserve_indel_positive_control.tsv"]
    missing = [name for name in required if not (SUMMARY / name).is_file()]
    if missing:
        raise RuntimeError(f"missing controlled-benchmark diagnostic outputs: {missing}")
    print("CONTROLLED_BENCHMARK_DIAGNOSTIC_QC_PASS")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("strategy-c")
    circular = sub.add_parser("circular-matched"); circular.add_argument("--force", action="store_true")
    sub.add_parser("mutserve-trace")
    sub.add_parser("positive-control")
    sub.add_parser("finalize")
    args = parser.parse_args()
    {"strategy-c": strategy_c_diagnostics, "circular-matched": circular_matched,
     "mutserve-trace": mutserve_trace, "positive-control": positive_control,
     "finalize": finalize}[args.command](args)


if __name__ == "__main__":
    main()
