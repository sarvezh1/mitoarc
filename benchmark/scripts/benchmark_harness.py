#!/usr/bin/env python3
"""Controlled benchmark harness for the frozen MitoArc workflow.

Only Python's standard library is used so contracts can be validated before
Docker or Nextflow is invoked.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import io
import math
import os
import random
import re
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path


MANIFEST_FIELDS = [
    "benchmark_id", "sample_id", "dataset_type", "reference",
    "input_bam_or_fastq", "truth_vcf", "depth_target", "vaf_target",
    "vaf_bin", "variant_class", "variant_position", "boundary_class",
    "numt_context", "caller", "numt_strategy", "circular_mode",
    "replicate", "seed", "numt_truth", "numt_design",
]
TRUTH_FIELDS = [
    "truth_id", "chrom", "pos", "ref", "alt", "variant_class",
    "expected_vaf", "boundary_class", "numt_context", "source", "notes",
]
CALL_FIELDS = [
    "benchmark_id", "caller", "chrom", "pos", "ref", "alt",
    "variant_key", "filter_status", "called_vaf", "called_depth",
    "truth_status",
]
RESOURCE_FIELDS = [
    "benchmark_id", "attempt", "wall_clock_seconds", "cpu_seconds",
    "cpu_hours", "peak_rss_mb", "max_container_memory_mb",
    "disk_input_bytes", "disk_work_bytes_peak", "disk_output_bytes",
    "process_count", "cached_process_count", "submitted_process_count",
    "nextflow_status",
]
DEPTHS = {100, 250, 500, 1000, 2000, 5000}
DATASET_TYPES = {"synthetic", "spike_in", "downsampled", "real_validation"}
VAF_BINS = {
    "LT_0_01", "VAF_0_01_0_02", "VAF_0_02_0_05",
    "VAF_0_05_0_10", "VAF_0_10_0_50", "GT_0_50", "MULTIPLEXED",
}
CALLERS = {"gatk", "mutserve2", "all_supported"}
NUMT_STRATEGIES = {"A", "B", "C", "comparison"}
CIRCULAR_MODES = {"standard_only", "standard_plus_shifted", "all_supported"}
VARIANT_CLASSES = {"SNV", "indel", "MIXED"}
BOUNDARY_CLASSES = {"BOUNDARY", "NON_BOUNDARY", "MIXED"}
NUMT_CONTEXTS = {"NO_NUMT_CHALLENGE", "NUMT_LIKE", "AMBIGUOUS_CONTEXT", "MIXED"}
ID_RE = re.compile(r"^[A-Za-z0-9_.-]+$")
DNA_RE = re.compile(r"^[ACGTN]+$")


class HarnessError(RuntimeError):
    pass


def fail(message: str) -> None:
    raise HarnessError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def directory_bytes(path: Path) -> int:
    if not path.exists():
        return 0
    return sum(p.stat().st_size for p in path.rglob("*") if p.is_file())


def disk_usage_bytes(path: Path) -> int:
    if not path.exists():
        return 0
    try:
        result = subprocess.run(
            ["du", "-sb", str(path)], text=True, stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL, check=True,
        )
        return int(result.stdout.split()[0])
    except (FileNotFoundError, subprocess.CalledProcessError, ValueError, IndexError):
        return directory_bytes(path)


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.exists():
        fail(f"file does not exist: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            fail(f"empty TSV: {path}")
        return list(reader.fieldnames), list(reader)


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "NA") for field in fields})


def parse_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    name = None
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                if not name or name in records:
                    fail(f"invalid or duplicate FASTA record at {path}:{line_number}")
                records[name] = []
            elif name is None:
                fail(f"sequence before FASTA header at {path}:{line_number}")
            else:
                records[name].append(line.upper())
    if not records:
        fail(f"no FASTA records: {path}")
    return {name: "".join(parts) for name, parts in records.items()}


def vaf_bin(value: float) -> str:
    if 0 <= value < 0.01:
        return "LT_0_01"
    if 0.01 <= value < 0.02:
        return "VAF_0_01_0_02"
    if 0.02 <= value < 0.05:
        return "VAF_0_02_0_05"
    if 0.05 <= value < 0.10:
        return "VAF_0_05_0_10"
    if 0.10 <= value <= 0.50:
        return "VAF_0_10_0_50"
    if 0.50 < value <= 1:
        return "GT_0_50"
    fail(f"VAF outside 0..1: {value}")


def resolve_path(value: str, base: Path) -> Path:
    path = Path(value)
    return (base / path).resolve() if not path.is_absolute() else path.resolve()


def generated_provenance(inputs: list[Path]) -> Path | None:
    if len(inputs) == 1 and inputs[0].suffix.lower() == ".bam":
        candidate = Path(str(inputs[0]) + ".provenance.tsv")
        return candidate if candidate.is_file() else None
    if len(inputs) == 2:
        names = [path.name for path in inputs]
        match1 = re.fullmatch(r"(.+)_R1\.f(?:ast)?q\.gz", names[0], re.I)
        match2 = re.fullmatch(r"(.+)_R2\.f(?:ast)?q\.gz", names[1], re.I)
        if match1 and match2 and match1.group(1) == match2.group(1):
            candidate = inputs[0].with_name(match1.group(1) + ".provenance.tsv")
            return candidate if candidate.is_file() else None
    return None


def validate_manifest(path: Path, require_files: bool = True) -> list[dict[str, str]]:
    fields, rows = read_tsv(path)
    missing = [field for field in MANIFEST_FIELDS if field not in fields]
    extra = [field for field in fields if field not in MANIFEST_FIELDS]
    if missing or extra:
        fail(f"manifest schema mismatch; missing={missing or 'none'} extra={extra or 'none'}")
    if not rows:
        fail("manifest contains no experiments")
    seen: set[str] = set()
    replicate_seeds: dict[tuple[str, ...], dict[int, int]] = {}
    base = path.parent
    for line_number, row in enumerate(rows, 2):
        prefix = f"manifest line {line_number}"
        if any(row.get(field, "") == "" for field in MANIFEST_FIELDS):
            fail(f"{prefix}: empty fields are not allowed; use NA where applicable")
        benchmark_id = row["benchmark_id"]
        if not ID_RE.fullmatch(benchmark_id):
            fail(f"{prefix}: invalid benchmark_id {benchmark_id!r}")
        if benchmark_id in seen:
            fail(f"{prefix}: duplicate benchmark_id {benchmark_id}")
        seen.add(benchmark_id)
        for field, allowed in [
            ("dataset_type", DATASET_TYPES), ("vaf_bin", VAF_BINS),
            ("variant_class", VARIANT_CLASSES), ("boundary_class", BOUNDARY_CLASSES),
            ("numt_context", NUMT_CONTEXTS), ("caller", CALLERS),
            ("numt_strategy", NUMT_STRATEGIES), ("circular_mode", CIRCULAR_MODES),
        ]:
            if row[field] not in allowed:
                fail(f"{prefix}: unsupported {field} {row[field]!r}; allowed={sorted(allowed)}")
        try:
            depth = int(row["depth_target"])
            target_vaf = None if row["vaf_target"] == "NA" else float(row["vaf_target"])
            replicate = int(row["replicate"])
            seed = int(row["seed"])
        except ValueError as error:
            fail(f"{prefix}: invalid numeric field: {error}")
        if depth not in DEPTHS:
            fail(f"{prefix}: unsupported depth_target {depth}")
        if target_vaf is None:
            if row["vaf_bin"] != "MULTIPLEXED":
                fail(f"{prefix}: NA vaf_target requires MULTIPLEXED vaf_bin")
        elif not 0 <= target_vaf <= 1 or vaf_bin(target_vaf) != row["vaf_bin"]:
            fail(f"{prefix}: vaf_target {target_vaf} does not match {row['vaf_bin']}")
        if replicate < 1 or not 0 <= seed <= 2147483647:
            fail(f"{prefix}: replicate must be >=1 and seed must be 0..2147483647")
        if row["variant_position"] != "NA":
            try:
                if int(row["variant_position"]) < 1:
                    raise ValueError
            except ValueError:
                fail(f"{prefix}: invalid variant_position {row['variant_position']!r}")
        if row["caller"] == "mutserve2" and row["circular_mode"] != "standard_only":
            fail(f"{prefix}: frozen mutserve2 supports only standard_only input")
        if row["caller"] == "all_supported" and not (
            row["numt_strategy"] == "comparison" and row["circular_mode"] == "all_supported"
        ):
            fail(f"{prefix}: all_supported executions require comparison/all_supported modes")
        if require_files:
            for value in [row["reference"], row["truth_vcf"]]:
                if value == "NA" or not resolve_path(value, base).is_file():
                    fail(f"{prefix}: required file does not exist: {value}")
            resolved_inputs = []
            for value in row["input_bam_or_fastq"].split(";"):
                if value == "NA" or not resolve_path(value, base).is_file():
                    fail(f"{prefix}: input file does not exist: {value}")
                resolved_inputs.append(resolve_path(value, base))
            if row["dataset_type"] in {"synthetic", "spike_in", "downsampled"} and generated_provenance(resolved_inputs) is None:
                fail(f"{prefix}: generated input is missing its required provenance TSV")
            for field in ("numt_truth", "numt_design"):
                if row[field] != "NA" and not resolve_path(row[field], base).is_file():
                    fail(f"{prefix}: {field} file does not exist: {row[field]}")
        signature = tuple(row[field] for field in MANIFEST_FIELDS if field not in {"benchmark_id", "replicate", "seed", "caller"})
        prior = replicate_seeds.setdefault(signature, {})
        if replicate in prior:
            fail(f"{prefix}: duplicate replicate {replicate} for the same condition")
        if seed in prior.values():
            fail(f"{prefix}: seed {seed} is reused across nominal replicates")
        prior[replicate] = seed
    return rows


def trim_alleles(pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        pos += 1
        ref, alt = ref[1:], alt[1:]
    return pos, ref, alt


def normalize_allele(chrom: str, pos: int, ref: str, alt: str, references: dict[str, str]) -> tuple[int, str, str]:
    sequence = references[chrom]
    pos, ref, alt = trim_alleles(pos, ref, alt)
    if len(ref) != len(alt):
        while pos > 1:
            previous = sequence[pos - 2]
            candidate_pos, candidate_ref, candidate_alt = trim_alleles(
                pos - 1, previous + ref, previous + alt
            )
            if candidate_pos >= pos:
                break
            pos, ref, alt = candidate_pos, candidate_ref, candidate_alt
    return pos, ref, alt


def validate_truth(path: Path, reference: Path, output: Path | None = None) -> list[dict[str, object]]:
    fields, rows = read_tsv(path)
    missing = [field for field in TRUTH_FIELDS if field not in fields]
    extra = [field for field in fields if field not in TRUTH_FIELDS and field != "variant_key"]
    if missing or extra:
        fail(f"truth schema mismatch; missing={missing or 'none'} extra={extra or 'none'}")
    refs = parse_fasta(reference)
    seen_ids: set[str] = set()
    seen_keys: set[str] = set()
    normalized: list[dict[str, object]] = []
    for line_number, row in enumerate(rows, 2):
        prefix = f"truth line {line_number}"
        truth_id = row["truth_id"]
        if not truth_id or truth_id in seen_ids:
            fail(f"{prefix}: duplicate or empty truth_id {truth_id!r}")
        seen_ids.add(truth_id)
        chrom = row["chrom"]
        if chrom not in refs:
            fail(f"{prefix}: contig {chrom!r} is absent from reference")
        try:
            pos = int(row["pos"])
            expected_vaf = float(row["expected_vaf"])
        except ValueError as error:
            fail(f"{prefix}: invalid numeric value: {error}")
        ref, alt = row["ref"].upper(), row["alt"].upper()
        if not DNA_RE.fullmatch(ref) or not DNA_RE.fullmatch(alt) or ref == alt:
            fail(f"{prefix}: invalid REF/ALT {ref}>{alt}")
        if pos < 1 or pos + len(ref) - 1 > len(refs[chrom]):
            fail(f"{prefix}: invalid mitochondrial coordinate {chrom}:{pos}")
        observed_ref = refs[chrom][pos - 1:pos - 1 + len(ref)]
        if observed_ref != ref:
            fail(f"{prefix}: REF mismatch at {chrom}:{pos}; truth={ref} reference={observed_ref}")
        if not 0 <= expected_vaf <= 1:
            fail(f"{prefix}: impossible expected_vaf {expected_vaf}")
        if row["variant_class"] not in VARIANT_CLASSES:
            fail(f"{prefix}: invalid variant_class {row['variant_class']!r}")
        inferred = "SNV" if len(ref) == len(alt) == 1 else "indel"
        if inferred != row["variant_class"]:
            fail(f"{prefix}: {row['variant_class']} is inconsistent with {ref}>{alt}")
        if row["boundary_class"] not in BOUNDARY_CLASSES or row["numt_context"] not in NUMT_CONTEXTS:
            fail(f"{prefix}: invalid boundary_class or numt_context")
        pos, ref, alt = normalize_allele(chrom, pos, ref, alt, refs)
        key = f"{chrom}:{pos}:{ref}:{alt}"
        if key in seen_keys:
            fail(f"{prefix}: duplicate normalized truth allele {key}")
        seen_keys.add(key)
        result = dict(row)
        result.update(pos=pos, ref=ref, alt=alt, expected_vaf=f"{expected_vaf:.10g}", variant_key=key)
        normalized.append(result)
    if output:
        write_tsv(output, TRUTH_FIELDS + ["variant_key"], normalized)
    return normalized


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def circular_slice(sequence: str, start: int, length: int) -> str:
    return "".join(sequence[(start + index) % len(sequence)] for index in range(length))


def allele_offset(start: int, read_length: int, locus: int, ref_length: int, genome_length: int) -> int | None:
    for offset in range(read_length):
        if (start + offset) % genome_length == locus and offset + ref_length <= read_length:
            return offset
    return None


def apply_allele_to_read(read: str, start: int, row: dict[str, object], sequence: str) -> str:
    locus = int(row["pos"]) - 1
    ref, alt = str(row["ref"]), str(row["alt"])
    offset = allele_offset(start, len(read), locus, len(ref), len(sequence))
    if offset is None:
        return read
    if read[offset:offset + len(ref)] != ref:
        fail(f"synthetic read REF mismatch for {row['truth_id']} at read offset {offset}")
    mutated = read[:offset] + alt + read[offset + len(ref):]
    if len(mutated) < len(read):
        continuation_start = start + len(read)
        mutated += circular_slice(sequence, continuation_start, len(read) - len(mutated))
    return mutated[:len(read)]


def deterministic_fastq(args: argparse.Namespace) -> None:
    reference = Path(args.reference).resolve()
    truth_path = Path(args.truth).resolve()
    output_prefix = Path(args.output_prefix).resolve()
    if output_prefix.with_name(output_prefix.name + "_R1.fastq.gz").exists() or output_prefix.with_name(output_prefix.name + "_R2.fastq.gz").exists():
        fail("refusing to overwrite existing generated FASTQ")
    refs = parse_fasta(reference)
    truth = validate_truth(truth_path, reference)
    if args.mt_contig not in refs:
        fail(f"mitochondrial contig {args.mt_contig!r} is absent")
    sequence = refs[args.mt_contig]
    if args.read_length < 20 or args.fragment_length < args.read_length or args.fragment_length >= len(sequence):
        fail("read length must be >=20 and fragment length must be >= read length and < contig length")
    fragments = max(1, math.ceil(args.depth * len(sequence) / (2 * args.read_length)))
    rng = random.Random(args.seed)
    starts = [int(index * len(sequence) / fragments) for index in range(fragments)]
    rng.shuffle(starts)
    selected_by_truth: dict[str, set[int]] = {}
    achieved: dict[str, float] = {}
    for row in truth:
        locus = int(row["pos"]) - 1
        eligible = []
        for index, start in enumerate(starts):
            mate_start = (start + args.fragment_length - args.read_length) % len(sequence)
            if (
                allele_offset(start, args.read_length, locus, len(str(row["ref"])), len(sequence)) is not None
                or allele_offset(mate_start, args.read_length, locus, len(str(row["ref"])), len(sequence)) is not None
            ):
                eligible.append(index)
        requested = float(row["expected_vaf"])
        count = round(requested * len(eligible))
        chosen = set(rng.sample(eligible, count)) if count else set()
        selected_by_truth[str(row["truth_id"])] = chosen
        achieved[str(row["truth_id"])] = count / len(eligible) if eligible else 0.0
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    r1_path = output_prefix.with_name(output_prefix.name + "_R1.fastq.gz")
    r2_path = output_prefix.with_name(output_prefix.name + "_R2.fastq.gz")
    quality = "I" * args.read_length
    r1_raw = r1_path.open("wb")
    r2_raw = r2_path.open("wb")
    r1_gzip = gzip.GzipFile(filename="", mode="wb", fileobj=r1_raw, mtime=0)
    r2_gzip = gzip.GzipFile(filename="", mode="wb", fileobj=r2_raw, mtime=0)
    with io.TextIOWrapper(r1_gzip, encoding="utf-8") as r1, io.TextIOWrapper(r2_gzip, encoding="utf-8") as r2:
        for index, start in enumerate(starts):
            read1 = circular_slice(sequence, start, args.read_length)
            mate_start = (start + args.fragment_length - args.read_length) % len(sequence)
            mate_forward = circular_slice(sequence, mate_start, args.read_length)
            for row in truth:
                if index not in selected_by_truth[str(row["truth_id"])]:
                    continue
                read1 = apply_allele_to_read(read1, start, row, sequence)
                mate_forward = apply_allele_to_read(mate_forward, mate_start, row, sequence)
            name = f"{args.sample}:{index + 1}:{start + 1}:{args.seed}"
            r1.write(f"@{name}/1\n{read1}\n+\n{quality}\n")
            r2_sequence = reverse_complement(mate_forward)
            r2.write(f"@{name}/2\n{r2_sequence}\n+\n{quality}\n")
    r1_raw.close()
    r2_raw.close()
    achieved_truth_path = output_prefix.with_suffix(".achieved_truth.tsv")
    achieved_truth = []
    for row in truth:
        result = dict(row)
        result["expected_vaf"] = f"{achieved[str(row['truth_id'])]:.10g}"
        result["source"] = f"{row['source']};achieved_by_deterministic_fastq_generation"
        achieved_truth.append(result)
    write_tsv(achieved_truth_path, TRUTH_FIELDS + ["variant_key"], achieved_truth)
    provenance = output_prefix.with_suffix(".provenance.tsv")
    rows = [{
        "benchmark_id": args.benchmark_id,
        "source_dataset": str(reference),
        "source_checksum": sha256(reference),
        "reference_checksum": sha256(reference),
        "generation_method": "deterministic_circular_paired_fastq",
        "generation_tool": "benchmark_harness.py generate-fastq",
        "tool_version": "1.0.0",
        "seed": args.seed,
        "requested_depth": args.depth,
        "achieved_depth": f"{2 * fragments * args.read_length / len(sequence):.6f}",
        "requested_vaf": ";".join(f"{row['truth_id']}={row['expected_vaf']}" for row in truth),
        "achieved_vaf": ";".join(f"{key}={value:.6f}" for key, value in achieved.items()),
        "variant_truth": str(achieved_truth_path),
        "generation_timestamp": dt.datetime.now(dt.timezone.utc).isoformat(),
        "output_checksum": f"{r1_path.name}={sha256(r1_path)};{r2_path.name}={sha256(r2_path)}",
    }]
    write_tsv(provenance, list(rows[0]), rows)
    print(json.dumps({"r1": str(r1_path), "r2": str(r2_path), "provenance": str(provenance), "achieved_truth": str(achieved_truth_path), "fragments": fragments, "achieved_vaf": achieved}, sort_keys=True))


def tool_version(command: list[str]) -> str:
    result = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=True)
    return result.stdout.strip().splitlines()[0]


def downsample_bam(args: argparse.Namespace) -> None:
    source = Path(args.source).resolve()
    reference = Path(args.reference).resolve()
    output = Path(args.output).resolve()
    if output.exists() or Path(str(output) + ".bai").exists():
        fail(f"refusing to overwrite output: {output}")
    if args.source_depth <= 0 or args.target_depth <= 0 or args.target_depth > args.source_depth:
        fail("target depth must be positive and no greater than source depth")
    fraction = args.target_depth / args.source_depth
    fraction_digits = f"{fraction:.8f}".split(".", 1)[1]
    seed_fraction = f"{args.seed}.{fraction_digits}"
    output.parent.mkdir(parents=True, exist_ok=True)
    command = [args.samtools, "view", "-@", str(args.threads), "-b"]
    if fraction < 1:
        command += ["-s", seed_fraction]
    command += ["-o", str(output), str(source)]
    subprocess.run(command, check=True)
    subprocess.run([args.samtools, "index", "-@", str(args.threads), str(output)], check=True)
    count = int(subprocess.run([args.samtools, "view", "-c", str(output)], text=True, stdout=subprocess.PIPE, check=True).stdout)
    depth_result = subprocess.run(
        [args.samtools, "depth", "-a", "-r", args.mt_contig, str(output)],
        text=True, stdout=subprocess.PIPE, check=True,
    )
    depth_values = [int(line.split("\t")[2]) for line in depth_result.stdout.splitlines() if line]
    if not depth_values:
        fail(f"no depth records produced for {args.mt_contig} in {output}")
    achieved_depth = sum(depth_values) / len(depth_values)
    provenance = Path(str(output) + ".provenance.tsv")
    row = {
        "benchmark_id": args.benchmark_id, "source_dataset": str(source),
        "source_checksum": sha256(source), "reference_checksum": sha256(reference),
        "generation_method": "deterministic_samtools_subsampling",
        "generation_tool": args.samtools, "tool_version": tool_version([args.samtools, "--version"]),
        "seed": args.seed, "requested_depth": args.target_depth,
        "achieved_depth": f"{achieved_depth:.6f}", "requested_vaf": "NA", "achieved_vaf": "NA",
        "variant_truth": args.truth, "generation_timestamp": dt.datetime.now(dt.timezone.utc).isoformat(),
        "output_checksum": sha256(output), "source_depth": args.source_depth,
        "sampling_fraction": f"{fraction:.10f}", "reads_retained": count,
        "command": shlex.join(command),
    }
    write_tsv(provenance, list(row), [row])


def pipeline_checksum(root: Path) -> str:
    digest = hashlib.sha256()
    paths = [root / "main.nf", root / "nextflow.config"]
    paths += sorted((root / "modules").glob("*.nf"))
    paths += sorted(path for path in (root / "bin").iterdir() if path.is_file())
    paths += sorted((root / "containers").glob("*/Dockerfile"))
    for path in paths:
        digest.update(str(path.relative_to(root)).encode())
        digest.update(bytes.fromhex(sha256(path)))
    return digest.hexdigest()


def latest_attempt(run_dir: Path) -> int:
    attempts = [int(path.name.split("-")[-1]) for path in (run_dir / "nxf").glob("attempt-*") if path.name.split("-")[-1].isdigit()]
    return max(attempts, default=0)


def nextflow_run_name(attempt_dir: Path) -> str | None:
    metadata_path = attempt_dir / "run_metadata.json"
    if metadata_path.exists():
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        if metadata.get("nextflow_run_name"):
            return str(metadata["nextflow_run_name"])
    stdout_path = attempt_dir / "command.stdout.log"
    if stdout_path.exists():
        match = re.search(r"Launching .* \[([^]]+)\]", stdout_path.read_text(encoding="utf-8", errors="replace"))
        if match:
            return match.group(1)
    return None


def prior_successful_run_name(run_dir: Path) -> str | None:
    attempts = sorted((run_dir / "nxf").glob("attempt-*"), reverse=True)
    for attempt_dir in attempts:
        metadata_path = attempt_dir / "run_metadata.json"
        if not metadata_path.exists():
            continue
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        if metadata.get("status") == "PASS":
            # Preserve the original cache lineage across repeated resumes.
            # A resumed invocation receives a new display run name, but its
            # resume_target remains the authoritative cache lineage root.
            root_target = metadata.get("resume_target")
            if root_target and root_target != "NA":
                return str(root_target)
            name = nextflow_run_name(attempt_dir)
            if name:
                return name
    return None


def run_benchmarks(args: argparse.Namespace) -> None:
    manifest = Path(args.manifest).resolve()
    rows = validate_manifest(manifest)
    selected = set(args.benchmark_id or [])
    if selected:
        unknown = selected - {row["benchmark_id"] for row in rows}
        if unknown:
            fail(f"unknown benchmark IDs: {sorted(unknown)}")
        rows = [row for row in rows if row["benchmark_id"] in selected]
    if args.resume_target and not args.resume:
        fail("--resume-target requires --resume")
    if args.resume_target and len(rows) != 1:
        fail("--resume-target requires exactly one selected benchmark row")
    root = Path(args.pipeline_root).resolve()
    run_root = Path(args.run_root).resolve()
    manifest_base = manifest.parent
    for row in rows:
        benchmark_id = row["benchmark_id"]
        run_dir = run_root / benchmark_id
        work = run_dir / "work"
        results = run_dir / "results"
        nxf = run_dir / "nxf"
        if run_dir.exists() and not args.resume:
            fail(f"run directory already exists; use --resume: {run_dir}")
        attempt = latest_attempt(run_dir) + 1
        attempt_dir = nxf / f"attempt-{attempt:03d}"
        attempt_dir.mkdir(parents=True, exist_ok=False)
        work.mkdir(parents=True, exist_ok=True)
        results.mkdir(parents=True, exist_ok=True)
        inputs = [resolve_path(value, manifest_base) for value in row["input_bam_or_fastq"].split(";")]
        sample_sheet = attempt_dir / "samplesheet.csv"
        with sample_sheet.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            if len(inputs) == 1 and inputs[0].suffix.lower() == ".bam":
                writer.writerow(["sample", "bam"])
                writer.writerow([row["sample_id"], str(inputs[0])])
            elif len(inputs) == 2:
                writer.writerow(["sample", "fastq_1", "fastq_2"])
                writer.writerow([row["sample_id"], str(inputs[0]), str(inputs[1])])
            else:
                fail(f"{benchmark_id}: expected one BAM or two FASTQ files")
        reference = resolve_path(row["reference"], manifest_base)
        nf_numt = "none" if row["numt_strategy"] == "A" else "comparison"
        command = [args.nextflow, "-log", str(attempt_dir / "nextflow.log")]
        if args.config:
            command += ["-c", str(Path(args.config).resolve())]
        command += [
            "run", str(root / "main.nf"),
            "-ansi-log", "false", "-profile", args.profile,
            "--input", str(sample_sheet), "--reference", str(reference),
            "--numt_strategy", nf_numt, "--outdir", str(results),
            "-work-dir", str(work), "-with-trace", str(attempt_dir / "trace.tsv"),
            "-with-report", str(attempt_dir / "report.html"),
            "-with-timeline", str(attempt_dir / "timeline.html"),
            "-with-dag", str(attempt_dir / "dag.dot"),
        ]
        if row["numt_truth"] != "NA":
            command += ["--numt_fixture_truth", str(resolve_path(row["numt_truth"], manifest_base))]
        if row["numt_design"] != "NA":
            command += ["--numt_fixture_design", str(resolve_path(row["numt_design"], manifest_base))]
        resume_target = args.resume_target or (prior_successful_run_name(run_dir) if args.resume else None)
        if args.resume:
            command += ["-resume"] + ([resume_target] if resume_target else [])
        metadata = {
            "benchmark_id": benchmark_id, "attempt": attempt,
            "started_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "command": command, "command_shell": shlex.join(command),
            "pipeline_root": str(root), "pipeline_checksum": pipeline_checksum(root),
            "main_nf_sha256": sha256(root / "main.nf"), "manifest_sha256": sha256(manifest),
            "nextflow_version": tool_version([args.nextflow, "-version"]),
            "resume_requested": bool(args.resume), "resume_target": resume_target or "NA", "row": row,
        }
        (attempt_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        start = time.monotonic()
        peak_work_bytes = disk_usage_bytes(work)
        with (attempt_dir / "command.stdout.log").open("w", encoding="utf-8") as stdout:
            process = subprocess.Popen(command, cwd=root, text=True, stdout=stdout, stderr=subprocess.STDOUT)
            while True:
                try:
                    process.wait(timeout=30)
                    break
                except subprocess.TimeoutExpired:
                    peak_work_bytes = max(peak_work_bytes, disk_usage_bytes(work))
            return_code = process.returncode
            peak_work_bytes = max(peak_work_bytes, disk_usage_bytes(work))
        metadata.update(
            completed_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
            wall_clock_seconds=round(time.monotonic() - start, 6), exit_code=return_code,
            status="PASS" if return_code == 0 else "FAIL",
            disk_work_bytes_peak_observed=peak_work_bytes,
        )
        metadata["nextflow_run_name"] = nextflow_run_name(attempt_dir) or "NA"
        (attempt_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        if return_code:
            fail(f"{benchmark_id}: Nextflow failed; see {attempt_dir / 'command.stdout.log'}")
        print(f"PASS\t{benchmark_id}\tattempt={attempt}\tresume={args.resume}")


def parse_info(text: str) -> dict[str, str]:
    result = {}
    for item in text.split(";"):
        key, _, value = item.partition("=")
        result[key] = value if value else "True"
    return result


def float_or_na(value: str) -> str:
    try:
        return f"{float(value):.10g}"
    except (TypeError, ValueError):
        return "NA"


def normalize_vcf(path: Path, benchmark_id: str, caller: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    sample_index = None
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip().split("\t")
                sample_index = 9 if len(header) > 9 else None
                continue
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip().split("\t")
            chrom, pos, ref = fields[0], int(fields[1]), fields[3].upper()
            alts = fields[4].upper().split(",")
            info = parse_info(fields[7])
            format_map: dict[str, str] = {}
            if sample_index is not None and len(fields) > sample_index:
                format_map = dict(zip(fields[8].split(":"), fields[sample_index].split(":")))
            af_values = format_map.get("AF", info.get("AF", "NA")).split(",")
            ad_values = format_map.get("AD", "").split(",")
            depth = format_map.get("DP", info.get("DP", "NA"))
            for index, alt in enumerate(alts):
                called_vaf = af_values[index] if index < len(af_values) else "NA"
                if called_vaf == "NA" and len(ad_values) > index + 1:
                    try:
                        called_vaf = str(int(ad_values[index + 1]) / sum(int(value) for value in ad_values))
                    except (ValueError, ZeroDivisionError):
                        pass
                rows.append({
                    "benchmark_id": benchmark_id, "caller": caller, "chrom": chrom,
                    "pos": pos, "ref": ref, "alt": alt,
                    "variant_key": f"{chrom}:{pos}:{ref}:{alt}",
                    "filter_status": fields[6], "called_vaf": float_or_na(called_vaf),
                    "called_depth": depth if depth not in {"", "."} else "NA", "truth_status": "UNMATCHED",
                })
    return rows


def normalize_mutserve(path: Path, benchmark_id: str) -> list[dict[str, object]]:
    fields, source = read_tsv(path)
    required = {"contig", "position", "ref", "alt", "depth", "vaf", "status"}
    if not required.issubset(fields):
        fail(f"mutserve2 call table lacks fields: {sorted(required - set(fields))}")
    return [{
        "benchmark_id": benchmark_id, "caller": "mutserve2", "chrom": row["contig"],
        "pos": int(row["position"]), "ref": row["ref"].upper(), "alt": row["alt"].upper(),
        "variant_key": f"{row['contig']}:{row['position']}:{row['ref'].upper()}:{row['alt'].upper()}",
        "filter_status": row["status"], "called_vaf": float_or_na(row["vaf"]),
        "called_depth": row["depth"] if row["depth"] else "NA", "truth_status": "UNMATCHED",
    } for row in source]


def find_call_source(row: dict[str, str], results: Path) -> Path:
    sample, caller = row["sample_id"], row["caller"]
    strategy, circular = row["numt_strategy"], row["circular_mode"]
    suffix = "_strategy_b" if strategy == "B" else ""
    if caller == "gatk":
        if strategy == "C":
            suffix = ""
        if circular == "standard_only":
            return results / "variants/gatk/intermediate" / f"{sample}{suffix}.standard.gatk.raw.vcf.gz"
        return results / "variants/gatk" / f"{sample}{suffix}.gatk.filtered.vcf.gz"
    if circular != "standard_only":
        fail("mutserve2 circular-aware output is unsupported by the frozen pipeline")
    return results / "variants/mutserve2" / f"{sample}{suffix}.mutserve2.calls.tsv"


def apply_numt_c(rows: list[dict[str, object]], evidence_path: Path, caller: str) -> None:
    _, evidence = read_tsv(evidence_path)
    statuses = {
        f"{row['contig']}:{row['position']}:{row['ref']}:{row['alt']}": row["strategy_c_status"]
        for row in evidence if row["caller"] == caller
    }
    for row in rows:
        row["filter_status"] = statuses.get(str(row["variant_key"]), "NOT_EVALUATED")


def normalize_calls(args: argparse.Namespace) -> None:
    manifest = Path(args.manifest).resolve()
    rows = validate_manifest(manifest)
    selected = [row for row in rows if not args.benchmark_id or row["benchmark_id"] in args.benchmark_id]
    run_root = Path(args.run_root).resolve()
    output_root = Path(args.output_root).resolve()
    for row in selected:
        results = run_root / row["benchmark_id"] / "results"
        source = find_call_source(row, results)
        if not source.exists():
            fail(f"call source does not exist for {row['benchmark_id']}: {source}")
        normalized = normalize_vcf(source, row["benchmark_id"], "gatk") if row["caller"] == "gatk" else normalize_mutserve(source, row["benchmark_id"])
        references = parse_fasta(resolve_path(row["reference"], manifest.parent))
        for call in normalized:
            chrom = str(call["chrom"])
            if chrom not in references:
                fail(f"call contig is absent from benchmark reference: {chrom}")
            pos, ref, alt = normalize_allele(
                chrom, int(call["pos"]), str(call["ref"]), str(call["alt"]), references
            )
            call.update(pos=pos, ref=ref, alt=alt, variant_key=f"{chrom}:{pos}:{ref}:{alt}")
        if row["numt_strategy"] == "C":
            evidence = results / "variants/numt" / f"{row['sample_id']}.numt_variant_evidence.tsv"
            apply_numt_c(normalized, evidence, row["caller"])
        output = output_root / row["benchmark_id"] / "normalized_calls.tsv"
        write_tsv(output, CALL_FIELDS, normalized)
        print(f"PASS\t{row['benchmark_id']}\t{len(normalized)} calls\t{output}")


def callable_record(row: dict[str, str]) -> bool:
    return row["filter_status"] in {"PASS", ".", "SITE", "CALLED", "RETAINED", "RAW_UNFILTERED"}


def safe_metric(numerator: int, denominator: int) -> str:
    return "NA" if denominator == 0 else f"{numerator / denominator:.10g}"


def match_one(row: dict[str, str], truth_path: Path, reference: Path, calls_path: Path) -> tuple[dict[str, object], list[dict[str, object]]]:
    truth = validate_truth(truth_path, reference)
    fields, calls = read_tsv(calls_path)
    if fields != CALL_FIELDS:
        fail(f"normalized call schema mismatch: expected {CALL_FIELDS}, observed {fields}")
    active = [call for call in calls if callable_record(call)]
    truth_map = {str(record["variant_key"]): record for record in truth}
    call_map: dict[str, dict[str, str]] = {}
    for call in active:
        key = call["variant_key"]
        if key in call_map:
            fail(f"duplicate normalized call allele: {key}")
        call_map[key] = call
    tp = len(set(truth_map) & set(call_map))
    fp = len(set(call_map) - set(truth_map))
    fn = len(set(truth_map) - set(call_map))
    precision = safe_metric(tp, tp + fp)
    recall = safe_metric(tp, tp + fn)
    if precision == "NA" or recall == "NA":
        f1 = "NA"
    elif float(precision) + float(recall) == 0:
        f1 = "0"
    else:
        f1 = f"{2 * float(precision) * float(recall) / (float(precision) + float(recall)):.10g}"
    common = {
        "benchmark_id": row["benchmark_id"], "sample_id": row["sample_id"],
        "depth_target": row["depth_target"], "vaf_bin": row["vaf_bin"],
        "caller": row["caller"], "numt_strategy": row["numt_strategy"],
        "circular_mode": row["circular_mode"], "variant_class": row["variant_class"],
        "boundary_class": row["boundary_class"], "numt_context": row["numt_context"],
        "replicate": row["replicate"], "seed": row["seed"],
    }
    summary = dict(common, tp=tp, fp=fp, fn=fn, precision=precision, recall=recall, f1=f1)
    detail: list[dict[str, object]] = []
    for key, record in truth_map.items():
        call = call_map.get(key)
        expected = float(record["expected_vaf"])
        called = float(call["called_vaf"]) if call and call["called_vaf"] != "NA" else None
        absolute = abs(called - expected) if called is not None else None
        relative = absolute / expected if absolute is not None and expected != 0 else None
        detail.append(dict(
            common, vaf_bin=vaf_bin(expected), variant_class=record["variant_class"],
            boundary_class=record["boundary_class"], numt_context=record["numt_context"],
            truth_id=record["truth_id"], variant_key=key,
            chrom=record["chrom"], pos=record["pos"], ref=record["ref"], alt=record["alt"],
            truth_status="TP" if call else "FN", expected_or_achieved_vaf=record["expected_vaf"],
            called_vaf=call["called_vaf"] if call else "NA",
            absolute_vaf_error=f"{absolute:.10g}" if absolute is not None else "NA",
            relative_vaf_error=f"{relative:.10g}" if relative is not None else "NA",
            called_depth=call["called_depth"] if call else "NA",
            filter_status=call["filter_status"] if call else "NA",
        ))
    for key in sorted(set(call_map) - set(truth_map)):
        call = call_map[key]
        called_for_bin = float(call["called_vaf"]) if call["called_vaf"] != "NA" else float(row["vaf_target"])
        detail.append(dict(
            common, vaf_bin=vaf_bin(called_for_bin),
            variant_class="SNV" if len(call["ref"]) == len(call["alt"]) == 1 else "indel",
            truth_id="NA", variant_key=key, chrom=call["chrom"], pos=call["pos"],
            ref=call["ref"], alt=call["alt"], truth_status="FP",
            expected_or_achieved_vaf="NA", called_vaf=call["called_vaf"],
            absolute_vaf_error="NA", relative_vaf_error="NA",
            called_depth=call["called_depth"], filter_status=call["filter_status"],
        ))
    return summary, detail


def match_calls(args: argparse.Namespace) -> None:
    manifest = Path(args.manifest).resolve()
    rows = validate_manifest(manifest)
    selected = [row for row in rows if not args.benchmark_id or row["benchmark_id"] in args.benchmark_id]
    metrics_root = Path(args.metrics_root).resolve()
    base = manifest.parent
    call_fields = [
        "benchmark_id", "sample_id", "depth_target", "vaf_bin", "caller", "numt_strategy",
        "circular_mode", "variant_class", "boundary_class", "numt_context", "replicate", "seed",
        "tp", "fp", "fn", "precision", "recall", "f1",
    ]
    variant_fields = [
        "benchmark_id", "sample_id", "depth_target", "vaf_bin", "caller", "numt_strategy",
        "circular_mode", "variant_class", "boundary_class", "numt_context", "replicate", "seed",
        "truth_id", "variant_key", "chrom", "pos", "ref", "alt", "truth_status",
        "expected_or_achieved_vaf", "called_vaf", "absolute_vaf_error", "relative_vaf_error",
        "called_depth", "filter_status",
    ]
    for row in selected:
        calls = metrics_root / row["benchmark_id"] / "normalized_calls.tsv"
        summary, details = match_one(row, resolve_path(row["truth_vcf"], base), resolve_path(row["reference"], base), calls)
        groups: dict[tuple[str, str, str, str], list[dict[str, object]]] = {}
        for detail in details:
            key = tuple(str(detail[field]) for field in ("vaf_bin", "variant_class", "boundary_class", "numt_context"))
            groups.setdefault(key, []).append(detail)
        grouped_summaries = []
        for key, records in sorted(groups.items()):
            tp = sum(record["truth_status"] == "TP" for record in records)
            fp = sum(record["truth_status"] == "FP" for record in records)
            fn = sum(record["truth_status"] == "FN" for record in records)
            precision = safe_metric(tp, tp + fp)
            recall = safe_metric(tp, tp + fn)
            if precision == "NA" or recall == "NA":
                f1 = "NA"
            elif float(precision) + float(recall) == 0:
                f1 = "0"
            else:
                f1 = f"{2 * float(precision) * float(recall) / (float(precision) + float(recall)):.10g}"
            grouped = {field: summary[field] for field in call_fields if field in summary}
            grouped.update(vaf_bin=key[0], variant_class=key[1], boundary_class=key[2], numt_context=key[3], tp=tp, fp=fp, fn=fn, precision=precision, recall=recall, f1=f1)
            grouped_summaries.append(grouped)
        target = metrics_root / row["benchmark_id"]
        write_tsv(target / "benchmark_call_metrics.tsv", call_fields, grouped_summaries)
        write_tsv(target / "benchmark_variant_metrics.tsv", variant_fields, details)
        print(f"PASS\t{row['benchmark_id']}\tTP={summary['tp']} FP={summary['fp']} FN={summary['fn']} strata={len(grouped_summaries)}")


def parse_duration(value: str) -> float:
    if not value or value in {"-", "NA"}:
        return 0.0
    total = 0.0
    for number, unit in re.findall(r"([0-9.]+)\s*(ms|us|ns|d|h|m|s)", value):
        factor = {"d": 86400, "h": 3600, "m": 60, "s": 1, "ms": 1e-3, "us": 1e-6, "ns": 1e-9}[unit]
        total += float(number) * factor
    return total


def parse_memory_mb(value: str) -> float:
    if not value or value in {"-", "NA"}:
        return 0.0
    match = re.fullmatch(r"([0-9.]+)\s*([KMGTP]?B)", value.strip(), re.I)
    if not match:
        return 0.0
    factor = {"B": 1 / 1e6, "KB": 1e3 / 1e6, "MB": 1, "GB": 1e3, "TB": 1e6, "PB": 1e9}[match.group(2).upper()]
    return float(match.group(1)) * factor


def resource_metrics(args: argparse.Namespace) -> None:
    manifest = Path(args.manifest).resolve()
    rows = validate_manifest(manifest)
    selected = [row for row in rows if not args.benchmark_id or row["benchmark_id"] in args.benchmark_id]
    run_root = Path(args.run_root).resolve()
    output_root = Path(args.output_root).resolve()
    base = manifest.parent
    for row in selected:
        run_dir = run_root / row["benchmark_id"]
        attempt = args.attempt or latest_attempt(run_dir)
        attempt_dir = run_dir / "nxf" / f"attempt-{attempt:03d}"
        metadata = json.loads((attempt_dir / "run_metadata.json").read_text(encoding="utf-8"))
        trace_fields, trace_rows = read_tsv(attempt_dir / "trace.tsv")
        process_count = len(trace_rows)
        cached = sum(record.get("status") == "CACHED" for record in trace_rows)
        submitted = process_count - cached
        cpu_seconds = 0.0
        peak_rss = 0.0
        max_container_memory = 0.0
        active_rows = [record for record in trace_rows if record.get("status") != "CACHED"]
        for record in active_rows:
            duration = parse_duration(record.get("realtime", record.get("duration", "")))
            cpu_value = record.get("%cpu", "").strip().rstrip("%")
            try:
                cpu_seconds += duration * float(cpu_value) / 100 if cpu_value else duration * float(record.get("cpus", "1") or 1)
            except ValueError:
                cpu_seconds += duration
            peak_rss = max(peak_rss, parse_memory_mb(record.get("peak_rss", record.get("rss", ""))))
            max_container_memory = max(max_container_memory, parse_memory_mb(record.get("memory", "")))
        inputs = [resolve_path(value, base) for value in row["input_bam_or_fastq"].split(";")]
        inputs += [resolve_path(row["reference"], base), resolve_path(row["truth_vcf"], base)]
        result = {
            "benchmark_id": row["benchmark_id"], "attempt": attempt,
            "wall_clock_seconds": metadata["wall_clock_seconds"],
            "cpu_seconds": f"{cpu_seconds:.6f}", "cpu_hours": f"{cpu_seconds / 3600:.10g}",
            "peak_rss_mb": f"{peak_rss:.6f}" if peak_rss else "NA",
            "max_container_memory_mb": f"{max_container_memory:.6f}" if max_container_memory else "NA",
            "disk_input_bytes": sum(path.stat().st_size for path in inputs),
            "disk_work_bytes_peak": metadata.get("disk_work_bytes_peak_observed", disk_usage_bytes(run_dir / "work")),
            "disk_output_bytes": disk_usage_bytes(run_dir / "results"),
            "process_count": process_count, "cached_process_count": cached,
            "submitted_process_count": submitted, "nextflow_status": metadata["status"],
        }
        output = output_root / row["benchmark_id"] / f"benchmark_resource_metrics.attempt-{attempt:03d}.tsv"
        write_tsv(output, RESOURCE_FIELDS, [result])
        print(f"PASS\t{row['benchmark_id']}\tattempt={attempt}\tcached={cached}\tsubmitted={submitted}")


def summarize(args: argparse.Namespace) -> None:
    manifest = Path(args.manifest).resolve()
    rows = validate_manifest(manifest)
    metrics_root = Path(args.metrics_root).resolve()
    output = Path(args.output).resolve()
    call_rows: list[dict[str, str]] = []
    variant_rows: list[dict[str, str]] = []
    resource_rows: list[dict[str, str]] = []
    for row in rows:
        target = metrics_root / row["benchmark_id"]
        _, calls = read_tsv(target / "benchmark_call_metrics.tsv")
        _, variants = read_tsv(target / "benchmark_variant_metrics.tsv")
        resources = sorted(target.glob("benchmark_resource_metrics.attempt-*.tsv"))
        if not resources:
            fail(f"missing resource metrics for {row['benchmark_id']}")
        resource_candidates = []
        for resource_path in resources:
            _, resource_rows_for_attempt = read_tsv(resource_path)
            resource_candidates.extend(resource_rows_for_attempt)
        resource = [max(resource_candidates, key=lambda item: (int(item["submitted_process_count"]), int(item["attempt"])))]
        call_rows.extend(calls)
        variant_rows.extend(variants)
        resource_rows.extend(resource)
    output.mkdir(parents=True, exist_ok=True)
    write_tsv(output / "benchmark_call_metrics.tsv", list(call_rows[0]), call_rows)
    write_tsv(output / "benchmark_variant_metrics.tsv", list(variant_rows[0]), variant_rows)
    write_tsv(output / "benchmark_resource_metrics.tsv", RESOURCE_FIELDS, resource_rows)
    resource_map = {row["benchmark_id"]: row for row in resource_rows}
    summaries = []
    for row in call_rows:
        resource = resource_map[row["benchmark_id"]]
        summaries.append({**row, **{field: resource[field] for field in RESOURCE_FIELDS if field != "benchmark_id"}})
    write_tsv(output / "benchmark_run_summary.tsv", list(summaries[0]), summaries)
    print(f"PASS\t{len(rows)} benchmark rows summarized\t{output}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    item = sub.add_parser("validate-manifest")
    item.add_argument("--manifest", required=True)
    item.add_argument("--allow-missing-files", action="store_true")

    item = sub.add_parser("validate-truth")
    item.add_argument("--truth", required=True)
    item.add_argument("--reference", required=True)
    item.add_argument("--output")

    item = sub.add_parser("generate-fastq")
    item.add_argument("--benchmark-id", required=True)
    item.add_argument("--sample", required=True)
    item.add_argument("--reference", required=True)
    item.add_argument("--truth", required=True)
    item.add_argument("--output-prefix", required=True)
    item.add_argument("--depth", type=int, required=True, choices=sorted(DEPTHS))
    item.add_argument("--seed", type=int, required=True)
    item.add_argument("--mt-contig", default="chrM")
    item.add_argument("--read-length", type=int, default=100)
    item.add_argument("--fragment-length", type=int, default=180)

    item = sub.add_parser("downsample-bam")
    item.add_argument("--benchmark-id", required=True)
    item.add_argument("--source", required=True)
    item.add_argument("--reference", required=True)
    item.add_argument("--truth", default="NA")
    item.add_argument("--output", required=True)
    item.add_argument("--source-depth", type=float, required=True)
    item.add_argument("--target-depth", type=float, required=True)
    item.add_argument("--seed", type=int, required=True)
    item.add_argument("--threads", type=int, default=2)
    item.add_argument("--mt-contig", default="chrM")
    item.add_argument("--samtools", default="samtools")

    item = sub.add_parser("run")
    item.add_argument("--manifest", required=True)
    item.add_argument("--benchmark-id", action="append")
    item.add_argument("--pipeline-root", default=str(Path(__file__).resolve().parents[2]))
    item.add_argument("--run-root", default=str(Path(__file__).resolve().parents[1] / "runs"))
    item.add_argument("--nextflow", default="nextflow")
    item.add_argument("--profile", default="docker")
    item.add_argument("--config")
    item.add_argument("--resume", action="store_true")
    item.add_argument("--resume-target", help="explicit benchmark-specific Nextflow run name for cache-lineage validation")

    item = sub.add_parser("normalize-calls")
    item.add_argument("--manifest", required=True)
    item.add_argument("--benchmark-id", action="append")
    item.add_argument("--run-root", default=str(Path(__file__).resolve().parents[1] / "runs"))
    item.add_argument("--output-root", default=str(Path(__file__).resolve().parents[1] / "metrics"))

    item = sub.add_parser("match")
    item.add_argument("--manifest", required=True)
    item.add_argument("--benchmark-id", action="append")
    item.add_argument("--metrics-root", default=str(Path(__file__).resolve().parents[1] / "metrics"))

    item = sub.add_parser("resources")
    item.add_argument("--manifest", required=True)
    item.add_argument("--benchmark-id", action="append")
    item.add_argument("--attempt", type=int)
    item.add_argument("--run-root", default=str(Path(__file__).resolve().parents[1] / "runs"))
    item.add_argument("--output-root", default=str(Path(__file__).resolve().parents[1] / "metrics"))

    item = sub.add_parser("summarize")
    item.add_argument("--manifest", required=True)
    item.add_argument("--metrics-root", default=str(Path(__file__).resolve().parents[1] / "metrics"))
    item.add_argument("--output", default=str(Path(__file__).resolve().parents[1] / "summaries"))
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    if args.command == "validate-manifest":
        rows = validate_manifest(Path(args.manifest).resolve(), not args.allow_missing_files)
        print(f"PASS\t{len(rows)} benchmark rows")
    elif args.command == "validate-truth":
        rows = validate_truth(Path(args.truth).resolve(), Path(args.reference).resolve(), Path(args.output).resolve() if args.output else None)
        print(f"PASS\t{len(rows)} normalized truth alleles")
    elif args.command == "generate-fastq":
        deterministic_fastq(args)
    elif args.command == "downsample-bam":
        downsample_bam(args)
    elif args.command == "run":
        run_benchmarks(args)
    elif args.command == "normalize-calls":
        normalize_calls(args)
    elif args.command == "match":
        match_calls(args)
    elif args.command == "resources":
        resource_metrics(args)
    elif args.command == "summarize":
        summarize(args)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except HarnessError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
