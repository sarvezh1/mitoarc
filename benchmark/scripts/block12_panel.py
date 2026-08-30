#!/usr/bin/env python3
"""Generate and QC the sparse controlled multiplexed benchmark matrix.

The generator creates deterministic circular chrM paired reads plus controlled
synthetic nuclear-homolog read groups. It never modifies its canonical source
reference and records achieved, rather than nominal, allele fractions.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import math
import random
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BENCH = ROOT / "benchmark"
SOURCE_REF = BENCH / "datasets/block12/reference/block12_reference.fa"
DEPTHS = (100, 250, 500, 1000, 2000, 5000)
TARGETS = (0.005, 0.015, 0.035, 0.075, 0.25, 0.75)
REPLICATE_SEEDS = {"R1": 120101, "R2": 120202, "R3": 120303}
READ_LENGTH = 150
FRAGMENT_LENGTH = 450
TRUTH_FIELDS = ["truth_id", "chrom", "pos", "ref", "alt", "variant_class",
                "expected_vaf", "boundary_class", "numt_context", "source", "notes"]
MANIFEST_FIELDS = ["benchmark_id", "sample_id", "dataset_type", "reference",
                   "input_bam_or_fastq", "truth_vcf", "depth_target", "vaf_target",
                   "vaf_bin", "variant_class", "variant_position", "boundary_class",
                   "numt_context", "caller", "numt_strategy", "circular_mode",
                   "replicate", "seed", "numt_truth", "numt_design"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_fasta(path: Path) -> dict[str, str]:
    records, name = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = []
        elif line.strip():
            records[name].append(line.strip().upper())
    return {key: "".join(value) for key, value in records.items()}


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as out:
        for name, sequence in records:
            out.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                out.write(sequence[start:start + 80] + "\n")


def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as out:
        writer = csv.DictWriter(out, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "NA") for field in fields})


def circular_slice(sequence: str, start: int, length: int) -> str:
    return "".join(sequence[(start + i) % len(sequence)] for i in range(length))


def revcomp(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def distance_to_boundary(pos: int, length: int) -> int:
    return min(pos - 1, length - pos)


def alt_base(ref: str, salt: int) -> str:
    bases = "ACGT"
    return bases[(bases.index(ref) + 1 + salt % 3) % 4]


def find_positions(sequence: str) -> tuple[list[int], list[int]]:
    """Select separated, non-homopolymer technical sites reproducibly."""
    nonboundary = []
    for target in range(1800, 15001, 540):
        for pos in range(target, target + 80):
            window = sequence[pos - 4:pos + 4]
            if len(window) == 8 and len(set(window)) >= 3 and sequence[pos - 2] != sequence[pos - 1]:
                if not nonboundary or pos - nonboundary[-1] > READ_LENGTH:
                    nonboundary.append(pos)
                    break
        if len(nonboundary) == 24:
            break
    candidates = list(range(35, 995, 80)) + list(range(len(sequence) - 955, len(sequence) + 1, 80))
    boundary = []
    for target in candidates:
        for delta in range(30):
            pos = target + delta
            if pos > len(sequence):
                continue
            window = sequence[pos - 4:pos + 4]
            if len(window) == 8 and len(set(window)) >= 3 and sequence[pos - 2] != sequence[pos - 1]:
                if all(min((pos - old) % len(sequence), (old - pos) % len(sequence)) > READ_LENGTH for old in boundary):
                    boundary.append(pos)
                    break
        if len(boundary) == 12:
            break
    if len(nonboundary) != 24 or len(boundary) != 12:
        raise RuntimeError(f"unable to select panel positions: nonboundary={len(nonboundary)} boundary={len(boundary)}")
    return nonboundary, boundary


def panel_template(sequence: str) -> list[dict]:
    nonboundary, boundary = find_positions(sequence)
    n_idx = b_idx = 0
    rows = []
    archetypes = [
        ("SNV", "NON_BOUNDARY", "NO_NUMT_CHALLENGE"),
        ("SNV", "BOUNDARY", "NO_NUMT_CHALLENGE"),
        ("SNV", "NON_BOUNDARY", "NUMT_LIKE"),
        ("SNV", "NON_BOUNDARY", "AMBIGUOUS_CONTEXT"),
        ("indel", "NON_BOUNDARY", "NO_NUMT_CHALLENGE"),
        ("indel", "BOUNDARY", "NO_NUMT_CHALLENGE"),
    ]
    for vi, target in enumerate(TARGETS, 1):
        for ai, (kind, boundary_class, context) in enumerate(archetypes, 1):
            if boundary_class == "BOUNDARY":
                pos, b_idx = boundary[b_idx], b_idx + 1
            else:
                pos, n_idx = nonboundary[n_idx], n_idx + 1
            ref = sequence[pos - 1]
            alt = alt_base(ref, vi + ai) if kind == "SNV" else ref + alt_base(ref, vi + ai)
            rows.append({
                "truth_id": f"V{vi}_A{ai}", "chrom": "chrM", "pos": pos,
                "ref": ref, "alt": alt, "variant_class": kind,
                "expected_vaf": target, "boundary_class": boundary_class,
                "numt_context": context, "source": "block12_controlled_multiplex_panel",
                "notes": (f"nominal_target={target};archetype={ai};"
                          f"distance_to_boundary={distance_to_boundary(pos, len(sequence))}"),
                "target_vaf": target, "archetype": ai,
                "distance_to_boundary": distance_to_boundary(pos, len(sequence)),
            })
    return rows


def challenge_template(sequence: str, panel: list[dict]) -> list[dict]:
    occupied = {int(row["pos"]) for row in panel}
    candidates = []
    # Scan, rather than use a second fixed grid, so the challenge sites fall
    # into available gaps between the 24 non-boundary truth contexts.
    for pos in range(1150, 15500):
        if all(abs(pos - old) > READ_LENGTH for old in occupied | set(candidates)):
            window = sequence[pos - 4:pos + 4]
            if len(window) == 8 and len(set(window)) >= 3:
                candidates.append(pos)
        if len(candidates) == 12:
            break
    if len(candidates) != 12:
        raise RuntimeError(f"unable to select 12 independent NUMT challenge sites; found {len(candidates)}")
    rows = []
    for vi, target in enumerate(TARGETS, 1):
        for label, truth_class, origin in (("NF", "NUMT_FALSE", "synthetic_nuclear_homolog"),
                                           ("AF", "AMBIGUOUS", "chrM_or_synthetic_nuclear")):
            pos = candidates[len(rows)]
            ref = sequence[pos - 1]
            alt = alt_base(ref, vi + len(rows))
            rows.append({"case_id": f"V{vi}_{label}", "pos": pos, "ref": ref, "alt": alt,
                         "target_vaf": target, "truth_class": truth_class,
                         "generating_origin": origin,
                         "mitochondrial_homologous_position": pos,
                         "expected_mitochondrial_allele_state": "REFERENCE",
                         "expected_false_signal_state": "ALT_PRESENT_IN_CHALLENGE_READS",
                         "mapping_characteristics": ("nuclear_unique_mate_plus_mt_homolog" if label == "NF"
                                                     else "equal_mt_nuclear_homology")})
    return rows


def overlapping_read(start: int, read_length: int, pos0: int, genome_length: int) -> int | None:
    for offset in range(read_length):
        if (start + offset) % genome_length == pos0:
            return offset
    return None


def mutate_read(read: str, offset: int, ref: str, alt: str, reference: str, continuation: int) -> str:
    if read[offset:offset + len(ref)] != ref:
        raise RuntimeError("read REF mismatch during multiplex generation")
    value = read[:offset] + alt + read[offset + len(ref):]
    if len(value) < len(read):
        value += circular_slice(reference, continuation, len(read) - len(value))
    return value[:len(read)]


def build_reference(sequence: str, challenges: list[dict], panel: list[dict], path: Path) -> None:
    records = [("chrM", sequence)]
    # The frozen FASTQ alignment contract enables NUMT comparison only when
    # all 22 canonical autosome names are present. Deterministic neutral
    # scaffolds satisfy that structural contract; they are explicitly not a
    # complete biological WGS background and are not used as accuracy truth.
    for chrom in range(1, 23):
        rng = random.Random(900000 + chrom)
        records.append((f"chr{chrom}", "".join(rng.choice("ACGT") for _ in range(4000))))
    # One controlled homolog per challenge. Unique deterministic flanks allow
    # NUMT-like pairs to retain known nuclear origin while interior pairs tie.
    for index, row in enumerate(challenges, 1):
        pos = int(row["pos"])
        homolog = circular_slice(sequence, pos - 301, 600)
        rng = random.Random(910000 + index)
        left = "".join(rng.choice("ACGT") for _ in range(240))
        right = "".join(rng.choice("ACGT") for _ in range(240))
        records.append((f"chrN_NUMT_{index:02d}", left + homolog + right))
    context_rows = [row for row in panel if row["numt_context"] != "NO_NUMT_CHALLENGE"]
    for index, row in enumerate(context_rows, 1):
        pos = int(row["pos"])
        homolog = circular_slice(sequence, pos - 301, 600)
        rng = random.Random(920000 + index)
        left = "".join(rng.choice("ACGT") for _ in range(240))
        right = "".join(rng.choice("ACGT") for _ in range(240))
        label = "LIKE" if row["numt_context"] == "NUMT_LIKE" else "AMB"
        records.append((f"chrN_TRUE_{label}_{index:02d}", left + homolog + right))
    write_fasta(path, records)


def generate_dataset(sequence: str, panel: list[dict], challenges: list[dict], replicate: str,
                     depth: int, seed: int, reference_path: Path, outdir: Path) -> dict:
    benchmark_id = f"block12_{replicate}_d{depth}"
    sample = f"B12{replicate}D{depth}"
    prefix = outdir / sample
    for suffix in ("_R1.fastq.gz", "_R2.fastq.gz", ".provenance.tsv", ".achieved_truth.tsv"):
        if prefix.with_name(prefix.name + suffix).exists() if suffix.startswith("_") else prefix.with_suffix(suffix).exists():
            raise RuntimeError(f"refusing to overwrite generated dataset {benchmark_id}")
    fragments = math.ceil(depth * len(sequence) / (2 * READ_LENGTH))
    rng = random.Random(seed)
    starts = [int(i * len(sequence) / fragments) for i in range(fragments)]
    rng.shuffle(starts)
    chosen_by_truth, achieved = {}, {}
    for row in panel:
        pos0 = int(row["pos"]) - 1
        eligible = []
        for index, start in enumerate(starts):
            mate = (start + FRAGMENT_LENGTH - READ_LENGTH) % len(sequence)
            if overlapping_read(start, READ_LENGTH, pos0, len(sequence)) is not None or overlapping_read(mate, READ_LENGTH, pos0, len(sequence)) is not None:
                eligible.append(index)
        requested = float(row["target_vaf"])
        count = min(len(eligible), max(1, math.floor(requested * len(eligible) + 0.5)))
        # Prefer junction-crossing fragments for boundary alleles, while still
        # filling high target fractions from all eligible reads.
        if row["boundary_class"] == "BOUNDARY":
            crossing = [i for i in eligible if starts[i] + READ_LENGTH > len(sequence)
                        or ((starts[i] + FRAGMENT_LENGTH - READ_LENGTH) % len(sequence)) + READ_LENGTH > len(sequence)]
            rest = [i for i in eligible if i not in set(crossing)]
            rng.shuffle(crossing); rng.shuffle(rest)
            selected = set((crossing + rest)[:count])
        else:
            selected = set(rng.sample(eligible, count))
        chosen_by_truth[row["truth_id"]] = selected
        achieved[row["truth_id"]] = count / len(eligible)

    outdir.mkdir(parents=True, exist_ok=True)
    r1_path = prefix.with_name(prefix.name + "_R1.fastq.gz")
    r2_path = prefix.with_name(prefix.name + "_R2.fastq.gz")
    quality = "I" * READ_LENGTH
    handles = []
    for path in (r1_path, r2_path):
        raw = path.open("wb")
        gz = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
        text = io.TextIOWrapper(gz, encoding="utf-8")
        handles.append((raw, text))
    r1, r2 = handles[0][1], handles[1][1]
    for index, start in enumerate(starts):
        mate_start = (start + FRAGMENT_LENGTH - READ_LENGTH) % len(sequence)
        read1 = circular_slice(sequence, start, READ_LENGTH)
        mate_forward = circular_slice(sequence, mate_start, READ_LENGTH)
        for row in panel:
            if index not in chosen_by_truth[row["truth_id"]]:
                continue
            for which, read, read_start in ((1, read1, start), (2, mate_forward, mate_start)):
                offset = overlapping_read(read_start, READ_LENGTH, int(row["pos"]) - 1, len(sequence))
                if offset is not None:
                    changed = mutate_read(read, offset, row["ref"], row["alt"], sequence, read_start + READ_LENGTH)
                    if which == 1: read1 = changed
                    else: mate_forward = changed
        name = f"MT{replicate}{depth}:{index + 1}:{start + 1}"
        r1.write(f"@{name}/1\n{read1}\n+\n{quality}\n")
        r2.write(f"@{name}/2\n{revcomp(mate_forward)}\n+\n{quality}\n")

    challenge_counts = {}
    for ci, row in enumerate(challenges, 1):
        count = max(1, math.floor(float(row["target_vaf"]) * depth + 0.5))
        challenge_counts[row["case_id"]] = count
        pos = int(row["pos"])
        homolog_start = pos - 76
        alt_read = circular_slice(sequence, homolog_start - 1, READ_LENGTH)
        alt_offset = pos - homolog_start
        alt_read = mutate_read(alt_read, alt_offset, row["ref"], row["alt"], sequence, homolog_start - 1 + READ_LENGTH)
        nuclear = read_fasta(reference_path)[f"chrN_NUMT_{ci:02d}"]
        unique_mate = nuclear[-READ_LENGTH:]
        for j in range(count):
            jitter = j % 7
            seq1 = circular_slice(sequence, homolog_start - 1 - jitter, READ_LENGTH)
            off = pos - (homolog_start - jitter)
            seq1 = mutate_read(seq1, off, row["ref"], row["alt"], sequence, homolog_start - 1 - jitter + READ_LENGTH)
            if row["truth_class"] == "NUMT_FALSE":
                seq2 = unique_mate
            else:
                seq2 = circular_slice(sequence, pos + 25 + jitter, READ_LENGTH)
            name = f"{row['case_id']}X{replicate}:{j + 1}"
            r1.write(f"@{name}/1\n{seq1}\n+\n{quality}\n")
            r2.write(f"@{name}/2\n{revcomp(seq2)}\n+\n{quality}\n")
    r1.close(); r2.close()
    for raw, _ in handles: raw.close()

    truth_rows = []
    for row in panel:
        result = {field: row.get(field, "NA") for field in TRUTH_FIELDS}
        result["expected_vaf"] = f"{achieved[row['truth_id']]:.10g}"
        result["source"] += ";achieved_by_block12_read_count"
        truth_rows.append(result)
    truth_path = prefix.with_suffix(".achieved_truth.tsv")
    write_tsv(truth_path, TRUTH_FIELDS, truth_rows)
    numt_rows = [{"contig": "chrM", "pos": row["pos"], "ref": row["ref"], "alt": row["alt"],
                  "truth_class": "TRUE_MT", "case_id": row["truth_id"], "generating_origin": "chrM"}
                 for row in panel]
    numt_rows += [{"contig": "chrM", "pos": row["pos"], "ref": row["ref"], "alt": row["alt"],
                   "truth_class": row["truth_class"], "case_id": row["case_id"],
                   "generating_origin": row["generating_origin"]} for row in challenges]
    numt_path = prefix.with_suffix(".numt_truth.tsv")
    write_tsv(numt_path, ["contig", "pos", "ref", "alt", "truth_class", "case_id", "generating_origin"], numt_rows)
    design_rows = []
    for row in challenges:
        design_rows.append({"read_prefix": f"{row['case_id']}X{replicate}:", "pairs": challenge_counts[row["case_id"]],
                            "purpose": "controlled_numt_false_signal", "canonical_variant": f"chrM:{row['pos']}:{row['ref']}>{row['alt']}",
                            "generating_origin": row["generating_origin"], "expected_fullref_context": row["mapping_characteristics"]})
    design_path = prefix.with_suffix(".numt_design.tsv")
    write_tsv(design_path, ["read_prefix", "pairs", "purpose", "canonical_variant", "generating_origin", "expected_fullref_context"], design_rows)
    achieved_depth = 2 * fragments * READ_LENGTH / len(sequence)
    provenance_path = prefix.with_suffix(".provenance.tsv")
    provenance = {"benchmark_id": benchmark_id, "source_dataset": str(SOURCE_REF),
                  "source_checksum": sha256(SOURCE_REF), "reference_checksum": sha256(reference_path),
                  "generation_method": "deterministic_multiplexed_circular_fastq_with_controlled_numt_challenges",
                  "generation_tool": "block12_panel.py", "tool_version": "1.0.0", "seed": seed,
                  "requested_depth": depth, "achieved_depth": f"{achieved_depth:.8f}",
                  "requested_vaf": ";".join(f"{r['truth_id']}={r['target_vaf']}" for r in panel),
                  "achieved_vaf": ";".join(f"{key}={value:.10g}" for key, value in achieved.items()),
                  "variant_truth": str(truth_path), "generation_timestamp": dt.datetime.now(dt.timezone.utc).isoformat(),
                  "output_checksum": f"{r1_path.name}={sha256(r1_path)};{r2_path.name}={sha256(r2_path)}",
                  "reads_retained": 2 * (fragments + sum(challenge_counts.values())),
                  "mt_fragments": fragments, "challenge_pairs": sum(challenge_counts.values()),
                  "read_length": READ_LENGTH, "fragment_length": FRAGMENT_LENGTH}
    write_tsv(provenance_path, list(provenance), [provenance])
    return {"benchmark_id": benchmark_id, "sample_id": sample, "depth": depth, "replicate": replicate,
            "seed": seed, "r1": r1_path, "r2": r2_path, "truth": truth_path,
            "numt_truth": numt_path, "numt_design": design_path, "provenance": provenance_path,
            "achieved_depth": achieved_depth, "fragments": fragments}


def dataset_seed(replicate: str, depth: int) -> int:
    return REPLICATE_SEEDS[replicate] + DEPTHS.index(depth) * 1009


def relative(path: Path, base: Path) -> str:
    return str(path.resolve().relative_to(base.resolve())) if path.resolve().is_relative_to(base.resolve()) else str(path.resolve())


def create(args) -> None:
    if not SOURCE_REF.is_file():
        raise RuntimeError(f"canonical reference not found: {SOURCE_REF}")
    sequence = read_fasta(SOURCE_REF).get("chrM")
    if not sequence or len(sequence) != 16569:
        raise RuntimeError("source must contain canonical chrM length 16569")
    panel = panel_template(sequence)
    challenges = challenge_template(sequence, panel)
    ref_path = BENCH / "datasets/block12/reference/block12_reference.fa"
    build_reference(sequence, challenges, panel, ref_path)
    design_dir = BENCH / "truth/block12"
    write_tsv(design_dir / "block12_truth_panel.nominal.tsv", TRUTH_FIELDS + ["target_vaf", "archetype", "distance_to_boundary"], panel)
    challenge_fields = ["case_id", "pos", "ref", "alt", "target_vaf", "truth_class", "generating_origin",
                        "mitochondrial_homologous_position", "mapping_characteristics", "expected_mitochondrial_allele_state",
                        "expected_false_signal_state"]
    write_tsv(design_dir / "block12_numt_challenges.tsv", challenge_fields, challenges)
    matrix_rows, records = [], []
    manifest_path = BENCH / "manifests/block12_primary_matrix.tsv"
    manifest_base = manifest_path.parent
    requested = set(args.only or [])
    for replicate in REPLICATE_SEEDS:
        for depth in DEPTHS:
            benchmark_id = f"block12_{replicate}_d{depth}"
            if requested and benchmark_id not in requested:
                continue
            outdir = BENCH / "generated/block12" / benchmark_id
            records.append(generate_dataset(sequence, panel, challenges, replicate, depth,
                                            dataset_seed(replicate, depth), ref_path, outdir))
    if requested:
        print(json.dumps({"generated": [r["benchmark_id"] for r in records]}, indent=2))
        return
    for rec in records:
        matrix_rows.append({"benchmark_id": rec["benchmark_id"], "sample_id": rec["sample_id"], "dataset_type": "synthetic",
                            "reference": relative(ref_path, manifest_base),
                            "input_bam_or_fastq": f"{relative(rec['r1'], manifest_base)};{relative(rec['r2'], manifest_base)}",
                            "truth_vcf": relative(rec["truth"], manifest_base), "depth_target": rec["depth"],
                            "vaf_target": "NA", "vaf_bin": "MULTIPLEXED", "variant_class": "MIXED", "variant_position": "NA",
                            "boundary_class": "MIXED", "numt_context": "MIXED", "caller": "all_supported",
                            "numt_strategy": "comparison", "circular_mode": "all_supported",
                            "replicate": int(rec["replicate"][1:]), "seed": rec["seed"],
                            "numt_truth": relative(rec["numt_truth"], manifest_base),
                            "numt_design": relative(rec["numt_design"], manifest_base)})
    write_tsv(manifest_path, MANIFEST_FIELDS, matrix_rows)
    summary_dir = BENCH / "summaries/block12"
    write_tsv(summary_dir / "block12_primary_matrix.tsv",
              ["benchmark_id", "replicate", "depth_target", "seed", "truth_alleles", "achieved_depth", "mt_fragments"],
              [{"benchmark_id": r["benchmark_id"], "replicate": r["replicate"], "depth_target": r["depth"],
                "seed": r["seed"], "truth_alleles": len(panel), "achieved_depth": f"{r['achieved_depth']:.8f}",
                "mt_fragments": r["fragments"]} for r in records])
    print(json.dumps({"generated_inputs": len(records), "unique_pipeline_executions": len(records),
                      "truth_alleles_per_input": len(panel), "underlying_allele_observations": len(records) * len(panel),
                      "reference": str(ref_path), "manifest": str(manifest_path)}, indent=2))


def qc(args) -> None:
    sequence = read_fasta(SOURCE_REF)["chrM"]
    panel = panel_template(sequence)
    errors = []
    if len(panel) != 36: errors.append("truth allele count is not 36")
    if len({row["truth_id"] for row in panel}) != 36: errors.append("duplicate truth IDs")
    if len({int(row["pos"]) for row in panel}) != 36: errors.append("duplicate positions")
    for row in panel:
        if sequence[int(row["pos"]) - 1:int(row["pos"]) - 1 + len(row["ref"])] != row["ref"]:
            errors.append(f"REF mismatch {row['truth_id']}")
    positions = [int(row["pos"]) for row in panel]
    for i, pos in enumerate(positions):
        for other in positions[i + 1:]:
            if min((pos - other) % len(sequence), (other - pos) % len(sequence)) <= READ_LENGTH:
                errors.append(f"overlapping read contexts {pos}/{other}")
    seeds = [dataset_seed(rep, depth) for rep in REPLICATE_SEEDS for depth in DEPTHS]
    if len(seeds) != len(set(seeds)): errors.append("dataset seeds are not unique")
    counts = {}
    for target in TARGETS:
        subset = [r for r in panel if r["target_vaf"] == target]
        counts[str(target)] = len(subset)
        if len(subset) != 6: errors.append(f"target {target} does not have six archetypes")
    if errors:
        raise RuntimeError("; ".join(errors))
    print(json.dumps({"status": "PASS", "alleles": len(panel), "target_counts": counts,
                      "unique_dataset_seeds": len(seeds), "minimum_circular_spacing": min(
                          min((a-b) % len(sequence), (b-a) % len(sequence))
                          for i, a in enumerate(positions) for b in positions[i+1:])}, indent=2))


def refresh_reference(args) -> None:
    sequence = read_fasta(SOURCE_REF)["chrM"]
    panel = panel_template(sequence)
    challenges = challenge_template(sequence, panel)
    ref_path = BENCH / "datasets/block12/reference/block12_reference.fa"
    build_reference(sequence, challenges, panel, ref_path)
    reference_checksum = sha256(ref_path)
    updated = 0
    for path in sorted((BENCH / "generated/block12").glob("*/*.provenance.tsv")):
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields, rows = list(reader.fieldnames or []), list(reader)
        for row in rows:
            row["reference_checksum"] = reference_checksum
        write_tsv(path, fields, rows)
        updated += 1
    print(json.dumps({"reference": str(ref_path), "reference_checksum": reference_checksum,
                      "reference_records": len(read_fasta(ref_path)), "provenance_rows_updated": updated}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    create_parser = sub.add_parser("create")
    create_parser.add_argument("--only", action="append", help="generate only the selected stable benchmark ID")
    sub.add_parser("qc")
    sub.add_parser("refresh-reference")
    args = parser.parse_args()
    {"create": create, "qc": qc, "refresh-reference": refresh_reference}[args.command](args)


if __name__ == "__main__":
    main()
