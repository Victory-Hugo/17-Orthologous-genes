#!/usr/bin/env python3
"""
依赖列表:
- clinker (CLI)
- clustermap.js (可选，用于交互式可视化)
- Python >=3.8
- biopython
- gffutils (可选，当前脚本未强制使用)
"""
from __future__ import annotations

import argparse
import concurrent.futures
import datetime as dt
import gzip
import logging
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord


@dataclass
class InputRecord:
    faa_path: str
    anno_path: str
    anchor_id: str
    line_no: int


@dataclass
class ProcessResult:
    key: str
    sample_id: str
    status: str
    message: str
    gbk_path: Optional[str] = None
    anchor_contig: Optional[str] = None
    anchor_start: Optional[int] = None
    anchor_end: Optional[int] = None
    match_count: int = 0
    flank: int = 0
    neighborhood_genes: int = 0


def setup_logging(log_path: str) -> logging.Logger:
    logger = logging.getLogger("clinker_pipeline")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    formatter = logging.Formatter("%(asctime)s\t%(levelname)s\t%(message)s")

    fh = logging.FileHandler(log_path)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(formatter)
    logger.addHandler(sh)
    return logger


def log_versions(logger: logging.Logger) -> None:
    logger.info("Dependencies:")
    logger.info("- python: %s", sys.version.replace("\n", " "))

    try:
        import Bio

        logger.info("- biopython: %s", Bio.__version__)
    except Exception as exc:  # pragma: no cover - best effort logging
        logger.info("- biopython: unavailable (%s)", exc)

    try:
        import gffutils  # noqa: F401

        logger.info("- gffutils: %s", gffutils.__version__)
    except Exception as exc:  # pragma: no cover
        logger.info("- gffutils: unavailable (%s)", exc)

    for cmd in ["clinker", "clustermap.js"]:
        path = shutil.which(cmd)
        if not path:
            logger.info("- %s: not found in PATH", cmd)
            continue
        try:
            completed = subprocess.run([cmd, "--version"], capture_output=True, text=True)
            version_line = (completed.stdout or completed.stderr).strip().split("\n")[0]
            logger.info("- %s: %s", cmd, version_line)
        except Exception as exc:  # pragma: no cover
            logger.info("- %s: version check failed (%s)", cmd, exc)


def read_input_tsv(conf_path: str) -> List[InputRecord]:
    records: List[InputRecord] = []
    with open(conf_path, "r", encoding="utf-8") as handle:
        for idx, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                records.append(
                    InputRecord(faa_path="", anno_path="", anchor_id="", line_no=idx)
                )
                continue
            faa_path, anno_path, anchor_id = fields[:3]
            records.append(
                InputRecord(
                    faa_path=faa_path.strip(),
                    anno_path=anno_path.strip(),
                    anchor_id=anchor_id.strip(),
                    line_no=idx,
                )
            )
    return records


def load_faa_data(faa_path: str) -> Tuple[Dict[str, int], Dict[str, Seq]]:
    lengths: Dict[str, int] = {}
    seqs: Dict[str, Seq] = {}
    if not faa_path or not os.path.exists(faa_path):
        return lengths, seqs
    with open_maybe_gzip(faa_path, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            prot_id = record.id
            seqs[prot_id] = record.seq
            lengths[prot_id] = len(record.seq)
    return lengths, seqs


def open_maybe_gzip(path: str, mode: str):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode, encoding="utf-8") if "t" in mode else open(path, mode)


def find_fna_for_gff(gff_path: str) -> Optional[str]:
    search_dirs = [os.path.dirname(gff_path), os.path.dirname(os.path.dirname(gff_path))]
    patterns = ["*.fna", "*.fa", "*.fasta", "*.fna.gz", "*.fa.gz", "*.fasta.gz"]
    candidates: List[str] = []
    for d in search_dirs:
        if not d:
            continue
        for pat in patterns:
            candidates.extend(sorted(str(p) for p in Path(d).glob(pat)))
    if not candidates:
        return None

    def score(p: str) -> Tuple[int, int]:
        basename = os.path.basename(p).lower()
        priority = 0
        if "genomic" in basename:
            priority -= 5
        if basename.endswith(".fna") or basename.endswith(".fna.gz"):
            priority -= 2
        if basename.endswith(".fa") or basename.endswith(".fasta") or basename.endswith(".fa.gz") or basename.endswith(".fasta.gz"):
            priority -= 1
        if "cds" in basename or "protein" in basename or "pep" in basename or "aa" in basename or "rna" in basename:
            priority += 5
        return (priority, len(basename))

    return sorted(set(candidates), key=score)[0]


def parse_gff_cds(gff_path: str) -> List[Dict]:
    cds_list: List[Dict] = []
    with open_maybe_gzip(gff_path, "rt") as handle:
        for raw in handle:
            if raw.startswith("#"):
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, _source, ftype, start, end, _score, strand, _phase, attrs = parts
            if ftype.lower() != "cds":
                continue
            attrs_dict = parse_gff_attributes(attrs)
            cds_list.append(
                {
                    "seqid": seqid,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "attrs": attrs_dict,
                }
            )
    return cds_list


def parse_gff_attributes(attr_str: str) -> Dict[str, List[str]]:
    attrs: Dict[str, List[str]] = {}
    for item in attr_str.split(";"):
        if not item:
            continue
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        key = key.strip()
        value = value.strip()
        if not key:
            continue
        values = [v for v in value.split(",") if v]
        if values:
            attrs.setdefault(key, []).extend(values)
    return attrs


def build_seq_lookup(fna_path: str) -> Dict[str, SeqRecord]:
    seqs: Dict[str, SeqRecord] = {}
    with open_maybe_gzip(fna_path, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            for key in {record.id, record.name, record.description.split()[0]}:
                if key:
                    seqs[key] = record
    return seqs


def matches_target(values: Iterable[str], target: str) -> bool:
    for val in values:
        if val == target:
            return True
        if ":" in val and val.split(":", 1)[1] == target:
            return True
    return False


def match_in_attributes(attrs: Dict[str, List[str]], target: str) -> bool:
    keys = ["protein_id", "ID", "Name", "Dbxref", "db_xref", "locus_tag", "gene"]
    for key in keys:
        values = attrs.get(key, [])
        if matches_target(values, target):
            return True
    return False


def match_in_qualifiers(qualifiers: Dict, target: str) -> bool:
    keys = ["protein_id", "ID", "Name", "db_xref", "locus_tag", "gene"]
    for key in keys:
        values = qualifiers.get(key, [])
        if matches_target(values, target):
            return True
    return False


def choose_anchor_index(
    candidates: List[Tuple[int, int, int]], protein_len: Optional[int]
) -> int:
    if not candidates:
        raise ValueError("No candidates to choose from")
    if protein_len is None:
        return candidates[0][0]
    scored = sorted(candidates, key=lambda x: (abs(x[1] - protein_len), x[2]))
    return scored[0][0]


def safe_sample_id(path: str) -> str:
    base = os.path.basename(path)
    return re.sub(r"[^A-Za-z0-9._-]+", "_", os.path.splitext(base)[0])


def write_neighborhood_gbk(
    out_path: str,
    record: SeqRecord,
    features: List[SeqFeature],
) -> None:
    record.features = features
    record.annotations["molecule_type"] = "DNA"
    with open(out_path, "w", encoding="utf-8") as handle:
        SeqIO.write(record, handle, "genbank")


def extract_neighborhood_from_gbk(
    gbk_path: str,
    anchor_id: str,
    flank: int,
    protein_len: Optional[int],
    protein_seqs: Dict[str, Seq],
) -> Tuple[SeqRecord, List[SeqFeature], Dict]:
    matches: List[Tuple[SeqRecord, int]] = []
    cds_lists: Dict[str, List[SeqFeature]] = {}

    for record in SeqIO.parse(gbk_path, "genbank"):
        cds_list = [f for f in record.features if f.type == "CDS"]
        if not cds_list:
            continue
        cds_lists[record.id] = cds_list
        for idx, feat in enumerate(cds_list):
            if match_in_qualifiers(feat.qualifiers, anchor_id):
                matches.append((record, idx))

    if not matches:
        raise ValueError("Anchor protein ID not found in GenBank")

    candidates: List[Tuple[int, int, int]] = []
    for record, idx in matches:
        feat = cds_lists[record.id][idx]
        aa_len = None
        if "translation" in feat.qualifiers:
            aa_len = len(feat.qualifiers["translation"][0])
        elif protein_len is not None:
            aa_len = protein_len
        if aa_len is None:
            aa_len = int((feat.location.end - feat.location.start) / 3)
        candidates.append((len(candidates), aa_len, int(feat.location.end - feat.location.start)))

    chosen_idx = choose_anchor_index(candidates, protein_len)
    chosen_record, chosen_feature_index = matches[chosen_idx]
    cds_list = cds_lists[chosen_record.id]

    anchor_feature = cds_list[chosen_feature_index]
    anchor_index = chosen_feature_index

    left = max(0, anchor_index - flank)
    right = min(len(cds_list) - 1, anchor_index + flank)
    subset = cds_list[left : right + 1]

    region_start = min(int(f.location.start) for f in subset)
    region_end = max(int(f.location.end) for f in subset)
    region_seq = chosen_record.seq[region_start:region_end]

    new_record = SeqRecord(
        region_seq,
        id=f"{chosen_record.id}_{region_start + 1}_{region_end}",
        description=f"{chosen_record.id}:{region_start + 1}-{region_end}",
    )

    new_features: List[SeqFeature] = []
    for feat in subset:
        start = int(feat.location.start) - region_start
        end = int(feat.location.end) - region_start
        location = FeatureLocation(start, end, strand=feat.location.strand)
        qualifiers = {k: list(v) for k, v in feat.qualifiers.items()}
        if feat is anchor_feature:
            qualifiers.setdefault("note", []).append("ANCHOR")
        prot_id = None
        for key in ["protein_id", "ID", "Name", "locus_tag"]:
            if key in qualifiers:
                prot_id = qualifiers[key][0]
                break
        if prot_id and prot_id in protein_seqs:
            qualifiers.setdefault("translation", []).append(str(protein_seqs[prot_id]))
        new_features.append(SeqFeature(location=location, type="CDS", qualifiers=qualifiers))

    return new_record, new_features, {
        "contig": chosen_record.id,
        "start": int(anchor_feature.location.start) + 1,
        "end": int(anchor_feature.location.end),
        "match_count": len(matches),
        "neighbor_count": len(subset),
    }


def extract_neighborhood_from_gff(
    gff_path: str,
    fna_path: str,
    anchor_id: str,
    flank: int,
    protein_len: Optional[int],
    protein_seqs: Dict[str, Seq],
) -> Tuple[SeqRecord, List[SeqFeature], Dict]:
    cds_list = parse_gff_cds(gff_path)
    if not cds_list:
        raise ValueError("No CDS features found in GFF")

    matches = []
    for idx, cds in enumerate(cds_list):
        if match_in_attributes(cds["attrs"], anchor_id):
            matches.append(idx)

    if not matches:
        raise ValueError("Anchor protein ID not found in GFF")

    candidates: List[Tuple[int, int, int]] = []
    for idx in matches:
        cds = cds_list[idx]
        cds_len_aa = int((cds["end"] - cds["start"] + 1) / 3)
        candidates.append((idx, cds_len_aa, cds["end"] - cds["start"]))

    chosen_index = choose_anchor_index(candidates, protein_len)
    anchor_cds = cds_list[chosen_index]

    seqs = build_seq_lookup(fna_path)
    seqid = anchor_cds["seqid"]
    if seqid not in seqs:
        raise ValueError(f"Contig {seqid} not found in FNA")
    contig = seqs[seqid]

    same_contig = [c for c in cds_list if c["seqid"] == seqid]
    same_contig.sort(key=lambda c: (c["start"], c["end"]))
    anchor_index = next(i for i, c in enumerate(same_contig) if c is anchor_cds)

    left = max(0, anchor_index - flank)
    right = min(len(same_contig) - 1, anchor_index + flank)
    subset = same_contig[left : right + 1]

    region_start = min(c["start"] for c in subset)
    region_end = max(c["end"] for c in subset)
    region_seq = contig.seq[region_start - 1 : region_end]

    new_record = SeqRecord(
        region_seq,
        id=f"{seqid}_{region_start}_{region_end}",
        description=f"{seqid}:{region_start}-{region_end}",
    )

    new_features: List[SeqFeature] = []
    for cds in subset:
        start = cds["start"] - region_start
        end = cds["end"] - region_start
        strand = 1 if cds["strand"] == "+" else -1
        location = FeatureLocation(start, end, strand=strand)
        qualifiers = {k: list(v) for k, v in cds["attrs"].items()}
        if cds is anchor_cds:
            qualifiers.setdefault("note", []).append("ANCHOR")
        prot_id = None
        for key in ["protein_id", "ID", "Name", "locus_tag"]:
            if key in qualifiers:
                prot_id = qualifiers[key][0]
                break
        if prot_id and prot_id in protein_seqs:
            qualifiers.setdefault("translation", []).append(str(protein_seqs[prot_id]))
        new_features.append(SeqFeature(location=location, type="CDS", qualifiers=qualifiers))

    return new_record, new_features, {
        "contig": seqid,
        "start": anchor_cds["start"],
        "end": anchor_cds["end"],
        "match_count": len(matches),
        "neighbor_count": len(subset),
    }


def process_record(
    record: InputRecord,
    outdir: str,
    flank: int,
    keep_tmp: bool,
) -> ProcessResult:
    key = f"{record.faa_path}\t{record.anno_path}\t{record.anchor_id}"
    sample_id = safe_sample_id(record.anno_path or record.faa_path or f"line{record.line_no}")

    if not record.faa_path or not record.anno_path or not record.anchor_id:
        return ProcessResult(
            key=key,
            sample_id=sample_id,
            status="fail",
            message=f"Line {record.line_no}: missing fields",
        )
    if not os.path.exists(record.faa_path):
        return ProcessResult(
            key=key,
            sample_id=sample_id,
            status="fail",
            message=f"FAA not found: {record.faa_path}",
        )
    if not os.path.exists(record.anno_path):
        return ProcessResult(
            key=key,
            sample_id=sample_id,
            status="fail",
            message=f"Annotation not found: {record.anno_path}",
        )

    protein_lengths, protein_seqs = load_faa_data(record.faa_path)
    protein_len = protein_lengths.get(record.anchor_id)

    anno_ext = Path(record.anno_path).suffix.lower()
    gbk_out = os.path.join(outdir, "gbk", f"{sample_id}.gbk")
    os.makedirs(os.path.dirname(gbk_out), exist_ok=True)

    try:
        if anno_ext in {".gbk", ".gbff", ".gb"}:
            new_record, new_features, meta = extract_neighborhood_from_gbk(
                record.anno_path,
                record.anchor_id,
                flank,
                protein_len,
                protein_seqs,
            )
        elif anno_ext in {".gff", ".gff3"}:
            fna_path = find_fna_for_gff(record.anno_path)
            if not fna_path:
                return ProcessResult(
                    key=key,
                    sample_id=sample_id,
                    status="fail",
                    message=f"No FNA/FA found for GFF: {record.anno_path}",
                )
            new_record, new_features, meta = extract_neighborhood_from_gff(
                record.anno_path,
                fna_path,
                record.anchor_id,
                flank,
                protein_len,
                protein_seqs,
            )
        else:
            return ProcessResult(
                key=key,
                sample_id=sample_id,
                status="fail",
                message=f"Unsupported annotation format: {record.anno_path}",
            )

        write_neighborhood_gbk(gbk_out, new_record, new_features)
        message = "success"
        if meta["match_count"] > 1:
            if protein_len is None:
                message = "multiple matches; protein length unavailable"
            else:
                message = "multiple matches; picked closest length"

        return ProcessResult(
            key=key,
            sample_id=sample_id,
            status="ok",
            message=message,
            gbk_path=gbk_out,
            anchor_contig=meta["contig"],
            anchor_start=meta["start"],
            anchor_end=meta["end"],
            match_count=meta["match_count"],
            flank=flank,
            neighborhood_genes=meta["neighbor_count"],
        )
    except Exception as exc:
        return ProcessResult(
            key=key,
            sample_id=sample_id,
            status="fail",
            message=str(exc),
        )
    finally:
        if not keep_tmp:
            pass




def _parse_processed_fields(fields: List[str]) -> Optional[Tuple[str, str, str, str, str]]:
    if len(fields) >= 8:
        faa_path = fields[1]
        anno_path = fields[2]
        anchor_id = fields[3]
        status = fields[5]
        gbk_path = fields[7]
        key = "\t".join([faa_path, anno_path, anchor_id])
        return key, status, faa_path, anno_path, gbk_path
    if len(fields) >= 5:
        key = fields[1]
        status = fields[3]
        gbk_path = fields[5] if len(fields) > 5 else ""
        parts = key.split("\t")
        faa_path = parts[0] if len(parts) > 0 else ""
        anno_path = parts[1] if len(parts) > 1 else ""
        return key, status, faa_path, anno_path, gbk_path
    return None


def read_processed(log_path: str) -> Dict[str, str]:
    processed: Dict[str, str] = {}
    if not os.path.exists(log_path):
        return processed
    with open(log_path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            parsed = _parse_processed_fields(fields)
            if not parsed:
                continue
            key, status, _faa, _anno, _gbk = parsed
            processed[key] = status
    return processed


def read_processed_success_gbks(log_path: str) -> List[str]:
    gbks: List[str] = []
    if not os.path.exists(log_path):
        return gbks
    with open(log_path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            parsed = _parse_processed_fields(fields)
            if not parsed:
                continue
            _key, status, _faa, _anno, gbk_path = parsed
            if status == "ok" and gbk_path and os.path.exists(gbk_path):
                gbks.append(gbk_path)
    return gbks


def append_processed(log_path: str, results: List[ProcessResult]) -> None:
    with open(log_path, "a", encoding="utf-8") as handle:
        for res in results:
            timestamp = dt.datetime.now().isoformat(timespec="seconds")
            handle.write(
                "\t".join(
                    [
                        timestamp,
                        res.key,
                        res.sample_id,
                        res.status,
                        res.message.replace("\t", " "),
                        res.gbk_path or "",
                    ]
                )
                + "\n"
            )


def run_clinker(gbk_paths: List[str], outdir: str, logger: logging.Logger) -> None:
    if not gbk_paths:
        logger.info("No GBK files to run clinker.")
        return
    if not shutil.which("clinker"):
        logger.error("clinker not found in PATH; skip running clinker.")
        return
    clinker_dir = os.path.join(outdir, "clinker")
    os.makedirs(clinker_dir, exist_ok=True)
    html_path = os.path.join(clinker_dir, "clinker.html")
    align_path = os.path.join(clinker_dir, "alignments.tsv")
    matrix_path = os.path.join(clinker_dir, "similarity_matrix.tsv")

    cmd = [
        "clinker",
        *gbk_paths,
        "-p",
        html_path,
        "-o",
        align_path,
        "-dl",
        "\t",
        "-dc",
        "4",
        "-mo",
        matrix_path,
        "-f",
    ]
    logger.info("Running: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    logger.info("clinker stdout: %s", result.stdout.strip())
    if result.returncode != 0:
        logger.error("clinker failed: %s", result.stderr.strip())
        return

    json_path = os.path.join(clinker_dir, "clinker.json")
    svg_path = os.path.join(clinker_dir, "clinker.svg")
    if os.path.exists(json_path):
        plot_cmd = ["clinker", "plot", json_path, "-o", svg_path]
        logger.info("Running: %s", " ".join(plot_cmd))
        plot_result = subprocess.run(plot_cmd, capture_output=True, text=True)
        if plot_result.returncode != 0:
            logger.warning("clinker plot failed: %s", plot_result.stderr.strip())


def main() -> int:
    parser = argparse.ArgumentParser(
        description="clinker + clustermap.js pipeline for gene neighborhood visualization"
    )
    parser.add_argument("--conf", required=True, help="Input TSV configuration file")
    parser.add_argument("--outdir", required=True, help="Output directory root")
    parser.add_argument("--logdir", required=True, help="Log directory root")
    parser.add_argument("--flank", type=int, default=None, help="Flanking genes on each side")
    parser.add_argument("--threads", type=int, default=None, help="Parallel threads")
    parser.add_argument("--keep_tmp", action="store_true", help="Keep intermediate files")
    parser.add_argument(
        "--skip_extract",
        action="store_true",
        help="Skip extraction and use existing processed gbk files",
    )

    args = parser.parse_args()

    if args.flank is None or args.threads is None:
        raise SystemExit("Missing --flank or --threads; please set in conf/clinker.conf or CLI.")

    run_id = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    run_outdir = os.path.join(args.outdir, f"run_{run_id}")
    run_logdir = os.path.join(args.logdir, f"run_{run_id}")
    os.makedirs(run_outdir, exist_ok=True)
    os.makedirs(run_logdir, exist_ok=True)

    log_path = os.path.join(run_logdir, "run.log")
    logger = setup_logging(log_path)
    log_versions(logger)

    symlink_log = os.path.join(run_outdir, "run.log")
    if not os.path.exists(symlink_log):
        try:
            os.symlink(log_path, symlink_log)
        except OSError:
            shutil.copy2(log_path, symlink_log)

    processed_log = os.path.join(args.logdir, "processed.tsv")
    processed = read_processed(processed_log)

    records = read_input_tsv(args.conf)
    logger.info("Total records in conf: %d", len(records))
    input_summary_path = os.path.join(run_outdir, "inputs.tsv")
    with open(input_summary_path, "w", encoding="utf-8") as handle:
        handle.write("line\tfaa_path\tanno_path\tanchor_id\n")
        for r in records:
            handle.write(
                f"{r.line_no}\t{r.faa_path}\t{r.anno_path}\t{r.anchor_id}\n"
            )

    todo = []
    skipped: List[ProcessResult] = []
    for r in records:
        key = f"{r.faa_path}\t{r.anno_path}\t{r.anchor_id}"
        if processed.get(key) == "ok":
            skipped.append(
                ProcessResult(
                    key=key,
                    sample_id=safe_sample_id(r.anno_path or r.faa_path or f"line{r.line_no}"),
                    status="skipped",
                    message="already processed",
                )
            )
            continue
        todo.append(r)

    results: List[ProcessResult] = []
    logger.info("Pending records to process: %d", len(todo))
    if not args.skip_extract:
        if args.threads <= 1:
            for r in todo:
                logger.info(
                    "Start: %s",
                    safe_sample_id(r.anno_path or r.faa_path or f"line{r.line_no}"),
                )
                results.append(process_record(r, run_outdir, args.flank, args.keep_tmp))
                logger.info("Done: %s (%s)", results[-1].sample_id, results[-1].status)
        else:
            with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as exc:
                futures = [
                    exc.submit(process_record, r, run_outdir, args.flank, args.keep_tmp)
                    for r in todo
                ]
                logger.info("Submitted %d tasks to process pool", len(futures))
                for fut in concurrent.futures.as_completed(futures):
                    try:
                        res = fut.result()
                    except Exception as exc:
                        res = ProcessResult(
                            key="",
                            sample_id="unknown",
                            status="fail",
                            message=str(exc),
                        )
                    results.append(res)
                    logger.info("Done: %s (%s)", res.sample_id, res.status)
    else:
        logger.info("Skip extraction enabled; using existing processed GBKs")

    all_results = skipped + results
    append_processed(processed_log, results)

    summary_path = os.path.join(run_outdir, "summary.tsv")
    failed_path = os.path.join(run_outdir, "failed.tsv")
    with open(summary_path, "w", encoding="utf-8") as sh, open(
        failed_path, "w", encoding="utf-8"
    ) as fh:
        sh.write(
            "sample_id\tstatus\tanchor_id\tmatch_count\tanchor_contig\tanchor_start\tanchor_end\t"
            "flank\tneighbor_genes\tgbk_path\tmessage\n"
        )
        fh.write("sample_id\tmessage\n")
        for res in all_results:
            sh.write(
                "\t".join(
                    [
                        res.sample_id,
                        res.status,
                        res.key.split("\t")[-1] if res.key else "",
                        str(res.match_count),
                        res.anchor_contig or "",
                        str(res.anchor_start or ""),
                        str(res.anchor_end or ""),
                        str(res.flank),
                        str(res.neighborhood_genes),
                        res.gbk_path or "",
                        res.message,
                    ]
                )
                + "\n"
            )
            if res.status == "fail":
                fh.write(f"{res.sample_id}\t{res.message}\n")

    gbk_paths = read_processed_success_gbks(processed_log)
    gbk_paths.extend([r.gbk_path for r in results if r.status == "ok" and r.gbk_path])
    gbk_paths = sorted(set(gbk_paths))
    run_clinker(gbk_paths, run_outdir, logger)

    logger.info("Run output: %s", run_outdir)
    logger.info("Run logs: %s", run_logdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
