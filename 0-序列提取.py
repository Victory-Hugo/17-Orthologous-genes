#!/usr/bin/env python3
"""
Extract CDS/nucleotide sequences that correspond to a given protein FASTA.

Workflow:
- Try to locate an existing *.fna (or *_CDS.fna) alongside the data. If found,
  pull only the entries whose IDs appear in the input .faa.
- If no CDS file is available but a genome FASTA + GFF are present, rebuild
  the CDS sequences from the annotations.
- When translation/verification is needed, rely on seqkit (per requirement).
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple, Any


def colorize(text: str, level: str = "info") -> str:
    """Add simple ANSI colors for log out."""
    colors = {
        "info": "\033[92m",   # green
        "warn": "\033[93m",   # yellow
        "error": "\033[91m",  # red
    }
    reset = "\033[0m"
    prefix = colors.get(level, "")
    return f"{prefix}{text}{reset}" if prefix else text


def log_line(message: str, log_path: Optional[Path]) -> None:
    """Append plain message to log file if provided."""
    if not log_path:
        return
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as fh:
        fh.write(message.rstrip("\n") + "\n")


def render_progress(current: int, total: int, width: int = 30) -> str:
    """生成简单的文本进度条。"""
    if total <= 0:
        return "[无任务]"
    ratio = min(max(current / total, 0), 1)
    filled = int(ratio * width)
    bar = "█" * filled + "-" * (width - filled)
    percent = int(ratio * 100)
    return f"[{bar}] {percent:3d}% ({current}/{total})"


def read_ids_from_fasta(path: Path) -> List[str]:
    """仅读取FASTA首列ID，节省内存。"""
    ids: List[str] = []
    try:
        with path.open() as fh:
            for line in fh:
                if line.startswith(">"):
                    ids.append(line[1:].strip().split()[0])
    except (FileNotFoundError, PermissionError) as e:
        sys.stderr.write(colorize(f"[error] 无法读取文件 {path}: {e}\n", "error"))
        raise
    return ids


def collect_alias_ids_from_fasta(path: Path) -> List[str]:
    """读取FASTA header，并收集首列ID + protein_id/locus_tag等别名。"""
    ids: List[str] = []
    seen: Set[str] = set()
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                header = line[1:].strip()
                for val in [header.split()[0], *header_aliases(header)]:
                    if val not in seen:
                        seen.add(val)
                        ids.append(val)
    return ids


def extract_ids_from_gff(gff_path: Path) -> List[str]:
    """从GFF/GFF3中提取ID/protein_id/Name/locus_tag列表（支持Prodigal/Prokka/NCBI）。"""
    ids: Set[str] = set()
    with gff_path.open() as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            feature = parts[2]
            # 只处理CDS特征
            if feature not in {"CDS", "gene"}:
                continue
            attrs = parts[8]
            # 解析属性，支持不同的分隔符和格式
            attr_dict = {}
            for chunk in attrs.split(";"):
                chunk = chunk.strip()
                if "=" in chunk:
                    k, v = chunk.split("=", 1)
                    attr_dict[k.strip()] = v.strip()
            
            # 提取各种可能的ID字段
            for key in ["ID", "protein_id", "Name", "locus_tag", "gene", "inference"]:
                if key in attr_dict:
                    value = attr_dict[key]
                    # 处理可能的引号
                    value = value.strip('"\'')
                    ids.add(value)
                    
            # Prodigal特殊处理：ID字段通常是主要标识符
            if "ID" in attr_dict:
                main_id = attr_dict["ID"].strip('"\'')
                ids.add(main_id)
                # Prodigal的ID格式可能是 contig_number，也加入别名
                for alias in expand_aliases(main_id):
                    ids.add(alias)
                    
    return sorted(ids)


def extract_ids_from_gtf(gtf_path: Path) -> List[str]:
    """从GTF中提取常见ID字段（gene_id/transcript_id/protein_id/ID/Name）。"""
    ids: Set[str] = set()
    with gtf_path.open() as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            attrs = parts[8]
            for chunk in attrs.split(";"):
                chunk = chunk.strip()
                if not chunk:
                    continue
                if " " in chunk:
                    k, v = chunk.split(" ", 1)
                    v = v.strip().strip('"')
                elif "=" in chunk:
                    k, v = chunk.split("=", 1)
                else:
                    continue
                if k in {"gene_id", "transcript_id", "protein_id", "ID", "Name", "locus_tag"}:
                    ids.add(v)
    return sorted(ids)


def extract_ids_from_table(
    table_path: Path, column: int = 1, sep: str = "auto"
) -> List[str]:
    """
    从table/csv/tsv中提取指定列的ID（1-based）。
    sep: auto | tab | comma | space
    """
    ids: List[str] = []
    delimiter = None
    if sep == "tab":
        delimiter = "\t"
    elif sep in {"comma", "csv"}:
        delimiter = ","
    elif sep == "space":
        delimiter = None
    with table_path.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("#"):
                continue
            if delimiter is None and sep == "auto":
                if "\t" in line:
                    delimiter = "\t"
                elif "," in line:
                    delimiter = ","
                else:
                    delimiter = None
            parts = line.rstrip("\n").split(delimiter) if delimiter is not None else line.split()
            if len(parts) < column:
                continue
            ids.append(parts[column - 1].strip())
    return ids


def collect_ids(table_path: Path, column: int, sep: str) -> List[str]:
    ext = table_path.suffix.lower()
    if ext in {".gff", ".gff3"}:
        return extract_ids_from_gff(table_path)
    if ext == ".gtf":
        return extract_ids_from_gtf(table_path)
    return extract_ids_from_table(table_path, column, sep)


def extract_sequences_by_ids(
    ids: List[str], seq_path: Path, output_path: Path
) -> int:
    """
    通用ID提取：优先用seqkit grep -f id.list，否则回退Python。
    """
    # 验证输入路径
    if not seq_path.exists():
        raise FileNotFoundError(f"输入序列文件不存在: {seq_path}")
    if not seq_path.is_file():
        raise ValueError(f"路径不是文件: {seq_path}")
    
    seqkit = shutil.which("seqkit")
    if seqkit:
        ids_file = None
        try:
            with tempfile.NamedTemporaryFile("w", delete=False) as tmp_ids:
                tmp_ids.write("\n".join(ids) + "\n")
                ids_file = tmp_ids.name
            
            output_path.parent.mkdir(parents=True, exist_ok=True)
            with output_path.open("w") as out_handle:
                proc = subprocess.run(
                    [seqkit, "grep", "-f", ids_file, str(seq_path)],
                    stdout=out_handle,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=False,
                )
            if proc.returncode == 0:
                with output_path.open() as fh:
                    return sum(1 for line in fh if line.startswith(">"))
            else:
                sys.stderr.write(colorize(f"[warn] seqkit grep 失败，回退到Python解析: {proc.stderr}\n", "warn"))
        except Exception as e:
            sys.stderr.write(colorize(f"[warn] seqkit 执行异常: {e}，回退到Python解析\n", "warn"))
        finally:
            if ids_file and os.path.exists(ids_file):
                os.unlink(ids_file)

    # 流式回退：逐行扫描FASTA，避免一次性载入内存。
    id_set = set(ids)
    source_tag = seq_path.stem
    count = 0
    with seq_path.open() as inp, output_path.open("w") as out:
        keep = False
        for line in inp:
            if line.startswith(">"):
                seq_id = line[1:].strip().split()[0]
                keep = seq_id in id_set
                if keep:
                    out.write(f">{seq_id}|src={source_tag}\n")
                    count += 1
            else:
                if keep:
                    out.write(line)
    return count


class ChineseArgumentParser(argparse.ArgumentParser):
    """Argparse子类，用于输出中文错误信息。"""

    def error(self, message: str) -> None:  # type: ignore[override]
        # 翻译常见的错误提示
        translated = message
        if "the following arguments are required:" in message:
            translated = translated.replace(
                "the following arguments are required:", "缺少必填参数:"
            )
        self.print_usage(sys.stderr)
        self.exit(2, colorize(f"参数错误: {translated}\n", "error"))


CODON_TABLE_NAMES = {
    "1": "标准核基因编码表 (Standard)",
    "2": "脊椎动物线粒体 (Vertebrate Mitochondrial)",
    "3": "酵母线粒体 (Yeast Mitochondrial)",
    "4": "霉菌/原生动物/共生体线粒体 (Mold/Protozoan/Coelenterate Mitochondrial)",
    "5": "无脊椎动物线粒体 (Invertebrate Mitochondrial)",
    "6": "纤毛虫/草履虫核 (Ciliate/Dasycladacean/Stylonichia Nuclear)",
    "9": "棘皮动物/扁形动物线粒体 (Echinoderm/Flatworm Mitochondrial)",
    "10": "纤毛虫Euplotid核 (Euplotid Nuclear)",
    "11": "细菌/古菌/植物质体 (Bacterial, Archaeal and Plant Plastid)",
}


def read_fasta(path: Path) -> Dict[str, str]:
    """Load sequences into a dict keyed by the first token of the header."""
    seqs: Dict[str, List[str]] = {}
    header = None
    with path.open() as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:].strip().split()[0]
                seqs[header] = []
            else:
                if header is None:
                    continue
                seqs[header].append(line.strip())
    return {h: "".join(parts) for h, parts in seqs.items()}


def expand_aliases(value: str) -> Set[str]:
    """Generate possible aliases from an identifier (handles various tool formats)."""
    aliases: Set[str] = set()
    aliases.add(value)
    
    # 处理管道符分隔的ID（NCBI格式）
    parts = [p for p in value.split("|") if p]
    aliases.update(parts)
    
    # 处理版本号（如WP_000000001.1 -> WP_000000001）
    for v in list(aliases):
        if "." in v:
            aliases.add(v.split(".", 1)[0])
    
    # 处理Prodigal的下划线格式（如1_1, 2_3）
    for v in list(aliases):
        if "_" in v and v.replace("_", "").isdigit():
            # 保留原始格式和不同的变体
            parts = v.split("_")
            if len(parts) == 2 and all(p.isdigit() for p in parts):
                aliases.add(f"{parts[0]}_{parts[1]}")
                aliases.add(f"{parts[0]}.{parts[1]}")
    
    # 清理空字符串和特殊字符
    return {a.strip() for a in aliases if a.strip()}


def read_ids_from_faa(path: Path) -> Tuple[List[str], Dict[str, str]]:
    """
    Return IDs in order plus an alias map (alias -> canonical FAA ID).
    The alias map allows matching against protein_id/locus_tag in CDS headers.
    """
    ids: List[str] = []
    alias_map: Dict[str, str] = {}
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                rid = line[1:].strip().split()[0]
                ids.append(rid)
                for alias in expand_aliases(rid):
                    alias_map.setdefault(alias, rid)
    return ids, alias_map


def write_fasta(
    records: Iterable[Tuple[str, str, Optional[str]]], out_path: Path, source_tag: Optional[str] = None
) -> int:
    count = 0
    with out_path.open("w") as out:
        for rid, seq, alias in records:
            if not seq:
                continue
            header_parts = [rid]
            if alias and alias != rid:
                header_parts.append(f"alias={alias}")
            if source_tag:
                header_parts.append(f"src={source_tag}")
            header = "|".join(header_parts)
            out.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")
            count += 1
    return count


def header_aliases(header: str) -> Set[str]:
    """Extract possible IDs from FASTA headers (NCBI, Prodigal, Prokka compatible)."""
    aliases = set()
    token = header.split()[0]
    aliases.update(expand_aliases(token))
    
    # NCBI/Prokka格式: [protein_id=xxx] [locus_tag=xxx] [gene=xxx]
    for match in re.finditer(r"\[(protein_id|locus_tag|gene)=([^\]]+)\]", header):
        aliases.update(expand_aliases(match.group(2)))
    
    # Prodigal/其他格式: protein_id=xxx;locus_tag=xxx
    for match in re.finditer(r"(protein_id|locus_tag|gene|ID)=([^;\s\]]+)", header):
        aliases.update(expand_aliases(match.group(2)))
    
    # Prokka格式: >ID_name description [gene=xxx]
    # 如果header中有空格分隔的描述，第一个token通常是主要ID
    parts = header.split()
    if len(parts) > 1:
        # 检查是否有基因名等信息
        for part in parts[1:]:
            if '=' not in part and len(part) > 2:  # 可能是基因名
                aliases.add(part)
    
    return {a for a in aliases if a}


def extract_from_cds_fna(
    ids: List[str], alias_map: Dict[str, str], fna_path: Path, out_path: Path
) -> int:
    """
    Extract matching IDs from an existing CDS .fna.

    Instead of relying purely on FASTA IDs, we parse protein_id/locus_tag entries
    in the header to match the provided FAA.
    """
    source_path = fna_path
    seqkit = shutil.which("seqkit")
    tmp_filtered: Optional[Path] = None
    if seqkit:
        with tempfile.NamedTemporaryFile("w", delete=False) as tmp_ids:
            tmp_ids.write("\n".join(sorted(set(alias_map.keys()))) + "\n")
            ids_file = tmp_ids.name
        tmp_filtered = Path(tempfile.mkstemp(suffix=".fna")[1])
        try:
            with tmp_filtered.open("w") as out_handle:
                proc = subprocess.run(
                    [seqkit, "grep", "-f", ids_file, str(fna_path)],
                    stdout=out_handle,
                    stderr=subprocess.PIPE,
                    text=True,
                )
            if proc.returncode == 0:
                source_path = tmp_filtered
            else:
                sys.stderr.write(colorize(f"[warn] seqkit grep 失败，回退到Python解析: {proc.stderr}\n", "warn"))
        finally:
            os.unlink(ids_file)

    def scan_file(path: Path) -> Dict[str, Tuple[str, Optional[str]]]:
        found: Dict[str, Tuple[str, Optional[str]]] = {}
        with path.open() as fh:
            header = None
            seq_parts: List[str] = []
            for line in fh:
                if line.startswith(">"):
                    if header:
                        aliases = header_aliases(header)
                        canonical = next((alias_map[a] for a in aliases if a in alias_map), None)
                        if canonical and canonical not in found:
                            matched_alias = next((a for a in aliases if a in alias_map), None)
                            found[canonical] = ("".join(seq_parts), matched_alias)
                    header = line[1:].strip()
                    seq_parts = []
                else:
                    seq_parts.append(line.strip())
            if header:
                aliases = header_aliases(header)
                canonical = next((alias_map[a] for a in aliases if a in alias_map), None)
                if canonical and canonical not in found:
                    matched_alias = next((a for a in aliases if a in alias_map), None)
                    found[canonical] = ("".join(seq_parts), matched_alias)
        return found

    matched = scan_file(source_path)
    if not matched and tmp_filtered:
        matched = scan_file(fna_path)

    if tmp_filtered:
        tmp_filtered.unlink(missing_ok=True)

    source_tag = fna_path.stem
    records = ((rid, matched[rid][0], matched[rid][1]) for rid in ids if rid in matched)
    return write_fasta(records, out_path, source_tag)


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(table)[::-1]


def parse_gff_cds(
    gff_path: Path, alias_map: Dict[str, str]
) -> Dict[str, List[Tuple[str, int, int, str]]]:
    """
    Return mapping canonical FAA ID -> list of (seqid, start, end, strand) segments.
    Start/end are 0-based, end-exclusive.
    """
    target_aliases = set(alias_map.keys())
    segments: Dict[str, List[Tuple[str, int, int, str]]] = defaultdict(list)
    with gff_path.open() as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seqid, _source, feature, start, end, _score, strand, _phase, attrs = parts
            if feature != "CDS":
                continue
            attr_map = {}
            for chunk in attrs.split(";"):
                if "=" in chunk:
                    k, v = chunk.split("=", 1)
                    attr_map[k] = v
            # 按优先级提取protein_id，兼容不同工具
            protein_id = (
                attr_map.get("protein_id")  # NCBI/Prokka优先
                or attr_map.get("locus_tag")  # 通用locus_tag
                or attr_map.get("ID")  # Prodigal/GFF3标准
                or attr_map.get("Name")  # 备选
                or attr_map.get("gene")  # Prokka基因名
            )
            # 清理可能的引号
            if protein_id:
                protein_id = protein_id.strip('"\'')
            if not protein_id:
                continue
            canonical = None
            for alias in expand_aliases(protein_id):
                if alias in target_aliases:
                    canonical = alias_map[alias]
                    break
            if canonical is None:
                continue
            try:
                s = int(start) - 1
                e = int(end)
            except ValueError:
                continue
            segments[canonical].append((seqid, s, e, strand))
    return segments


def rebuild_cds_with_gff(
    ids: List[str],
    alias_map: Dict[str, str],
    genome_fasta: Path,
    gff_path: Path,
    out_path: Path,
) -> int:
    """备用：纯Python根据GFF重建（当gffread不可用时使用）。"""
    genome = read_fasta(genome_fasta)
    cds_segments = parse_gff_cds(gff_path, alias_map)

    records = []
    for rid in ids:
        if rid not in cds_segments:
            continue
        segs = cds_segments[rid]
        if not segs:
            continue
        # Assume all segments share the same strand.
        strand = segs[0][3]
        # Sort by start to ensure correct order before reverse complement if needed.
        segs_sorted = sorted(segs, key=lambda x: x[1])
        seq_parts: List[str] = []
        for seqid, start, end, _ in segs_sorted:
            contig = genome.get(seqid)
            if contig is None:
                continue
            seq_parts.append(contig[start:end])
        cds_seq = "".join(seq_parts)
        if strand == "-":
            cds_seq = revcomp(cds_seq)
        records.append((rid, cds_seq, None))

    return write_fasta(records, out_path, genome_fasta.stem)


def rebuild_cds_with_gffread(
    ids: List[str],
    alias_map: Dict[str, str],
    genome_fasta: Path,
    gff_path: Path,
    out_path: Path,
) -> int:
    """
    使用gffread生成CDS再抽取；自动过滤伪基因行和非编码特征。
    支持Prodigal、Prokka和NCBI格式。
    """
    gffread = shutil.which("gffread")
    if not gffread:
        return 0

    # 过滤伪基因和非编码特征到临时GFF
    filtered_gff = Path(tempfile.mkstemp(suffix=".gff")[1])
    with gff_path.open() as src, filtered_gff.open("w") as dst:
        for line in src:
            # 保留注释行
            if line.startswith("#"):
                dst.write(line)
                continue
            # 过滤伪基因和非CDS特征
            line_lower = line.lower()
            if ("pseudo" in line_lower or 
                "pseudogene" in line_lower or
                ("\t" in line and line.split("\t")[2] not in {"CDS", "gene"})):
                continue
            dst.write(line)

    tmp_cds = Path(tempfile.mkstemp(suffix=".fna")[1])
    try:
        proc = subprocess.run(
            [gffread, "-g", str(genome_fasta), "-x", str(tmp_cds), str(filtered_gff)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if proc.returncode != 0:
            sys.stderr.write(colorize(f"[warn] gffread 生成CDS失败: {proc.stderr}\n", "warn"))
            return 0
        return extract_from_cds_fna(ids, alias_map, tmp_cds, out_path)
    finally:
        filtered_gff.unlink(missing_ok=True)
        tmp_cds.unlink(missing_ok=True)
def translate_with_seqkit(
    fna_path: Path, table: str = "11", output_path: Optional[Path] = None
) -> Path | None:
    seqkit = shutil.which("seqkit")
    if not seqkit:
        sys.stderr.write(colorize("[warn] 未找到seqkit，跳过翻译。\n", "warn"))
        return None
    tmp_created = False
    if output_path is None:
        output_path = Path(tempfile.mkstemp(suffix=".faa")[1])
        tmp_created = True
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [seqkit, "translate", "-T", table, "-o", str(output_path), str(fna_path)]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        sys.stderr.write(colorize(f"[warn] seqkit translate 失败: {proc.stderr}\n", "warn"))
        if tmp_created:
            output_path.unlink(missing_ok=True)
        return None
    return output_path


def check_consistency(
    translated_faa: Path, original_faa: Path
) -> Tuple[int, int, int]:
    """
    比对翻译后的CDS与原始FAA，返回 (匹配数, 长度不一致, 序列不一致)。
    """
    trans = read_fasta(translated_faa)
    orig = read_fasta(original_faa)
    match = length_diff = seq_diff = 0
    for rid, t_seq in trans.items():
        o_seq = orig.get(rid)
        if o_seq is None:
            continue
        if len(t_seq) != len(o_seq):
            length_diff += 1
            continue
        if t_seq != o_seq:
            seq_diff += 1
        else:
            match += 1
    return match, length_diff, seq_diff


def parse_faa_map_file(map_path: Path) -> List[Dict[str, Any]]:
    """
    解析FAA批量map文件：三列为 faa路径、data目录、out目录。
    """
    entries: List[Dict[str, Any]] = []
    base_dir = map_path.parent
    with map_path.open() as fh:
        for line_no, line in enumerate(fh, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if len(parts) < 3:
                sys.stderr.write(
                    colorize(f"[warn] map文件第{line_no}行列数不足(需3列)，已跳过: {stripped}\n", "warn")
                )
                continue
            col1, col2, col3 = parts[:3]
            faa_path = (base_dir / col1).resolve()
            data_path = (base_dir / col2).resolve()
            out_path = (base_dir / col3).resolve()
            entries.append({
                "faa": faa_path,
                "data": data_path,
                "out": out_path
            })
    return entries


def parse_fna_map_file(map_path: Path) -> List[Dict[str, Any]]:
    """
    解析FNA批量map文件：
    - 三列格式：fna路径、out目录、outfaa（可为文件或目录）
    - 四列格式：fna路径、data目录、out目录、outfaa（可为文件或目录）
    """
    entries: List[Dict[str, Any]] = []
    base_dir = map_path.parent
    with map_path.open() as fh:
        for line_no, line in enumerate(fh, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if len(parts) < 3:
                sys.stderr.write(
                    colorize(f"[warn] map文件第{line_no}行列数不足(需3列)，已跳过: {stripped}\n", "warn")
                )
                continue
            
            col1 = parts[0]
            input_path = (base_dir / col1).resolve()
            
            data_path = None
            out_dir = None
            target = None
            
            if len(parts) >= 4:
                # 四列格式：fna data out outfaa
                data_path = (base_dir / parts[1]).resolve()
                out_dir = (base_dir / parts[2]).resolve()
                target = (base_dir / parts[3]).resolve()
            else:
                # 三列格式：fna out outfaa
                out_dir = (base_dir / parts[1]).resolve()
                target = (base_dir / parts[2]).resolve()
                
            out_dir.mkdir(parents=True, exist_ok=True)
            if target.suffix.lower() == ".faa":
                out_faa = target
            else:
                out_faa = target if target.is_dir() else out_dir / target.name
                if out_faa.is_dir():
                    out_faa = out_faa / f"{input_path.stem}.faa"
                    
            entries.append({
                "input": input_path,
                "outdir": out_dir,
                "outfaa": out_faa,
                "data": data_path,
            })
    return entries


def process_faa_entry(args: Tuple[Dict[str, Any], str, bool, Optional[Path]]) -> Tuple[bool, str, str]:
    """
    并行处理单个FAA条目的worker函数
    """
    entry, translation_table, check_consistency_flag, log_path = args
    try:
        run(
            entry["faa"],
            entry["data"],
            entry["out"],
            translation_table,
            check_consistency_flag,
            log_path,
        )
        return True, str(entry["faa"]), ""
    except KeyboardInterrupt:
        return False, str(entry["faa"]), "用户中断"
    except Exception as exc:
        import traceback
        error_detail = f"{exc}\n{traceback.format_exc()}"
        return False, str(entry["faa"]), error_detail


def process_fna_entry(args: Tuple[Dict[str, Any], str, Optional[Path]]) -> Tuple[bool, str, str]:
    """
    并行处理单个FNA条目的worker函数
    """
    entry, translation_table, log_path = args
    try:
        run_fna_to_faa(
            entry["input"],
            entry["outfaa"],
            entry.get("data"),
            translation_table,
            log_path,
        )
        return True, str(entry["input"]), ""
    except Exception as exc:
        return False, str(entry["input"]), str(exc)


def run_batch_faa_to_fna(
    map_path: Path,
    translation_table: str = "11",
    check_consistency_flag: bool = False,
    log_path: Optional[Path] = None,
    processes: int = 1,
) -> None:
    """批量执行FAA到FNA转换（支持并行处理）。"""
    entries = parse_faa_map_file(map_path)
    if not entries:
        raise ValueError("map文件中没有有效的任务行。")
    
    start_msg = f"[info] 开始批量FAA->FNA处理: {len(entries)} 个任务，使用 {processes} 个进程"
    sys.stderr.write(colorize(start_msg + "\n", "info"))
    log_line(start_msg, log_path)
    
    if processes == 1:
        # 串行处理
        successes = []
        failures = []
        for idx, entry in enumerate(entries, 1):
            prog = render_progress(idx - 1, len(entries))
            sys.stderr.write(colorize(f"\r{prog} 处理: {entry['faa'].name}", "info"))
            sys.stderr.flush()
            
            success, name, error = process_faa_entry((entry, translation_table, check_consistency_flag, log_path))
            if success:
                successes.append(name)
            else:
                failures.append(f"{name} -> {error}")
    else:
        # 并行处理
        args_list = [(entry, translation_table, check_consistency_flag, log_path) for entry in entries]
        
        with Pool(processes=processes) as pool:
            results = []
            for i, result in enumerate(pool.imap(process_faa_entry, args_list)):
                prog = render_progress(i + 1, len(entries))
                sys.stderr.write(colorize(f"\r{prog}", "info"))
                sys.stderr.flush()
                results.append(result)
        
        successes = [name for success, name, _ in results if success]
        failures = [f"{name} -> {error}" for success, name, error in results if not success]
    
    sys.stderr.write("\n")
    summary = f"[info] FAA->FNA批量完成: 成功 {len(successes)} 条, 失败 {len(failures)} 条"
    sys.stderr.write(colorize(summary + "\n", "info"))
    log_line(summary, log_path)
    
    if failures:
        for item in failures:
            warn_msg = f"[warn] 失败记录: {item}"
            sys.stderr.write(colorize(warn_msg + "\n", "warn"))
            log_line(warn_msg, log_path)


def run_batch_fna_to_faa(
    map_path: Path,
    translation_table: str = "11",
    log_path: Optional[Path] = None,
    processes: int = 1,
) -> None:
    """批量执行FNA到FAA转换（支持并行处理）。"""
    entries = parse_fna_map_file(map_path)
    if not entries:
        raise ValueError("map文件中没有有效的任务行。")
    
    start_msg = f"[info] 开始批量FNA->FAA处理: {len(entries)} 个任务，使用 {processes} 个进程"
    sys.stderr.write(colorize(start_msg + "\n", "info"))
    log_line(start_msg, log_path)
    
    if processes == 1:
        # 串行处理
        successes = []
        failures = []
        for idx, entry in enumerate(entries, 1):
            prog = render_progress(idx - 1, len(entries))
            sys.stderr.write(colorize(f"\r{prog} 处理: {entry['input'].name}", "info"))
            sys.stderr.flush()
            
            success, name, error = process_fna_entry((entry, translation_table, log_path))
            if success:
                successes.append(name)
            else:
                failures.append(f"{name} -> {error}")
    else:
        # 并行处理
        args_list = [(entry, translation_table, log_path) for entry in entries]
        
        with Pool(processes=processes) as pool:
            results = []
            for i, result in enumerate(pool.imap(process_fna_entry, args_list)):
                prog = render_progress(i + 1, len(entries))
                sys.stderr.write(colorize(f"\r{prog}", "info"))
                sys.stderr.flush()
                results.append(result)
        
        successes = [name for success, name, _ in results if success]
        failures = [f"{name} -> {error}" for success, name, error in results if not success]
    
    sys.stderr.write("\n")
    summary = f"[info] FNA->FAA批量完成: 成功 {len(successes)} 条, 失败 {len(failures)} 条"
    sys.stderr.write(colorize(summary + "\n", "info"))
    log_line(summary, log_path)
    
    if failures:
        for item in failures:
            warn_msg = f"[warn] 失败记录: {item}"
            sys.stderr.write(colorize(warn_msg + "\n", "warn"))
            log_line(warn_msg, log_path)


def run(
    faa_path: Path,
    data_dir: Path,
    out_dir: Path,
    translation_table: str = "11",
    check_consistency_flag: bool = False,
    log_path: Optional[Path] = None,
) -> Tuple[Path, int, Optional[Path]]:
    """
    Extract nucleotide CDS sequences matching the given protein FASTA.

    Returns:
        out_path: Path to the resulting .fna file.
        extracted: Number of sequences written.
        translated_path: Path to translated FAA (if generated), else None.
    """
    start_msg = f"[info] 开始处理蛋白文件: {faa_path}"
    sys.stderr.write(colorize(start_msg + "\n", "info"))
    log_line(start_msg, log_path)
    out_dir.mkdir(parents=True, exist_ok=True)

    if not faa_path.exists():
        raise FileNotFoundError(f"未找到输入FAA文件: {faa_path}")

    ids, alias_map = read_ids_from_faa(faa_path)
    if not ids:
        raise ValueError("在输入FAA中未找到任何序列ID。")

    base = faa_path.stem
    out_path = out_dir / f"{base}.fna"

    candidates = [
        data_dir / f"{base}_CDS.fna",
        data_dir / f"{base}_CDS.fasta",
        data_dir / f"{base}.fna",
        data_dir / f"{base}.fa",
    ]

    extracted = 0
    for cand in candidates:
        if cand.exists():
            msg = f"[info] 使用已存在的核酸文件: {cand}"
            sys.stderr.write(colorize(msg + "\n", "info"))
            log_line(msg, log_path)
            extracted = extract_from_cds_fna(ids, alias_map, cand, out_path)
            if extracted:
                break

    if extracted == 0:
        genome_fa = data_dir / f"{base}.fasta"
        gff = data_dir / f"{base}.gff"
        if genome_fa.exists() and gff.exists():
            msg = "[info] 未找到CDS文件，改用基因组FASTA+GFF重建CDS。"
            sys.stderr.write(colorize(msg + "\n", "info"))
            log_line(msg, log_path)
            extracted = rebuild_cds_with_gffread(ids, alias_map, genome_fa, gff, out_path)
            if extracted == 0:
                extracted = rebuild_cds_with_gff(ids, alias_map, genome_fa, gff, out_path)
        else:
            msg = "[warn] 无可用CDS或基因组/GFF文件，无法重建序列。"
            sys.stderr.write(colorize(msg + "\n", "warn"))
            log_line(msg, log_path)

    if extracted == 0:
        raise RuntimeError("未能提取任何序列，请检查输入文件。")

    ok_msg = f"[info] 成功提取 {extracted} 条序列到: {out_path}"
    sys.stderr.write(colorize(ok_msg + "\n", "info"))
    log_line(ok_msg, log_path)

    translated = translate_with_seqkit(out_path, translation_table)
    if translated:
        trans_msg = f"[info] 已用seqkit翻译，结果(仅供核查): {translated}"
        sys.stderr.write(colorize(trans_msg + "\n", "info"))
        log_line(trans_msg, log_path)
        if check_consistency_flag:
            match, length_diff, seq_diff = check_consistency(translated, faa_path)
            cons_msg = (
                f"[info] 一致性校验: 完全一致 {match} 条, 长度不一致 {length_diff} 条, 序列不一致 {seq_diff} 条"
            )
            sys.stderr.write(colorize(cons_msg + "\n", "info"))
            log_line(cons_msg, log_path)

    return out_path, extracted, translated


def run_fna_to_faa(
    fna_path: Path,
    out_faa: Path,
    data_dir: Optional[Path] = None,
    translation_table: str = "11",
    log_path: Optional[Path] = None,
) -> Path:
    """从FNA生成FAA：优先使用data目录中已有FAA文件按ID抽取，缺失则翻译。"""
    start_msg = f"[info] 开始处理核酸文件: {fna_path}"
    sys.stderr.write(colorize(start_msg + "\n", "info"))
    log_line(start_msg, log_path)
    if not fna_path.exists():
        raise FileNotFoundError(f"未找到输入FNA文件: {fna_path}")
    out_faa.parent.mkdir(parents=True, exist_ok=True)
    ids = collect_alias_ids_from_fasta(fna_path)

    extracted = 0
    if data_dir:
        msg = f"[info] 在data目录中尝试按ID抽取FAA: {data_dir}"
        sys.stderr.write(colorize(msg + "\n", "info"))
        log_line(msg, log_path)
        candidates = sorted(
            list(data_dir.glob("*.faa"))
            + list(data_dir.glob("*.fa"))
            + list(data_dir.glob("*.fasta"))
        )
        for cand in candidates:
            extracted = extract_sequences_by_ids(ids, cand, out_faa)
            if extracted:
                ok_msg = f"[info] 已从 {cand} 抽取 {extracted} 条序列到 {out_faa}"
                sys.stderr.write(colorize(ok_msg + "\n", "info"))
                log_line(ok_msg, log_path)
                break

    if extracted == 0:
        warn_msg = "[warn] data目录未找到可抽取的FAA，改为翻译。"
        sys.stderr.write(colorize(warn_msg + "\n", "warn"))
        log_line(warn_msg, log_path)
        translated = translate_with_seqkit(fna_path, translation_table, out_faa)
        if translated is None or not translated.exists():
            raise RuntimeError("翻译失败，未生成FAA文件。")
        ok_msg = f"[info] 翻译完成，输出: {translated}"
        sys.stderr.write(colorize(ok_msg + "\n", "info"))
        log_line(ok_msg, log_path)
        return translated

    ok_msg = f"[info] 从已有FAA提取完成，输出: {out_faa}，序列数: {extracted}"
    sys.stderr.write(colorize(ok_msg + "\n", "info"))
    log_line(ok_msg, log_path)
    return out_faa


def run_table_extract(
    table_path: Path,
    seq_path: Path,
    out_path: Path,
    column: int = 1,
    sep: str = "auto",
    log_path: Optional[Path] = None,
) -> int:
    """
    根据table/gff/csv/tsv中的ID列表，从给定FASTA中提取序列（FNA或FAA均可）。
    """
    start_msg = f"[info] 读取ID列表: {table_path} (列{column}, 分隔符={sep})"
    sys.stderr.write(colorize(start_msg + "\n", "info"))
    log_line(start_msg, log_path)
    ids = collect_ids(table_path, column, sep)
    if not ids:
        raise ValueError("未能从表格中提取到任何ID。")
    seq_msg = f"[info] 从 {seq_path} 按ID提取序列 -> {out_path}"
    sys.stderr.write(colorize(seq_msg + "\n", "info"))
    log_line(seq_msg, log_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    count = extract_sequences_by_ids(ids, seq_path, out_path)
    done_msg = f"[info] 提取完成，共 {count} 条"
    sys.stderr.write(colorize(done_msg + "\n", "info"))
    log_line(done_msg, log_path)
    return count


def build_parser() -> argparse.ArgumentParser:
    epilog_lines = [
        "子命令：",
        "  faa-to-fna       FAA抽取/重建对应CDS (FNA)，可选一致性校验",
        "  fna-to-faa       FNA生成FAA（优先用data中已有FAA提取，缺失则翻译）", 
        "  subset           表格/CSV/TSV/GFF/GTF + FASTA 按ID提取序列",
        "  batch-faa-to-fna 批量FAA->FNA转换（map格式：faa data out）",
        "  batch-fna-to-faa 批量FNA->FAA转换（map格式：fna [data] out outfaa）",
        "",
        "常见遗传密码表编号对应：",
    ]
    for num, name in sorted(CODON_TABLE_NAMES.items(), key=lambda kv: int(kv[0])):
        epilog_lines.append(f"  {num}: {name}")

    parser = ChineseArgumentParser(
        description="子命令式：faa-to-fna | fna-to-faa | subset | batch-faa-to-fna | batch-fna-to-faa",
        epilog="\n".join(epilog_lines),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    subparsers = parser.add_subparsers(dest="cmd", required=True)

    # FAA -> FNA
    p_faa = subparsers.add_parser("faa-to-fna", help="FAA抽取/重建对应CDS(FNA)")
    p_faa.add_argument("--faa", type=Path, required=True, help="输入FAA")
    p_faa.add_argument("--data", type=Path, required=True, help="包含CDS.fna或genome+gff的目录")
    p_faa.add_argument("--outdir", type=Path, required=True, help="输出目录")
    p_faa.add_argument("--translation-table", default="11", help="遗传密码表编号")
    p_faa.add_argument("--check-consistency", action="store_true", help="翻译结果与原FAA一致性检查")
    p_faa.add_argument("--log", type=Path, help="日志文件")

    # FNA -> FAA
    p_fna = subparsers.add_parser("fna-to-faa", help="FNA生成FAA（先查已有FAA，再翻译）")
    p_fna.add_argument("--fna", type=Path, required=True, help="输入FNA")
    p_fna.add_argument("--data", type=Path, required=True, help="同名FAA所在目录，用于优先提取")
    p_fna.add_argument("--outdir", type=Path, required=True, help="输出目录")
    p_fna.add_argument("--outfaa", type=Path, help="输出FAA路径，默认 outdir/<name>.faa")
    p_fna.add_argument("--translation-table", default="11", help="遗传密码表编号")
    p_fna.add_argument("--log", type=Path, help="日志文件")

    # Subset
    p_sub = subparsers.add_parser("subset", help="表格/GFF/GTF+FASTA按ID提取")
    p_sub.add_argument("--table", type=Path, required=True, help="ID表格/GFF/GTF/CSV/TSV")
    p_sub.add_argument("--seq", type=Path, required=True, help="输入FASTA（FNA或FAA）")
    p_sub.add_argument("--table-col", type=int, default=1, help="ID列(1-based)")
    p_sub.add_argument(
        "--table-sep",
        choices=["auto", "tab", "comma", "space"],
        default="auto",
        help="表格分隔符",
    )
    p_sub.add_argument("--seq-out", type=Path, help="输出文件，未给则默认 outdir/<name>.subset.*")
    p_sub.add_argument("--outdir", type=Path, help="输出目录（用于默认输出路径）")
    p_sub.add_argument("--log", type=Path, help="日志文件")

    # Batch FAA -> FNA
    p_batch_faa = subparsers.add_parser(
        "batch-faa-to-fna",
        help="批量FAA->FNA转换（map格式：faa data out）",
    )
    p_batch_faa.add_argument("--map", type=Path, required=True, help="map文件")
    p_batch_faa.add_argument("--translation-table", default="11", help="遗传密码表编号")
    p_batch_faa.add_argument("--check-consistency", action="store_true", help="一致性检查")
    p_batch_faa.add_argument("--log", type=Path, help="日志文件")
    p_batch_faa.add_argument(
        "--processes", 
        type=int, 
        default=1, 
        help=f"并行进程数（默认1，最大{cpu_count()}）"
    )

    # Batch FNA -> FAA
    p_batch_fna = subparsers.add_parser(
        "batch-fna-to-faa",
        help="批量FNA->FAA转换（map格式：fna [data] out outfaa）",
    )
    p_batch_fna.add_argument("--map", type=Path, required=True, help="map文件")
    p_batch_fna.add_argument("--translation-table", default="11", help="遗传密码表编号")
    p_batch_fna.add_argument("--log", type=Path, help="日志文件")
    p_batch_fna.add_argument(
        "--processes", 
        type=int, 
        default=1, 
        help=f"并行进程数（默认1，最大{cpu_count()}）"
    )

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    try:
        if args.cmd == "faa-to-fna":
            run(
                args.faa,
                args.data,
                args.outdir,
                args.translation_table,
                args.check_consistency,
                args.log,
            )
        elif args.cmd == "fna-to-faa":
            out_faa = args.outfaa or (args.outdir / f"{args.fna.stem}.faa")
            run_fna_to_faa(args.fna, out_faa, args.data, args.translation_table, args.log)
        elif args.cmd == "subset":
            if args.seq_out:
                out_path = args.seq_out
            else:
                if args.outdir is None:
                    parser.error("subset模式需提供 --seq-out 或 --outdir")
                seq_suffix = Path(args.seq).suffix or ".fasta"
                out_path = Path(args.outdir) / f"{Path(args.seq).stem}.subset{seq_suffix}"
            run_table_extract(
                args.table,
                args.seq,
                out_path,
                args.table_col,
                args.table_sep,
                args.log,
            )
        elif args.cmd == "batch-faa-to-fna":
            processes = min(max(args.processes, 1), cpu_count())
            run_batch_faa_to_fna(
                args.map, 
                args.translation_table, 
                args.check_consistency, 
                args.log,
                processes
            )
        elif args.cmd == "batch-fna-to-faa":
            processes = min(max(args.processes, 1), cpu_count())
            run_batch_fna_to_faa(
                args.map, 
                args.translation_table, 
                args.log,
                processes
            )
    except Exception as exc:  # pylint: disable=broad-except
        sys.exit(colorize(str(exc), "error"))


if __name__ == "__main__":
    main()
