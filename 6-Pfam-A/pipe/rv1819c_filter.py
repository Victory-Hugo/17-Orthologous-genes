#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import os
import re
import sys


NBD_DOMAINS = {"ABC_tran", "ABC_N", "ABC_ATPase"} #? NBD 域集合
TMD_DOMAINS = {
    "ABC2_membrane",
    "ABC2_membrane_2",
    "ABC2_membrane_3",
    "ABC2_membrane_4",
    "ABC2_membrane_5",
    "ABC2_membrane_6",
    "ABC2_membrane_7",
    "ABC_membrane",
    "ABC_membrane_2",
    "ABC_membrane_3",
    "ABC_export",
}                                                   #? TMD 域集合
EXCLUDE_DOMAINS = {"ABC_sub_bind", "BCA_ABC_TP_C"}  #? 排除域集合


def _safe_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _safe_int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def load_pfam_hits(pfam_csv, evalue_cutoff): #? 加载 Pfam 命中结果
    info = {}
    with open(pfam_csv, "r", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            seq_id = (row.get("target_name") or "").strip()
            domain = (row.get("query_name") or "").strip()
            if not seq_id or not domain:
                continue
            # 兼容不同 CSV 列名：i_E_value vs i_Evalue
            evalue = _safe_float(row.get("i_E_value") or row.get("i_Evalue"))
            ali_from = _safe_int(row.get("ali_from"))
            ali_to = _safe_int(row.get("ali_to"))
            rec = info.setdefault(
                seq_id,
                {
                    "nbd_start": None,
                    "tmd_end": None,
                    "nbd_domains": set(),
                    "tmd_domains": set(),
                    "exclude_domains": set(),
                    "best_nbd_any": None,
                    "best_tmd_any": None,
                    "best_exclude_any": None,
                },
            )

            # 记录所有命中（不考虑阈值）用于解释“有命中但未过阈值”的情况
            if domain in NBD_DOMAINS and evalue is not None:
                rec["best_nbd_any"] = (
                    evalue
                    if rec["best_nbd_any"] is None
                    else min(rec["best_nbd_any"], evalue)
                )
            if domain in TMD_DOMAINS and evalue is not None:
                rec["best_tmd_any"] = (
                    evalue
                    if rec["best_tmd_any"] is None
                    else min(rec["best_tmd_any"], evalue)
                )
            if domain in EXCLUDE_DOMAINS and evalue is not None:
                rec["best_exclude_any"] = (
                    evalue
                    if rec["best_exclude_any"] is None
                    else min(rec["best_exclude_any"], evalue)
                )

            # 仅在 i-Evalue 达标时计入筛选依据
            if evalue is None or evalue > evalue_cutoff:
                continue

            if domain in NBD_DOMAINS and ali_from is not None:
                rec["nbd_domains"].add(domain)
                rec["nbd_start"] = (
                    ali_from
                    if rec["nbd_start"] is None
                    else min(rec["nbd_start"], ali_from)
                )
            if domain in TMD_DOMAINS and ali_to is not None:
                rec["tmd_domains"].add(domain)
                rec["tmd_end"] = (
                    ali_to if rec["tmd_end"] is None else max(rec["tmd_end"], ali_to)
                )
            if domain in EXCLUDE_DOMAINS:
                rec["exclude_domains"].add(domain)

    return info


def load_tmhmm(tmhmm_merge):                    #? 加载 TMHMM 结果
    info = {}
    with open(tmhmm_merge, "r") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            seq_id = parts[0].strip()
            predhel = None
            topology = None
            for item in parts[1:]:
                if item.startswith("PredHel="):
                    predhel = _safe_int(item.split("=", 1)[1])
                elif item.startswith("Topology="):
                    topology = item.split("=", 1)[1]
            info[seq_id] = {"predhel": predhel, "topology": topology}
    return info


def load_signalp(signal_merge):
    info = {}
    with open(signal_merge, "r") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 10:
                continue
            seq_id = parts[0].strip()
            signal_flag = parts[9].strip()
            info[seq_id] = signal_flag == "Y"
    return info


def classify_sequences(pfam_info, tmhmm_info, signal_info, args):
    all_seq_ids = set(pfam_info) | set(tmhmm_info) | set(signal_info)
    results = []

    for seq_id in sorted(all_seq_ids):
        pfam = pfam_info.get(seq_id, {})
        tmhmm = tmhmm_info.get(seq_id, {})
        signal = signal_info.get(seq_id)

        reasons = []
        nbd_start = pfam.get("nbd_start")
        tmd_end = pfam.get("tmd_end")
        nbd_domains = pfam.get("nbd_domains", set())
        tmd_domains = pfam.get("tmd_domains", set())
        exclude_domains = pfam.get("exclude_domains", set())

        has_nbd = bool(nbd_domains)
        has_tmd = bool(tmd_domains)

        if not has_nbd:
            if pfam.get("best_nbd_any") is not None:
                reasons.append(
                    f"NBD 命中但 i-Evalue={pfam['best_nbd_any']:.2e} > {args.evalue_cutoff}"
                )
            else:
                reasons.append("未命中 NBD 域（ABC_tran/ABC_N/ABC_ATPase）")
        if not has_tmd:
            if pfam.get("best_tmd_any") is not None:
                reasons.append(
                    f"TMD 命中但 i-Evalue={pfam['best_tmd_any']:.2e} > {args.evalue_cutoff}"
                )
            else:
                reasons.append("未命中 TMD 域（ABC2_membrane/ABC_membrane/ABC_export 等）")

        predhel = tmhmm.get("predhel")
        if predhel is None:
            reasons.append("缺少 TMHMM PredHel 信息")
        elif predhel < args.predhel_min:
            reasons.append(f"PredHel={predhel} < {args.predhel_min}")

        gap = None
        if has_nbd and has_tmd and nbd_start is not None and tmd_end is not None: #! 计算 gap
            gap = nbd_start - tmd_end   #? gap = NBD 起始位置 - TMD 结束位置
            if gap < 0:                 #? TMD 在 NBD 之后
                reasons.append("TMD 在 NBD 之后（gap < 0）") 
            elif gap > args.gap_max:    #? gap 超过最大值
                reasons.append(f"gap={gap} > {args.gap_max}")
        else:
            if has_nbd and has_tmd:
                reasons.append("缺少 NBD/TMD 边界，无法计算 gap")

        # 排除条件：命中排除域 或 (PredHel<=2 且 SignalP 阳性)
        exclude_condition = False
        if exclude_domains:
            exclude_condition = True
            reasons.append(f"命中排除域：{','.join(sorted(exclude_domains))}")
        if predhel is not None and predhel <= args.predhel_exclude_max and signal is True:
            exclude_condition = True
            reasons.append("PredHel<=2 且 SignalP 阳性")

        # 必须条件：NBD+TMD 达标，PredHel>=4，TMD 在前且 gap>=0
        pass_mandatory = (
            has_nbd
            and has_tmd
            and predhel is not None
            and predhel >= args.predhel_min
            and gap is not None
            and gap >= 0
        )

        if exclude_condition:
            clazz = "C"
        elif pass_mandatory and gap is not None and gap <= args.gap_max:
            clazz = "A"
            reasons.append("满足全部必须条件")
        elif (
            has_nbd
            and has_tmd
            and predhel is not None
            and gap is not None
            and gap >= 0
        ):
            clazz = "B"
            reasons.append("满足 NBD+TMD，但 gap 或 TM 数异常")
        else:
            clazz = "C"

        results.append(
            {
                "seq_id": seq_id,
                "class": clazz,
                "reasons": "；".join(reasons) if reasons else "无",
                "nbd_domains": ",".join(sorted(nbd_domains)),
                "tmd_domains": ",".join(sorted(tmd_domains)),
                "exclude_domains": ",".join(sorted(exclude_domains)),
                "best_nbd_evalue": pfam.get("best_nbd_any"),
                "best_tmd_evalue": pfam.get("best_tmd_any"),
                "best_exclude_evalue": pfam.get("best_exclude_any"),
                "nbd_start": nbd_start,
                "tmd_end": tmd_end,
                "gap": gap,
                "predhel": predhel,
                "signalp_positive": signal,
            }
        )

    return results


def write_results(results, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    all_path = os.path.join(output_dir, "rv1819c_filter_all.tsv")
    a_path = os.path.join(output_dir, "rv1819c_class_A.tsv")
    b_path = os.path.join(output_dir, "rv1819c_class_B.tsv")
    c_path = os.path.join(output_dir, "rv1819c_class_C.tsv")

    headers = [
        "seq_id",
        "class",
        "reasons",
        "nbd_domains",
        "tmd_domains",
        "exclude_domains",
        "best_nbd_evalue",
        "best_tmd_evalue",
        "best_exclude_evalue",
        "nbd_start",
        "tmd_end",
        "gap",
        "predhel",
        "signalp_positive",
    ]

    def _open_writer(path):
        handle = open(path, "w", newline="")
        writer = csv.DictWriter(handle, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        return handle, writer

    all_handle, all_writer = _open_writer(all_path)
    a_handle, a_writer = _open_writer(a_path)
    b_handle, b_writer = _open_writer(b_path)
    c_handle, c_writer = _open_writer(c_path)

    for row in results:
        all_writer.writerow(row)
        if row["class"] == "A":
            a_writer.writerow(row)
        elif row["class"] == "B":
            b_writer.writerow(row)
        else:
            c_writer.writerow(row)

    all_handle.close()
    a_handle.close()
    b_handle.close()
    c_handle.close()


def main():
    parser = argparse.ArgumentParser(
        description="按 Rv1819c-like 过滤标准筛选 Pfam/TMHMM/SignalP 结果"
    )
    parser.add_argument("--pfam-csv", required=True, help="pfam_combined.csv 路径")
    parser.add_argument("--tmhmm-merge", required=True, help="tmhmm_merge.tsv 路径")
    parser.add_argument("--signal-merge", required=True, help="signal_merge.tsv 路径")
    parser.add_argument("--output-dir", required=True, help="输出目录")
    parser.add_argument("--evalue-cutoff", type=float, default=1e-10)
    parser.add_argument("--gap-max", type=int, default=250)
    parser.add_argument("--predhel-min", type=int, default=4)
    parser.add_argument("--predhel-exclude-max", type=int, default=2)
    args = parser.parse_args()

    pfam_info = load_pfam_hits(args.pfam_csv, args.evalue_cutoff)
    tmhmm_info = load_tmhmm(args.tmhmm_merge)
    signal_info = load_signalp(args.signal_merge)

    results = classify_sequences(pfam_info, tmhmm_info, signal_info, args)
    write_results(results, args.output_dir)


if __name__ == "__main__":
    main()
