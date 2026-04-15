#!/usr/bin/env python3
from __future__ import annotations
import argparse
import csv
import math
import re
import sys
from pathlib import Path

REQUIRED_SOURCE_PATTERNS = {
    "two_layer_main": [
        r"_central_mpc",
        r"_anc_mpc",
        r"_actuator_sat",
        r"set_reference_horizon_control_preview",
        r"_dhat_hor",
        r"_zref_hor",
        r"_vref_hor",
        r"predict_horizon",
        r"_hopf_preview_ok",
    ],
    "guarded_mpc": [
        r"bridge_owner_match",
        r"last_export_ok",
        r"rodas_stub_touch_count",
        r"mark_rodas_stub_touched",
        r"activate_bridge\(",
    ],
    "profile_bundle": [
        r"make_runtime_profile_bundle",
        r"engineering_default",
        r"paper_reproduction",
    ],
}

DEFAULT_THRESHOLDS = {
    "resid": None,
    "prvupd": None,
    "tubee": None,
    "tubev": None,
    "tubedv": None,
    "zviol": 1e-5,
    "xviol": 1e-5,
    "cvu": 1e-5,
    "avu": 1e-5,
    "clipn": None,
    "cexok": 0.5,
    "cbrdg": 0.5,
    "abrdg": 0.5,
    "hprev": 0.5,
    "prvspn": None,
}

def check_patterns(name: str, text: str, patterns: list[str]) -> list[str]:
    missing = []
    for pat in patterns:
        if not re.search(pat, text):
            missing.append(pat)
    return missing

def load_metric_csv(path: Path) -> dict[str, list[float]]:
    data: dict[str, list[float]] = {}
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            key = row.get("key")
            value = row.get("value")
            if key is None or value is None:
                continue
            try:
                fv = float(value)
            except ValueError:
                continue
            data.setdefault(key, []).append(fv)
    return data

def summarize_metrics(data: dict[str, list[float]]) -> list[str]:
    lines: list[str] = []
    for key, threshold in DEFAULT_THRESHOLDS.items():
        vals = data.get(key, [])
        if not vals:
            lines.append(f"[WARN] metric {key} not found")
            continue

        vmax = max(vals)
        vmin = min(vals)
        vlast = vals[-1]
        if threshold is None:
            lines.append(f"[INFO] {key}: min={vmin:.6g} max={vmax:.6g} last={vlast:.6g}")
        else:
            ok = vmax <= threshold if key not in {"cexok", "cbrdg", "abrdg", "hprev"} else vmin >= threshold
            verdict = "PASS" if ok else "FAIL"
            lines.append(
                f"[{verdict}] {key}: min={vmin:.6g} max={vmax:.6g} last={vlast:.6g} threshold={threshold:.6g}"
            )
    return lines

def main() -> int:
    ap = argparse.ArgumentParser(description="Validate Hopf + Tube-MPC整改后的结构与日志指标")
    ap.add_argument("--main", type=Path, required=True, help="AuvHopf 主文件")
    ap.add_argument("--mpc", type=Path, required=True, help="guarded mpc_optim.cpp")
    ap.add_argument("--common", type=Path, required=True, help="auv_common.hpp")
    ap.add_argument("--hopf", type=Path, help="可选：hopf_observer.cpp，用于检查 autonomous preview rollout")
    ap.add_argument("--metrics-csv", type=Path, help="可选：debug_key_value 导出 CSV，至少包含 key,value 两列")
    args = ap.parse_args()

    main_text = args.main.read_text(encoding="utf-8")
    mpc_text = args.mpc.read_text(encoding="utf-8")
    common_text = args.common.read_text(encoding="utf-8")
    hopf_text = args.hopf.read_text(encoding="utf-8") if args.hopf else ""

    failures = 0
    source_texts = {
        "two_layer_main": main_text + "\n" + hopf_text,
        "guarded_mpc": mpc_text,
        "profile_bundle": common_text,
    }

    for name, pats in REQUIRED_SOURCE_PATTERNS.items():
        missing = check_patterns(name, source_texts[name], pats)
        if missing:
            failures += 1
            print(f"[FAIL] {name} 缺少关键模式:")
            for pat in missing:
                print(f"  - {pat}")
        else:
            print(f"[PASS] {name} 结构关键模式齐全")

    if args.metrics_csv:
        if not args.metrics_csv.exists():
            print(f"[FAIL] metrics csv not found: {args.metrics_csv}")
            return 1
        print("\n[Metrics]")
        data = load_metric_csv(args.metrics_csv)
        for line in summarize_metrics(data):
            print(line)
            if line.startswith("[FAIL]"):
                failures += 1

    return 1 if failures else 0

if __name__ == "__main__":
    raise SystemExit(main())
