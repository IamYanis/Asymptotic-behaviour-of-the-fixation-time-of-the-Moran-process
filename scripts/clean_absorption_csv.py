#!/usr/bin/env python3
"""
Clean 'absorption_times.csv' files that may contain broken lines (extra commas, merged rows, or truncations).
- Reconstructs rows based on the expected 7 columns.
- Coerces numeric types.
- Drops incomplete/invalid rows with a detailed report.
Usage:
    python clean_absorption_csv.py --in data/absorption_times.csv --out data/absorption_times.cleaned.csv
"""

import argparse
import csv
import os
import sys
from typing import List
import math

EXPECTED_HEADER = ["N","gamma","r","Sim_ID","AbsorptionTime","AbsorbingState","Absorbed"]
NCOLS = len(EXPECTED_HEADER)

def is_float(x: str) -> bool:
    try:
        float(x)
        return True
    except:
        return False

def clean_file(in_path: str, out_path: str) -> int:
    # Read header
    with open(in_path, "r", encoding="utf-8", errors="replace", newline="") as f:
        header_line = f.readline().strip()
        header = [h.strip() for h in header_line.split(",")] if header_line else []

    if header != EXPECTED_HEADER:
        print(f"[WARN] Header mismatch. Found: {header}. Expected: {EXPECTED_HEADER}")
        # If file has no header, assume tokens start immediately and use expected header
        if header and all(is_float(h) for h in header):  # file likely has no header, it started with numbers
            # rewind and treat first line as data
            with open(in_path, "r", encoding="utf-8", errors="replace", newline="") as f:
                raw = f.read()
            data_text = raw
            header = EXPECTED_HEADER[:]
        else:
            with open(in_path, "r", encoding="utf-8", errors="replace", newline="") as f:
                data_text = f.read()
            if header_line:
                data_text = data_text[len(header_line):]  # remove the mismatched header line
            header = EXPECTED_HEADER[:]
    else:
        with open(in_path, "r", encoding="utf-8", errors="replace", newline="") as f:
            _ = f.readline()  # skip header
            data_text = f.read()

    # Tokenize on commas or newlines
    tokens = []
    field = []
    for ch in data_text:
        if ch in (',', '\n', '\r'):
            token = ''.join(field).strip()
            if token != "":
                tokens.append(token)
            field = []
        else:
            field.append(ch)
    # last token
    if field:
        token = ''.join(field).strip()
        if token != "":
            tokens.append(token)

    # Build rows of 7, dropping leftovers
    rows = [tokens[i:i+NCOLS] for i in range(0, len(tokens) - len(tokens) % NCOLS, NCOLS)]
    leftovers = len(tokens) % NCOLS

    # Validate and coerce
    def try_int(x: str):
        try:
            # some integers may arrive as '1000.0'
            as_float = float(x)
            if math.isfinite(as_float) and abs(as_float - round(as_float)) < 1e-9:
                return int(round(as_float))
            return int(x)
        except:
            return None

    def try_float(x: str):
        try:
            v = float(x)
            return v
        except:
            return None

    cleaned = []
    bad = 0
    for rec in rows:
        if len(rec) != NCOLS:
            bad += 1
            continue
        N, gamma, r, Sim_ID, AbsorptionTime, AbsorbingState, Absorbed = rec
        N = try_int(N)
        Sim_ID = try_int(Sim_ID)
        AbsorbingState = try_int(AbsorbingState)
        Absorbed = try_int(Absorbed)
        gamma = try_float(gamma)
        r = try_float(r)
        AbsorptionTime = try_float(AbsorptionTime)

        if None in (N, gamma, r, Sim_ID, AbsorptionTime, AbsorbingState, Absorbed):
            bad += 1
            continue
        cleaned.append([N, gamma, r, Sim_ID, AbsorptionTime, AbsorbingState, Absorbed])

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(EXPECTED_HEADER)
        w.writerows(cleaned)

    print(f"[OK] Wrote {len(cleaned)} valid rows to: {out_path}")
    print(f"[INFO] Dropped {bad} malformed rows. Leftover tokens (ignored): {leftovers}")
    return 0

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--in", dest="in_path", required=True, help="Input CSV path (possibly corrupted)")
    p.add_argument("--out", dest="out_path", required=True, help="Output cleaned CSV path")
    args = p.parse_args()
    return clean_file(args.in_path, args.out_path)

if __name__ == "__main__":
    sys.exit(main())
