#!/usr/bin/env python3
"""Compare two BaMM!motif output directories with a numeric tolerance.

Every file present in the reference directory must also exist in the
candidate directory. Files are compared line by line and token by token:
numeric tokens must agree within ``--rtol``/``--atol``; all other tokens
must match exactly. Runtime lines (``Runtime ...``) are ignored.
Reference files may be gzip-compressed (``name.gz`` is compared with
``name`` in the candidate directory).

Usage:
    compare_outputs.py REFERENCE_DIR CANDIDATE_DIR [--rtol 1e-3] [--atol 1e-5]

Exit status is 0 when all files agree and 1 otherwise.
"""

import argparse
import gzip
import math
import sys
from pathlib import Path

IGNORED_SUFFIXES = {".log"}


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("reference", type=Path, help="directory with reference outputs")
    parser.add_argument("candidate", type=Path, help="directory with new outputs")
    parser.add_argument("--rtol", type=float, default=1e-3, help="relative tolerance (default: 1e-3)")
    parser.add_argument("--atol", type=float, default=1e-5, help="absolute tolerance (default: 1e-5)")
    parser.add_argument("--max-report", type=int, default=5, help="differences to print per file (default: 5)")
    return parser.parse_args(argv)


def to_float(token):
    try:
        return float(token)
    except ValueError:
        return None


def read_lines(path):
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as handle:
            return handle.read().splitlines()
    return path.read_text().splitlines()


def compare_file(ref_path, cand_path, rtol, atol, max_report):
    """Return (n_mismatches, max_abs_diff, messages) for one file pair."""
    ref_lines = read_lines(ref_path)
    cand_lines = read_lines(cand_path)
    messages = []
    mismatches = 0
    max_diff = 0.0

    if len(ref_lines) != len(cand_lines):
        messages.append(f"line count differs: {len(ref_lines)} vs {len(cand_lines)}")
        mismatches += 1

    for lineno, (ref_line, cand_line) in enumerate(zip(ref_lines, cand_lines), start=1):
        if "Runtime" in ref_line:
            continue
        ref_tokens, cand_tokens = ref_line.split(), cand_line.split()
        if len(ref_tokens) != len(cand_tokens):
            mismatches += 1
            if len(messages) < max_report:
                messages.append(f"line {lineno}: token count differs")
            continue
        for ref_tok, cand_tok in zip(ref_tokens, cand_tokens):
            if ref_tok == cand_tok:
                continue
            ref_val, cand_val = to_float(ref_tok), to_float(cand_tok)
            if ref_val is None or cand_val is None:
                ok = False
            elif math.isinf(ref_val) or math.isinf(cand_val) or math.isnan(ref_val):
                ok = ref_tok == cand_tok
            else:
                diff = abs(ref_val - cand_val)
                max_diff = max(max_diff, diff)
                ok = diff <= atol + rtol * abs(ref_val)
            if not ok:
                mismatches += 1
                if len(messages) < max_report:
                    messages.append(f"line {lineno}: {ref_tok!r} vs {cand_tok!r}")
    return mismatches, max_diff, messages


def main(argv=None):
    args = parse_args(argv)
    ref_files = sorted(p for p in args.reference.rglob("*") if p.is_file() and p.suffix not in IGNORED_SUFFIXES)
    if not ref_files:
        print(f"no reference files found in {args.reference}", file=sys.stderr)
        return 1

    failed = False
    for ref_path in ref_files:
        rel = ref_path.relative_to(args.reference)
        if rel.suffix == ".gz":
            rel = rel.with_suffix("")
        cand_path = args.candidate / rel
        if not cand_path.exists():
            print(f"MISSING  {rel}")
            failed = True
            continue
        n_bad, max_diff, messages = compare_file(ref_path, cand_path, args.rtol, args.atol, args.max_report)
        status = "OK      " if n_bad == 0 else "DIFFERS "
        print(f"{status} {rel}  (max |diff| = {max_diff:.3g}, mismatches = {n_bad})")
        for msg in messages:
            print(f"           {msg}")
        failed |= n_bad > 0
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
