#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Mapping file (header optional; tab/comma/space OK):
    chr1_H1   ch1
    chr2_H1   ch2
or
    key   replacement
    chr1_H1   ch1
"""

import argparse
import re
import sys
from typing import Dict, Optional

def parse_args():
    ap = argparse.ArgumentParser(
        description="Rewrite GFF3 IDs/Parents: chrom replacement via seqname map, genename padding, optional species/cultivar override. Seqname (col1) is NOT renamed."
    )
    ap.add_argument("-i", "--input",  required=True, help="Input GFF3 ('-' for stdin)")
    ap.add_argument("-o", "--output", default="-", help="Output GFF3 ('-' for stdout)")
    ap.add_argument("-m", "--chrom-map", required=True,
                    help="Seqname→replacement mapping file (2 columns; header optional).")
    ap.add_argument("-p", "--pad", type=int, default=5,
                    help="Zero-pad width for genename numbers [default: 5]")
    ap.add_argument("-s", "--species", type=str, default=None,
                    help="Override species token in IDs/Parents (if present)")
    ap.add_argument("-c", "--cultivar", type=str, default=None,
                    help="Override cultivar token in IDs/Parents (if present)")
    return ap.parse_args()

# Robust mapping loader (tab/comma/space; header optional)
_WS_SPLIT = re.compile(r"[,\t ]+")

def load_mapping(path: str) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with open(path, "r", encoding="utf-8") as fh:
        saw_first = False
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = _WS_SPLIT.split(line)
            if len(parts) < 2:
                continue
            k0, v0 = parts[0], parts[1]
            if not saw_first:
                l0, l1 = k0.lower(), v0.lower()
                if (l0 in ("key", "seqname")) and (l1 in ("replacement", "value", "target")):
                    saw_first = True
                    continue
                saw_first = True
            mapping[k0] = v0
    if not mapping:
        raise ValueError("Empty mapping. Check mapping file format.")
    return mapping

# Rewriting utilities
LEADING_TOKENS_RE = re.compile(
    r'^(?P<species>[^\.]+)\.(?P<cultivar>[^\.]+)\.(?P<haplo>H[0-9]+)\.(?P<chrom>chrom)\.(?P<rest>.+)$'
)
GENE_NUMBER_RE = re.compile(r'(?P<prefix>\b)genename(?P<num>\d+)\b')
ATTR_SPLIT_RE = re.compile(r';(?=(?:[^"]*"[^"]*")*[^"]*$)')  # split on ';' not within quotes

def pad_genename_only(s: str, pad_width: int) -> str:
    def _sub(m: re.Match) -> str:
        return m.group("prefix") + m.group("num").zfill(pad_width)
    return GENE_NUMBER_RE.sub(_sub, s)

def rewrite_token_block(token: str,
                        seqname_for_map: str,
                        chrom_map: Dict[str, str],
                        pad_width: int,
                        species_override: Optional[str],
                        cultivar_override: Optional[str]) -> str:
    # If token matches the structured prefix ... .chrom. ...
    m = LEADING_TOKENS_RE.match(token)
    if not m:
        # Still apply genename padding even if pattern doesn't match
        return pad_genename_only(token, pad_width)

    species  = species_override if species_override else m.group("species")
    cultivar = cultivar_override if cultivar_override else m.group("cultivar")
    haplo    = m.group("haplo")
    rest     = m.group("rest")

    chrom_repl = chrom_map.get(seqname_for_map, "chrom")
    rest2 = pad_genename_only(rest, pad_width)
    return f"{species}.{cultivar}.{haplo}.{chrom_repl}.{rest2}"

def rewrite_attr_value(val: str,
                       seqname_for_map: str,
                       chrom_map: Dict[str, str],
                       pad_width: int,
                       species_override: Optional[str],
                       cultivar_override: Optional[str]) -> str:
    parts = [p.strip() for p in val.split(",")]
    return ",".join(
        rewrite_token_block(p, seqname_for_map, chrom_map, pad_width, species_override, cultivar_override)
        for p in parts if p
    )

def process_attributes(attr_field: str,
                       seqname_for_map: str,
                       chrom_map: Dict[str, str],
                       pad_width: int,
                       species_override: Optional[str],
                       cultivar_override: Optional[str]) -> str:
    if attr_field == "." or not attr_field.strip():
        return attr_field
    items = ATTR_SPLIT_RE.split(attr_field.strip())
    out = []
    for item in items:
        if not item:
            continue
        if "=" not in item:
            out.append(item)
            continue
        k, v = item.split("=", 1)
        k = k.strip(); v = v.strip()
        if k in ("ID", "Parent"):
            v = rewrite_attr_value(v, seqname_for_map, chrom_map, pad_width, species_override, cultivar_override)
        else:
            v = pad_genename_only(v, pad_width)
        out.append(f"{k}={v}")
    return ";".join(out)

def main():
    args = parse_args()
    chrom_map = load_mapping(args.chrom_map)

    infh = sys.stdin if args.input == "-" else open(args.input, "r", encoding="utf-8", newline="")
    outfh = sys.stdout if args.output == "-" else open(args.output, "w", encoding="utf-8", newline="")

    try:
        for line in infh:
            if not line.strip() or line.startswith("#"):
                outfh.write(line)
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                outfh.write(line)
                continue

            seqname = cols[0]  # kept unchanged

            # Attributes rewrite (ID/Parent & genename padding; chrom token via mapping)
            cols[8] = process_attributes(
                cols[8], seqname, chrom_map, args.pad, args.species, args.cultivar
            )

            # Write back (seqname untouched)
            outfh.write("\t".join(cols) + "\n")
    finally:
        if infh is not sys.stdin: infh.close()
        if outfh is not sys.stdout: outfh.close()

if __name__ == "__main__":
    main()
