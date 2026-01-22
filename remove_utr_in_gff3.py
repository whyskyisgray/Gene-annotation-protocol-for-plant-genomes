#!/usr/bin/env python3
import sys, argparse, re
from collections import defaultdict

def parse_attrs(s):
    d = {}
    for part in s.split(";"):
        if not part or "=" not in part: 
            continue
        k, v = part.split("=", 1)
        d[k] = v
    return d

def attrs_to_str(d):
    items = []
    if "ID" in d:
        items.append(("ID", d["ID"]))
    for k in sorted(k for k in d.keys() if k != "ID"):
        items.append((k, d[k]))
    return ";".join(f"{k}={v}" for k,v in items)

def main():
    ap = argparse.ArgumentParser(
        description="Redefine gene ranges based on union of CDS or exon (UTRs excluded)."
    )
    ap.add_argument("-i", "--input", help="GFF3 input (default: stdin)")
    ap.add_argument("-o", "--output", help="GFF3 output (default: stdout)")
    ap.add_argument("--mode", choices=["cds","exon"], default="cds",
                    help="Basis for redefining gene range: cds or exon (default: cds)")
    ap.add_argument("--fallback-to-exon-if-no-cds", action="store_true",
                    help="When mode=cds and a gene has no CDS, fall back to exon union.")
    ap.add_argument("--drop-genes-without-basis", action="store_true",
                    help="If no basis features (CDS/exon) exist, drop the gene line. Default: keep original gene range.")
    ap.add_argument("--touch-transcripts", action="store_true",
                    help="Also redefine transcript ranges (min/max of chosen basis under each transcript).")
    args = ap.parse_args()

    fin  = sys.stdin if not args.input  else open(args.input, "r")
    fout = sys.stdout if not args.output else open(args.output, "w")


    lines = []  # list of (is_feature, parsed_fields or raw line)
    # feature: (seqid, source, type, start, end, score, strand, phase, attrs_dict)
    features = []
    id2idx = {}  # feature ID -> index in features

    gene_children = defaultdict(set)       # gene_id -> transcript_ids
    tx_parent = {}                         # tx_id -> gene_id
    tx_children = defaultdict(lambda: {"exon":[], "CDS":[]})  # tx_id -> dict


    for raw in fin:
        if raw.startswith("#"):
            lines.append((False, raw.rstrip("\n")))
            continue
        raw = raw.rstrip("\n")
        if not raw:
            lines.append((False, raw))
            continue
        parts = raw.split("\t")
        if len(parts) != 9:
            lines.append((False, raw))
            continue

        seqid, source, ftype, start, end, score, strand, phase, attr = parts
        try:
            start_i = int(start); end_i = int(end)
        except:
            lines.append((False, raw))
            continue
        attrs = parse_attrs(attr)
        feat = [seqid, source, ftype, start_i, end_i, score, strand, phase, attrs]
        idx = len(features)
        features.append(feat)

        fid = attrs.get("ID")
        if fid:
            id2idx[fid] = idx

        parents = []
        if "Parent" in attrs:
            parents = attrs["Parent"].split(",")
        if ftype.lower() in ("mrna","transcript","ncrna","lncrna","snorna","snrna","trna"):

            for p in parents:
                if p:
                    tx_parent[attrs.get("ID","")] = p
                    gene_children[p].add(attrs.get("ID",""))
        elif ftype == "exon" or ftype == "CDS":
            for p in parents:
                if p:
                    tx_children[p][ftype].append(idx)


        lines.append((True, idx))

    tx_new_range = {}
    for tx_id, gene_id in tx_parent.items():
        child = tx_children.get(tx_id)
        if not child: 
            continue

        basis_type = "CDS" if args.mode == "cds" else "exon"
        idxs = child[basis_type]
        if args.mode == "cds" and not idxs and args.fallback_to_exon_if_no_cds:
            idxs = child["exon"]
        if not idxs:
            continue
        smin = min(features[i][3] for i in idxs)
        smax = max(features[i][4] for i in idxs)
        tx_new_range[tx_id] = (smin, smax)

    gene_new_range = {}
    for gene_id, tx_ids in gene_children.items():

        spans = []
        for tx_id in tx_ids:
            child = tx_children.get(tx_id)
            if not child:
                continue
            basis_type = "CDS" if args.mode == "cds" else "exon"
            idxs = child[basis_type]
            if args.mode == "cds" and not idxs and args.fallback_to_exon_if_no_cds:
                idxs = child["exon"]
            if idxs:
                spans.extend((features[i][3], features[i][4]) for i in idxs)
        if spans:
            smin = min(s for s,_ in spans)
            smax = max(e for _,e in spans)
            gene_new_range[gene_id] = (smin, smax)

    for is_feat, payload in lines:
        if not is_feat:
            print(payload, file=fout)
            continue
        idx = payload
        seqid, source, ftype, start_i, end_i, score, strand, phase, attrs = features[idx]
        out_start, out_end = start_i, end_i

        if ftype == "gene":
            gid = attrs.get("ID")
            if gid in gene_new_range:
                out_start, out_end = gene_new_range[gid]
            else:
                if args.drop_genes_without_basis:
                    continue

        elif args.touch_transcripts and ftype.lower() in ("mrna","transcript","ncrna","lncrna","snorna","snrna","trna"):
            tid = attrs.get("ID")
            if tid in tx_new_range:
                out_start, out_end = tx_new_range[tid]

        out = [
            seqid, source, ftype,
            str(out_start), str(out_end),
            score, strand, phase,
            attrs_to_str(attrs)
        ]
        print("\t".join(out), file=fout)

    if fout is not sys.stdout:
        fout.close()
    if fin is not sys.stdin:
        fin.close()

if __name__ == "__main__":
    main()

