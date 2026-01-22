#!/usr/bin/env python3


import argparse
import os
import re
from collections import defaultdict
from typing import Dict, List, Tuple, Optional

##############################################################################
# 1) Attribute parsing (GTF/GFF3)
##############################################################################

_attr_split_re = re.compile(r';\s*')
_key_val_gff3  = re.compile(r'^(?P<k>[^=]+)=(?P<v>.*)$')
_key_val_gtf   = re.compile(r'^(?P<k>\S+)\s+(?P<v>.*)$')

def parse_attr(attr_field: str, feature_type: str) -> Dict[str, str]:
    attr_field = attr_field.strip().strip(';')
    if not attr_field:
        return {}
    attrs: Dict[str, str] = {}
    for token in _attr_split_re.split(attr_field):
        token = token.strip()
        if not token:
            continue
        if '=' in token:  # GFF3
            m = _key_val_gff3.match(token)
            if m:
                attrs[m['k'].strip()] = m['v'].strip().strip('"')
        elif ' ' in token:  # GTF
            m = _key_val_gtf.match(token)
            if m:
                attrs[m['k']] = m['v'].strip().strip('"')
        else:  # bare value fallback
            if feature_type.lower() == 'gene':
                attrs['gene_id'] = token
            else:
                attrs.setdefault('gene_id', token)
    return attrs

##############################################################################
# 2) Optional prefix handling
##############################################################################

_id_key_gtf   = ('gene_id', 'transcript_id', 'protein_id')

def _handle_value(val: str, prefix: str) -> str:
    parts = [p for p in val.split(',') if p]
    return ','.join(prefix + p if not p.startswith(prefix) else p for p in parts)

def apply_prefix_to_attr_string(attr_field: str, prefix: str, is_gff3: bool) -> str:
    if not prefix:
        return attr_field
    prefix = prefix + '_'
    if is_gff3:
        attrs = _attr_split_re.split(attr_field.strip())
        new_fields = []
        for a in attrs:
            if not a.strip():
                continue
            if a.startswith('ID=') or a.startswith('Parent='):
                k, v = a.split('=', 1)
                new_fields.append(f"{k}={_handle_value(v, prefix)}")
            else:
                new_fields.append(a)
        return ';'.join(new_fields)
    else:  # GTF
        attrs = _attr_split_re.split(attr_field.strip())
        new_fields = []
        for a in attrs:
            if not a.strip():
                continue
            for k in _id_key_gtf:
                if a.startswith(k):
                    key, rest = a.split(' ', 1)
                    val = rest.strip().strip('"')
                    new_fields.append(f'{key} "{_handle_value(val, prefix)}"')
                    break
            else:
                new_fields.append(a)
        return '; '.join(new_fields) + ';'

def maybe_prefix_line(line: str, prefix: Optional[str]) -> str:
    if not prefix or line.startswith('#'):
        return line
    cols = line.rstrip('\n').split('\t')
    if len(cols) != 9:
        return line
    is_gff3 = '=' in cols[8] and ('"' not in cols[8])
    cols[8] = apply_prefix_to_attr_string(cols[8], prefix, is_gff3)
    return '\t'.join(cols) + '\n'

##############################################################################
# 3) Data structures & loaders
##############################################################################

class GeneRec:
    def __init__(self, gid: str, chr_: str, strand: str, start: int, end: int, origin: str):
        self.gid = gid
        self.chr = chr_
        self.strand = strand
        self.start = start
        self.end = end
        self.origin = origin  # 'short' or 'long'
        self.lines: List[str] = []
        self.exons: List[Tuple[int, int]] = []
        self.cds:   List[Tuple[int, int]] = []

def _ensure_gene(genes: Dict[str, 'GeneRec'], gid: str,
                 chr_: str, strand: str, start: int, end: int, origin_label: str,
                 order: List['GeneRec']) -> GeneRec:
    if gid not in genes:
        genes[gid] = GeneRec(gid, chr_, strand, int(start), int(end), origin_label)
        order.append(genes[gid])
    else:
        g = genes[gid]
        g.start = min(g.start, int(start))
        g.end   = max(g.end,   int(end))
    return genes[gid]

def load_gtf(path: str, origin_label: str) -> Tuple[List[GeneRec], Dict[str, GeneRec]]:
    genes: Dict[str, GeneRec] = {}
    transcript_to_gene: Dict[str, str] = {}
    order: List[GeneRec] = []
    other_feature_lines = []

    with open(path) as fh:
        for ln in fh:
            if ln.startswith('#') or not ln.strip():
                continue
            cols = ln.rstrip('\n').split('\t')
            if len(cols) < 9:
                continue
            seqid, _, feature, start, end, _, strand, _, attr = cols
            attrs = parse_attr(attr, feature)
            ftype = feature.lower()

            if ftype == 'gene':
                gid = attrs.get('ID') or attrs.get('gene_id')
                if gid:
                    g = _ensure_gene(genes, gid, seqid, strand, int(start), int(end), origin_label, order)
                    g.lines.append(ln)

            elif ftype in ('transcript', 'mrna'):
                tid = attrs.get('ID') or attrs.get('transcript_id')
                parent_gid = attrs.get('Parent') or attrs.get('gene_id')
                if not (tid and parent_gid) and attrs.get('transcript_id') and attrs.get('gene_id'):
                    tid = attrs['transcript_id']
                    parent_gid = attrs['gene_id']
                if tid and parent_gid:
                    transcript_to_gene[tid] = parent_gid
                    g = _ensure_gene(genes, parent_gid, seqid, strand, int(start), int(end), origin_label, order)
                    g.lines.append(ln)
                else:
                    other_feature_lines.append((cols, ln))
            else:
                other_feature_lines.append((cols, ln))

    # Pass 2
    for cols, ln in other_feature_lines:
        seqid, _, feature, start, end, _, strand, _, attr = cols
        attrs = parse_attr(attr, feature)
        ftype = feature.lower()

        parent_tid = attrs.get('Parent') or attrs.get('transcript_id')
        parent_gid = attrs.get('gene_id')
        if not parent_gid and parent_tid:
            parent_gid = transcript_to_gene.get(parent_tid)
        if not parent_gid:
            continue

        g = _ensure_gene(genes, parent_gid, seqid, strand, int(start), int(end), origin_label, order)
        g.lines.append(ln)

        if ftype == 'exon':
            g.exons.append((int(start), int(end)))
        elif ftype == 'cds':
            g.cds.append((int(start), int(end)))

    # Fallback spans for overlap math
    for g in genes.values():
        if not g.cds and not g.exons:
            g.exons.append((g.start, g.end))

    return order, genes

##############################################################################
# 4) Clustering
##############################################################################

def cluster_genes(genes: List[GeneRec]) -> List[List[GeneRec]]:
    clusters: List[List[GeneRec]] = []
    by_chr_strand: Dict[Tuple[str, str], List[GeneRec]] = defaultdict(list)
    for g in genes:
        by_chr_strand[(g.chr, g.strand)].append(g)

    for (chr_, strand), lst in by_chr_strand.items():
        lst.sort(key=lambda x: x.start)
        cur: List[GeneRec] = []
        cur_end = -1
        for g in lst:
            if not cur:
                cur = [g]
                cur_end = g.end
            else:
                if g.start <= cur_end:
                    cur.append(g)
                    cur_end = max(cur_end, g.end)
                else:
                    clusters.append(cur)
                    cur = [g]
                    cur_end = g.end
        if cur:
            clusters.append(cur)
    return clusters

##############################################################################
# 5) Case classification & CDS-overlap
##############################################################################

def classify(cluster: List[GeneRec]) -> int:
    short_cnt = sum(1 for g in cluster if g.origin == 'short')
    long_cnt  = sum(1 for g in cluster if g.origin == 'long')

    if short_cnt == 0 and long_cnt > 0:
        return 4
    if short_cnt > 0 and long_cnt == 0:
        return 3
    if short_cnt == 1 and long_cnt == 1:
        s = next(g for g in cluster if g.origin == 'short')
        l = next(g for g in cluster if g.origin == 'long')
        if s.start == l.start and s.end == l.end:
            return 6
    if short_cnt == 1 and long_cnt >= 2:
        return 1
    if short_cnt >= 2 and long_cnt == 1:
        return 2
    if short_cnt >= 1 and long_cnt >= 1:
        return 5
    return 5

def _span_len(spans: List[Tuple[int, int]]) -> int:
    return sum(e - s + 1 for s, e in spans)

def _calc_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)

def _get_overlap_basis(g: GeneRec) -> List[Tuple[int, int]]:
    """Prefer CDS; fallback to exons; fallback to gene span."""
    if g.cds:
        return g.cds
    if g.exons:
        return g.exons
    return [(g.start, g.end)]

def find_cds_overlaps(cluster: List[GeneRec]) -> List[Tuple[str, str, int, float]]:
    """
    Return (short_gid, long_gid, overlap_bp, overlap_fraction_wrt_short_basis),
    where basis is CDS (preferred), else exon, else gene span.
    """
    pairs: List[Tuple[str, str, int, float]] = []
    shorts = [g for g in cluster if g.origin == 'short']
    longs  = [g for g in cluster if g.origin == 'long']
    if not shorts or not longs:
        return pairs

    short_basis_len: Dict[str, int] = {g.gid: _span_len(_get_overlap_basis(g)) for g in shorts}

    for s in shorts:
        s_basis = _get_overlap_basis(s)
        s_bbox = (s.start, s.end)
        for l in longs:
            # quick bbox check
            if _calc_overlap(s_bbox[0], s_bbox[1], l.start, l.end) == 0:
                continue
            l_basis = _get_overlap_basis(l)
            ov_bp = 0
            for ss, se in s_basis:
                for ls, le in l_basis:
                    ov_bp += _calc_overlap(ss, se, ls, le)
            denom = short_basis_len[s.gid] or 1
            frac = ov_bp / denom
            pairs.append((s.gid, l.gid, ov_bp, frac))
    return pairs

##############################################################################
# 6) Writers & merger
##############################################################################

def write_cluster_file(cluster: List[GeneRec], cluster_id: str, out_dir: str, prefix: Optional[str]):
    path = os.path.join(out_dir, f'{cluster_id}.gff3')
    with open(path, 'w') as fh:
        fh.write('##gff-version 3\n')
        for g in sorted(cluster, key=lambda x: x.start):
            for ln in g.lines:
                fh.write(maybe_prefix_line(ln, prefix))

def merge_annotations(short_recs: Dict[str, GeneRec], long_recs: Dict[str, GeneRec],
                      to_remove: set, to_include_long: set, merged_path: str, prefix: Optional[str]):
    with open(merged_path, 'w') as out:
        out.write('##gff-version 3\n')
        retained = [g for g in short_recs.values() if g.gid not in to_remove]
        added    = [long_recs[gid] for gid in to_include_long if gid in long_recs]
        all_genes = sorted(retained + added, key=lambda x: (x.chr, x.start))
        for g in all_genes:
            for ln in g.lines:
                out.write(maybe_prefix_line(ln, prefix))

##############################################################################
# 7) Main
##############################################################################

def main():
    parser = argparse.ArgumentParser(
        description="Merge short/long annotations; decide replacement by CDS-overlap."
    )
    parser.add_argument('-s', '--short', required=True, help='Short-read annotation (GTF/GFF3)')
    parser.add_argument('-l', '--long',  required=True, help='Long-read annotation (GTF/GFF3)')
    parser.add_argument('-o', '--output_prefix', required=True, help='Output prefix')
    parser.add_argument('--prefix', default=None, help='Optional ID/Parent prefix in outputs')
    parser.add_argument(
        '--min_short_overlap', type=float, default=0.7,
        help='[0–1] Short CDS-length (fallback exon/gene). '
             '0 → >0. 1 → 100%% overlapped. '
             'default 0.7'
    )
    args = parser.parse_args()

    short_list, short_map = load_gtf(args.short, 'short')
    long_list,  long_map  = load_gtf(args.long,  'long')

    clusters = cluster_genes(short_list + long_list)

    dir_case1 = f'{args.output_prefix}_one_short_multiple_long'
    dir_case2 = f'{args.output_prefix}_multiple_short_one_long'
    dir_case4 = f'{args.output_prefix}_long_read_only'
    for d in (dir_case1, dir_case2, dir_case4):
        os.makedirs(d, exist_ok=True)

    meta_path      = f'{args.output_prefix}_clusters_metadata.tsv'
    overlap_path   = f'{args.output_prefix}_cds_overlap.tsv'
    decisions_path = f'{args.output_prefix}_replacement_decisions.tsv'
    merged_path    = f'{args.output_prefix}_merged.gff3'

    remove_short_ids: set = set()
    include_long_ids: set = set()
    case_idx = defaultdict(int)

    def _passed(best_f: float, thr: float) -> bool:
        if thr == 0.0:
            return best_f > 0.0     
        if thr == 1.0:
            return best_f >= 1.0   
        return best_f >= thr

    with open(meta_path, 'w') as meta_f, \
         open(overlap_path, 'w') as ov_f, \
         open(decisions_path, 'w') as decisions_f:

        meta_f.write("Cluster_ID\tCase\tChr\tShort_cnt\tLong_cnt\tShort_IDs\tLong_IDs\tStart\tEnd\n")
        ov_f.write("Cluster_ID\tShort_ID\tLong_ID\tCDS_Overlap_bp\tOverlap_fraction\n")
        decisions_f.write("Cluster_ID\tShort_ID\tBest_Long_ID\tBest_Overlap_Fraction\tPass\tThreshold\n")

        for cl in clusters:
            case = classify(cl)
            if case not in (1, 2, 4):
                continue

            case_idx[case] += 1
            out_dir = {1: dir_case1, 2: dir_case2, 4: dir_case4}[case]
            cl_prefix = {1: 'one_short_multiple_long', 2: 'multiple_short_one_long', 4: 'long_read_only'}[case]
            cluster_id = f"{cl_prefix}_{case_idx[case]}"

            short_ids = [g.gid for g in cl if g.origin == 'short']
            long_ids  = [g.gid for g in cl if g.origin == 'long']

            pairs = find_cds_overlaps(cl)
            for s_id, l_id, ov_bp, ov_frac in pairs:
                ov_f.write(f"{cluster_id}\t{s_id}\t{l_id}\t{ov_bp}\t{ov_frac:.6f}\n")

            if case in (1, 2):

                include_long_ids.update(long_ids)
                best_for_short: Dict[str, Tuple[str, float]] = {s: ('-', 0.0) for s in short_ids}
                for s_id, l_id, _ov_bp, ov_frac in pairs:
                    if ov_frac > best_for_short[s_id][1]:
                        best_for_short[s_id] = (l_id, ov_frac)

                for s_id in short_ids:
                    best_l, best_f = best_for_short.get(s_id, ('-', 0.0))
                    pass_flag = int(_passed(best_f, args.min_short_overlap) and best_l != '-')
                    decisions_f.write(
                        f"{cluster_id}\t{s_id}\t{best_l}\t{best_f:.6f}\t{pass_flag}\t{args.min_short_overlap}\n"
                    )
                    if pass_flag:
                        remove_short_ids.add(s_id)

            elif case == 4:
                include_long_ids.update(long_ids)
                # no short decisions to write

            meta_f.write(
                f"{cluster_id}\t{case}\t{cl[0].chr}\t{len(short_ids)}\t{len(long_ids)}\t"
                f"{','.join(short_ids) if short_ids else '-'}\t{','.join(long_ids) if long_ids else '-'}\t"
                f"{min(g.start for g in cl)}\t{max(g.end for g in cl)}\n"
            )
            write_cluster_file(cl, cluster_id, out_dir, args.prefix)

    merge_annotations(short_map, long_map, remove_short_ids, include_long_ids, merged_path, args.prefix)

    print("[Summary]")
    print(f"  Total short genes     : {len(short_map)}")
    print(f"  Total long genes      : {len(long_map)}")
    print(f"  Removed short genes   : {len(remove_short_ids)}")
    print(f"  Included long genes   : {len(include_long_ids)}")
    print("Outputs:")
    print(f"  {merged_path}")
    print(f"  {meta_path}")
    print(f"  {overlap_path}")
    print(f"  {decisions_path}")
    print(f"  {dir_case1}/  {dir_case2}/  {dir_case4}/")

if __name__ == '__main__':
    main()

