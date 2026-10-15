#!/usr/bin/env python3
"""
Trim every sequence in a cluster fasta to the local-alignment window it
shares (directly or transitively) with an anchor sequence from the same
cluster, using coordinates from a single mmseqs2 all-vs-all search pass.

Anchor choice is unchanged from the previous version: the longest member
in the cluster fasta. Any sequence header in the fasta is guaranteed to
also appear as a query/target in alignments.m8, since only cluster
*labels* get remapped by remap_cluster_ids.py, not per-sequence
accessions.

Unlike the previous version, a member with no *direct* hit to the anchor
is no longer kept at full length. Its window is instead propagated
through the alignment graph: if C only hit B, and B's window (relative to
the anchor) is already known, C's window is estimated by mapping B's
known window through the B-C alignment coordinates. Members with no path
at all back to the anchor are dropped -- there's no basis for picking a
window for them.

Propagation is an approximation: each hop assumes a roughly linear
(ungapped) correspondence between the two aligned spans, so coordinate
error can compound across multiple hops. Direct hits are exact; expanding
the steepest-bitscore edge first keeps hop counts low where possible, and
--max-hops lets you cut off long, low-confidence chains entirely.
"""

import argparse
import heapq
import itertools
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass

from DarwinsRNAHunt.genome_analysis import parse_fasta, write_fasta_output
from DarwinsRNAHunt.sequence import NucleotideSequence


@dataclass
class Hit:
    other: str
    self_lo: int
    self_hi: int
    self_is_minus: bool
    other_lo: int
    other_hi: int
    other_is_minus: bool
    bits: float


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--fasta", required=True, help="cluster_{id}.fasta (full-length members)")
    p.add_argument("--alignments", required=True, help="flanking_cds_alignments.m8 (whole dataset)")
    p.add_argument("--padding", type=int, default=5)
    p.add_argument("--min-length", type=int, default=30)
    p.add_argument("--max-length", type=int, default=10000)
    p.add_argument("--max-hops", type=int, default=0,
                    help="Drop members only reachable via more than this many hops "
                         "from the anchor (0 = unlimited).")
    p.add_argument("--diagnostics", action="store_true",
                    help="Print a hop-count histogram and drop-reason breakdown to stderr.")
    p.add_argument("--output", required=True)
    return p.parse_args()

def print_diagnostics(windows, kept_ids, dropped_out_of_range_ids, dropped_unreachable_count, anchor_id):
    hops_kept = Counter()
    hops_dropped_length = Counter()

    for seq_id in kept_ids:
        _, _, _, hops = windows[seq_id]
        hops_kept[hops] += 1
    for seq_id in dropped_out_of_range_ids:
        _, _, _, hops = windows[seq_id]
        hops_dropped_length[hops] += 1

    all_hops = sorted(set(hops_kept) | set(hops_dropped_length))
    if not all_hops:
        print("[diagnostics] no members reachable from anchor besides itself", file=sys.stderr)
        return

    max_count = max(hops_kept[h] + hops_dropped_length[h] for h in all_hops)
    scale = 40 / max_count if max_count > 0 else 1

    print(f"[diagnostics] anchor={anchor_id}", file=sys.stderr)
    print(f"{'hops':>5} {'kept':>6} {'dropped(len)':>13}  bar", file=sys.stderr)
    for h in all_hops:
        k, d = hops_kept[h], hops_dropped_length[h]
        bar = "#" * round(k * scale) + "." * round(d * scale)
        label = "anchor" if h == 0 else str(h)
        print(f"{label:>5} {k:>6} {d:>13}  {bar}", file=sys.stderr)

    print(f"[diagnostics] unreachable from anchor (dropped): {dropped_unreachable_count}", file=sys.stderr)


def normalise_span(start, end):
    """mmseqs coords: start > end means the hit sits on the minus strand
    of *that* sequence. Returns (lo, hi, is_minus) with lo < hi."""
    if start > end:
        return end, start, True
    return start, end, False


def build_graph(m8_path, member_ids):
    """
    Undirected adjacency: seq_id -> list[Hit]. Keeps only the best-bitscore
    hit per unordered pair, since mmseqs2 typically reports both A-vs-B
    and B-vs-A rows for the same pair.
    """
    best = {}
    with open(m8_path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            query, target, qstart, qend, tstart, tend, qlen, tlen, fident, bits = fields[:10]
            if query == target or query not in member_ids or target not in member_ids:
                continue
            bits = float(bits)
            key = frozenset((query, target))
            if key not in best or bits > best[key][0]:
                best[key] = (bits, query, target, int(qstart), int(qend), int(tstart), int(tend))

    graph = defaultdict(list)
    for bits, query, target, qs, qe, ts, te in best.values():
        q_lo, q_hi, q_minus = normalise_span(qs, qe)
        t_lo, t_hi, t_minus = normalise_span(ts, te)
        graph[query].append(Hit(target, q_lo, q_hi, q_minus, t_lo, t_hi, t_minus, bits))
        graph[target].append(Hit(query, t_lo, t_hi, t_minus, q_lo, q_hi, q_minus, bits))
    return graph


def propagate_windows(graph, anchor_id, anchor_lo, anchor_hi, max_hops):
    """
    Best-first propagation of the anchor's core window outward through the
    alignment graph. Returns {seq_id: (lo, hi, is_minus, hops)}, where
    is_minus means "reverse-complement this slice to match the anchor's
    reading direction" and hops is 0 for the anchor / direct hits.
    """
    resolved = {anchor_id: (anchor_lo, anchor_hi, False, 0)}
    heap = []
    counter = itertools.count()  # tie-breaker so ties never compare Hit objects
    for hit in graph[anchor_id]:
        heapq.heappush(heap, (-hit.bits, next(counter), anchor_id, hit))

    while heap:
        _, _, from_id, hit = heapq.heappop(heap)
        if hit.other in resolved:
            continue
        from_lo, from_hi, from_is_minus, from_hops = resolved[from_id]
        if max_hops and from_hops >= max_hops:
            continue

        lo = max(from_lo, hit.self_lo)
        hi = min(from_hi, hit.self_hi)
        if lo >= hi:
            continue

        self_span = hit.self_hi - hit.self_lo
        other_span = hit.other_hi - hit.other_lo
        if self_span <= 0 or other_span <= 0:
            continue

        frac_lo = (lo - hit.self_lo) / self_span
        frac_hi = (hi - hit.self_lo) / self_span

        hit_flip = hit.self_is_minus != hit.other_is_minus
        if hit_flip:
            other_lo = hit.other_lo + (1 - frac_hi) * other_span
            other_hi = hit.other_lo + (1 - frac_lo) * other_span
        else:
            other_lo = hit.other_lo + frac_lo * other_span
            other_hi = hit.other_lo + frac_hi * other_span

        other_is_minus = from_is_minus != hit_flip
        resolved[hit.other] = (round(other_lo), round(other_hi), other_is_minus, from_hops + 1)

        for next_hit in graph[hit.other]:
            if next_hit.other not in resolved:
                heapq.heappush(heap, (-next_hit.bits, next(counter), hit.other, next_hit))

    return resolved


def main():
    args = parse_args()
    fasta_records = parse_fasta(args.fasta)
    records = {}
    for sequence_id, record in fasta_records.items():
        sequence = NucleotideSequence.from_header(sequence_id, sequence=record.seq)
        if sequence is None:
            sys.exit(f"Could not parse sequence header: {sequence_id}")
        records[sequence_id] = sequence

    if not records:
        sys.exit(f"No sequences found in {args.fasta}")

    member_ids = set(records.keys())
    anchor_id = max(records, key=lambda sequence_id: len(records[sequence_id]))
    graph = build_graph(args.alignments, member_ids)

    # Anchor's own core window = union of every span other members aligned to.
    anchor_spans = [(h.self_lo, h.self_hi) for h in graph[anchor_id]]
    if anchor_spans:
        a_lo = min(s[0] for s in anchor_spans)
        a_hi = max(s[1] for s in anchor_spans)
    else:
        a_lo, a_hi = 1, len(records[anchor_id].seq)

    windows = propagate_windows(graph, anchor_id, a_lo, a_hi, args.max_hops)

    kept = []
    kept_direct = kept_propagated = 0
    dropped_unreachable = []
    dropped_out_of_range = []

    for seq_id, rec in records.items():
        window = windows.get(seq_id)
        if window is None:
            dropped_unreachable.append(seq_id)
            continue

        w_lo, w_hi, is_minus, hops = window
        seqlen = len(rec)
        lo = max(0, w_lo - 1 - args.padding)  # mmseqs coords are 1-based
        hi = min(seqlen, w_hi + args.padding)

        if hi - lo < args.min_length or hi - lo > args.max_length:
            dropped_out_of_range.append(seq_id)
            continue

        window_seq = rec.seq[lo:hi]
        if is_minus:
            window_seq = window_seq.reverse_complement()

        trimmed_rec = rec.with_sequence(window_seq)
        kept.append(trimmed_rec)
        if hops == 0:
            kept_direct += 1
        else:
            kept_propagated += 1

    write_fasta_output(kept, args.output)

    print(f"[trim] anchor={anchor_id} kept_direct={kept_direct} kept_propagated={kept_propagated} "
          f"dropped_unreachable={len(dropped_unreachable)} dropped_out_of_range={len(dropped_out_of_range)}",
          file=sys.stderr)

    if args.diagnostics:
        kept_ids = [rec.id for rec in kept]
        print_diagnostics(windows, kept_ids, dropped_out_of_range, len(dropped_unreachable), anchor_id)


if __name__ == "__main__":
    main()