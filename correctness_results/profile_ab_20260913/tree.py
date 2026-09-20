#!/usr/bin/env python3
"""Nested (inclusive) breakdown of the new build's CPE cycles.

The LWPF regions are hierarchical: a region's counter includes any region opened
inside it. The tree was read off the sources:

  WORKER_ALIGNMENT                      src/slave/slave.c: worker12_s_pre_fast
    MEM_CHAIN                           src/slave/bwamem.c: mem_align1_core_impl
      MEM_CHAIN_COLLECT                 bwamem.c: mem_chain
        MEM_COLLECT_{FIRST,SPLIT,LAST,SORT}   bwamem.c: mem_collect_intv
      MEM_CHAIN_BUILD                   bwamem.c: mem_chain
        CHAIN_BUILD_{REPETITIVE,SA,RID,TREE_SEARCH,MERGE,INSERT,FINALIZE}
    CHAIN_FILTER                        mem_align1_core_impl
    CHAIN_EXTENSION                     mem_align1_core_impl
      CHAIN_EXTENSION_DP                bwamem.c: mem_chain2aln
    ALIGNMENT_FINALIZE                  mem_align1_core_impl
    SAM_FORMAT                          bwamem.c: worker12_pre_fast
      MATE_RESCUE                       bwamem_pair.c: mem_sam_pe   (PE only)
        MATE_REF_FETCH                  bwamem_pair.c: swbwa_matesw_prepare
        MATE_KSW_ALIGN                  bwamem_pair.c: swbwa_matesw_run_one/_pair
          KSW_QUERY_INIT_{FORWARD,REVERSE}, KSW_DP_{FORWARD,REVERSE}
        MATE_DEDUP                      bwamem_pair.c: swbwa_matesw_finish
          DEDUP_{SORT_END,REDUNDANCY,SORT_SCORE}
      PAIRING                           bwamem_pair.c: mem_sam_pe   (PE only)
  SAM_COPY                              slave.c: worker12_s_fast   (true sibling)

A child is a SUBSET of its parent, so rows must NOT be added together. Only
WORKER_ALIGNMENT and SAM_COPY are disjoint.
"""
import sys

TREE = [
    ("WORKER_ALIGNMENT", 0),
    ("MEM_CHAIN", 1),
    ("MEM_CHAIN_COLLECT", 2),
    ("MEM_COLLECT_FIRST", 3),
    ("MEM_COLLECT_SPLIT", 3),
    ("MEM_COLLECT_LAST", 3),
    ("MEM_COLLECT_SORT", 3),
    ("MEM_CHAIN_BUILD", 2),
    ("CHAIN_BUILD_REPETITIVE", 3),
    ("CHAIN_BUILD_SA", 3),
    ("CHAIN_BUILD_RID", 3),
    ("CHAIN_BUILD_TREE_SEARCH", 3),
    ("CHAIN_BUILD_MERGE", 3),
    ("CHAIN_BUILD_INSERT", 3),
    ("CHAIN_BUILD_FINALIZE", 3),
    ("CHAIN_FILTER", 1),
    ("CHAIN_EXTENSION", 1),
    ("CHAIN_EXTENSION_DP", 2),
    ("ALIGNMENT_FINALIZE", 1),
    ("SAM_FORMAT", 1),
    ("MATE_RESCUE", 2),
    ("MATE_REF_FETCH", 3),
    ("MATE_KSW_ALIGN", 3),
    ("KSW_QUERY_INIT_FORWARD", 4),
    ("KSW_DP_FORWARD", 4),
    ("KSW_QUERY_INIT_REVERSE", 4),
    ("KSW_DP_REVERSE", 4),
    ("MATE_DEDUP", 3),
    ("DEDUP_SORT_END", 4),
    ("DEDUP_REDUNDANCY", 4),
    ("DEDUP_SORT_SCORE", 4),
    ("PAIRING", 2),
    ("SAM_COPY", 0),
]
TOP = ("WORKER_ALIGNMENT", "SAM_COPY")


def rows(path):
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(head, l.rstrip("\n").split("\t")))
                for l in fh if l.strip()]


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "out/kernel_cycles.tsv"
    variant = sys.argv[2] if len(sys.argv) > 2 else "new_p"
    data = {}
    for r in rows(path):
        if r["variant"] == variant:
            data[(r["cfg"], r["dataset"], r["mode"], r["kernel"])] = \
                float(r["cycle_avg"])

    cfg = sys.argv[3] if len(sys.argv) > 3 else None
    ds = sys.argv[4] if len(sys.argv) > 4 else None
    mode = sys.argv[5] if len(sys.argv) > 5 else None
    conds = sorted({k[:3] for k in data})
    if cfg:
        conds = [c for c in conds if c[0] == cfg]
    if ds:
        conds = [c for c in conds if c[1] == ds]
    if mode:
        conds = [c for c in conds if c[2] == mode]

    print("variant=%s   denominator = WORKER_ALIGNMENT + SAM_COPY" % variant)
    print("child regions are SUBSETS of their parent; do not sum rows\n")
    print("| region | " + " | ".join("%s/%s/%s" % c for c in conds) + " |")
    print("| --- | " + " | ".join("---:" for _ in conds) + " |")
    for node, depth in TREE:
        cells = []
        for (c, d, m) in conds:
            tot = sum(data.get((c, d, m, t), 0.0) for t in TOP)
            own = data.get((c, d, m, node), 0.0)
            cells.append("%.1f%%" % (100.0 * own / tot if tot else 0))
        indent = "\u00a0" * 3 * depth
        print("| %s%s | %s |" % (indent, node, " | ".join(cells)))


if __name__ == "__main__":
    main()
