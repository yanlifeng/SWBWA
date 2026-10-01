"""Compare actual CPE rescue scheduling with the host BWA direction loop.

KSW and reference fetching are stubbed to exhaustively exercise failed rescue,
successful insertion, skipped directions and rejected windows. Dedup is an
order-sensitive stand-in: both its input arrays and the final output must match.
This tests orchestration, not SIMD arithmetic or target-platform performance.
"""
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
HARNESS = r'''
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bwamem.h"
#include "kvec.h"
#include "ksw.h"
#define swbwa_cpe_profile_start(...) ((void)0)
#define swbwa_cpe_profile_stop(...) ((void)0)
#define ksw_align2_matesw ksw_align2

static int modes[2][4], lengths[2], calls[2], event_count[2];
static unsigned long long events[2][16];
static int mem_infer_dir(int64_t lp, int64_t a, int64_t b, int64_t *dist)
{
    (void)lp; (void)a; (void)b;
    *dist = 1000000;
    return 0;
}

uint8_t *bns_fetch_seq(const bntseq_t *bns, const uint8_t *pac,
        int64_t *beg, int64_t mid, int64_t *end, int *rid)
{
    int owner = *beg > 50000, r;
    int64_t anchor = 10000 + owner * 100000;
    int64_t starts[4] = {anchor + 20, anchor + 20 - lengths[owner],
                        anchor - 200 - lengths[owner], anchor - 200};
    uint8_t *ref = calloc((size_t)(*end - *beg), 1);
    (void)bns; (void)pac; (void)mid;
    assert(ref && *end - *beg > 2);
    for (r = 0; r < 4 && *beg != starts[r]; ++r) {}
    assert(r < 4 && modes[owner][r] != 0);
    *rid = modes[owner][r] == 1 ? -1 : owner;
    if (modes[owner][r] == 2) *end = *beg + 5;
    ref[0] = modes[owner][r];
    ref[1] = owner;
    return ref;
}

kswr_t ksw_align2(int qlen, uint8_t *query, int tlen, uint8_t *target,
        int m, const int8_t *mat, int od, int ed, int oi, int ei,
        int xtra, kswq_t **qry)
{
    kswr_t result = {0};
    int owner = target[1];
    (void)query; (void)tlen; (void)m; (void)mat;
    (void)od; (void)ed; (void)oi; (void)ei; (void)xtra; (void)qry;
    assert(qlen == lengths[owner] && target[0] >= 3);
    assert(event_count[owner] < 16);
    events[owner][event_count[owner]++] = target[0];
    ++calls[owner];
    result.score = target[0] == 3 ? 18 : 50;
    result.qb = target[0] == 4 ? -1 : 0;
    result.qe = result.te = 29;
    return result;
}

static void ksw_align2_matesw_dual_forward(
        int ql0, uint8_t *q0, int tl0, uint8_t *t0,
        int ql1, uint8_t *q1, int tl1, uint8_t *t1,
        int m, const int8_t *mat, int od, int ed, int oi, int ei,
        int xtra, kswr_t results[2])
{
    results[0] = ksw_align2(ql0, q0, tl0, t0, m, mat, od, ed, oi, ei, xtra, 0);
    results[1] = ksw_align2(ql1, q1, tl1, t1, m, mat, od, ed, oi, ei, xtra, 0);
}

int mem_sort_dedup_patch(const mem_opt_t *opt, const bntseq_t *bns,
        const uint8_t *pac, uint8_t *query, int n, mem_alnreg_t *a)
{
    unsigned long long hash = 14695981039346656037ULL;
    int owner = a[0].rid, i;
    (void)opt;
    assert(!bns && !pac && !query && n > 0);
    for (i = 0; i < n; ++i) {
        hash = (hash ^ a[i].rb) * 1099511628211ULL;
        hash = (hash ^ a[i].re) * 1099511628211ULL;
        hash = (hash ^ a[i].qb) * 1099511628211ULL;
        hash = (hash ^ a[i].qe) * 1099511628211ULL;
        hash = (hash ^ a[i].score) * 1099511628211ULL;
    }
    assert(event_count[owner] < 16);
    events[owner][event_count[owner]++] = hash;
    /* Greedy dedup need not be idempotent; make omission/order observable. */
    return n > 1 ? n - 1 : n;
}

/* RESCUE_IMPLEMENTATIONS */

static void initialize(mem_alnreg_v *v, int owner)
{
    int j;
    v->n = v->m = 13;
    v->a = calloc(v->n, sizeof(*v->a));
    assert(v->a);
    for (j = 0; j < (int)v->n; ++j) {
        v->a[j].rid = owner;
        v->a[j].rb = 1000000 + j * 1000;
        v->a[j].re = v->a[j].rb + 100;
        v->a[j].qe = 100;
        v->a[j].score = 70 - j;
    }
}

int main(void)
{
    mem_opt_t opt = {0};
    bntseq_t bns = {0};
    mem_alnreg_t anchors[2] = {{0}, {0}};
    uint8_t query[2][160] = {{0}, {0}}, pac = 0;
    const uint8_t *queries[2] = {query[0], query[1]};
    int test, variant, owner, r;
    opt.min_seed_len = 19;
    opt.a = 1;
    bns.l_pac = 10000000;
    for (owner = 0; owner < 2; ++owner) {
        anchors[owner].rb = 10000 + owner * 100000;
        anchors[owner].rid = owner;
    }
    for (test = 0; test < 1296; ++test) {
        mem_pestat_t pes[2][4] = {{{0}}};
        for (owner = 0; owner < 2; ++owner) {
            int code = owner ? (test * 587 + 23) % 1296 : test;
            for (r = 0; r < 4; ++r, code /= 6) {
                modes[owner][r] = code % 6;
                pes[owner][r].failed = modes[owner][r] == 0;
                pes[owner][r].low = 20;
                pes[owner][r].high = 200;
            }
        }
        for (variant = 0; variant < 3; ++variant) {
            mem_alnreg_v expected[2], actual[2];
            int expected_calls[2], expected_count[2];
            unsigned long long expected_events[2][16];
            lengths[0] = 150;
            lengths[1] = variant == 2 ? 151 : 150;
            memset(calls, 0, sizeof(calls));
            memset(event_count, 0, sizeof(event_count));
            for (owner = 0; owner < 2; ++owner) {
                initialize(&expected[owner], owner);
                initialize(&actual[owner], owner);
                expected_calls[owner] = reference_matesw(&opt, &bns, &pac,
                    pes[owner], &anchors[owner], lengths[owner], query[owner],
                    &expected[owner]);
                assert(expected_calls[owner] == calls[owner]);
            }
            memcpy(expected_count, event_count, sizeof(event_count));
            memcpy(expected_events, events, sizeof(events));
            memset(calls, 0, sizeof(calls));
            memset(event_count, 0, sizeof(event_count));
            if (variant == 0) {
                for (owner = 0; owner < 2; ++owner)
                    assert(mem_matesw(&opt, &bns, &pac, pes[owner],
                        &anchors[owner], lengths[owner], query[owner],
                        &actual[owner]) == expected_calls[owner]);
            } else {
                /* Prepare independently to exercise unequal candidate counts. */
                swbwa_matesw_task_t tasks[2];
                int paired, i;
                for (owner = 0; owner < 2; ++owner)
                    swbwa_matesw_prepare(&tasks[owner], &opt, &bns, &pac,
                        pes[owner], &anchors[owner], lengths[owner], queries[owner],
                        &actual[owner], owner);
                paired = tasks[0].candidate_count < tasks[1].candidate_count
                       ? tasks[0].candidate_count : tasks[1].candidate_count;
                for (i = 0; i < paired; ++i) swbwa_matesw_run_pair(tasks, i);
                for (owner = 0; owner < 2; ++owner) {
                    for (i = paired; i < tasks[owner].candidate_count; ++i)
                        swbwa_matesw_run_one(&tasks[owner], i);
                    swbwa_matesw_finish(&tasks[owner]);
                }
            }
            for (owner = 0; owner < 2; ++owner) {
                assert(calls[owner] == expected_calls[owner]);
                assert(event_count[owner] == expected_count[owner]);
                assert(memcmp(events[owner], expected_events[owner],
                              expected_count[owner] * sizeof(events[0][0])) == 0);
                assert(actual[owner].n == expected[owner].n);
                assert(memcmp(actual[owner].a, expected[owner].a,
                              actual[owner].n * sizeof(mem_alnreg_t)) == 0);
                free(actual[owner].a);
                free(expected[owner].a);
            }
        }
    }
    puts("PASS: 1296 direction patterns x serial/dual/equal-and-unequal lengths");
    return 0;
}
'''


def main():
    host = (ROOT / "src/host/bwamem_pair.c").read_text()
    slave = (ROOT / "src/slave/bwamem_pair.c").read_text()
    reference = host[host.index("int mem_matesw("):host.index("\nint mem_pair(")]
    reference = reference.replace("int mem_matesw(", "static int reference_matesw(", 1)
    start = slave.index("typedef struct {\n\tuint8_t *seq;")
    actual = slave[start:slave.index("\nint mem_pair(", start)]
    source = HARNESS.replace("/* RESCUE_IMPLEMENTATIONS */", reference + actual)
    with tempfile.TemporaryDirectory(prefix="swbwa-matesw-") as temp:
        path = Path(temp)
        (path / "test.c").write_text(source)
        command = shlex.split(os.environ.get("CC", "cc")) + [
            "-std=gnu11", "-O1", "-g", "-Wall", "-Wextra", "-Wno-unused-function",
            "-fsanitize=address,undefined", "-fno-omit-frame-pointer", "-no-pie",
            "-I" + str(ROOT / "include"), str(path / "test.c"),
            "-o", str(path / "test")]
        subprocess.run(command, check=True)
        subprocess.run([str(path / "test")], check=True)


if __name__ == "__main__":
    main()
