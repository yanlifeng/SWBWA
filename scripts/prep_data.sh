#!/usr/bin/env bash
# Selling-point B experiment: prepare the input data.
#   tiny_*  : original small dataset, for smoke tests and startup overhead
#             measurement (10000 reads).
#   big4_*  : 4x amplification, for performance comparison (so that stage2
#             dominates the total).
# Must run on a compute node (submitted via bsub), not on the login node.
#
# Concurrency hardening (after the P1 incident):
#   1. flock single-instance lock: an instance that cannot grab the lock exits
#      immediately without writing.
#   2. output is first written to a unique temp file in the same directory,
#      then moved with an atomic rename: even with concurrent instances, the
#      file that lands is always a fully-written copy from one instance, never
#      interleaved.
set -euo pipefail

D=${PREP_DATA_DIR:?Set PREP_DATA_DIR to the explicitly approved data directory}
BASE=${WORK_ROOT:?Set WORK_ROOT to the experiment directory}
gfsquota
[[ "${PREP_SPACE_CONFIRMED:-0}" == 1 ]] || {
    echo "Check gfsquota and temporary/output space, then set PREP_SPACE_CONFIRMED=1" >&2
    exit 1
}
for name in tiny_1 tiny_2 big4_SRR7963242_1 big4_SRR7963242_2; do
    [[ ! -e "$D/$name.fastq" ]] || { echo "Preserving $D/$name.fastq" >&2; exit 1; }
done
LOCK="$BASE/logs/prep_data.lock"

# ---- single-instance protection ----
mkdir -p "$BASE/logs"
exec 9>"$LOCK"
if ! flock -n 9; then
    echo "SKIP: 已有另一个 prep_data 实例持有锁 ($LOCK)，本实例退出以免并发写坏数据"
    echo "PREP_RESULT=SKIPPED_DUPLICATE"
    exit 0
fi
echo "已获得单实例锁: $LOCK"
echo "运行于节点: $(hostname)  PID=$$"
echo "开始时间: $(date '+%F %T')"

# Temp files live in the target directory: the same filesystem is required
# for an atomic mv.
TMPPFX="$D/.prep_tmp_$$"
cleanup() { rm -f "$D"/.prep_tmp_"$$"_* 2>/dev/null || true; }
trap cleanup EXIT

echo
echo "=== 输入源 ==="
ls -l "$D/small_SRR7963242_1.fastq" "$D/small_SRR7963242_2.fastq"

# ---- generate tiny ----
echo
echo "=== 生成 tiny（前 40000 行 = 10000 条 read）==="
head -n 40000 "$D/small_SRR7963242_1.fastq" > "${TMPPFX}_tiny_1.fastq"
head -n 40000 "$D/small_SRR7963242_2.fastq" > "${TMPPFX}_tiny_2.fastq"
mv -f "${TMPPFX}_tiny_1.fastq" "$D/tiny_1.fastq"
mv -f "${TMPPFX}_tiny_2.fastq" "$D/tiny_2.fastq"
echo "tiny 完成: $(date '+%F %T')"

# ---- generate big4 (R1/R2 concatenated the same number of times, in the
#      same order, keeping pairs aligned) ----
echo
echo "=== 生成 big4（4 倍放大）==="
for i in 1 2; do
    src="$D/small_SRR7963242_${i}.fastq"
    t0=$(date +%s)
    cat "$src" "$src" "$src" "$src" > "${TMPPFX}_big4_SRR7963242_${i}.fastq"
    mv -f "${TMPPFX}_big4_SRR7963242_${i}.fastq" "$D/big4_SRR7963242_${i}.fastq"
    t1=$(date +%s)
    echo "big4 R${i} 完成: 用时 $((t1 - t0))s, 字节 $(stat -c %s "$D/big4_SRR7963242_${i}.fastq")  ($(date '+%F %T'))"
done

echo
echo "=== 校验表 ==="
printf '%-22s %14s %14s %10s\n' NAME BYTES LINES LINES_MOD4
for f in tiny_1 tiny_2 big4_SRR7963242_1 big4_SRR7963242_2; do
    p="$D/$f.fastq"
    b=$(stat -c %s "$p")
    l=$(wc -l < "$p")
    printf '%-22s %14s %14s %10s\n' "$f" "$b" "$l" "$((l % 4))"
done

echo
echo "=== 断言 ==="
s1=$(stat -c %s "$D/small_SRR7963242_1.fastq")
s2=$(stat -c %s "$D/small_SRR7963242_2.fastq")
b1=$(stat -c %s "$D/big4_SRR7963242_1.fastq")
b2=$(stat -c %s "$D/big4_SRR7963242_2.fastq")
t1=$(stat -c %s "$D/tiny_1.fastq")
t2=$(stat -c %s "$D/tiny_2.fastq")

fail=0
[[ "$b1" == "$b2" ]]           && echo "OK   big4 R1/R2 字节数相等 ($b1)"      || { echo "FAIL big4 R1/R2 字节数不等: $b1 vs $b2"; fail=1; }
[[ "$b1" -eq $((s1 * 4)) ]]    && echo "OK   big4 R1 字节数 = 原始 × 4"        || { echo "FAIL big4 R1 字节数 != 原始×4"; fail=1; }
[[ "$b2" -eq $((s2 * 4)) ]]    && echo "OK   big4 R2 字节数 = 原始 × 4"        || { echo "FAIL big4 R2 字节数 != 原始×4"; fail=1; }
[[ "$t1" == "$t2" ]]           && echo "OK   tiny R1/R2 字节数相等 ($t1)"      || { echo "WARN tiny R1/R2 字节数不等: $t1 vs $t2"; }

# Segment boundary spot check: the first and last 1 MiB of each of big4's
# segments 1..4 must exactly match the corresponding position in the source.
# Note: compute nodes run BusyBox, whose cmp lacks -n/-i; use
# `tail -c +N | head -c W | md5sum` instead (BusyBox tail seeks on regular
# files, measured <10 ms per call).
W=1048576
for i in 1 2; do
    p="$D/big4_SRR7963242_${i}.fastq"
    src="$D/small_SRR7963242_${i}.fastq"
    seg=$(stat -c %s "$src")
    h_head=$(head -c "$W" "$src" | md5sum | awk '{print $1}')
    h_tail=$(tail -c "$W" "$src" | md5sum | awk '{print $1}')
    for k in 0 1 2 3; do
        off=$((k * seg))
        # Bounded dd avoids tail receiving SIGPIPE under set -o pipefail.
        wh=$(dd if="$p" bs=1 skip="$off" count="$W" 2>/dev/null | md5sum | awk '{print $1}')
        we=$(dd if="$p" bs=1 skip=$((off + seg - W)) count="$W" 2>/dev/null | md5sum | awk '{print $1}')
        if [[ "$wh" == "$h_head" && "$we" == "$h_tail" ]]; then
            echo "OK   big4 R${i} 第$((k+1))段 首/尾1MiB 与原始一致 (offset=$off)"
        else
            echo "FAIL big4 R${i} 第$((k+1))段 首/尾1MiB 不一致 (offset=$off)"
            fail=1
        fi
    done
done

echo
echo "结束时间: $(date '+%F %T')"
if [[ "$fail" == 0 ]]; then
    echo "PREP_RESULT=PASS"
else
    echo "PREP_RESULT=FAIL"
    exit 1
fi
