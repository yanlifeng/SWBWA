#!/usr/bin/env python3
"""Spot-check big4_* sizes and segment ends (not a full byte-wise check)."""
import sys

if len(sys.argv) != 2:
    sys.exit("usage: verify_big4.py DATA_DIRECTORY")
D = sys.argv[1]
W = 1 << 20  # 1 MiB spot-check window

fail = 0
for i in (1, 2):
    big = f"{D}/big4_SRR7963242_{i}.fastq"
    src = f"{D}/small_SRR7963242_{i}.fastq"

    with open(src, "rb") as f:
        f.seek(0, 2)
        seg = f.tell()
        f.seek(0)
        head = f.read(W)

    with open(big, "rb") as f:
        f.seek(0, 2)
        bsz = f.tell()

    print(f"--- R{i} ---")
    print(f"    small={seg}  big4={bsz}  倍数={bsz / seg:.6f}")

    if bsz != seg * 4:
        print(f"    FAIL 字节数不是 4 倍")
        fail = 1

    with open(big, "rb") as f:
        for k in range(4):
            f.seek(k * seg)
            w = f.read(W)
            if w == head:
                print(f"    OK   第{k+1}段起始 1MiB 与原始一致 (offset={k * seg})")
            else:
                print(f"    FAIL 第{k+1}段起始 1MiB 不匹配 (offset={k * seg})")
                fail = 1

    # The last 1 MiB of each segment must also match, to catch a wrong-length
    # but correct-beginning misplacement.
    with open(src, "rb") as f:
        f.seek(-W, 2)
        tail = f.read(W)
    with open(big, "rb") as f:
        for k in range(4):
            end = (k + 1) * seg
            f.seek(end - W)
            w = f.read(W)
            if w == tail:
                print(f"    OK   第{k+1}段末尾 1MiB 与原始一致")
            else:
                print(f"    FAIL 第{k+1}段末尾 1MiB 不匹配")
                fail = 1

print()
print("BOUNDARY_CHECK=" + ("PASS" if fail == 0 else "FAIL"))
sys.exit(fail)
