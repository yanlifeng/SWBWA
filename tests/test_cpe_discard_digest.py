#!/usr/bin/env python3
"""Native regression of the current checkout's CPE discard digest path.

Extract the actual host slice loop, SE/PE copy worker and Stage3 output loop;
link the actual output writer. An independent scalar reference is in the C
harness. No Git history, archived source, target SDK, or input data is needed.
An independent preparation harness checks aligned/unaligned metadata and the
unknown-endian fallback against the current header's actual prepare function.
The full alignment engine and target memory-visibility contract are not tested.
Before digest production code is present this test explicitly skips.
"""
import argparse
import itertools
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
HASH = re.compile(r"calls=(\d+) bytes=(\d+) sum=0x([0-9a-f]+) xor=0x([0-9a-f]+) enabled=([01]) hash_prefix_bytes=0")


def between(source, start, end):
    if source.count(start) != 1:
        raise ValueError(f"source anchor changed or is ambiguous: {start!r}")
    begin = source.index(start)
    return source[begin:source.index(end, begin)]


def expected_sam():
    records = []
    for paired in (False, True):
        for batch in range(1, 9):
            for i in range(12 if paired else 11):
                length = 0 if i == 0 else (i * 31 + batch * 7) % 200
                records.append(bytes(10 if j % 43 == 42 else 33 + (i + j + batch) % 90 for j in range(length)))
    return b"".join(records)


def run(command, **kwargs):
    result = subprocess.run(command, text=True, capture_output=True, **kwargs)
    if result.returncode:
        raise RuntimeError(f"{shlex.join(map(str, command))}\n{result.stdout}\n{result.stderr}")
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=ROOT, help="another checkout to test; default is this checkout")
    args = parser.parse_args(argv)
    source_root = args.source_root.resolve()
    if not (source_root / "include/swbwa_discard_digest.h").is_file():
        print("SKIP: CPE discard digest production header is not present in", source_root)
        return 0
    cc = shlex.split(os.environ.get("CC", "cc"))
    common = ["-std=gnu11", "-O2", "-g", "-Wall", "-Wextra", "-Wno-unused-function",
              "-Wno-unused-parameter", "-Wno-unused-variable", "-Wno-clobbered",
              "-DSWBWA_USE_MPI=0", "-DSWBWA_EXEC_MODE=SWBWA_EXEC_CGS_CROSS",
              "-DSWBWA_CPE_ALLOC_MODE=SWBWA_CPE_ALLOC_POOL", "-DSWBWA_DISCARD_HASH_BYTES=0",
              "-DSWBWA_ENABLE_HOST_MALLOC_WRAPPER=0", "-DSWBWA_CPE_FORMAT_BUFFER_BYTES=4096"]
    sanitizer = ["-fsanitize=address,undefined", "-fno-sanitize-recover=all",
                 "-fno-omit-frame-pointer", "-fno-pie", "-no-pie"]
    with tempfile.TemporaryDirectory(prefix="swbwa-digest-native-") as tmp:
        work = Path(tmp)
        slave = (source_root / "src/slave/bwamem.c").read_text()
        host = (source_root / "src/host/bwamem.c").read_text()
        output = (source_root / "src/host/fastmap.c").read_text()
        (work / "actual_worker.inc").write_text(between(slave, "void worker12_fast(", "\nvoid mem_process_seqs("))
        merge = host[host.index("void mem_process_seqs_merge2("):]
        (work / "actual_slices.inc").write_text(between(merge, "    size_t now_f_block_pos = 0;", "    t_work1_4 +="))
        (work / "actual_output_loop.inc").write_text(between(output, "\t\tfor (i = 0; i < data->n_seqs; ++i)", "\t\tt_step3_1 +="))
        includes = ["-I"+str(work), "-I"+str(source_root / "include"), "-I"+str(source_root / "src/host")]
        (work / "athread.h").write_text("#pragma once\n#define SLAVE_FUN(x) slave_##x\n#define __uncached\n"
            "void athread_init(void); void athread_init_cgs(void);\n"
            "void athread_join(void); void athread_join_cgs(void);\n"
            "void __real_athread_spawn(void *, void *, int);\n"
            "void __real_athread_spawn_cgs(void *, void *, int);\n")
        expected = expected_sam()
        assert len(expected) == 17300
        for digest, reuse, terminators, mode in itertools.product((0, 1), (0, 1), (0, 1), ("DISCARD", "SPLIT")):
            defines = [f"-DSWBWA_CPE_DISCARD_DIGEST={digest}", f"-DSWBWA_HOST_PREP_REUSE={reuse}",
                       f"-DSWBWA_HOST_PREP_CPE_TERMINATORS={terminators}", f"-DSWBWA_OUTPUT_MODE=SWBWA_OUTPUT_{mode}"]
            executable = work / "test"
            run(cc + common + sanitizer + defines + includes +
                [str(ROOT / "tests/cpe_discard_digest_harness.c"), str(source_root / "src/host/swbwa_output.c"), "-o", str(executable)])
            for policy in (None, "", "0", "1") if mode == "DISCARD" else ("1",):
                env = dict(os.environ, ASAN_OPTIONS="detect_leaks=1:halt_on_error=1")
                env.pop("SWBWA_DISCARD_HASH", None)
                if policy is not None:
                    env["SWBWA_DISCARD_HASH"] = policy
                sam = work / "output.sam"
                result = run([str(executable), str(sam)], env=env)
                if mode == "DISCARD":
                    actual, reference = HASH.search(result.stderr), HASH.search(result.stdout)
                    assert actual and reference and actual.groups() == reference.groups(), result.stderr
                    assert not sam.exists()
                    assert ("hash_scope=generated_sam" in result.stderr) == bool(digest)
                else:
                    assert sam.read_bytes() == expected
                    sam.unlink()
            if mode == "DISCARD":
                invalid = subprocess.run([str(executable), str(work / "invalid.sam")],
                    env=dict(os.environ, SWBWA_DISCARD_HASH="invalid"), text=True, capture_output=True)
                assert invalid.returncode == 9 and not (work / "invalid.sam").exists()
            for name in ("bwamem.c", "fastmap.c"):
                run(cc + common + defines + includes + ["-D_GNU_SOURCE", "-fsyntax-only", str(source_root / "src/host" / name)])
            print(f"PASS digest={digest} reuse={reuse} terminators={terminators} mode={mode}: SE/PE, scalar reference, actual output loop, native host syntax")
        guards = [([], True), (["-DSWBWA_CPE_DISCARD_DIGEST=1"], True),
                  (["-DSWBWA_CPE_DISCARD_DIGEST=2"], False),
                  (["-DSWBWA_CPE_DISCARD_DIGEST=1", "-USWBWA_DISCARD_HASH_BYTES", "-DSWBWA_DISCARD_HASH_BYTES=64"], False),
                  (["-DSWBWA_CPE_DISCARD_DIGEST=1", "-USWBWA_USE_MPI", "-DSWBWA_USE_MPI=1"], False),
                  (["-DSWBWA_CPE_DISCARD_DIGEST=1", "-USWBWA_EXEC_MODE", "-DSWBWA_EXEC_MODE=SWBWA_EXEC_SINGLE_CG"], False),
                  (["-DSWBWA_CPE_DISCARD_DIGEST=1", "-USWBWA_CPE_ALLOC_MODE", "-DSWBWA_CPE_ALLOC_MODE=SWBWA_CPE_ALLOC_SYSTEM"], False)]
        for flags, success in guards:
            result = subprocess.run(cc + common + includes + ["-DSWBWA_OUTPUT_MODE=SWBWA_OUTPUT_DISCARD", *flags,
                "-fsyntax-only", "-x", "c", "-"], input='#include "swbwa_discard_digest.h"\n', text=True, capture_output=True)
            assert (result.returncode == 0) == success, result.stderr
        run(cc + common + includes + ["-DSWBWA_OUTPUT_MODE=SWBWA_OUTPUT_DISCARD", "-fsyntax-only", "-x", "c", "-"],
            input='#include "swbwa_discard_digest.h"\n_Static_assert(SWBWA_CPE_DISCARD_DIGEST == 0, "standalone header default must be off");\n')
        print("PASS standalone-header default-off, FULL-only, execution/allocator guards")
        for unknown_endian in (False, True):
            flags = ["-DSWBWA_CPE_DISCARD_DIGEST=1", "-DSWBWA_OUTPUT_MODE=SWBWA_OUTPUT_DISCARD"]
            if unknown_endian:
                flags += ["-U__BYTE_ORDER__", "-DSWBWA_TEST_UNKNOWN_ENDIAN=1"]
            executable = work / "prepare"
            run(cc + common + sanitizer + includes + flags +
                [str(ROOT / "tests/cpe_discard_prepare_harness.c"), "-o", str(executable)])
            run([str(executable)], env=dict(os.environ, ASAN_OPTIONS="detect_leaks=1:halt_on_error=1"))
            print(f"PASS actual prepare: offsets 0..31, 6 lengths x 4 policies, byte format/canaries; unknown_endian={int(unknown_endian)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
