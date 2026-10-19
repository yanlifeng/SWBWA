#!/usr/bin/env python3
"""Test the real utils.c fatal functions with a terminating MPI-abort stub.

No MPI library, launcher, SDK or remote execution is required. MPI stderr is
fully buffered and the stub uses _Exit without flushing, so missing flushes
are observable. Non-MPI core paths must retain SIGABRT; core dumps are disabled.
"""
import os
from pathlib import Path
import shlex
import signal
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]

HARNESS = r"""
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <unistd.h>
#include "utils.c"

static void record_local_exit(void) { puts("LOCAL_EXIT"); }

void swbwa_mpi_abort(const char *message)
{
    if (message != NULL) _Exit(90);
    if (write(STDOUT_FILENO, "MPI_ABORT\n", 10) != 10) _Exit(91);
    _Exit(86);
}

int main(int argc, char **argv)
{
    const struct rlimit no_core = {0, 0};
    static char log_buffer[65536];
    assert(argc == 3);
    assert(setrlimit(RLIMIT_CORE, &no_core) == 0);
    assert(atexit(record_local_exit) == 0);
#if SWBWA_USE_MPI
    assert(setvbuf(stderr, log_buffer, _IOFBF, sizeof(log_buffer)) == 0);
#else
    (void)log_buffer;
#endif
    switch (atoi(argv[1])) {
        case 0: err_fatal("probe", "number=%d text=%s %%", 17, argv[2]); break;
        case 1: err_fatal_core("probe", "number=%d text=%s %%", 17, argv[2]); break;
        case 2: _err_fatal_simple("probe", argv[2]); break;
        case 3: _err_fatal_simple_core("probe", argv[2]); break;
        default: return 92;
    }
    return 93;
}
"""


class UtilsFatalTest(unittest.TestCase):
    def test_logs_flush_and_termination(self):
        compiler = shlex.split(os.environ.get("CC", "cc"))
        with tempfile.TemporaryDirectory(prefix="swbwa-utils-fatal-") as temporary:
            work = Path(temporary)
            harness = work / "fatal.c"
            harness.write_text(HARNESS)
            for mpi in (0, 1):
                with self.subTest(build_mpi=mpi):
                    executable = work / f"fatal-{mpi}"
                    command = compiler + [
                        "-std=gnu11", "-O1", "-g", "-Wall", "-Wextra",
                        "-Werror=implicit-function-declaration", "-Werror=return-type",
                        "-ffunction-sections", "-fdata-sections",
                        "-I" + str(ROOT / "include"), f"-DSWBWA_USE_MPI={mpi}",
                        "-I" + str(ROOT / "src/host"),
                        "-DSWBWA_ENABLE_HOST_MALLOC_WRAPPER=0",
                        str(harness),
                        "-Wl,--gc-sections", "-o", str(executable),
                    ]
                    build = subprocess.run(command, text=True, capture_output=True, timeout=60)
                    self.assertEqual(build.returncode, 0, build.stdout + build.stderr)
                    for function in range(4):
                        for message in ("detail 100% %s", "x" * 16000):
                            with self.subTest(mpi=mpi, function=function, length=len(message)):
                                result = subprocess.run(
                                    [str(executable), str(function), message],
                                    text=True, capture_output=True, timeout=10)
                                core = function in (1, 3)
                                expected_code = 86 if mpi else -signal.SIGABRT if core else 1
                                self.assertEqual(result.returncode, expected_code, result.stderr)
                                text = f"number=17 text={message} %" if function < 2 else message
                                self.assertEqual(result.stderr, "[probe] " + text + (" Abort!\n" if core else "\n"))
                                self.assertEqual(result.stdout, "MPI_ABORT\n" if mpi else "" if core else "LOCAL_EXIT\n")


if __name__ == "__main__":
    unittest.main(verbosity=2)
