"""Test actual host policy parsing and shared C/Python production defaults."""
import os
from pathlib import Path
import subprocess
import tempfile
import shlex

from test_ldm_policy import p

ROOT = Path(__file__).resolve().parents[1]


def main():
    source = (ROOT / "src/host/bwamem.c").read_text()
    start = source.index("static unsigned ldm_policy_flag(")
    stop = source.index("\n#endif", start)
    harness = '''
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "swbwa_ldm_alloc.h"
int bwa_verbose = 3;
''' + source[start:stop] + '''
int main(void) {
    swbwa_ldm_policy_t policy;
    ldm_policy_init(&policy);
    for (int i=0; i<SWBWA_LDM_SITE_COUNT; ++i)
        printf("%s%u", i ? "," : "", policy.cap[i]);
    puts("");
    for (int i=0; i<SWBWA_LDM_SITE_COUNT; ++i)
        printf("%s%u", i ? "," : "", policy.reserve[i]);
    puts("");
    printf("%u %u %u %u %u\\n", policy.profile, policy.fast_realloc,
           policy.scratch_hints, policy.heap_cache, policy.pool_bitmap);
    return 0;
}
'''
    env = {k: v for k, v in os.environ.items() if not k.startswith("SWBWA_LDM_")}
    with tempfile.TemporaryDirectory(prefix="ldm-config-") as temp:
        temp = Path(temp)
        code, exe = temp / "policy.c", temp / "policy"
        code.write_text(harness)
        subprocess.run(shlex.split(os.environ.get("CC", "cc")) + [
            "-std=c11", "-O2", "-Wall", "-Wextra", "-Werror",
            "-I", str(ROOT / "include"), "-DSWBWA_EXEC_MODE=3",
            "-DSWBWA_CPE_ALLOC_MODE=2", "-DSWBWA_USE_MPI=0",
            str(code), "-o", str(exe)], check=True)
        policy = p.Policy()
        expected = [",".join(map(str, policy.caps)), ",".join(map(str, policy.reserves)),
                    "0 1 255 1 2"]
        for overrides in ({}, policy.env()):
            out = subprocess.run([str(exe)], env=dict(env, **overrides),
                                 capture_output=True, text=True, check=True)
            assert out.stdout.splitlines() == expected, out.stdout
        for name, value in (("SWBWA_LDM_CAPS", "1,0,0,0,0,0,0,0,0,0,0"),
                            ("SWBWA_LDM_CAPS", "0,1"), ("SWBWA_LDM_RESERVES", ""),
                            ("SWBWA_LDM_SCRATCH_HINTS", "256"),
                            ("SWBWA_LDM_POOL_BITMAP", "3"),
                            ("SWBWA_LDM_PROFILE", "yes"),
                            ("SWBWA_LDM_GROW_IN_PLACE", "2"),
                            ("SWBWA_LDM_POOL_CACHE", "")):
            result = subprocess.run([str(exe)], env=dict(env, **{name: value}),
                                    capture_output=True)
            assert result.returncode != 0, (name, value)
    print("PASS production defaults, offline policy agreement, invalid env rejection")


if __name__ == "__main__":
    main()
