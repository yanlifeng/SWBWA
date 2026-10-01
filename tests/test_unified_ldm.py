"""Reuse the real allocator/DP differential harness for the unified policy."""
import os
from pathlib import Path
import runpy
import shutil
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    # Generate only test scaffolding; the allocator under test is not rewritten.
    with tempfile.TemporaryDirectory(prefix="swbwa-unified-test-") as name:
        root = Path(name)
        for directory in ("src", "include", "tests"):
            shutil.copytree(ROOT / directory, root / directory)
        test = root / "tests/test_ldm_allocator.c"
        text = test.read_text()
        text = text.replace("set_big_buffer(pool, per_core);", """
        set_big_buffer(pool, per_core);
        swbwa_ldm_policy_t policy = {0};
        for (unsigned s = 1; s < SWBWA_LDM_SITE_COUNT; ++s) {
            policy.cap[s] = swbwa_ldm_alloc_tier(s) == 2
                          ? SWBWA_CPE_LDM_BYTES / 8 : SWBWA_CPE_LDM_BYTES / 2;
            policy.reserve[s] = swbwa_ldm_alloc_tier(s) == 2
                              ? SWBWA_CPE_LDM_BYTES / 4 : 0;
        }
        policy.profile = getenv("TEST_PROFILE") != NULL;
        policy.fast_realloc = getenv("TEST_INPLACE") != NULL;
        policy.heap_cache = 1;
        policy.scratch_hints = 255;
        if (SWBWA_CPE_LDM_BYTES == 65536) {
            policy.heap_cache = 0;
            policy.cap[SWBWA_LDM_SITE_CHAIN] = 65536;
            policy.reserve[SWBWA_LDM_SITE_CHAIN] = 0;
        }
        swbwa_ldm_set_policy(&policy);
        """)
        text = text.replace(
            "assert(stats.site[SWBWA_LDM_SITE_GLOBAL_DP].requests == 2400);",
            "assert(stats.site[SWBWA_LDM_SITE_GLOBAL_DP].requests == "
            "(policy.profile ? 2400 : 0));")
        text = text.replace("if (!refuse) assert(stats.site[SWBWA_LDM_SITE_GLOBAL_DP].placed > 0);",
                            "if (!refuse && policy.profile) assert(stats.site[SWBWA_LDM_SITE_GLOBAL_DP].placed > 0);")
        text = text.replace("is_ldm(low[i]) == (i < 6)",
                            "is_ldm(low[i]) == (i < (policy.heap_cache ? 5 : 6))")
        text = text.replace("swbwa_auto_free(NULL);", """
            swbwa_auto_free(NULL);
            unsigned char *invalid_site = swbwa_auto_malloc(32, SWBWA_LDM_SITE_CHAIN);
            memset(invalid_site, 0x3b, 32);
            invalid_site = swbwa_auto_realloc(invalid_site, 128, (unsigned)-1);
            for (int i = 0; i < 32; ++i) assert(invalid_site[i] == 0x3b);
            swbwa_auto_free(invalid_site);
            if (SWBWA_CPE_LDM_BYTES == 65536) {
                unsigned char *full = swbwa_auto_malloc(65536, SWBWA_LDM_SITE_CHAIN);
                assert(is_ldm(full) == !refuse);
                memset(full, 0x67, 65536);
                full = swbwa_auto_realloc(full, 65537, SWBWA_LDM_SITE_CHAIN);
                assert(!is_ldm(full));
                for (unsigned i = 0; i < 65536; ++i) assert(full[i] == 0x67);
                swbwa_auto_free(full);
            }
        """)
        text = text.replace("void *cached_query = swbwa_test_cached_query();", """
            swbwa_ldm_policy_t saved_policy = policy;
            /* Policy is set before begin; exercise hints in a fresh epoch. */
            swbwa_ldm_allocator_end();
            policy.scratch_hints = 255;
            swbwa_ldm_set_policy(&policy);
            swbwa_ldm_allocator_begin();
            unsigned char *smem = swbwa_ldm_alloc(16384, 1);
            assert((smem != NULL) == !refuse);
            if (smem) memset(smem, 0x4d, 16384);
            void *query_dp = swbwa_ldm_alloc(4096, 6);
            assert((query_dp != NULL) == !refuse);
            if (query_dp) memset(query_dp, 0x62, 4096);
            assert(swbwa_ldm_alloc(SWBWA_CPE_LDM_BYTES + 1, 5) == NULL);
            if (smem) {
                for (unsigned j = 0; j < 16384; ++j) assert(smem[j] == 0x4d);
                swbwa_ldm_release(smem, 16384);
            }
            if (query_dp) swbwa_ldm_release(query_dp, 4096);
            for (int s = 1; s <= 8; ++s) {
                void *scratch = swbwa_ldm_alloc(73, s);
                assert((scratch != NULL) == !refuse);
                if (scratch) {
                    memset(scratch, s, 73);
                    swbwa_ldm_release(scratch, 73);
                }
            }
            swbwa_ldm_allocator_end();
            policy = saved_policy;
            swbwa_ldm_set_policy(&policy);
            swbwa_ldm_allocator_begin();
            void *cached_query = swbwa_test_cached_query();
        """)
        test.write_text(text)
        driver = root / "tests/test_ldm_allocator.py"
        text = driver.read_text()
        a = text.index("        for mode, cap, manual, budget in (")
        b = text.index("        ):", a) + len("        ):")
        text = text[:a] + "        for mode, cap, manual, budget in ((4, 32768, 0, 40960), (4, 57344, 0, 73728), (4, 65536, 0, 81920)):" + text[b:]
        text = text.replace("#define SWBWA_CPE_LDM_MODE 1", "#define SWBWA_CPE_LDM_MODE 3")
        driver.write_text(text)
        for profile in (False, True):
            for inplace in (False, True):
                for key, on in (("TEST_PROFILE", profile), ("TEST_INPLACE", inplace)):
                    if on:
                        os.environ[key] = "1"
                    else:
                        os.environ.pop(key, None)
                print(f"unified profile={profile} inplace={inplace}", flush=True)
                runpy.run_path(str(driver), run_name="__main__")


if __name__ == "__main__":
    main()
