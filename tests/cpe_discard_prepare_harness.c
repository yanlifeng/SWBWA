#include <assert.h>
#include <stdint.h>
#include <string.h>
#include "swbwa_discard_digest.h"

#if !SWBWA_CPE_DISCARD_DIGEST_ACTIVE
#error "Preparation test requires the active discard digest path"
#endif
#if defined(SWBWA_TEST_UNKNOWN_ENDIAN) && defined(__BYTE_ORDER__)
#error "Unknown-endian test must undefine __BYTE_ORDER__"
#endif

int main(void)
{
    const uint64_t lengths[] = {0, 1, 255, 256, UINT64_C(0x0123456789abcdef), UINT64_MAX};
    const int policies[] = {0, 1, -1, 2};
    _Alignas(8) unsigned char bytes[80];
    assert(((uintptr_t)bytes & 7) == 0);
    for (unsigned int offset = 0; offset < 32; ++offset) {
        for (unsigned int n = 0; n < sizeof(lengths) / sizeof(lengths[0]); ++n) {
            for (unsigned int policy = 0; policy < sizeof(policies) / sizeof(policies[0]); ++policy) {
                /* Independent 24-byte wire format, not the production accessors. */
                unsigned char expected[24] = {0};
                memset(bytes, 0xa5, sizeof(bytes));
                for (unsigned int i = 0; i < 8; ++i)
                    expected[i] = (unsigned char)(lengths[n] >> (i * 8));
                expected[16] = policies[policy] ? 1 : 2;
                swbwa_digest_prepare(bytes + offset, lengths[n], policies[policy]);
                assert(memcmp(bytes + offset, expected, sizeof(expected)) == 0);
                for (unsigned int i = 0; i < offset; ++i) assert(bytes[i] == 0xa5);
                for (unsigned int i = offset + sizeof(expected); i < sizeof(bytes); ++i)
                    assert(bytes[i] == 0xa5);
            }
        }
    }
    return 0;
}
