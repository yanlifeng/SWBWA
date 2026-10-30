#include "swbwa_cross_layout.h"

#include <assert.h>
#include <stdio.h>

int main(void)
{
    const uintptr_t linked = UINT64_C(0x4ffff0410000);
    const uintptr_t bad = UINT64_C(0x5803bffe2000);
    const size_t text = 0x3d1861;
    uintptr_t state = 1;
    unsigned i;

    assert(!swbwa_cross_text_fits(bad, text));
    assert(swbwa_cross_text_fits(UINT64_C(0x5803c40c0000), text));
    assert(swbwa_cross_text_fits(UINT64_C(0x5803f6f84000), text));
    assert(!swbwa_cross_text_fits(0, 0));
    assert(swbwa_cross_text_fits(SWBWA_CROSS_PC_BYTES - text, text));
    assert(!swbwa_cross_text_fits(SWBWA_CROSS_PC_BYTES - text + 1, text));
    assert(!swbwa_cross_text_fits(0, SIZE_MAX));
    for (i = 0; i < 100000; ++i) {
        uintptr_t raw, padding, placed;
        state = state * UINT64_C(6364136223846793005) + 1;
        raw = UINT64_C(0x580000000000) + (state & UINT64_C(0xffffffff));
        padding = swbwa_cross_relocation_padding(raw, linked);
        placed = raw + padding;
        assert(padding < SWBWA_CROSS_PC_BYTES);
        assert(((placed - linked) & (SWBWA_CROSS_PC_BYTES - 1)) == 0);
        assert(placed >= raw && placed + text <= raw + text + SWBWA_CROSS_PC_BYTES - 1);
        assert(swbwa_cross_text_fits(placed, text));
    }
    puts("CROSS_LAYOUT_PASS");
    return 0;
}
