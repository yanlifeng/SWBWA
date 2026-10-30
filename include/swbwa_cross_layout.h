#ifndef SWBWA_CROSS_LAYOUT_H
#define SWBWA_CROSS_LAYOUT_H

#include <stddef.h>
#include <stdint.h>

/* Short-PC instruction addresses are word-addressed (24 bits x 4 bytes). */
#define SWBWA_CROSS_PC_BYTES ((uintptr_t)1 << 26)

static inline int swbwa_cross_text_fits(uintptr_t address, size_t bytes)
{
    return bytes != 0 && bytes <= SWBWA_CROSS_PC_BYTES &&
           (address & (SWBWA_CROSS_PC_BYTES - 1)) <=
               SWBWA_CROSS_PC_BYTES - bytes;
}

static inline uintptr_t swbwa_cross_relocation_padding(uintptr_t allocation,
                                                       uintptr_t linked_text)
{
    return (linked_text - allocation) & (SWBWA_CROSS_PC_BYTES - 1);
}

#endif
