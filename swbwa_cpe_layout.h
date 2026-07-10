#ifndef SWBWA_CPE_LAYOUT_H
#define SWBWA_CPE_LAYOUT_H

/*
 * CPE ELF layout used by the cross-segment runtime. Regenerate this file with
 * build_cross.sh after changes that affect the linked CPE image.
 */
#define SWBWA_CPE_TEXT_START_ADDRESS 0x00004ffff0410000UL
#define SWBWA_CPE_TEXT_SEGMENT_BYTES 0x00000000002b1a38UL
#define SWBWA_CPE_DATA_START_ADDRESS 0x0000500000004040UL
#define SWBWA_CPE_DATA_SEGMENT_BYTES 0x000000000bc43fc8UL

#endif /* SWBWA_CPE_LAYOUT_H */
