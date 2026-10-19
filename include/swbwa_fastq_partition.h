#ifndef SWBWA_FASTQ_PARTITION_H
#define SWBWA_FASTQ_PARTITION_H

#include <stddef.h>
#include <stdint.h>

/* Boundaries describe equal record ordinals, not necessarily equal byte offsets. */
typedef struct {
    int64_t begin[2], end[2];
    int64_t first_record, records;
} swbwa_fastq_partition_t;

typedef struct {
    const char *name, *seq, *qual;
    size_t name_size, length;
} swbwa_fastq_record_view;

/* Full validation belongs in the CPE formatting pass, not in the MPE reader.
 * Return 1 for a record, 0 at the end, -1 for invalid four-line FASTQ. */
static inline int swbwa_fastq_record_next(const char *data, size_t size,
                                         size_t *position,
                                         swbwa_fastq_record_view *record)
{
    size_t begin = *position, p = begin, ends[4], name_end;
    unsigned line;
    if (begin == size) return 0;
    if (begin > size || data == NULL || data[begin] != '@') return -1;
    for (line = 0; line < 4; ++line) {
        while (p < size && data[p] != '\n') {
            if (data[p] == '\0' || data[p] == '\r' || p - begin >= (1U << 20))
                return -1;
            ++p;
        }
        if (p == size || p - begin >= (1U << 20)) return -1;
        ends[line] = ++p;
    }
    if (ends[0] - begin < 3 || ends[2] - ends[1] < 2 ||
        data[ends[1]] != '+' || ends[1] - ends[0] != ends[3] - ends[2])
        return -1;
    name_end = begin + 1;
    while (name_end < ends[0] - 1 && data[name_end] != ' ' && data[name_end] != '\t')
        ++name_end;
    if (name_end == begin + 1) return -1;
    record->name = data + begin + 1;
    record->name_size = name_end - begin - 1;
    record->seq = data + ends[0];
    record->qual = data + ends[2];
    record->length = ends[1] - ends[0] - 1;
    *position = p;
    return 1;
}
#endif
