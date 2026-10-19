#ifndef SWBWA_DECODED_VIEW_H
#define SWBWA_DECODED_VIEW_H
#include "line_index.h"
#include <memory>

// Byte ownership may outlive the reader and its next decode batch.
struct SwbwaDecodedView {
    std::shared_ptr<unsigned char> owner;
    unsigned char *data = nullptr;
    size_t size = 0;
    SwbwaLineIndex lines;
};
#endif
