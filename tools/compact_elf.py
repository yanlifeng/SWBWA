#!/usr/bin/env python3
"""Strip nonallocated metadata after CPE relocation, retaining a debug image."""
import os
import pathlib
import subprocess
import sys

from check_elf_layout import Image, check_pair


def main():
    source, destination = map(pathlib.Path, sys.argv[1:])
    before = Image(source.read_bytes())
    subprocess.check_call([os.environ.get("OBJCOPY", "swobjcopy"), "--strip-all",
                           "--keep-symbol=main", str(source), str(destination)])
    check_pair(before, Image(destination.read_bytes()))
    print("ELF metadata compacted; runtime mappings and allocated bytes unchanged")


if __name__ == "__main__":
    main()
