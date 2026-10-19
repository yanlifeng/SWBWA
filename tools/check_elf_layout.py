#!/usr/bin/env python3
"""Read-only, conservative ELF64 layout check before/after metadata stripping.

File offsets may move; runtime segment descriptors and allocated sections must
not. This is not a CPE relocation, private-CSR mapping, or execution validator.
Only standard little-endian ELF64 is supported; unfamiliar layouts fail closed.
No symbol table, target SDK, or third-party Python modules are required.
"""
import argparse
from collections import Counter, namedtuple
from pathlib import Path
import re
import struct
import sys

from xlink import ELF64_Hdr, ELF64_Shdr, ELF64_SHDR_FMT

HEADER = struct.Struct("<16sHHIQQQIHHHHHH")
PHDR = struct.Struct("<IIQQQQQQ")
SHDR = struct.Struct("<" + ELF64_SHDR_FMT)
Program = namedtuple("Program", "p_type p_flags p_offset p_vaddr p_paddr p_filesz p_memsz p_align")
PT_LOAD, PT_PHDR, PT_TLS = 1, 6, 7
SHT_NOBITS, SHF_ALLOC, SHF_TLS = 8, 2, 0x400


def require(condition, message):
    if not condition:
        raise ValueError(message)


def contains(start, size, inner, length):
    return start <= inner and length <= size and inner - start <= size - length


def segment_key(segment):
    return tuple(value for field, value in zip(Program._fields, segment) if field != "p_offset")


class Image:
    def __init__(self, data):
        self.data = data
        require(len(data) >= HEADER.size, "truncated ELF header")
        self.header = h = ELF64_Hdr(*HEADER.unpack_from(data))
        require(h.e_ident[:7] == b"\x7fELF\x02\x01\x01", "expected little-endian ELF64 version 1")
        require(h.e_type in (2, 3) and h.e_version == 1, "expected executable/shared ELF version 1")
        require(h.e_ehsize == HEADER.size and h.e_phentsize == PHDR.size and
                h.e_shentsize == SHDR.size, "unexpected ELF table entry size")
        require(0 < h.e_phnum < 0xffff and 0 < h.e_shnum and
                0 < h.e_shstrndx < h.e_shnum, "missing/extended ELF numbering is unsupported")
        self.programs = [Program(*p) for p in PHDR.iter_unpack(
            self.bytes(h.e_phoff, h.e_phnum * PHDR.size))]
        sections = [ELF64_Shdr(*s) for s in SHDR.iter_unpack(
            self.bytes(h.e_shoff, h.e_shnum * SHDR.size))]
        names_header = sections[h.e_shstrndx]
        require(names_header.sh_type == 3, "invalid section-name string table")
        names = self.bytes(names_header.sh_offset, names_header.sh_size)
        self.sections = {}
        for section in sections[1:]:
            require(section.sh_name < len(names), "invalid section-name offset")
            end = names.find(b"\0", section.sh_name)
            require(end >= 0, "unterminated section name")
            name = names[section.sh_name:end].decode("ascii")
            require(name not in self.sections, "duplicate section name: " + name)
            self.sections[name] = section
            if section.sh_type != SHT_NOBITS:
                self.bytes(section.sh_offset, section.sh_size)

    def bytes(self, offset, size):
        require(contains(0, len(self.data), offset, size), "ELF file range is out of bounds")
        return self.data[offset:offset + size]

    def mapping(self, section, kind):
        result = []
        for p in self.programs:
            if p.p_type != kind or not contains(p.p_vaddr, p.p_memsz, section.sh_addr, section.sh_size):
                continue
            if section.sh_type == SHT_NOBITS:
                if section.sh_addr < p.p_vaddr + p.p_filesz:
                    continue
            elif not (contains(p.p_offset, p.p_filesz, section.sh_offset, section.sh_size) and
                      section.sh_addr - p.p_vaddr == section.sh_offset - p.p_offset):
                continue
            require(not section.sh_flags & 1 or p.p_flags & 2 or kind == PT_TLS,
                    "writable section mapped by non-writable LOAD")
            require(not section.sh_flags & 4 or p.p_flags & 1,
                    "executable section mapped by non-executable segment")
            result.append(segment_key(p))
        return sorted(result)

    def validate(self):
        loads = [p for p in self.programs if p.p_type == PT_LOAD]
        require(loads, "no PT_LOAD segments")
        require([p.p_vaddr for p in loads] == sorted(p.p_vaddr for p in loads),
                "PT_LOAD virtual addresses are not ordered")
        for p in self.programs:
            self.bytes(p.p_offset, p.p_filesz)
            if p.p_type not in (PT_LOAD, PT_TLS):
                continue
            require(p.p_filesz <= p.p_memsz, "segment filesz exceeds memsz")
            require(p.p_vaddr + p.p_memsz <= 1 << 64, "segment virtual address overflow")
            require(p.p_align in (0, 1) or p.p_align & (p.p_align - 1) == 0,
                    "segment alignment is not a power of two")
            if p.p_type == PT_LOAD and p.p_align > 1:
                require(p.p_vaddr % p.p_align == p.p_offset % p.p_align,
                        "PT_LOAD offset/virtual-address alignment mismatch")
        h = self.header
        for p in self.programs:
            if p.p_type == PT_PHDR:
                require(p.p_offset == h.e_phoff and
                        p.p_filesz == p.p_memsz == h.e_phnum * h.e_phentsize,
                        "PT_PHDR does not describe the actual program-header table")
                require(any(contains(q.p_offset, q.p_filesz, p.p_offset, p.p_filesz) and
                            contains(q.p_vaddr, q.p_memsz, p.p_vaddr, p.p_memsz) and
                            p.p_vaddr - q.p_vaddr == p.p_offset - q.p_offset for q in loads),
                        "PT_PHDR is not mapped by PT_LOAD")
        if h.e_entry:
            require(any(p.p_flags & 1 and contains(p.p_vaddr, p.p_filesz, h.e_entry, 1)
                        for p in loads), "ELF entry is outside executable file-backed LOAD")
        mappings = {}
        for name, s in self.sections.items():
            if not s.sh_flags & SHF_ALLOC or not s.sh_size:
                continue
            require(s.sh_addr + s.sh_size <= 1 << 64, "section virtual address overflow: " + name)
            kinds = [PT_TLS] if s.sh_flags & SHF_TLS else []
            # TLS NOBITS has template addresses, not ordinary process BSS.
            if not (s.sh_flags & SHF_TLS and s.sh_type == SHT_NOBITS):
                kinds.append(PT_LOAD)
            mappings[name] = []
            for kind in kinds:
                mapping = self.mapping(s, kind)
                require(mapping, "allocated section lacks valid " +
                        ("TLS" if kind == PT_TLS else "LOAD") + " mapping: " + name)
                mappings[name].append((kind, mapping))
        return mappings


def check_pair(before, after):
    old_map, new_map = before.validate(), after.validate()
    fields = ("e_type", "e_machine", "e_version", "e_entry", "e_flags")
    require(before.header.e_ident[:9] == after.header.e_ident[:9] and
            all(getattr(before.header, f) == getattr(after.header, f) for f in fields),
            "ELF identity, ABI, entry, or flags changed")
    # NULL/NOTE describe metadata; PHDR is checked against each image's own table.
    runtime = lambda image: Counter(segment_key(p) for p in image.programs if p.p_type not in (0, 4, PT_PHDR))
    require(runtime(before) == runtime(after),
            "runtime PHDR changed (type/flags/vaddr/paddr/filesz/memsz/align)")
    allocated = lambda image: {n: s for n, s in image.sections.items() if s.sh_flags & SHF_ALLOC}
    old, new = allocated(before), allocated(after)
    require(old.keys() == new.keys(), "allocated section set changed")
    fields = ("sh_type", "sh_flags", "sh_addr", "sh_size", "sh_addralign", "sh_entsize")
    for name, section in old.items():
        other = new[name]
        require(all(getattr(section, f) == getattr(other, f) for f in fields),
                "allocated section metadata changed: " + name)
        if section.sh_type != SHT_NOBITS:
            require(before.bytes(section.sh_offset, section.sh_size) ==
                    after.bytes(other.sh_offset, other.sh_size), "allocated section bytes changed: " + name)
    require(old_map == new_map, "allocated section-to-segment mapping changed")


def check_cross_header(image, text):
    values = []
    for name in ("TEXT_START_ADDRESS", "TEXT_SEGMENT_BYTES", "DATA_START_ADDRESS", "DATA_SEGMENT_BYTES"):
        matches = re.findall(r"^[ \t]*#[ \t]*define[ \t]+SWBWA_CPE_" + name +
                             r"[ \t]+(0[xX][0-9a-fA-F]+|[0-9]+)[uUlL]*[ \t]*$", text, re.MULTILINE)
        require(len(matches) == 1, "missing/ambiguous literal layout macro: " + name)
        values.append(int(matches[0], 16 if matches[0].lower().startswith("0x") else 10))
    tx, ts, da, ds = values
    require(ts > 0 and ds > 0 and tx + ts <= da and da + ds <= 1 << 64, "invalid cross copy ranges")
    loads = [p for p in image.programs if p.p_type == PT_LOAD]
    require(sum(p.p_vaddr == tx and p.p_filesz == ts and bool(p.p_flags & 1) for p in loads) == 1,
            "cross text constants do not match one executable LOAD")
    require(sum(p.p_vaddr == da and p.p_memsz == ds and bool(p.p_flags & 2) for p in loads) == 1,
            "cross data constants do not match one writable LOAD")
    for name, start, size in ((".text1", tx, ts), (".got", da, ds)):
        s = image.sections.get(name)
        require(s is not None and s.sh_flags & SHF_ALLOC and contains(start, size, s.sh_addr, s.sh_size),
                "required CPE section is outside its copy range: " + name)
    for name, s in image.sections.items():
        if s.sh_flags & SHF_ALLOC and s.sh_size and s.sh_addr < da + ds and s.sh_addr + s.sh_size > tx:
            require(contains(tx, ts, s.sh_addr, s.sh_size) or contains(da, ds, s.sh_addr, s.sh_size),
                    "allocated section intersects an uncopied cross range: " + name)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("before", type=Path)
    parser.add_argument("after", type=Path)
    parser.add_argument("--layout-header", type=Path, help="also check generated copy constants against BOTH images")
    args = parser.parse_args(argv)
    try:
        before, after = Image(args.before.read_bytes()), Image(args.after.read_bytes())
        check_pair(before, after)
        if args.layout_header:
            text = args.layout_header.read_text()
            check_cross_header(before, text)
            check_cross_header(after, text)
        for p in after.programs:
            if p.p_type in (PT_LOAD, PT_TLS):
                print("{} va={:#x} filesz={:#x} memsz={:#x} flags={:#x} align={:#x}".format(
                    "LOAD" if p.p_type == PT_LOAD else "TLS", p.p_vaddr, p.p_filesz, p.p_memsz, p.p_flags, p.p_align))
        for name in (".slave_tls", ".tdata", ".tbss"):
            if name in after.sections:
                s = after.sections[name]
                print("{} va={:#x} size={:#x} flags={:#x}".format(name, s.sh_addr, s.sh_size, s.sh_flags))
        print("PASS: runtime descriptors, allocated bytes and mappings preserved; file offsets may differ")
        print("NOTE: does not validate CPE relocation completeness, private CSR coverage or dynamic LDM/heap use")
    except (OSError, ValueError, struct.error) as error:
        print("FAIL: " + str(error), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
