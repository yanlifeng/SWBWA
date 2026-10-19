"""Read-only ELF comparison: offset movement, TLS templates and corruption."""
from pathlib import Path
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from check_elf_layout import HEADER, PHDR, SHDR, Image, check_pair, check_cross_header


def fixture(shift=0):
    text, data, names_offset = 0x1000 + shift, 0x2000 + shift, 0x3000 + shift
    shoff = names_offset + 0x100
    names = b"\0.text1\0.got\0.tdata\0.tbss\0.bss\0.shstrtab\0"
    sections = [
        (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        (names.index(b".text1"), 1, 6, 0x100010, text + 0x10, 8, 0, 0, 8, 0),
        (names.index(b".got"), 1, 3, 0x200000, data, 8, 0, 0, 8, 0),
        (names.index(b".tdata"), 1, 0x403, 0x200010, data + 0x10, 8, 0, 0, 8, 0),
        # .tbss deliberately overlaps file-backed LOAD bytes, as real ELF can.
        (names.index(b".tbss"), 8, 0x403, 0x200018, data + 0x18, 0x10, 0, 0, 8, 0),
        (names.index(b".bss"), 8, 3, 0x200040, data + 0x40, 0x20, 0, 0, 8, 0),
        (names.index(b".shstrtab"), 3, 0, 0, names_offset, len(names), 0, 0, 1, 0),
    ]
    programs = [
        (1, 5, text, 0x100000, 0x100000, 0x40, 0x40, 0x1000),
        (1, 6, data, 0x200000, 0x200000, 0x40, 0x100, 0x1000),
        (7, 4, data + 0x10, 0x200010, 0x200010, 8, 0x18, 8),
    ]
    result = bytearray(shoff + len(sections) * SHDR.size)
    HEADER.pack_into(result, 0, b"\x7fELF\x02\x01\x01" + bytes(9), 2, 62, 1,
                     0x100010, HEADER.size, shoff, 0, HEADER.size, PHDR.size,
                     len(programs), SHDR.size, len(sections), len(sections) - 1)
    for i, program in enumerate(programs):
        PHDR.pack_into(result, HEADER.size + i * PHDR.size, *program)
    for i, section in enumerate(sections):
        SHDR.pack_into(result, shoff + i * SHDR.size, *section)
    result[text:text + 0x40] = bytes(range(0x40))
    result[data:data + 0x40] = bytes(range(0x40, 0x80))
    result[names_offset:names_offset + len(names)] = names
    return result


def patch(data, fmt, offset, field, value):
    result = bytearray(data)
    values = list(fmt.unpack_from(result, offset))
    values[field] = value
    fmt.pack_into(result, offset, *values)
    return result


class LayoutTest(unittest.TestCase):
    def test_offsets_may_move_without_symbols(self):
        check_pair(Image(fixture()), Image(fixture(0x1000)))

    def test_loaded_section_wrong_file_mapping(self):
        original = fixture()
        shoff = Image(original).header.e_shoff
        changed = patch(original, SHDR, shoff + SHDR.size, 4, 0x1020)
        changed[0x1020:0x1028] = original[0x1010:0x1018]
        with self.assertRaisesRegex(ValueError, "LOAD mapping"):
            check_pair(Image(original), Image(changed))

    def test_runtime_descriptor_changes(self):
        original = fixture()
        for index, field, value in ((0, 4, 0x110000), (0, 7, 0x100),
                                    (1, 5, 0x30), (1, 6, 0x200), (1, 1, 7),
                                    (2, 6, 0x20)):
            with self.subTest(index=index, field=field):
                changed = patch(original, PHDR, HEADER.size + index * PHDR.size, field, value)
                with self.assertRaisesRegex(ValueError, "runtime PHDR changed"):
                    check_pair(Image(original), Image(changed))

    def test_invalid_load_alignment(self):
        original = fixture()
        changed = patch(original, PHDR, HEADER.size, 2, 0x1001)
        with self.assertRaisesRegex(ValueError, "alignment mismatch"):
            check_pair(Image(original), Image(changed))

    def test_entry_and_identity(self):
        original = fixture()
        for field, value in ((4, 0x100011), (2, 0x9026), (7, 1)):
            with self.subTest(field=field):
                with self.assertRaisesRegex(ValueError, "identity, ABI, entry, or flags"):
                    check_pair(Image(original), Image(patch(original, HEADER, 0, field, value)))

    def test_allocated_bytes(self):
        original = fixture()
        changed = bytearray(original)
        changed[0x1010] ^= 1
        with self.assertRaisesRegex(ValueError, "section bytes changed"):
            check_pair(Image(original), Image(changed))

    def test_allocated_set_and_flags(self):
        original = fixture()
        shoff = Image(original).header.e_shoff
        for value, error in ((0, "section set changed"), (2, "section metadata changed")):
            changed = patch(original, SHDR, shoff + SHDR.size, 2, value)
            with self.assertRaisesRegex(ValueError, error):
                check_pair(Image(original), Image(changed))

    def test_bss_is_zero_fill_but_tls_bss_may_overlap_load_bytes(self):
        original = fixture()
        check_pair(Image(original), Image(original))
        shoff = Image(original).header.e_shoff
        changed = patch(original, SHDR, shoff + 5 * SHDR.size, 3, 0x200020)
        with self.assertRaisesRegex(ValueError, "LOAD mapping"):
            check_pair(Image(original), Image(changed))

    def test_truncated_or_unknown_format(self):
        original = fixture()
        for data in (original[:30], original[:-1], b"NOTELF" + original[6:]):
            with self.assertRaises(ValueError):
                Image(data)

    def test_cross_header(self):
        header = "\n".join("#define SWBWA_CPE_" + name + " " + value + "UL" for name, value in (
            ("TEXT_START_ADDRESS", "0x100000"), ("TEXT_SEGMENT_BYTES", "0x40"),
            ("DATA_START_ADDRESS", "0x200000"), ("DATA_SEGMENT_BYTES", "0x100")))
        check_cross_header(Image(fixture()), header)
        with self.assertRaisesRegex(ValueError, "text constants"):
            check_cross_header(Image(fixture()), header.replace("0x40UL", "0x41UL"))
        with self.assertRaisesRegex(ValueError, "data constants"):
            check_cross_header(Image(fixture()), header.replace("0x100UL", "0x101UL"))

    @unittest.skipUnless(shutil.which("cc") and shutil.which("objcopy"), "native compiler/binutils unavailable")
    def test_native_compaction_and_cli(self):
        with tempfile.TemporaryDirectory(prefix="swbwa-elf-check-") as tmp:
            original, compact = Path(tmp) / "before", Path(tmp) / "after"
            subprocess.run(shlex.split(os.environ.get("CC", "cc")) + [
                "-g", "-O0", "-fno-pie", "-no-pie", "-x", "c", "-", "-o", str(original)],
                input="__thread int x = 3; __thread char y[64]; char z[128];\n"
                      "int main(void) { return x + y[0] + z[0]; }\n", text=True, check=True)
            subprocess.run([sys.executable, "-B", str(ROOT / "tools/compact_elf.py"),
                            str(original), str(compact)], check=True,
                           env=dict(os.environ, OBJCOPY=shutil.which("objcopy")))
            check_pair(Image(original.read_bytes()), Image(compact.read_bytes()))
            result = subprocess.run([sys.executable, "-B", str(ROOT / "tools/check_elf_layout.py"),
                                     str(original), str(compact)], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("PASS:", result.stdout)


if __name__ == "__main__":
    unittest.main()
