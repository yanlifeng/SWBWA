"""Exact bytes, arbitrary chunk boundaries, writer integration and >4 GiB count."""
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SAM = b"r1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n"


def build_md5_object(directory):
    obj = str(Path(directory) / "sam_md5.o")
    subprocess.run(shlex.split(os.getenv("CXX", "c++")) + [
        "-std=c++11", "-O2", "-Iinclude", "-c", "src/host/swbwa_sam_md5.cpp",
        "-o", obj], cwd=ROOT, check=True)
    return [obj, "-lstdc++"]


def report(text):
    return json.loads(re.search(rb"SAM_STREAM_MD5 (\{[^\n]+\})", text)[1])


def main():
    with tempfile.TemporaryDirectory(prefix="sam-stream-md5-") as tmp:
        tmp = Path(tmp)
        harness = tmp / "harness.cpp"
        harness.write_text(r'''
#include "swbwa_sam_md5.h"
#include <cstdio>
#include <cstdlib>
#include <vector>
int main(int argc, char **argv) {
    if (argc != 3) return 2;
    std::vector<unsigned char> buffer(std::strtoul(argv[1], 0, 10));
    unsigned repeat = std::strtoul(argv[2], 0, 10);
    if (swbwa_sam_md5_open() || swbwa_sam_md5_update(NULL, 0)) return 1;
    size_t n;
    while ((n = std::fread(buffer.data(), 1, buffer.size(), stdin)))
        for (unsigned i = 0; i < repeat; ++i)
            if (swbwa_sam_md5_update(buffer.data(), n)) return 1;
    return std::ferror(stdin) || swbwa_sam_md5_close(stdout) ? 1 : 0;
}
''')
        common = ["-g", "-O2", "-fsanitize=address,undefined", "-fno-sanitize-recover=all",
                  "-fno-pie", "-no-pie", "-I" + str(ROOT / "include")]
        helper = str(ROOT / "src/host/swbwa_sam_md5.cpp")
        binary = tmp / "stream"
        subprocess.run(["c++", "-std=c++11"] + common + [str(harness), helper, "-o", str(binary)], check=True)
        headers = b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n@CO\t" + b"x"*70000 + b"\n"
        for prefix in (b"", headers):
            for body in (SAM, SAM * 10, SAM + b"@SQ\tnot-a-leading-header\n", b"x\n"):
                for block in (1, 3, 4, 7, 63, 64, 65, 4096, 1048576):
                    result = subprocess.run([str(binary), str(block), "1"], input=prefix+body,
                                            capture_output=True, check=True)
                    value = report(result.stdout)
                    assert value["md5"] == hashlib.md5(body).hexdigest(), (prefix[:4], block)
                    assert value["bytes"] == len(body) and value["records"] == body.count(b"\n")
        for bad in (b"", headers, b"@SQ\tunterminated", SAM[:-1], SAM+b"\0", b"@CO\tbad\0\n"+SAM):
            assert subprocess.run([str(binary), "3", "1"], input=bad, capture_output=True).returncode != 0
        print("PASS arbitrary header/record splits, exact bytes, malformed/truncated input", flush=True)

        ctest = tmp / "writer.c"
        ctest.write_text(r'''
#include "swbwa_output.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
int swbwa_mpi_rank(void) { return 0; }
int swbwa_mpi_size(void) { return 1; }
int swbwa_mpi_is_root(void) { return 1; }
void swbwa_mpi_print_rank_ordered(void (*f)(void)) { f(); }
int main(int argc, char **argv) {
    const char *s="r1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n";
    assert(argc==2 && access(argv[1], F_OK)!=0);
    setenv("SWBWA_OUTPUT_MD5", "1", 1);
    assert(swbwa_output_open(argv[1], 0)==0);
    for (int i=0; i<100; ++i) {
        assert(swbwa_output_write(s, 5)==0);
        assert(swbwa_output_flush()==0);
        assert(swbwa_output_write(s+5, strlen(s)-5)==0);
    }
    assert(swbwa_output_close()==0 && access(argv[1], F_OK)!=0);
    setenv("SWBWA_OUTPUT_MD5", "0", 1);
    assert(swbwa_output_open(argv[1], 0)==0);
    assert(swbwa_output_write(s, strlen(s))==0);
    assert(swbwa_output_close()==0 && access(argv[1], F_OK)==0);
    return 0;
}
''')
        objects = []
        for source in (ctest, ROOT / "src/host/swbwa_output.c"):
            obj = tmp / (source.stem + ".o")
            subprocess.run(["cc", "-std=gnu11"] + common + ["-DSWBWA_USE_MPI=0",
                "-DSWBWA_OUTPUT_MODE=SWBWA_OUTPUT_SPLIT", "-DSWBWA_ENABLE_HOST_MALLOC_WRAPPER=0",
                "-c", str(source), "-o", str(obj)], check=True)
            objects.append(str(obj))
        writer = tmp / "writer"
        subprocess.run(["c++", "-std=c++11"] + common + objects + [helper, "-o", str(writer)], check=True)
        output = tmp / "actual.sam"
        result = subprocess.run([str(writer), str(output)], capture_output=True, check=True)
        assert report(result.stderr)["md5"] == hashlib.md5(SAM*100).hexdigest()
        assert output.read_bytes() == SAM
        print("PASS real output writer: no file in MD5 mode; ordinary file mode unchanged", flush=True)

        # Repeated in memory, without creating a multi-GB fixture on disk.
        block = (b"A"*63+b"\n") * (1048576//64)
        expected = hashlib.md5()
        for _ in range(4097):
            expected.update(block)
        result = subprocess.run([str(binary), "1048576", "4097"], input=block,
                                capture_output=True, check=True)
        value = report(result.stdout)
        assert value["md5"] == expected.hexdigest()
        assert value["bytes"] == 1048576*4097 and value["records"] == 16384*4097
        print("PASS >4 GiB stream:", value, flush=True)


if __name__ == "__main__":
    main()
