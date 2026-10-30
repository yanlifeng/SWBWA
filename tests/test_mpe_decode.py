#!/usr/bin/env python3
"""ASan/UBSan/LSan check of reused MPE batches with malloc-backed embedded state."""
import os
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
GZIP = ROOT / 'src/gzip'


def main():
    with tempfile.TemporaryDirectory(prefix='swbwa-mpe-decode-') as tmp:
        tmp = Path(tmp)
        for workers in (1, 6):
            flags = ['-O1', '-g', '-pthread', '-fsanitize=address,undefined',
                     '-fno-sanitize-recover=all', '-fno-pie',
                     '-DSWBWA_CODEC_MPE=1', '-DSWBWA_CODEC_LDM=0',
                     '-I' + str(ROOT / 'include'), '-I' + str(GZIP),
                     '-DSWBWA_HOST_MPE_THREADS=' + str(workers)]
            aliases = ['-Dswbwa_codec_execute=swbwa_codec_execute_mpe',
                       '-Dswbwa_codec_gzip_worker=swbwa_codec_gzip_worker_mpe',
                       '-Dswbwa_rapid_discover=swbwa_rapid_discover_mpe']
            objects = []
            for source in (ROOT / 'src/host/swbwa_host_workers.c',
                           GZIP / 'mpe_decode.c', GZIP / 'codec_job.c'):
                obj = tmp / (source.stem + '.o')
                special = ['-DSWBWA_CODEC_RAPID=1', *aliases] if source.name == 'codec_job.c' else []
                subprocess.run([os.getenv('CC', 'cc'), '-std=gnu11', *flags, *special,
                                '-c', str(source), '-o', str(obj)], check=True)
                objects.append(str(obj))
            rapid = tmp / 'rapid.o'
            subprocess.run([os.getenv('CXX', 'c++'), '-std=c++17', *flags, *aliases,
                            '-I' + str(GZIP / 'vendor/rapidgzip'), '-c',
                            str(GZIP / 'rapid_firstpass.cpp'), '-o', str(rapid)], check=True)
            binary = tmp / 'test'
            subprocess.run([os.getenv('CXX', 'c++'), '-std=c++17', *flags, '-no-pie',
                            str(ROOT / 'tests/test_mpe_decode.cpp'), *objects, str(rapid),
                            '-lz', '-o', str(binary)], check=True)
            env = dict(os.environ, ASAN_OPTIONS='detect_leaks=1:halt_on_error=1')
            subprocess.run([str(binary)], check=True, timeout=90, env=env)
            print(f'PASS MPE decode workers={workers} SWBWA_CODEC_MPE=1 SWBWA_CODEC_LDM=0 ASan/UBSan/LSan',
                  flush=True)


if __name__ == '__main__':
    main()
