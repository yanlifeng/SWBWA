#!/usr/bin/env bash
# Build the old/new A/B + LWPF-profile binaries for the stage2 kernel comparison,
# and lay out one self-contained run directory per (variant, config, profile).
#
# Prerequisites (set up manually before first use):
#   $AB_ROOT/old   = checkout of the baseline commit (fc658b8)
#   $AB_ROOT/new   = checkout of the current HEAD (a3e7beb)
#
#   variant : old = fc658b8        | new = a3e7beb (HEAD)
#   config  : single (EXEC_MODE=single,    CPE_ALLOCATOR=system)
#             cross  (EXEC_MODE=cgs_cross, CPE_ALLOCATOR=pool)
#   profile : 0 = no LWPF, 1 = CPE_PROFILE=1 (LWPF kernel counters)
#
# Each run directory holds SWBWA + data.bin + data.bin2, because SWBWA reads
# data.bin / data.bin2 from the CURRENT DIRECTORY at runtime (see
# src/host/bwamem.c load_cpe_relocation). Running from a directory that lacks
# them makes the binary die at startup.
#
# NOTE the build-entry difference: fc658b8 only ships build_cross.sh (it does the
# two-pass cross build), while HEAD ships build.sh.
set -euo pipefail

AB_ROOT=${AB_ROOT:?Set AB_ROOT to the experiment directory containing old/ and new/}
JOBS=${JOBS:-8}
mkdir -p "$AB_ROOT/logs" "$AB_ROOT/runs" "$AB_ROOT/scratch"

build_cross()
{
    local variant=$1 prof=$2
    local driver=./build_cross.sh
    [[ "$variant" == new ]] && driver=./build.sh
    ( cd "$AB_ROOT/$variant" && "$driver" "$JOBS" \
        EXEC_MODE=cgs_cross CPE_ALLOCATOR=pool \
        CPE_PROFILE="$prof" HOST_MALLOC_WRAPPER=1 HOST_MALLOC_STATS=0 USE_MPI=0 )
}

build_single()
{
    local variant=$1 prof=$2
    ( cd "$AB_ROOT/$variant" && make clean && make -j "$JOBS" \
        EXEC_MODE=single CPE_ALLOCATOR=system \
        CPE_PROFILE="$prof" HOST_MALLOC_WRAPPER=1 HOST_MALLOC_STATS=0 USE_MPI=0 )
}

for variant in old new; do
    for prof in 0 1; do
        cdir="$AB_ROOT/runs/${variant}_cross_p${prof}"
        sdir="$AB_ROOT/runs/${variant}_single_p${prof}"
        [[ ! -e "$cdir" && ! -e "$sdir" ]] || {
            echo "Refusing to overwrite existing run directories" >&2; exit 1;
        }
        mkdir -p "$cdir" "$sdir"

        # Cross first: this is the build that generates data.bin / data.bin2.
        echo "[BUILD] ${variant}_cross_p${prof}"
        if build_cross "$variant" "$prof" \
               > "$AB_ROOT/logs/build_${variant}_cross_p${prof}.log" 2>&1; then
            cp -f "$AB_ROOT/$variant/SWBWA" "$cdir/SWBWA"
            cp -f "$AB_ROOT/$variant/data.bin"  "$cdir/" 2>/dev/null
            cp -f "$AB_ROOT/$variant/data.bin2" "$cdir/" 2>/dev/null
            echo "[BUILD]   ok: $(ls "$cdir" | tr '\n' ' ')"
        else
            echo "[BUILD]   FAILED (see logs/build_${variant}_cross_p${prof}.log)"
            exit 1
        fi

        echo "[BUILD] ${variant}_single_p${prof}"
        if build_single "$variant" "$prof" \
               > "$AB_ROOT/logs/build_${variant}_single_p${prof}.log" 2>&1; then
            cp -f "$AB_ROOT/$variant/SWBWA" "$sdir/SWBWA"
            echo "[BUILD]   ok: $(ls "$sdir" | tr '\n' ' ')"
        else
            echo "[BUILD]   FAILED (see logs/build_${variant}_single_p${prof}.log)"
            exit 1
        fi
    done
done

echo "[BUILD] ---- run directories ----"
for d in "$AB_ROOT"/runs/*/; do
    printf '%-46s %s\n' "$(basename "$d")" "$(ls "$d" | tr '\n' ' ')"
done
echo "[BUILD] ---- effective CPE_PROFILE / layout per build ----"
for f in "$AB_ROOT"/logs/build_*.log; do
    printf '%-30s profile=%s cg=%s\n' "$(basename "$f")" \
        "$(grep -o 'DSWBWA_ENABLE_CPE_PROFILE=[0-9]*' "$f" | head -1 | cut -d= -f2)" \
        "$(grep -o 'DSWBWA_CPE_PROFILE_CG=[0-9]*' "$f" | head -1 | cut -d= -f2)"
done
