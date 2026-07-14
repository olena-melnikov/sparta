#!/usr/bin/env python3
"""
ad_convert.py -- mechanical double -> sfloat conversion of the SPARTA
physics core (Phase 2 of AD_PLAN.md).

Rules:
  - token `double`     -> `sfloat`   (type only; comments may be touched,
                                      which is harmless)
  - token `MPI_DOUBLE` -> `MPI_SFLOAT` (macro resolves per build)
  - ubuf bit-punning unions are reverted to `double d;` afterward
    (they reinterpret bits; the pun target must stay 8 bytes)

Passive files (never converted): RNG (random_*), timers, memory
allocator, library C API, STUBS, KOKKOS/FFT packages, math_extra
(already templated by hand), sfloat.h itself.

Idempotent: running twice changes nothing. Reapply after pulling
upstream SPARTA updates. A marker comment is prepended to converted
files so they can be identified (and the set audited) later.
"""

import os, re, sys

SRC = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "..", "src"))

MARKER = "/* AD-CONVERTED: double->sfloat by ad_convert.py (see sfloat.h) */\n"

EXCLUDE_FILES = {
    "sfloat.h", "spatype.h",
    "math_extra.h", "math_extra.cpp",        # templated by hand
    "random_mars.h", "random_mars.cpp",      # RNG: passive by design
    "random_knuth.h", "random_knuth.cpp",
    "timer.h", "timer.cpp",                  # wallclock: passive
    "memory.h", "memory.cpp",                # byte-level allocator
    "library.h", "library.cpp",              # C API stays double-facing
    "main.cpp",
    "spawindows.h",
}
EXCLUDE_DIRS = {"STUBS", "KOKKOS", "FFT", "MAKE", "Obj_serial",
                "Obj_serial_ad", "PYTHON"}

UBUF_RE = re.compile(r"(union\s+ubuf\s*\{[^}]*\})", re.S)

def convert_text(text):
    text = re.sub(r"\bMPI_DOUBLE\b", "MPI_SFLOAT", text)
    text = re.sub(r"\bdouble\b", "sfloat", text)
    # revert ubuf unions: the pun target must remain a real double
    def fix_ubuf(m):
        return m.group(1).replace("sfloat", "double")
    text = UBUF_RE.sub(fix_ubuf, text)
    return text

def main():
    converted = skipped = 0
    for fname in sorted(os.listdir(SRC)):
        path = os.path.join(SRC, fname)
        if os.path.isdir(path):
            continue
        if not (fname.endswith(".cpp") or fname.endswith(".h")):
            continue
        if fname in EXCLUDE_FILES:
            skipped += 1
            continue
        with open(path) as f:
            text = f.read()
        if text.startswith(MARKER):        # idempotency
            continue
        new = convert_text(text)
        if new != text or True:            # mark even if no doubles
            with open(path, "w") as f:
                f.write(MARKER + new)
            converted += 1
    print(f"converted {converted} files, {skipped} passive files skipped")

if __name__ == "__main__":
    main()
