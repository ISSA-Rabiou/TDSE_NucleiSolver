#!/usr/bin/env python3
"""Apply the TD_SCHROD re-entrancy fix to Tnum-Tana's Bunch/Poly reader."""
from pathlib import Path
import re
import sys

if len(sys.argv) != 2:
    raise SystemExit(f"Usage: {sys.argv[0]} path/to/BunchPolyTransfo.f90")
path = Path(sys.argv[1])
if not path.is_file():
    raise SystemExit(f"Tnum-Tana source file not found: {path}")
text = path.read_text()
marker = "The counters below have SAVE semantics because this routine is"
if marker in text:
    print("Tnum-Tana Bunch/Poly fix already present")
    raise SystemExit(0)

sub = re.search(r"(?im)^\s*(?:recursive\s+)?subroutine\s+RecRead_BFTransfo\b", text)
if not sub:
    raise SystemExit(f"RecRead_BFTransfo not found in {path}")
anchor = re.search(r"(?im)^\s*IF\s*\(\s*iv_tot\s*==\s*0\s*\)\s*THEN\s*$", text[sub.start():])
if not anchor:
    raise SystemExit(f"Initialization anchor not found in {path}")
pos = sub.start() + anchor.start()
fix = (
    "    ! The counters below have SAVE semantics because this routine is\n"
    "    ! recursive. Reset them only on entry to a new top-level transformation;\n"
    "    ! recursive calls see rec_level > 0 and must retain the shared state.\n"
    "    IF (rec_level == 0) THEN\n"
    "      iv_tot            = 0\n"
    "      i_Qpoly           = 0\n"
    "      num_Frame_in_BF   = 0\n"
    "      num_vect_in_BF    = 0\n"
    "      zmat_order_save   = .FALSE.\n"
    "    END IF\n\n"
)
text = text[:pos] + fix + text[pos:]
path.write_text(text)
print(f"Applied Tnum-Tana Bunch/Poly re-entrancy fix: {path}")
