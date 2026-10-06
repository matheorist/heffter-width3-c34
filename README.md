# Reproducibility package for the width-three Heffter array sets

This repository accompanies the manuscript
`Existence of integer Heffter array sets IHS(m,3;3) and IHS(m,3;4)`.
It contains only the reproducibility material for the complete fixed slices
`c=3` and `c=4`.

## Contents

- `paper/`: the standalone LaTeX source, the c=3 certificate files, the c=4 support/seed
  certificate, and the appendix containing the 192 c=4 balance patterns;
- `code/c3/`: the period-24 support and balance rules, JSON certificates, and c=3 verifiers;
- `code/c4/`: the period-48 support and balance rules, finite seeds, and c=4 verifiers.

The repository intentionally excludes the unfinished `c=5`, `c=7`, width-five experiments,
temporary logs, generated PDFs, and local runtime directories.

## Requirements

Python 3.10 or later is recommended. The c=3 support and transport checks use only the Python
standard library. The c=4 generation and affine-balance regeneration scripts additionally use
the packages listed in `requirements.txt` (`ortools` and `python-sat`).

Install them with:

```text
python -m pip install -r requirements.txt
```

## Verification

From the corresponding code directory, run:

```text
python verify_c3_plus4_transport.py --minimum 67 --maximum 200
python prove_c3_support_affine.py
python prove_c3_balance_affine.py
python prove_c4_support_affine.py
python verify_c4_all_range.py --maximum 100
python verify_c4_balance_patterns.py
```

The c=4 balance prover can regenerate the 192 affine sign patterns with:

```text
python prove_c4_balance_affine.py
```

The verification scripts are deterministic checkers for the canonical certificate data.  The c=3
transport checker covers the canonical periodic regime (the displayed command starts at (m=67));
the six smaller c=3 preperiod orientations are printed in the manuscript and are not regenerated
by `generate_c3_all.py`.  The c=4 balance prover is an optional solver-backed regeneration tool:
it may return a different valid sign certificate, depending on solver version and search settings,
whereas `verify_c4_balance_patterns.py` checks the canonical JSON shipped in this repository.
None of these programs is used as a logical oracle in the manuscript: the finite words, tables,
sign patterns, and seed certificates are included explicitly in the paper source package.  The
long 192-row c=4 balance table is printed in the appendix, while the other certificate data are
kept in the main text.

The environment used for the checked release was Python 3.12, OR-Tools 9.15.6755, and
python-sat 1.9.dev5.  Other compatible versions may also work, but can produce a different valid
solver certificate.

## Paper compilation

Open `paper/IHS_m_3_3_standalone.tex` in a LaTeX environment and compile it with BibTeX in the
usual sequence. The four `\\input{...}` files must remain in the same `paper/` directory.
