# Reproducibility package for the width-three Heffter array sets

This repository accompanies the manuscript
`Existence of integer Heffter array sets IHS(m,3;3) and IHS(m,3;4)`.
It contains only the reproducibility material for the complete fixed slices
`c=3` and `c=4`.

## Contents

- `paper/`: the standalone LaTeX source and the three finite-certificate source files;
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

The scripts are deterministic certificate checkers or generators. Their output is not used as a
logical oracle in the manuscript: the finite words, tables, sign patterns, and seed certificates
are included explicitly in the paper source package.

## Paper compilation

Open `paper/IHS_m_3_3_standalone.tex` in a LaTeX environment and compile it with BibTeX in the
usual sequence. The three `\\input{...}` files must remain in the same `paper/` directory.
