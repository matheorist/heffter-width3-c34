"""Generate the periodic c=3 construction from the ancillary rule files.

This implements the periodic construction for admissible m >= 35.
It is deterministic: no SAT/CP search is performed.
"""

import argparse
import ast
import json
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def zhou_value(n, i, j):
    r = 2 * n + 1
    parity = -1 if (j + n) % 2 else 1
    return int(
        Fraction(1 - parity * r, 4) * i
        + Fraction(2 - 3 * abs(i), 2) * j
    )


def load_json(name):
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def rule_data(n):
    residue = n % 24
    if residue == 9:
        words = load_json("joint_words_P24.json")
        boundary_name = "boundary_105.json" if n % 48 == 9 else "boundary_129.json"
        phases = (6, 12)
    elif residue == 21:
        words = load_json("words21.json")
        boundary_name = (
            "boundary21_117.json" if n % 48 == 21 else "boundary21_141.json"
        )
        phases = (10, 13)
    else:
        raise ValueError("the c=3 rule requires n = 9 or 21 (mod 24)")
    return words, load_json(boundary_name), phases


def support_rows(n):
    words, boundary, phases = rule_data(n)
    kappa = (n + 1) // 2
    rows = []
    for j in range(1, n + 1):
        if j <= 20:
            delta = boundary["TL"][str(j)]
        elif j <= kappa - 21:
            delta = words["W1"][(j - phases[0]) % 24]
        elif j <= kappa + 20:
            delta = boundary["TM"][str(j - kappa)]
        elif j <= n - 21:
            delta = words["W2"][(j - phases[1]) % 24]
        else:
            delta = boundary["TU"][str(j - n)]
        row = [
            3 * zhou_value(n, -1, j) + delta[0],
            3 * j + delta[1],
            3 * zhou_value(n, 1, j) + delta[2],
        ]
        if sum(row) != 0:
            raise AssertionError((n, j, row))
        rows.append(row)
    return rows


def category_key(j, m, n):
    kappa = (n + 1) // 2
    if j <= 20:
        return ("bnd", "low", j)
    if abs(j - m) <= 12:
        return ("bnd", "m", j - m)
    if abs(j - kappa) <= 20:
        return ("bnd", "k", j - kappa)
    if j >= n - 20:
        return ("bnd", "top", j - n)
    if j < m - 12:
        zone = 1
    elif j < kappa - 20:
        zone = 2
    else:
        zone = 3
    return ("bulk", zone, j % 24)


def pattern_for(n, group, patterns):
    cls = n % 24
    t = (n - cls) // 24
    state = 15 + ((t - 15) % 6)
    return patterns[f"cls{cls}_t{state:03d}_r{group}"]


def orient(mu_sign, diff_sign, triple, mode):
    a, b, c = sorted(abs(value) for value in triple)
    if a + b != c:
        raise AssertionError((a, b, c))
    if mode == "c":
        mu, diff, third = c, b - a, c
    elif mode == "a":
        mu, diff, third = a, b + c, a
    elif mode == "b":
        mu, diff, third = b, a + c, b
    else:
        raise ValueError(mode)
    first_numerator = mu_sign * mu + diff_sign * diff
    second_numerator = mu_sign * mu - diff_sign * diff
    if first_numerator % 2 or second_numerator % 2:
        raise AssertionError((triple, mode, mu_sign, diff_sign))
    row = [
        first_numerator // 2,
        second_numerator // 2,
        -mu_sign * third,
    ]
    if sum(row) != 0 or sorted(abs(value) for value in row) != [a, b, c]:
        raise AssertionError((triple, mode, row))
    return row


def generate(m):
    if m < 35 or m % 4 != 3 or m % 3 == 0:
        raise ValueError("periodic c=3 generator requires m>=35, m=3 (mod 4), gcd(m,3)=1")
    n = 3 * m
    triples = support_rows(n)
    patterns = load_json("lemma2_patterns.json")
    arrays = []
    for group in range(3):
        pattern = pattern_for(n, group, patterns)
        c_signs = {ast.literal_eval(key): value for key, value in pattern["c"].items()}
        d_signs = {ast.literal_eval(key): value for key, value in pattern["d"].items()}
        switch_key, switch_mode = pattern["switch"]
        switch_key = None if switch_key == "None" else ast.literal_eval(switch_key)
        array = []
        for j in range(1, n + 1):
            if j % 3 != group:
                continue
            key = category_key(j, m, n)
            mode = switch_mode if key == switch_key else "c"
            array.append(orient(c_signs[key], d_signs[key], triples[j - 1], mode))
        arrays.append(array)
    return arrays


def verify(m, arrays):
    assert len(arrays) == 3
    assert all(len(array) == m for array in arrays)
    assert all(sum(row) == 0 for array in arrays for row in array)
    assert all(
        [sum(row[column] for row in array) for column in range(3)] == [0, 0, 0]
        for array in arrays
    )
    support = [abs(value) for array in arrays for row in array for value in row]
    assert sorted(support) == list(range(1, 9 * m + 1))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("m", type=int)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    arrays = generate(args.m)
    verify(args.m, arrays)
    print(f"VERIFIED IHS({args.m},3;3): support [1,{9 * args.m}]")
    if args.output:
        data = {"m": args.m, "c": 3, "arrays": arrays}
        args.output.write_text(json.dumps(data, indent=2), encoding="utf-8")
        print(f"WROTE {args.output}")


if __name__ == "__main__":
    main()
