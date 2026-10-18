"""Read-only audit of the finite constructions printed in the c=3,4 paper.

Standard library only; no JSON, generators, solvers, or output files are used.
Run with python -B.  Optional --compare-root checks exact decoded row equality
against an earlier directory containing the same TeX sources.
"""

import argparse
import math
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
if (ROOT / "paper").is_dir():
    ROOT = ROOT / "paper"
HEX_ROW = r"&\\texttt\{([0-9A-F]+)\}&\\texttt\{([0-9A-F]+)\}"


def alphabet(source, command, symbols):
    headers = re.findall(r"\\text\{symbol\}&([^\\\r\n]+)", source)
    assert "".join(item.strip() for row in headers for item in row.split("&")) == symbols
    rows = re.findall(re.escape(command) + r"&([^\r\n]+)", source)
    values = [tuple(map(int, match)) for row in rows
              for match in re.findall(r"\((-?\d+),(-?\d+),(-?\d+)\)", row)]
    assert len(values) == len(symbols), (command, len(values))
    assert len(set(values)) == len(values)
    assert all(sum(d) == 0 and max(map(abs, d)) <= 2 for d in values)
    return dict(zip(symbols, values))


def small_words(source, command, expected, c):
    pairs = re.findall(r"^(\d+)&" + re.escape(command)
                       + r"\{([0-9A-I]+)\}", source, re.M)
    result = {int(m): word for m, word in pairs}
    assert len(result) == len(pairs) and set(result) == set(expected)
    assert all(len(word) == c * m for m, word in result.items())
    return result


def triple(n, j, delta):
    """Apply the manuscript's Kotzig formula and corrected lift exactly."""
    parity = (-1) ** (j + n)
    row = []
    for i, d in zip((-1, 0, 1), delta):
        numerator = (1 - parity * (2 * n + 1)) * i
        numerator += 2 * (2 - 3 * abs(i)) * j
        assert numerator % 4 == 0
        row.append(3 * (numerator // 4) + d)
    assert sum(row) == 0
    a, b, c = sorted(map(abs, row))
    assert 0 < a < b < c and a + b == c, (n, j, row)
    return a, b, c


def signs(word, length):
    bits = "".join(f"{int(digit, 16):04b}" for digit in word)
    assert len(bits) == 4 * ((length + 3) // 4)
    return [1 if bit == "1" else -1 for bit in bits[:length]]


def orient(values, mode, epsilon, eta):
    a, b, c = values
    mu, nu = {"c": (c, b - a), "a": (a, b + c), "b": (b, a + c)}[mode]
    assert (epsilon * mu + eta * nu) % 2 == 0
    return ((epsilon * mu + eta * nu) // 2,
            (epsilon * mu - eta * nu) // 2, -epsilon * mu)


def check_arrays(m, c, arrays):
    assert len(arrays) == c and all(len(array) == m for array in arrays)
    assert all(len(row) == 3 and sum(row) == 0 for array in arrays for row in array)
    assert all(tuple(map(sum, zip(*array))) == (0, 0, 0) for array in arrays)
    support = [abs(x) for array in arrays for row in array for x in row]
    assert sorted(support) == list(range(1, 3 * c * m + 1))


def category3(m, j):
    n, k = 3 * m, (3 * m + 1) // 2
    if j <= 20:
        return ("bnd", "low", j)
    if abs(j - m) <= 12:
        return ("bnd", "m", j - m)
    if abs(j - k) <= 20:
        return ("bnd", "k", j - k)
    if j >= n - 20:
        return ("bnd", "top", j - n)
    if j < m - 12:
        return ("bulk", 1, j % 24)
    if m + 12 < j < k - 20:
        return ("bulk", 2, j % 24)
    return ("bulk", 3, j % 24)


def literal_key(key):
    # Literal-string order, NOT numeric tuple order; no inserted whitespace.
    return "(" + ",".join(repr(item) for item in key) + ")"


def audit(root):
    def read(name):
        return (root / name).read_text(encoding="utf-8")

    c3 = read("c3_finite_certificate_body.tex")
    s3 = read("c3_exceptional_constructions_body.tex")
    c4 = read("c4_main_certificate_body.tex")
    main = read("IHS_m_3_3_standalone.tex")
    result = {}
    a3 = alphabet(c3, r"\text{offset}", "0123456789ABC")
    a_small3 = alphabet(s3, r"\delta", "0123456789ABCDEFGHI")
    small4 = c4.split(r"\subsection{The six bounded constructions}", 1)[1]
    a_small4 = alphabet(small4, r"\delta", "0123456789ABCDEFGHI")
    assert a_small3 == a_small4
    master_pairs = re.findall(r"W_([12])&\\mathtt\{([0-9A-C]+)\}", c3)
    master = {int(i): word for i, word in master_pairs}
    assert len(master_pairs) == 2 and set(master) == {1, 2}
    assert all(len(word) == 24 for word in master.values())
    boundary_pairs = re.findall(
        r"^(9|33|21|45)&([LMUKT])&[^&]+&\\mathtt\{([0-9A-C]+)\}", c3, re.M)
    # M/U support labels are accepted only for comparing pre-revision sources.
    boundary = {(int(e), {"M": "K", "U": "T"}.get(name, name)): word
                for e, name, word in boundary_pairs}
    assert len(boundary_pairs) == len(boundary) == 12
    assert set(boundary) == {(e, name) for e in (9, 21, 33, 45) for name in "LKT"}
    assert all(len(word) == {"L": 20, "K": 41, "T": 21}[name]
               for (_, name), word in boundary.items())

    def offset3(m, j):
        n, k = 3 * m, (3 * m + 1) // 2
        if j <= 20:
            return a3[boundary[n % 48, "L"][j - 1]]
        if abs(j - k) <= 20:
            return a3[boundary[n % 48, "K"][j - k + 20]]
        if j >= n - 20:
            return a3[boundary[n % 48, "T"][j - n + 20]]
        if n % 24 == 9:
            i, phase, shift = (1, 6, 0) if j < k else (2, 12, 0)
        else:
            assert n % 24 == 21
            i, phase, shift = (1, 10, 16) if j < k else (2, 13, 13)
        return a3[master[i][(j - phase + shift) % 24]]

    words3 = small_words(s3, r"\mathtt", (11, 19, 23), 3)
    words4 = small_words(small4, r"\texttt", (5, 7, 11, 13, 17, 19), 4)
    parameters3 = (11, 19, 23, 31, 35, 43, 47, 55, 59)
    for c, source, expected, words, alpha in (
            (3, s3, parameters3, words3, a_small3),
            (4, small4, tuple(words4), words4, a_small4)):
        records = re.findall(r"^(\d+)&(\d+)&(.*?)" + HEX_ROW, source, re.M)
        assert len(records) == c * len(expected)
        seen, arrays = set(), {m: [None] * c for m in expected}
        for mm, rr, switch_text, E, H in records:
            m, r = int(mm), int(rr)
            assert m in arrays and 0 <= r < c and (m, r) not in seen
            seen.add((m, r))
            entries = re.findall(r"(\d+):([ab])", switch_text)
            switches = {int(t): mode for t, mode in entries}
            assert len(switches) == len(entries)
            assert all(0 <= t < m for t in switches)
            assert entries or switch_text.strip() == "--"
            epsilon, eta = signs(E, m), signs(H, m)
            rows = []
            for t in range(m):
                j = c * t + r + 1
                d = alpha[words[m][j - 1]] if m in words else offset3(m, j)
                rows.append(orient(triple(c * m, j, d), switches.get(t, "c"),
                                   epsilon[t], eta[t]))
            arrays[m][r] = rows
        for m, rows in arrays.items():
            check_arrays(m, c, rows)
            result[(f"c{c}-bounded", m)] = rows

    patterns = {}
    records = re.findall(r"^(\d+)&(\d+)&(\d+)&(\d+)&(.*?)" + HEX_ROW, c3, re.M)
    assert len(records) == 24
    for ee, ss, rr, ll, switch_text, C, D in records:
        e, s, r, length = map(int, (ee, ss, rr, ll))
        assert (e, s, r) not in patterns
        category_sets = []
        for tau in (s, s + 6):
            keys = set()
            for Q in range(7):
                m = (e + 24 * (tau + 12 * Q)) // 3
                keys.update(category3(m, j) for j in range(1, 3 * m + 1) if j % 3 == r)
            category_sets.append(keys)
        assert category_sets[0] == category_sets[1]
        keys = sorted(category_sets[0], key=literal_key)
        assert len(keys) == length
        switches = re.findall(r"K:(-?\d+);([ab])", switch_text)
        assert len(switches) <= 1
        assert switches or switch_text.strip() == "--"
        switch = ("bnd", "k", int(switches[0][0])) if switches else None
        mode = switches[0][1] if switches else "c"
        patterns[e, s, r] = (dict(zip(keys, signs(C, length))),
                              dict(zip(keys, signs(D, length))), switch, mode)
    assert set(patterns) == {(e, s, r) for e, states in
                            ((9, (16, 17, 19, 20)), (21, (15, 17, 18, 20)))
                            for s in states for r in range(3)}
    bridge = (67, 71, 79, 83, 91, 95, 103, 107, 115, 119)
    assert list(bridge) == [m for m in range(67, 120) if m % 4 == 3 and math.gcd(m, 3) == 1]
    for m in bridge:
        n, e = 3 * m, (3 * m) % 24
        s = 15 + (((n - e) // 24 - 15) % 6)
        arrays = []
        for r in range(3):
            epsilon, eta, switch, mode = patterns[e, s, r]
            rows = []
            for j in range(1, n + 1):
                if j % 3 == r:
                    key = category3(m, j)
                    rows.append(orient(triple(n, j, offset3(m, j)),
                                       mode if key == switch else "c", epsilon[key], eta[key]))
            arrays.append(rows)
        check_arrays(m, 3, arrays)
        result[("c3-bridge", m)] = arrays

    for c, marker in ((3, "The smallest admissible instance"),
                      (4, "The compact data in Section")):
        example = main.split(marker, 1)[1].split(r"\end{center}", 1)[0]
        arrays = []
        for body in re.findall(r"\\begin\{array\}\{rrr\}(.*?)\\end\{array\}", example, re.S):
            arrays.append([tuple(map(int, row.strip().split("&")))
                           for row in body.strip().split(r"\\")])
        check_arrays(7, c, arrays)
        result[(f"c{c}-example", 7)] = arrays
    assert result["c4-example", 7] == result["c4-bounded", 7]
    assert [7] + list(parameters3) == [m for m in range(3, 67)
                                      if m % 4 == 3 and math.gcd(m, 3) == 1]
    assert min(m for m in range(120, 128) if m % 4 == 3 and math.gcd(m, 3) == 1) == 127
    starts = {1: 25, 5: 29, 7: 31, 11: 23}
    assert sorted(words4) == [m for m in range(5, 32, 2)
                              if m % 3 and m < starts[m % 12]]
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--compare-root", type=Path)
    args = parser.parse_args()
    decoded = audit(args.root)
    print("PASS: c=3 nine bounded parameters / 27 records / 969 rows")
    print("PASS: c=4 six bounded parameters / 24 records / 288 rows")
    print("PASS: c=3 ten finite-bridge parameters; fixed category lists and literal-string ordering")
    print("PASS: both printed m=7 examples; c=4 example equals its compact decoding")
    print("PASS: word lengths, exact supports, integer entries, row/column sums, small-parameter coverage")
    if args.compare_root:
        assert decoded == audit(args.compare_root), "Decoded signed entries changed"
        print("PASS: exact decoded signed entries unchanged against comparison TeX sources")


if __name__ == "__main__":
    main()
