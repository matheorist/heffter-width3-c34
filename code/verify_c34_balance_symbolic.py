"""Exact symbolic summation of the printed c=3 and c=4 balance patterns.

This is not a finite-difference test.  Each bulk category is an arithmetic
progression with affine endpoints.  Its sum is expanded in Q with Fraction
arithmetic, and all three coefficients of each signed total must be zero.
For c=4 the w cusp is split at h=2Q+floor(-B/432); Q=0 is checked separately.
Only the Python standard library is required.  No solver is imported.
"""

from __future__ import annotations

import ast
import json
import re
from fractions import Fraction as F
from pathlib import Path


HERE = Path(__file__).resolve().parent
PAPER = HERE.parent / "paper" if (HERE.parent / "paper").is_dir() else HERE.parent
ZERO = (F(0), F(0), F(0))


def load(name):
    candidates = [folder / name for folder in (HERE, HERE / "c3", HERE / "c4")
                  if (folder / name).is_file()]
    assert len(candidates) == 1, (name, candidates)
    return json.loads(candidates[0].read_text(encoding="utf-8"))


def poly(c=0, q=0, qq=0):
    return (F(c), F(q), F(qq))


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def scale(a, s):
    return tuple(x * s for x in a)


def multiply(a, b):
    result = [F(0)] * 5
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i + j] += x * y
    assert result[3:] == [0, 0], result
    return tuple(result[:3])


def value(a, q):
    return sum(x * q**i for i, x in enumerate(a))


def nonnegative(a, qmin):
    assert a[2] == 0
    assert a[1] >= 0 and value(a, qmin) >= 0, (a, qmin)


def positive(a, qmin):
    assert a[2] == 0
    assert a[1] >= 0 and value(a, qmin) > 0, (a, qmin)


def sum_family(aq, bh, constant, low, high, qmin):
    """Sum aq*Q+bh*h+constant for low(Q)<=h<=high(Q)."""
    count = add(add(high, scale(low, -1)), poly(1))
    nonnegative(count, qmin)  # Zero-length progressions are allowed.
    return add(
        multiply(poly(constant, aq), count),
        scale(multiply(add(low, high), count), F(bh, 2)),
    )


def endpoint(aq, bh, constant, bound):
    return add(poly(constant, aq), scale(bound, bh))


def prove_family_positive(aq, bh, constant, low, high, qmin):
    positive(endpoint(aq, bh, constant, low), qmin)
    positive(endpoint(aq, bh, constant, high), qmin)


def raw_row(n0, nstep, j0, jstep, delta):
    n, j = poly(n0, nstep), poly(j0, jstep)
    parity = -1 if (n0 + j0) % 2 else 1
    result = []
    for i, d in zip((-1, 0, 1), delta):
        term = scale(add(poly(1 - parity), scale(n, -2 * parity)), F(3 * i, 4))
        term = add(term, scale(j, F(3 * (2 - 3 * abs(i)), 2)))
        result.append(add(term, poly(d)))
    assert tuple(sum(row[k] for row in result) for k in range(3)) == ZERO
    return result


def boundary_values(n0, nstep, j0, jstep, delta, mode, qmin):
    absolute = []
    for row in raw_row(n0, nstep, j0, jstep, delta):
        s = 1 if value(row, qmin) > 0 else -1
        row = scale(row, s)
        positive(row, qmin)
        absolute.append(row)
    a, b, c = sorted(absolute, key=lambda row: value(row, qmin))
    positive(add(b, scale(a, -1)), qmin)
    positive(add(c, scale(b, -1)), qmin)
    assert add(a, b) == c
    if mode == "c":
        return c, add(b, scale(a, -1))
    if mode == "a":
        return a, add(b, c)
    assert mode == "b"
    return b, add(a, c)


def bulk_coefficients(n0, nstep, jstar, period, delta):
    """Return H and D=M-L as affine forms (Q,h,constant)."""
    assert sum(delta) == 0
    if (n0 + jstar) % 2 == 0:
        low_delta, high_delta, nshift = delta[0], -delta[2], 0
    else:
        low_delta, high_delta, nshift = delta[2], -delta[0], 1
    high = (F(3 * nstep, 2), F(3 * period, 2),
            F(3 * (n0 + nshift + jstar), 2) + high_delta)
    diff = (-F(3 * nstep, 2), F(9 * period, 2),
            F(9 * jstar - 3 * (n0 + nshift), 2) + delta[1] - low_delta)
    # Both summands stay positive in these bulk zones.
    small = (F(3 * nstep, 2), -F(3 * period, 2),
             F(3 * (n0 + nshift - jstar), 2) + low_delta)
    middle = (F(0), F(3 * period), F(3 * jstar + delta[1]))
    assert high == tuple(x + y for x, y in zip(small, middle))
    return high, diff, small, middle


def pack(signs, keys):
    bits = "".join("1" if signs[key] == 1 else "0" for key in keys)
    bits += "0" * (-len(bits) % 4)
    return "".join(f"{int(bits[i:i + 4], 2):X}" for i in range(0, len(bits), 4))


def unpack(encoded, keys):
    bits = "".join(f"{int(digit, 16):04b}" for digit in encoded)
    assert len(bits) >= len(keys)
    return {key: 1 if bit == "1" else -1 for key, bit in zip(keys, bits)}


def verify_printed_signs():
    """Check that the symbols actually printed are the patterns being proved."""
    row_re = re.compile(
        r"^(\d+)&(\d+)&(\d+)&(\d+)&(.*?)&"
        r"\\texttt\{([0-9A-F]+)\}&\\texttt\{([0-9A-F]+)\}", re.M
    )
    patterns = load("lemma2_patterns.json")
    rows = row_re.findall((PAPER / "c3_finite_certificate_body.tex").read_text(encoding="utf-8"))
    assert len(rows) == 24
    for e, s, group, length, switch, c_hex, d_hex in rows:
        pattern = patterns[f"cls{e}_t{int(s):03d}_r{group}"]
        keys = sorted(pattern["c"])
        assert len(keys) == int(length)
        assert c_hex == pack(pattern["c"], keys)
        assert d_hex == pack(pattern["d"], keys)
        switch_key, mode = pattern["switch"]
        expected = "--" if switch_key == "None" else rf"\(K:{ast.literal_eval(switch_key)[2]};{mode}\)"
        assert switch.strip() == expected, (switch, expected)
    patterns = {(p["parameter_residue_mod_12"], p["step_state_mod_12"], p["group"]): p
                for p in load("c4_balance_affine_patterns.json")}
    rows = row_re.findall((PAPER / "c4_balance_patterns_appendix.tex").read_text(encoding="utf-8"))
    assert len(rows) == 192
    family_order = {"low": 0, "middle": 1, "top": 2, "word1": 3, "word2": 4}
    for rho, state, group, length, switch, c_hex, d_hex in rows:
        pattern = patterns[(int(rho), int(state), int(group))]
        def order(key):
            key = ast.literal_eval(key)
            return (family_order[key[0]],) + key[1:]
        keys = sorted(pattern["mu_signs"], key=order)
        assert len(keys) == int(length)
        pattern["mu_signs"] = unpack(c_hex, keys)
        pattern["w_signs"] = unpack(d_hex, keys)
        if switch.strip() == "--":
            pattern["switch_category"], pattern["switch_mode"] = None, 0
        else:
            sw = re.fullmatch(r"\$L_\{(-?\d+)\}:([ab])\$", switch.strip())
            assert sw, switch
            pattern["switch_category"] = ["low", int(sw[1])]
            pattern["switch_mode"] = {"a": 1, "b": 2}[sw[2]]
    print("Printed data: decoded all 24 c=3 and 192 c=4 sign rows and mode switches")
    return list(patterns.values())


def affine_index_bound(qcoef, constant, jstar, period, ceil=False):
    assert qcoef % period == 0
    numerator = constant - jstar
    c = -((-numerator) // period) if ceil else numerator // period
    return poly(c, qcoef // period)


def c3_rule(n0):
    if n0 % 24 == 9:
        return (load("joint_words_P24.json"),
                load("boundary_105.json" if n0 % 48 == 9 else "boundary_129.json"),
                (6, 12))
    return (load("words21.json"),
            load("boundary21_117.json" if n0 % 48 == 21 else "boundary21_141.json"),
            (10, 13))


def verify_printed_c3_support():
    """Link the printed c=3 words and boundaries to the symbolic inputs."""
    source = (PAPER / "c3_finite_certificate_body.tex").read_text(encoding="utf-8")
    rows = re.findall(r"\\text\{offset\}&([^\\]+)\\\\", source)
    offsets = [tuple(map(int, match)) for row in rows
               for match in re.findall(r"\((-?\d+),(-?\d+),(-?\d+)\)", row)]
    assert len(offsets) == 13 and all(sum(delta) == 0 for delta in offsets)
    alphabet = dict(zip("0123456789ABC", offsets))
    master_rows = re.findall(r"W_([12])&\\mathtt\{([0-9A-C]{24})\}", source)
    assert len(master_rows) == 2
    master = {int(i): word for i, word in master_rows}
    assert set(master) == {1, 2}
    rows = re.findall(
        r"^(9|21|33|45)&([LKT])&[^&]+&\\mathtt\{([0-9A-C]+)\}", source, re.M
    )
    assert len(rows) == 12
    printed_boundary = {(int(residue), name): word for residue, name, word in rows}
    assert len(printed_boundary) == 12
    word_checks, boundary_checks = 0, 0
    for e in (9, 21):
        shifts = (0, 0) if e == 9 else (16, 13)
        for tau in range(15, 27):
            n0 = e + 24 * tau
            if (n0 // 3) % 3 == 0:
                continue
            words, boundary, phases = c3_rule(n0)
            assert phases == ((6, 12) if e == 9 else (10, 13))
            for i in (1, 2):
                for p in range(24):
                    printed = alphabet[master[i][(p + shifts[i - 1]) % 24]]
                    assert printed == tuple(words[f"W{i}"][p]), (e, tau, i, p)
                    word_checks += 1
            for name, field, low, length in (("L", "TL", 1, 20),
                                             ("K", "TM", -20, 41),
                                             ("T", "TU", -20, 21)):
                word = printed_boundary[n0 % 48, name]
                assert len(word) == length
                for position, letter in enumerate(word, low):
                    assert alphabet[letter] == tuple(boundary[field][str(position)])
                    boundary_checks += 1
    assert (word_checks, boundary_checks) == (768, 1312)
    print("Printed c=3 support: all 768 word offsets and 1312 boundary offsets verified")


def decode_printed_polynomial(expression):
    expression = expression.replace(" ", "")
    terms = re.findall(r"([+-]?\d+)(Q\^2|Q)?", expression)
    assert "".join(a + b for a, b in terms) == expression, expression
    coefficients = [F(0), F(0), F(0)]
    for coefficient, power in terms:
        coefficients[{"": 0, "Q": 1, "Q^2": 2}[power]] += F(coefficient)
    return tuple(coefficients)


def printed_worked_c3():
    """Decode the polynomial table in the manuscript's worked packet example."""
    source = (PAPER / "IHS_m_3_3_standalone.tex").read_text(encoding="utf-8")
    example = source.split(r"\paragraph{A packet-sum calculation.}", 1)[1]
    example = example.split(r"\begin{lemma}", 1)[0]

    rows = re.findall(r"^(L,M,K,T|B_[123])&([^&]+)&([^\\]+)\\\\", example, re.M)
    assert len(rows) == 4
    totals = {name: [decode_printed_polynomial(mu), decode_printed_polynomial(nu)]
              for name, mu, nu in rows}
    packet = {}
    for line in example.splitlines():
        for chain in ("mu", "nu"):
            if line.strip().startswith(rf"\sum\{chain}_j&="):
                packet[chain] = decode_printed_polynomial(
                    line.rsplit("=", 1)[1].rstrip("\\, ."))
    assert set(packet) == {"mu", "nu"}
    return totals, packet


def printed_worked_c4():
    """Decode the printed branch sums, not a numerical fit to sample levels."""
    source = (PAPER / "IHS_m_3_3_standalone.tex").read_text(encoding="utf-8")
    example = source.split(r"\paragraph{A base-level branch change.}", 1)[1]
    example = example.split(r"\begin{proof}", 1)[0]
    rows = re.findall(r"\\sum\\(mu|nu)_j&=([^\\]+)", example)
    assert len(rows) == 2
    totals = {name: decode_printed_polynomial(expression.rstrip(", .\n\r"))
              for name, expression in rows}
    branches = re.search(r"give \$\\nu\$-totals \$([^$]+)\$ and \$([^$]+)\$", example)
    assert branches
    return totals, tuple(decode_printed_polynomial(x) for x in branches.groups())


def verify_c3():
    patterns = load("lemma2_patterns.json")
    printed_totals, printed_packet = printed_worked_c3()
    checked = 0
    for e in (9, 21):
        for tau in range(15, 27):
            n0 = e + 24 * tau
            if (n0 // 3) % 3 == 0:
                continue
            m0, k0 = n0 // 3, (n0 + 1) // 2
            words, boundary, phases = c3_rule(n0)
            s = 15 + (tau - 15) % 6
            for group in range(3):
                pattern = patterns[f"cls{e}_t{s:03d}_r{group}"]
                signs = [{ast.literal_eval(k): v for k, v in pattern[name].items()}
                         for name in ("c", "d")]
                switch, mode = pattern["switch"]
                switch = None if switch == "None" else ast.literal_eval(switch)
                sums, seen = [ZERO, ZERO], set()
                pieces = {name: [ZERO, ZERO] for name in ("L,M,K,T", "B_1", "B_2", "B_3")}
                for name, offsets, jbase, jstep in (
                    ("low", range(1, 21), 0, 0),
                    ("m", range(-12, 13), m0, 96),
                    ("k", range(-20, 21), k0, 144),
                    ("top", range(-20, 1), n0, 288),
                ):
                    for offset in offsets:
                        j = jbase + offset
                        if j % 3 != group:
                            continue
                        key = ("bnd", name, offset)
                        if name == "m":
                            delta = words["W1"][(j - phases[0]) % 24]
                        else:
                            field = {"low": "TL", "k": "TM", "top": "TU"}[name]
                            delta = boundary[field][str(offset)]
                        vals = boundary_values(n0, 288, j, jstep, delta,
                                               mode if key == switch else "c", 0)
                        for chain in range(2):
                            contribution = scale(vals[chain], signs[chain][key])
                            sums[chain] = add(sums[chain], contribution)
                            pieces["L,M,K,T"][chain] = add(pieces["L,M,K,T"][chain], contribution)
                        seen.add(key)
                for zone in (1, 2, 3):
                    for p in range(24):
                        jstar = p or 24
                        if jstar % 3 != group:
                            continue
                        key = ("bulk", zone, p)
                        word = 0 if zone < 3 else 1
                        delta = words[f"W{word + 1}"][(jstar - phases[word]) % 24]
                        if zone == 1:
                            lo = affine_index_bound(0, 21, jstar, 24, True)
                            hi = affine_index_bound(96, m0 - 13, jstar, 24)
                        elif zone == 2:
                            lo = affine_index_bound(96, m0 + 13, jstar, 24, True)
                            hi = affine_index_bound(144, k0 - 21, jstar, 24)
                        else:
                            lo = affine_index_bound(144, k0 + 21, jstar, 24, True)
                            hi = affine_index_bound(288, n0 - 21, jstar, 24)
                        high, diff, small, middle = bulk_coefficients(n0, 288, jstar, 24, delta)
                        prove_family_positive(*small, lo, hi, 0)
                        prove_family_positive(*middle, lo, hi, 0)
                        ds = -1 if zone == 1 else 1
                        prove_family_positive(*(ds * x for x in diff), lo, hi, 0)
                        vals = (sum_family(*high, lo, hi, 0),
                                scale(sum_family(*diff, lo, hi, 0), ds))
                        if (e, tau, group, zone, p) == (21, 15, 0, 2, 0):
                            assert vals == (printed_packet["mu"], printed_packet["nu"])
                            assert signs[0][key] == signs[1][key] == 1
                        for chain in range(2):
                            contribution = scale(vals[chain], signs[chain][key])
                            sums[chain] = add(sums[chain], contribution)
                            pieces[f"B_{zone}"][chain] = add(pieces[f"B_{zone}"][chain], contribution)
                        seen.add(key)
                assert set(signs[0]) == set(signs[1]) == seen
                assert sums == [ZERO, ZERO], (e, tau, group, sums)
                if (e, tau, group) == (21, 15, 0):
                    assert pieces == printed_totals, (pieces, printed_totals)
                checked += 1
    assert checked == 48
    print("c=3: 48 affine group states; all 288 signed coefficients are exactly zero")
    print("c=3: printed worked packet and all four polynomial-table rows verified exactly")


def c4_rule(rho, state):
    minimum = {1: 25, 5: 29, 7: 31, 11: 23}[rho]
    m0 = minimum + 12 * state
    if rho == 1:
        rule = load("c4_period48_class4_rule.json")
        odd = load("c4_period48_class4_boundary_odd.json")
        base, ph = 37, (24, 36)
    else:
        nbase = {5: 356, 7: 364, 11: 380}[rho]
        rule = load(f"c4_period48_base{nbase}_rule.json")
        odd = load(f"c4_period48_base{nbase}_boundary_odd.json")
        base, ph = nbase // 4, (0, 0)
    step = (m0 - base) // 12
    phase = tuple((x + 16 * step) % 48 for x in ph)
    boundary = rule if step % 2 == 0 else odd
    return m0, phase, rule, boundary


def verify_printed_c4_support():
    """Link the two printed short words to all JSON support refinements."""
    source = (PAPER / "c4_main_certificate_body.tex").read_text(encoding="utf-8")
    alphabet_rows = re.findall(r"\\delta&([^\\]+)\\\\", source.split(r"For $\rho", 1)[0])
    offsets = [tuple(map(int, match)) for row in alphabet_rows
               for match in re.findall(r"\((-?\d+),(-?\d+),(-?\d+)\)", row)]
    assert len(offsets) == 12 and all(sum(delta) == 0 for delta in offsets)
    alphabet = dict(zip("0123456789AB", offsets))
    words = {int(i): word for i, word in
             re.findall(r"P_([12])&\\texttt\{([0-9A-B]{16})\}", source)}
    assert set(words) == {1, 2}
    phase_rows = re.findall(
        r"^(1|5|7|11)&(\d+)&\((\d+),(\d+)\)&\((\d+),(\d+)\)\\\\", source, re.M
    )
    assert len(phase_rows) == 4
    boundary_rows = re.findall(
        r"^(1|5|7|11)&([01])&([LMT])&[^&]+&\\texttt\{([0-9A-B]+)\}", source, re.M
    )
    assert len(boundary_rows) == 24
    printed_boundaries = {(int(rho), int(parity), name): word
                          for rho, parity, name, word in boundary_rows}
    checks = 0
    for rho, minimum, psi1, psi2, phi1, phi2 in phase_rows:
        rho, minimum = int(rho), int(minimum)
        psi, phi = (int(psi1), int(psi2)), (int(phi1), int(phi2))
        for state in range(12):
            m0, phases, rule, boundary = c4_rule(rho, state)
            assert m0 == minimum + 12 * state
            assert phases == tuple((p + 16 * state) % 48 for p in phi)
            for zone in (1, 2):
                for j in range(1, 49):
                    printed = alphabet[words[zone][(j - psi[zone - 1]) % 16]]
                    stored = tuple(rule[f"word{zone}"][str((j - phases[zone - 1]) % 48)])
                    assert printed == stored, (rho, state, zone, j)
                    checks += 1
            for name, key, low in (("L", "low", 1), ("M", "middle", -20),
                                   ("T", "top", -20)):
                word = printed_boundaries[(rho, state % 2, name)]
                expected_length = {"L": 20, "M": 41, "T": 21}[name]
                assert len(word) == expected_length
                for position, letter in enumerate(word, low):
                    assert alphabet[letter] == tuple(boundary[key][str(position)])
    assert checks == 4608
    print("Printed c=4 support: two 16-letter words, phases and 24 boundary strings verified")
    print("c=4: all 4608 short-word offsets match the 48-position refinements")


def verify_c4(patterns):
    printed_totals, printed_branches = printed_worked_c4()
    checked, base_checks = 0, 0
    example_checked = False
    for pattern in patterns:
        rho, state, group = (pattern[name] for name in
                              ("parameter_residue_mod_12", "step_state_mod_12", "group"))
        m0, phases, rule, boundary = c4_rule(rho, state)
        assert m0 == pattern["base_m"] and list(phases) == pattern["phase"]
        n0, k0 = 4 * m0, 2 * m0
        signs = [{ast.literal_eval(k): v for k, v in pattern[f"{name}_signs"].items()}
                 for name in ("mu", "w")]
        switch = pattern["switch_category"]
        switch = None if switch is None else tuple(switch)
        mode = ("c", "a", "b")[pattern["switch_mode"]]
        sums, seen = [ZERO, ZERO], set()
        for name, offsets, jbase, jstep in (
            ("low", range(1, 21), 0, 0),
            ("middle", range(-20, 21), k0, 288),
            ("top", range(-20, 1), n0, 576),
        ):
            for offset in offsets:
                j = jbase + offset
                if (j - 1) % 4 != group:
                    continue
                key = (name, offset)
                vals = boundary_values(n0, 576, j, jstep, boundary[name][str(offset)],
                                       mode if key == switch else "c", 1)
                for chain in range(2):
                    sums[chain] = add(sums[chain], scale(vals[chain], signs[chain][key]))
                seen.add(key)
        for zone in (1, 2):
            for p in range(48):
                j0 = (phases[zone - 1] + p - 1) % 48 + 1
                if (j0 - 1) % 4 != group:
                    continue
                for layer in (0, 1):
                    key, jstar = (f"word{zone}", p, layer), j0 + 48 * layer
                    delta = rule[f"word{zone}"][str(p)]
                    if zone == 1:
                        lo = affine_index_bound(0, 21, jstar, 96, True)
                        hi = affine_index_bound(288, k0 - 21, jstar, 96)
                    else:
                        lo = affine_index_bound(288, k0 + 21, jstar, 96, True)
                        hi = affine_index_bound(576, n0 - 21, jstar, 96)
                    high, diff, small, middle = bulk_coefficients(n0, 576, jstar, 96, delta)
                    prove_family_positive(*small, lo, hi, 1)
                    prove_family_positive(*middle, lo, hi, 1)
                    mu = sum_family(*high, lo, hi, 1)
                    if zone == 1:
                        assert diff[:2] == (-864, 432)
                        cutoff = poly((-diff[2]) // 432, 2)
                        w = add(scale(sum_family(*diff, lo, cutoff, 1), -1),
                                sum_family(*diff, add(cutoff, poly(1)), hi, 1))
                    else:
                        prove_family_positive(*diff, lo, hi, 1)
                        w = sum_family(*diff, lo, hi, 1)
                    if (rho, state, group, zone, p, layer) == (11, 0, 1, 1, 26, 0):
                        assert jstar == 26 and tuple(delta) == (0, -1, 1)
                        assert lo == poly() and hi == poly(-1, 3)
                        assert high == (864, 144, 176)
                        assert diff == (-864, 432, -22)
                        assert cutoff == poly(0, 2)
                        branches = (
                            scale(sum_family(*diff, lo, cutoff, 1), -1),
                            sum_family(*diff, add(cutoff, poly(1)), hi, 1),
                        )
                        assert branches == printed_branches
                        assert (mu, w) == (printed_totals["mu"], printed_totals["nu"])
                        assert signs[0][key] == signs[1][key] == 1
                        # At Q=0 the original category is empty, while the
                        # second refined branch has count -1 and cannot be summed.
                        assert value(lo, 0) == 0 and value(hi, 0) == -1
                        assert value(add(hi, scale(cutoff, -1)), 0) == -1
                        assert value(w, 0) == 44
                        example_checked = True
                    for chain, result in enumerate((mu, w)):
                        sums[chain] = add(sums[chain], scale(result, signs[chain][key]))
                    seen.add(key)
        assert set(signs[0]) == set(signs[1]) == seen
        assert len(seen) in (68, 69)
        assert sums == [ZERO, ZERO], (rho, state, group, sums)
        # Separate base equations, with no polynomial extrapolation through Q=0.
        base_sums = [F(0), F(0)]
        for j in range(group + 1, n0 + 1, 4):
            if j <= 20:
                key, delta = ("low", j), boundary["low"][str(j)]
            elif j <= k0 - 21:
                zone = 1
                p = (j - phases[0]) % 48
                j0 = (j - 1) % 48 + 1
                key, delta = ("word1", p, ((j - j0) // 48) % 2), rule["word1"][str(p)]
            elif j <= k0 + 20:
                key, delta = ("middle", j - k0), boundary["middle"][str(j - k0)]
            elif j <= n0 - 21:
                p = (j - phases[1]) % 48
                j0 = (j - 1) % 48 + 1
                key, delta = ("word2", p, ((j - j0) // 48) % 2), rule["word2"][str(p)]
            else:
                key, delta = ("top", j - n0), boundary["top"][str(j - n0)]
            a, b, c = sorted(abs(row[0]) for row in raw_row(n0, 0, j, 0, delta))
            assert a > 0 and a + b == c and a < b
            md = mode if key == switch else "c"
            vals = {"c": (c, b - a), "a": (a, b + c), "b": (b, a + c)}[md]
            for chain in range(2):
                base_sums[chain] += vals[chain] * signs[chain][key]
        assert base_sums == [0, 0], (rho, state, group, base_sums)
        checked += 1
        base_checks += 2
    assert checked == 192
    assert example_checked
    print("c=4: 192 affine group states; all 1152 signed coefficients are exactly zero")
    print(f"c=4: {base_checks} separate signed base equations at Q=0 are exactly zero")
    print("c=4: printed branch sums verified exactly; empty Q=0 packet is not extrapolated")


if __name__ == "__main__":
    verify_printed_c3_support()
    verify_printed_c4_support()
    c4_patterns = verify_printed_signs()
    verify_c3()
    verify_c4(c4_patterns)
    print("PASS: exact symbolic coefficient proof; no finite-difference inference")
