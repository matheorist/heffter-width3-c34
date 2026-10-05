"""Prove the recovered c=4 support rules by finite affine interval tilings.

For each of 4 parameter classes and 12 step residues, write step=s0+12Q.
All phases, boundary choices, and parities are then fixed.  Every bulk output,
after splitting a 72-step family by layer parity, is an arithmetic progression
with difference 144.  For each residue modulo 144 we verify symbolically that
the quotient intervals tile the target interval with adjacent affine endpoints.
"""

from fractions import Fraction

from generate_c4_all import PERIODIC_MINIMUM, load_periodic_rule


def ceil_div(a, b):
    return -((-a) // b)


def zhou(n, i, j):
    modulus = 2 * n + 1
    parity = -1 if (j + n) % 2 else 1
    return int(
        Fraction(1 - parity * modulus, 4) * i
        + Fraction(2 - 3 * abs(i), 2) * j
    )


def output(n, j, delta, coordinate):
    bases = (zhou(n, -1, j), j, zhou(n, 1, j))
    return abs(3 * bases[coordinate] + delta[coordinate])


def affine_bulk(n0, j0, h_sample, delta, coordinate):
    # n=n0+576Q and j=j0+48h.  Sample inside a stable bulk zone.
    q = 2
    h = h_sample
    value = output(n0 + 576 * q, j0 + 48 * h, delta, coordinate)
    aq = output(n0 + 576 * (q + 1), j0 + 48 * h, delta, coordinate) - value
    bh = output(n0 + 576 * q, j0 + 48 * (h + 1), delta, coordinate) - value
    constant = value - aq * q - bh * h
    return aq, bh, constant


def affine_boundary(n0, j_function, delta, coordinate):
    values = [
        output(n0 + 576 * q, j_function(q), delta, coordinate) for q in (2, 3)
    ]
    aq = values[1] - values[0]
    return aq, values[0] - 2 * aq


def residue_and_quotient(aq, constant):
    residue = ((constant - 1) % 144) + 1
    assert aq % 144 == 0 and (constant - residue) % 144 == 0
    return residue, aq // 144, (constant - residue) // 144


def add_bulk_intervals(intervals, n0, phase, word, zone):
    kappa0 = n0 // 2
    for word_position, delta in word.items():
        j0 = ((phase + word_position - 1) % 48) + 1
        if zone == 1:
            h_lower = (0, ceil_div(21 - j0, 48))
            h_upper = (6, (kappa0 - 21 - j0) // 48)
        else:
            h_lower = (6, ceil_div(kappa0 + 21 - j0, 48))
            h_upper = (12, (n0 - 21 - j0) // 48)

        h_sample = h_lower[0] * 2 + h_lower[1] + 2
        for coordinate in range(3):
            aq, bh, constant = affine_bulk(
                n0, j0, h_sample, delta, coordinate
            )
            assert abs(bh) in (72, 144)
            parities = (None,) if abs(bh) == 144 else (0, 1)
            for parity in parities:
                if parity is None:
                    low = h_lower
                    high = h_upper
                    step_multiplier = 1
                else:
                    # Q coefficients are even (0, 6, or 12), so endpoint parity
                    # corrections depend only on the constant terms.
                    low = (
                        h_lower[0],
                        h_lower[1] + ((parity - h_lower[1]) % 2),
                    )
                    high = (
                        h_upper[0],
                        h_upper[1] - ((h_upper[1] - parity) % 2),
                    )
                    step_multiplier = 2

                def x_at(endpoint):
                    return (
                        aq + bh * endpoint[0],
                        constant + bh * endpoint[1],
                    )

                first = x_at(low)
                last = x_at(high)
                if bh < 0:
                    first, last = last, first
                # Splitting a 72-family by parity, or retaining a 144-family,
                # always gives step +144 in increasing order.
                assert abs(bh * step_multiplier) == 144
                residue, low_a, low_c = residue_and_quotient(*first)
                residue2, high_a, high_c = residue_and_quotient(*last)
                assert residue == residue2
                intervals.setdefault(residue, []).append(
                    ((low_a, low_c), (high_a, high_c), f"z{zone}:p{word_position}:c{coordinate}:e{parity}")
                )


def add_boundary_intervals(intervals, n0, boundary):
    kappa0 = n0 // 2
    for section, table in enumerate(boundary):
        for relative, delta in table.items():
            if section == 0:
                j_function = lambda q, relative=relative: relative
            elif section == 1:
                j_function = lambda q, relative=relative: kappa0 + 288 * q + relative
            else:
                j_function = lambda q, relative=relative: n0 + 576 * q + relative
            for coordinate in range(3):
                aq, constant = affine_boundary(
                    n0, j_function, delta, coordinate
                )
                residue, qa, qc = residue_and_quotient(aq, constant)
                intervals.setdefault(residue, []).append(
                    ((qa, qc), (qa, qc), f"b{section}:{relative}:c{coordinate}")
                )


def evaluate(endpoint, q):
    return endpoint[0] * q + endpoint[1]


def prove_state(parameter_residue, state):
    minimum = PERIODIC_MINIMUM[parameter_residue]
    m0 = minimum + 12 * state
    n0 = 4 * m0
    base_m, base_phase, words, even_boundary, odd_boundary = load_periodic_rule(
        parameter_residue
    )
    original_step = (m0 - base_m) // 12
    phase = (
        (base_phase[0] + 16 * original_step) % 48,
        (base_phase[1] + 16 * original_step) % 48,
    )
    boundary = even_boundary if original_step % 2 == 0 else odd_boundary

    intervals = {}
    add_bulk_intervals(intervals, n0, phase[0], words[0], 1)
    add_bulk_intervals(intervals, n0, phase[1], words[1], 2)
    add_boundary_intervals(intervals, n0, boundary)

    target_length_constant = 3 * n0
    for residue in range(1, 145):
        pieces = intervals.get(residue, [])
        # Q=2 lies safely beyond any initially empty short bulk family.
        pieces.sort(key=lambda item: evaluate(item[0], 2))
        assert pieces, (parameter_residue, state, residue, "empty")
        expected_low = (0, 0)
        for low, high, label in pieces:
            assert low == expected_low, (
                parameter_residue,
                state,
                residue,
                "gap/overlap",
                expected_low,
                low,
                label,
            )
            # The next interval must begin one quotient after this one ends.
            expected_low = (high[0], high[1] + 1)
        target_high = (
            12,
            (target_length_constant - residue) // 144,
        )
        assert (expected_low[0], expected_low[1] - 1) == target_high, (
            parameter_residue,
            state,
            residue,
            "wrong final endpoint",
            (expected_low[0], expected_low[1] - 1),
            target_high,
        )

    # Directly check Q=0, which also covers any family empty only at the base.
    from generate_c4_all import corrected_triples

    triples, _, _ = corrected_triples(m0)
    assert sorted(value for triple in triples for value in triple) == list(
        range(1, 12 * m0 + 1)
    )


def main():
    count = 0
    for parameter_residue in (1, 5, 7, 11):
        for state in range(12):
            prove_state(parameter_residue, state)
            count += 1
            print(f"PROVED residue={parameter_residue}, step_state={state}")
    print(f"AFFINE SUPPORT PROOF VERIFIED: {count} states")


if __name__ == "__main__":
    main()
