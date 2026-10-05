"""Prove the c=3 periodic support rules by finite affine interval tilings."""

from generate_c3_all import rule_data, support_rows, zhou_value


MODULUS = 72


def admissible_state(cls, state):
    """Return whether n=cls+24*state has gcd(n/3,3)=1."""
    return state % 3 != (0 if cls == 9 else 1)


def ceil_div(a, b):
    return -((-a) // b)


def output(n, j, delta, coordinate):
    bases = (zhou_value(n, -1, j), j, zhou_value(n, 1, j))
    return abs(3 * bases[coordinate] + delta[coordinate])


def affine_bulk(n0, j0, h_sample, delta, coordinate):
    # A fixed state has n=n0+288Q and j=j0+24h.
    q = 2
    h = h_sample
    value = output(n0 + 288 * q, j0 + 24 * h, delta, coordinate)
    aq = output(n0 + 288 * (q + 1), j0 + 24 * h, delta, coordinate) - value
    bh = output(n0 + 288 * q, j0 + 24 * (h + 1), delta, coordinate) - value
    return aq, bh, value - aq * q - bh * h


def affine_boundary(n0, j_function, delta, coordinate):
    values = [output(n0 + 288 * q, j_function(q), delta, coordinate) for q in (2, 3)]
    aq = values[1] - values[0]
    return aq, values[0] - 2 * aq


def residue_and_quotient(aq, constant):
    residue = ((constant - 1) % MODULUS) + 1
    assert aq % MODULUS == 0 and (constant - residue) % MODULUS == 0
    return residue, aq // MODULUS, (constant - residue) // MODULUS


def add_bulk_intervals(intervals, n0, phase, word, zone):
    kappa0 = (n0 + 1) // 2
    for word_position, delta in enumerate(word):
        j0 = ((phase + word_position - 1) % 24) + 1
        if zone == 1:
            h_lower = (0, ceil_div(21 - j0, 24))
            h_upper = (6, (kappa0 - 21 - j0) // 24)
        else:
            h_lower = (6, ceil_div(kappa0 + 21 - j0, 24))
            h_upper = (12, (n0 - 21 - j0) // 24)

        h_sample = h_lower[0] * 2 + h_lower[1] + 2
        for coordinate in range(3):
            aq, bh, constant = affine_bulk(n0, j0, h_sample, delta, coordinate)
            assert abs(bh) in (36, 72)
            parities = (None,) if abs(bh) == 72 else (0, 1)
            for parity in parities:
                if parity is None:
                    low, high, step_multiplier = h_lower, h_upper, 1
                else:
                    low = (h_lower[0], h_lower[1] + ((parity - h_lower[1]) % 2))
                    high = (h_upper[0], h_upper[1] - ((h_upper[1] - parity) % 2))
                    step_multiplier = 2

                def x_at(endpoint):
                    return aq + bh * endpoint[0], constant + bh * endpoint[1]

                first, last = x_at(low), x_at(high)
                if bh < 0:
                    first, last = last, first
                assert abs(bh * step_multiplier) == MODULUS
                residue, low_a, low_c = residue_and_quotient(*first)
                residue2, high_a, high_c = residue_and_quotient(*last)
                assert residue == residue2
                intervals.setdefault(residue, []).append(
                    ((low_a, low_c), (high_a, high_c),
                     f"z{zone}:p{word_position}:c{coordinate}:e{parity}")
                )


def add_boundary_intervals(intervals, n0, boundary):
    kappa0 = (n0 + 1) // 2
    sections = (boundary["TL"], boundary["TM"], boundary["TU"])
    for section, table in enumerate(sections):
        for relative_text, delta in table.items():
            relative = int(relative_text)
            if section == 0:
                j_function = lambda q, relative=relative: relative
            elif section == 1:
                j_function = lambda q, relative=relative: kappa0 + 144 * q + relative
            else:
                j_function = lambda q, relative=relative: n0 + 288 * q + relative
            for coordinate in range(3):
                aq, constant = affine_boundary(n0, j_function, delta, coordinate)
                residue, qa, qc = residue_and_quotient(aq, constant)
                intervals.setdefault(residue, []).append(
                    ((qa, qc), (qa, qc), f"b{section}:{relative}:c{coordinate}")
                )


def evaluate(endpoint, q):
    return endpoint[0] * q + endpoint[1]


def prove_state(cls, state):
    # state is an actual t in 15..26; fixing t mod 12 fixes all parities.
    n0 = 24 * state + cls
    words, boundary, phases = rule_data(n0)
    intervals = {}
    add_bulk_intervals(intervals, n0, phases[0], words["W1"], 1)
    add_bulk_intervals(intervals, n0, phases[1], words["W2"], 2)
    add_boundary_intervals(intervals, n0, boundary)

    for residue in range(1, MODULUS + 1):
        pieces = intervals.get(residue, [])
        pieces.sort(key=lambda item: evaluate(item[0], 2))
        assert pieces, (cls, state, residue, "empty")
        expected_low = (0, 0)
        for low, high, label in pieces:
            assert low == expected_low, (
                cls, state, residue, "gap/overlap", expected_low, low, label
            )
            expected_low = (high[0], high[1] + 1)
        target_high = (12, (3 * n0 - residue) // MODULUS)
        assert (expected_low[0], expected_low[1] - 1) == target_high, (
            cls, state, residue, "wrong final endpoint",
            (expected_low[0], expected_low[1] - 1), target_high
        )

    support = sorted(abs(value) for row in support_rows(n0) for value in row)
    assert support == list(range(1, 3 * n0 + 1))


def main():
    count = 0
    for cls in (9, 21):
        for state in range(15, 27):
            if not admissible_state(cls, state):
                continue
            prove_state(cls, state)
            count += 1
            print(f"PROVED class={cls}, t_state={state}")
    print(f"AFFINE C3 SUPPORT PROOF VERIFIED: {count} states")


if __name__ == "__main__":
    main()
