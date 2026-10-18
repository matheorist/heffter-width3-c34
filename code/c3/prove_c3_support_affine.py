"""Prove the c=3 support rules by exact affine expansion and interval tilings.

Coordinate signs are certified throughout Q>=0 before absolute values are
removed. No coefficients are inferred from sampled absolute-value outputs.
"""

from generate_c3_all import rule_data, support_rows


MODULUS = 72
N_STEP, H_STEP = 288, 24


def admissible_state(cls, state):
    """Return whether n=cls+24*state has gcd(n/3,3)=1."""
    return state % 3 != (0 if cls == 9 else 1)


def ceil_div(a, b):
    return -((-a) // b)


def raw_affine(n0, j0, jq, jh, delta, coordinate):
    """Expand the signed Kotzig lift for j=j0+jq*Q+jh*h exactly."""
    assert N_STEP % 2 == jq % 2 == jh % 2 == 0
    if coordinate == 1:
        return 3 * jq, 3 * jh, 3 * j0 + delta[1]
    i = (-1, 0, 1)[coordinate]
    parity = -1 if (n0 + j0) % 2 else 1
    numerators = (
        -6 * i * parity * N_STEP - 6 * jq,
        -6 * jh,
        3 * i * (1 - parity * (2 * n0 + 1)) - 6 * j0 + 4 * delta[coordinate],
    )
    assert all(x % 4 == 0 for x in numerators)
    return tuple(x // 4 for x in numerators)


def affine_bulk(n0, j0, h_lower, h_upper, delta, coordinate):
    aq, bh, constant = raw_affine(n0, j0, 0, H_STEP, delta, coordinate)
    parity = -1 if (n0 + j0) % 2 else 1
    sign = (parity, 1, -parity)[coordinate]
    # The smaller signed endpoint is a lower bound whenever the packet
    # is nonempty. An empty packet imposes no further endpoint condition.
    endpoint = h_lower if sign * bh >= 0 else h_upper
    min_a = sign * (aq + bh * endpoint[0])
    min_c = sign * (constant + bh * endpoint[1])
    assert min_a >= 0 and min_c > 0, (n0, j0, coordinate, min_a, min_c)
    return sign * aq, sign * bh, sign * constant


def affine_boundary(n0, j_function, delta, coordinate):
    # Each caller explicitly supplies an affine index, not sampled output.
    j0 = j_function(0)
    jq = j_function(1) - j0
    aq, bh, constant = raw_affine(n0, j0, jq, 0, delta, coordinate)
    assert bh == 0 and constant != 0
    sign = 1 if constant > 0 else -1
    assert sign * aq >= 0
    return sign * aq, sign * constant


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

        for coordinate in range(3):
            aq, bh, constant = affine_bulk(
                n0, j0, h_lower, h_upper, delta, coordinate
            )
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
            # An affine interval may be empty at Q=0, but its length
            # must never be negative anywhere on the Q>=0 tail.
            assert high[0] - low[0] >= 0 and high[1] - low[1] + 1 >= 0, (
                cls, state, residue, "negative packet length", low, high, label
            )
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
