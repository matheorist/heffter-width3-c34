"""Audit the canonical c=3 transport m -> m+4.

The step is admissible precisely for m = 7 (mod 12).  This script checks
the half-period relation between the two printed period-24 bulk words, the
four boundary-state transitions, the exact inherited-row identity, and the
complete IHS support/balance conditions over a user-specified finite range.
No search is performed.
"""

import argparse

from generate_c3_all import generate, rule_data, support_rows, verify, zhou_value


def correction(n, j):
    words, boundary, phases = rule_data(n)
    kappa = (n + 1) // 2
    if j <= 20:
        return tuple(boundary["TL"][str(j)])
    if j <= kappa - 21:
        return tuple(words["W1"][(j - phases[0]) % 24])
    if j <= kappa + 20:
        return tuple(boundary["TM"][str(j - kappa)])
    if j <= n - 21:
        return tuple(words["W2"][(j - phases[1]) % 24])
    return tuple(boundary["TU"][str(j - n)])


def check_word_relation():
    words9, _, phases9 = rule_data(9)
    words21, _, phases21 = rule_data(21)
    assert phases9 == (6, 12)
    assert phases21 == (10, 13)
    assert all(
        tuple(words21["W1"][p]) == tuple(words9["W1"][(p + 16) % 24])
        for p in range(24)
    )
    assert all(
        tuple(words21["W2"][p]) == tuple(words9["W2"][(p + 13) % 24])
        for p in range(24)
    )


def check_pair(m):
    if m % 12 != 7:
        raise ValueError("the +4 target is admissible only for m = 7 (mod 12)")
    source_n = 3 * m
    target_n = 3 * (m + 4)
    source_rows = support_rows(source_n)
    target_rows = support_rows(target_n)
    for j in range(1, source_n + 1):
        source_delta = correction(source_n, j)
        target_delta = correction(target_n, j)
        sign = -1 if (j + source_n) % 2 else 1
        predicted = [
            source_rows[j - 1][0] + 18 * sign
            + target_delta[0] - source_delta[0],
            source_rows[j - 1][1]
            + target_delta[1] - source_delta[1],
            source_rows[j - 1][2] - 18 * sign
            + target_delta[2] - source_delta[2],
        ]
        assert predicted == target_rows[j - 1], (m, j, predicted, target_rows[j - 1])

    source_arrays = generate(m)
    target_arrays = generate(m + 4)
    verify(m, source_arrays)
    verify(m + 4, target_arrays)

    new_counts = [0, 0, 0]
    for j in range(source_n + 1, target_n + 1):
        new_counts[j % 3] += 1
    assert new_counts == [4, 4, 4]
    return source_n % 48, target_n % 48


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--minimum", type=int, default=67)
    parser.add_argument("--maximum", type=int, default=1003)
    args = parser.parse_args()
    check_word_relation()
    pairs = []
    transitions = set()
    for m in range(args.minimum, args.maximum + 1):
        if m % 12 != 7:
            continue
        transitions.add(check_pair(m))
        pairs.append((m, m + 4))
    assert transitions == {(9, 21), (21, 33), (33, 45), (45, 9)}
    print(
        "VERIFIED canonical c=3 +4 transport for "
        f"{len(pairs)} pairs from {pairs[0]} through {pairs[-1]}; "
        "bulk shift=12; boundary cycle=9->21->33->45->9; "
        "new rows per group=(4,4,4)"
    )


if __name__ == "__main__":
    main()
