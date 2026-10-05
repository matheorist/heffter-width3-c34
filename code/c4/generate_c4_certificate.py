"""Generate a certificate for IHS(m,3;4) from a corrected Zhou lift.

For m <= 7 all 3^m mode assignments are tested.  For larger m the search is
ordered by the number of rows whose mode is not c; this reflects the observed
bulk-c structure and is deterministic.
"""

import argparse
import itertools
import json
from pathlib import Path

from sat_batch import zhou_partition_sat


def signed_subset(values):
    """Return indices of a half-sum subset, or None if no signed zero-sum exists."""
    total = sum(values)
    if total % 2:
        return None
    target = total // 2
    bits = 1
    history = []
    for value in values:
        history.append(bits)
        bits |= bits << value
    if not ((bits >> target) & 1):
        return None

    result = []
    for index in range(len(values) - 1, -1, -1):
        value = values[index]
        if target >= value and ((history[index] >> (target - value)) & 1):
            result.append(index)
            target -= value
    assert target == 0
    return result


def mode_values(triple, mode):
    a, b, c = triple
    if mode == 0:  # c in column 3
        return c, b - a
    if mode == 1:  # a in column 3
        return a, b + c
    return b, a + c  # b in column 3


def test_modes(group, modes):
    pairs = [mode_values(triple, mode) for triple, mode in zip(group, modes)]
    mu_subset = signed_subset([pair[0] for pair in pairs])
    if mu_subset is None:
        return None
    w_subset = signed_subset([pair[1] for pair in pairs])
    if w_subset is None:
        return None
    return list(modes), mu_subset, w_subset


def find_modes(group, max_switches):
    m = len(group)
    if m <= 7:
        for modes in itertools.product(range(3), repeat=m):
            certificate = test_modes(group, modes)
            if certificate is not None:
                return certificate
        return None

    for weight in range(max_switches + 1):
        for indices in itertools.combinations(range(m), weight):
            for choices in itertools.product((1, 2), repeat=weight):
                modes = [0] * m
                for index, mode in zip(indices, choices):
                    modes[index] = mode
                certificate = test_modes(group, modes)
                if certificate is not None:
                    return certificate
    return None


def oriented_row(triple, mode, mu_positive, w_positive):
    mu, w = mode_values(triple, mode)
    row_sum_12 = mu if mu_positive else -mu
    row_diff_12 = w if w_positive else -w
    assert (row_sum_12 + row_diff_12) % 2 == 0
    row = [
        (row_sum_12 + row_diff_12) // 2,
        (row_sum_12 - row_diff_12) // 2,
        -row_sum_12,
    ]
    assert sum(row) == 0
    assert sorted(map(abs, row)) == list(triple)
    return row


def build_certificate(m, radius, max_switches):
    n = 4 * m
    triples = zhou_partition_sat(n, rad=radius)
    if not isinstance(triples, list):
        raise RuntimeError(f"no radius-{radius} corrected-lift partition for m={m}")

    arrays = []
    balance_certificates = []
    for residue in range(4):
        group = [triples[j] for j in range(n) if j % 4 == residue]
        certificate = find_modes(group, max_switches)
        if certificate is None:
            raise RuntimeError(
                f"group {residue} has no certificate with <= {max_switches} switches"
            )
        modes, mu_subset, w_subset = certificate
        array = [
            oriented_row(
                triple,
                mode,
                index in mu_subset,
                index in w_subset,
            )
            for index, (triple, mode) in enumerate(zip(group, modes))
        ]
        arrays.append(array)
        balance_certificates.append(
            {
                "residue": residue,
                "modes": modes,
                "mu_half_subset": sorted(mu_subset),
                "w_half_subset": sorted(w_subset),
                "switch_count": sum(mode != 0 for mode in modes),
            }
        )

    return {
        "type": "IHS(m,3;c)",
        "construction": "radius-corrected Zhou half-column lift, grouped by j mod 4",
        "m": m,
        "c": 4,
        "radius": radius,
        "arrays": arrays,
        "balance_certificates": balance_certificates,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("m", type=int)
    parser.add_argument("--radius", type=int, default=2)
    parser.add_argument("--max-switches", type=int, default=3)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    data = build_certificate(args.m, args.radius, args.max_switches)
    output = args.output or Path(f"c4_m{args.m}_certificate.json")
    output.write_text(json.dumps(data, indent=2), encoding="utf-8")
    counts = [item["switch_count"] for item in data["balance_certificates"]]
    print(f"WROTE {output}: switch counts {counts}")


if __name__ == "__main__":
    main()
