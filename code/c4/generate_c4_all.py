"""Generate a candidate certificate for every admissible IHS(m,3;4).

Six small values use separately verified finite certificates.  All remaining
values use one of four recovered period-48 corrected-lift rules, according to
m modulo 12, followed by an exact <=1-switch balance certificate search.
"""

import argparse
import json
from fractions import Fraction
from pathlib import Path

from generate_c4_certificate import find_modes, oriented_row


HERE = Path(__file__).resolve().parent
SMALL = {5, 7, 11, 13, 17, 19}
PERIODIC_MINIMUM = {1: 25, 5: 29, 7: 31, 11: 23}


def integer_keys(data):
    return {int(key): tuple(value) for key, value in data.items()}


def load_periodic_rule(residue):
    if residue == 1:
        even = json.loads((HERE / "c4_period48_class4_rule.json").read_text())
        odd = json.loads(
            (HERE / "c4_period48_class4_boundary_odd.json").read_text()
        )
        base_m = 37
        base_phase = (24, 36)
    else:
        base_n = {5: 356, 7: 364, 11: 380}[residue]
        even = json.loads((HERE / f"c4_period48_base{base_n}_rule.json").read_text())
        odd = json.loads(
            (HERE / f"c4_period48_base{base_n}_boundary_odd.json").read_text()
        )
        base_m = base_n // 4
        base_phase = (0, 0)

    words = [integer_keys(even["word1"]), integer_keys(even["word2"])]
    even_boundary = [
        integer_keys(even[name]) for name in ("low", "middle", "top")
    ]
    odd_boundary = [
        integer_keys(odd[name]) for name in ("low", "middle", "top")
    ]
    return base_m, base_phase, words, even_boundary, odd_boundary


def corrected_triples(m):
    residue = m % 12
    if residue not in PERIODIC_MINIMUM or m < PERIODIC_MINIMUM[residue]:
        raise ValueError("m is below the periodic range")
    base_m, base_phase, words, even_boundary, odd_boundary = load_periodic_rule(
        residue
    )
    step = (m - base_m) // 12
    phase = (
        (base_phase[0] + 16 * step) % 48,
        (base_phase[1] + 16 * step) % 48,
    )
    boundary = even_boundary if step % 2 == 0 else odd_boundary
    n = 4 * m
    kappa = n // 2
    modulus = 2 * n + 1

    def zhou(i, j):
        parity = -1 if (j + n) % 2 else 1
        return int(
            Fraction(1 - parity * modulus, 4) * i
            + Fraction(2 - 3 * abs(i), 2) * j
        )

    triples = []
    for j in range(1, n + 1):
        if j <= 20:
            delta = boundary[0][j]
        elif j <= kappa - 21:
            delta = words[0][(j - phase[0]) % 48]
        elif j <= kappa + 20:
            delta = boundary[1][j - kappa]
        elif j <= n - 21:
            delta = words[1][(j - phase[1]) % 48]
        else:
            delta = boundary[2][j - n]
        values = (
            abs(3 * zhou(-1, j) + delta[0]),
            abs(3 * j + delta[1]),
            abs(3 * zhou(1, j) + delta[2]),
        )
        triples.append(tuple(sorted(values)))

    assert sorted(value for triple in triples for value in triple) == list(
        range(1, 3 * n + 1)
    )
    return triples, step, phase


def build_periodic(m):
    triples, step, phase = corrected_triples(m)
    n = 4 * m
    arrays = []
    certificates = []
    for residue in range(4):
        group = [triples[j] for j in range(n) if j % 4 == residue]
        certificate = find_modes(group, 1)
        if certificate is None:
            raise RuntimeError(f"no <=1-switch certificate for group {residue}")
        modes, mu_subset, w_subset = certificate
        arrays.append(
            [
                oriented_row(
                    triple,
                    mode,
                    index in mu_subset,
                    index in w_subset,
                )
                for index, (triple, mode) in enumerate(zip(group, modes))
            ]
        )
        certificates.append(
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
        "construction": f"period-48 rule for m={m % 12} mod 12",
        "m": m,
        "c": 4,
        "step": step,
        "phase": list(phase),
        "arrays": arrays,
        "balance_certificates": certificates,
    }


def build(m):
    if m < 5 or m % 2 == 0 or m % 3 == 0:
        raise ValueError("admissible m must be odd, >=5, and coprime to 3")
    if m in SMALL:
        path = HERE / f"c4_m{m}_certificate.json"
        data = json.loads(path.read_text())
        data["construction"] = "finite base certificate"
        return data
    return build_periodic(m)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("m", type=int)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = build(args.m)
    output = args.output or Path(f"c4_m{args.m}_all_certificate.json")
    output.write_text(json.dumps(data, indent=2), encoding="utf-8")
    switches = [
        item["switch_count"] for item in data.get("balance_certificates", [])
    ]
    print(f"WROTE {output}: construction={data['construction']}, switches={switches}")


if __name__ == "__main__":
    main()
