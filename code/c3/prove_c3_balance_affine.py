"""Verify the stored c=3 sign patterns by quadratic finite differences."""

import ast

from generate_c3_all import category_key, load_json, pattern_for, support_rows


def admissible_state(cls, state):
    """Return whether n=cls+24*state has gcd(n/3,3)=1."""
    return state % 3 != (0 if cls == 9 else 1)


def mode_values(triple, mode):
    a, b, c = sorted(abs(value) for value in triple)
    assert a + b == c
    if mode == "c":
        return c, b - a
    if mode == "a":
        return a, b + c
    if mode == "b":
        return b, a + c
    raise ValueError(mode)


def third_differences(values):
    current = list(values)
    for _ in range(3):
        current = [current[index + 1] - current[index] for index in range(len(current) - 1)]
    return current


def category_sums(n, group, value_kind, patterns):
    m = n // 3
    triples = support_rows(n)
    pattern = pattern_for(n, group, patterns)
    switch_key, switch_mode = pattern["switch"]
    switch_key = None if switch_key == "None" else ast.literal_eval(switch_key)
    result = {}
    for j in range(1, n + 1):
        if j % 3 != group:
            continue
        key = category_key(j, m, n)
        mode = switch_mode if key == switch_key else "c"
        mu, w = mode_values(triples[j - 1], mode)
        value = mu if value_kind == "mu" else w
        result[key] = result.get(key, 0) + value
    return result


def prove_one(cls, state, group, patterns):
    n0 = 24 * state + cls
    pattern = pattern_for(n0, group, patterns)
    for value_kind, sign_name in (("mu", "c"), ("w", "d")):
        series = [
            category_sums(n0 + 288 * q, group, value_kind, patterns)
            for q in range(7)
        ]
        categories = sorted({key for row in series for key in row}, key=str)
        signs = {ast.literal_eval(key): value for key, value in pattern[sign_name].items()}
        assert set(categories) == set(signs), (
            cls, state, group, value_kind, set(categories) - set(signs), set(signs) - set(categories)
        )
        for key in categories:
            values = [row.get(key, 0) for row in series]
            assert all(value == 0 for value in third_differences(values)), (
                cls, state, group, value_kind, key, values
            )
        totals = [sum(row[key] * signs[key] for key in categories) for row in series]
        assert totals == [0] * len(series), (cls, state, group, value_kind, totals)


def main():
    patterns = load_json("lemma2_patterns.json")
    count = 0
    for cls in (9, 21):
        for state in range(15, 27):
            if not admissible_state(cls, state):
                continue
            for group in range(3):
                prove_one(cls, state, group, patterns)
                count += 1
            print(f"PROVED balance class={cls}, t_state={state}")
    print(f"AFFINE C3 BALANCE PROOF VERIFIED: {count} group states")


if __name__ == "__main__":
    main()
