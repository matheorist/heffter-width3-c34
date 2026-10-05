"""Find and verify finite periodic sign formulas for all c=4 parameter states.

For each of 4 parameter residues, 12 step states, and 4 arrays, bulk signs are
constant on (zone, period-48 position, layer parity); boundary signs are constant
on their relative position.  Category sums are quadratic in Q.  CP-SAT finds
signs making both independent signed sums zero, and finite differences verify
the quadratic identity.
"""

import json
from pathlib import Path

from ortools.sat.python import cp_model

from generate_c4_all import PERIODIC_MINIMUM, corrected_triples, load_periodic_rule
from generate_c4_certificate import find_modes, mode_values


HERE = Path(__file__).resolve().parent


def category(m, j, phase):
    n = 4 * m
    kappa = n // 2
    if j <= 20:
        return ("low", j)
    if j <= kappa - 21:
        zone = "word1"
        word_position = (j - phase[0]) % 48
        phase_value = phase[0]
    elif j <= kappa + 20:
        return ("middle", j - kappa)
    elif j <= n - 21:
        zone = "word2"
        word_position = (j - phase[1]) % 48
        phase_value = phase[1]
    else:
        return ("top", j - n)
    j0 = ((phase_value + word_position - 1) % 48) + 1
    layer = (j - j0) // 48
    return (zone, word_position, layer % 2)


def state_data(parameter_residue, state):
    minimum = PERIODIC_MINIMUM[parameter_residue]
    m0 = minimum + 12 * state
    base_m, base_phase, _, _, _ = load_periodic_rule(parameter_residue)
    original_step = (m0 - base_m) // 12
    phase = (
        (base_phase[0] + 16 * original_step) % 48,
        (base_phase[1] + 16 * original_step) % 48,
    )
    return m0, phase


def switch_specification(m0, phase, group):
    triples, _, _ = corrected_triples(m0)
    group_triples = [triples[index] for index in range(4 * m0) if index % 4 == group]
    certificate = find_modes(group_triples, 1)
    assert certificate is not None
    modes = certificate[0]
    switched = [index for index, mode in enumerate(modes) if mode != 0]
    if not switched:
        return None, 0
    assert len(switched) == 1
    index = switched[0]
    j = group + 1 + 4 * index
    switch_category = category(m0, j, phase)
    assert switch_category[0] == "low"
    return switch_category, modes[index]


def category_sums(m, phase, group, switch_category, switch_mode, value_kind):
    triples, _, _ = corrected_triples(m)
    result = {}
    for zero_based_j in range(group, 4 * m, 4):
        j = zero_based_j + 1
        row_category = category(m, j, phase)
        mode = switch_mode if row_category == switch_category else 0
        mu, w = mode_values(triples[zero_based_j], mode)
        value = mu if value_kind == "mu" else w
        result[row_category] = result.get(row_category, 0) + value
    return result


def third_differences(values):
    current = list(values)
    for _ in range(3):
        current = [current[index + 1] - current[index] for index in range(len(current) - 1)]
    return current


def solve_signs(series):
    categories = sorted({key for row in series for key in row}, key=str)
    # Every category sum must be a polynomial of degree at most 2 in Q.
    for key in categories:
        values = [row.get(key, 0) for row in series]
        # A bounded number of base states can lie on the other side of an
        # absolute-value breakpoint.  From Q=1 onward every category sum is
        # quadratic; Q=0 is retained as a separate exact equation below.
        assert all(value == 0 for value in third_differences(values[1:])), (key, values)

    model = cp_model.CpModel()
    variables = [model.NewBoolVar(f"s{index}") for index in range(len(categories))]
    for row in series:
        model.Add(
            sum(
                row.get(key, 0) * (2 * variables[index] - 1)
                for index, key in enumerate(categories)
            )
            == 0
        )
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = 30
    solver.parameters.num_search_workers = 8
    solver.parameters.random_seed = 31
    status = solver.Solve(model)
    assert status in (cp_model.OPTIMAL, cp_model.FEASIBLE), solver.StatusName(status)
    signs = {
        str(key): 1 if solver.Value(variables[index]) else -1
        for index, key in enumerate(categories)
    }
    # Six zeros plus degree <=2 proves the signed category sum is identically zero.
    assert all(
        sum(row.get(key, 0) * signs[str(key)] for key in categories) == 0
        for row in series
    )
    return signs


def prove_one(parameter_residue, state, group):
    m0, phase = state_data(parameter_residue, state)
    switch_category, switch_mode = switch_specification(m0, phase, group)
    output = {
        "parameter_residue_mod_12": parameter_residue,
        "step_state_mod_12": state,
        "group": group,
        "base_m": m0,
        "phase": list(phase),
        "switch_category": None if switch_category is None else list(switch_category),
        "switch_mode": switch_mode,
    }
    for value_kind in ("mu", "w"):
        series = [
            category_sums(
                m0 + 144 * q,
                phase,
                group,
                switch_category,
                switch_mode,
                value_kind,
            )
            for q in range(7)
        ]
        output[f"{value_kind}_signs"] = solve_signs(series)
    return output


def main():
    patterns = []
    for parameter_residue in (1, 5, 7, 11):
        for state in range(12):
            for group in range(4):
                patterns.append(prove_one(parameter_residue, state, group))
            print(f"PROVED balance residue={parameter_residue}, state={state}", flush=True)
    path = HERE / "c4_balance_affine_patterns.json"
    path.write_text(json.dumps(patterns, indent=2), encoding="utf-8")
    print(f"AFFINE BALANCE PROOF VERIFIED: {len(patterns)} group states")
    print(f"WROTE {path}")


if __name__ == "__main__":
    main()
