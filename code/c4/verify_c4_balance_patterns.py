"""Independently verify saved c=4 affine balance patterns without a solver."""

import ast
import json
from pathlib import Path

from prove_c4_balance_affine import category_sums, state_data, third_differences


HERE = Path(__file__).resolve().parent


def main():
    patterns = json.loads((HERE / "c4_balance_affine_patterns.json").read_text())
    assert len(patterns) == 192
    for pattern in patterns:
        parameter_residue = pattern["parameter_residue_mod_12"]
        state = pattern["step_state_mod_12"]
        group = pattern["group"]
        m0, phase = state_data(parameter_residue, state)
        assert m0 == pattern["base_m"]
        assert list(phase) == pattern["phase"]
        switch_category = pattern["switch_category"]
        if switch_category is not None:
            switch_category = tuple(switch_category)
            assert switch_category[0] == "low"
        switch_mode = pattern["switch_mode"]

        for value_kind in ("mu", "w"):
            signs = pattern[f"{value_kind}_signs"]
            assert set(signs.values()) <= {-1, 1}
            decoded_signs = {ast.literal_eval(key): value for key, value in signs.items()}
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
            categories = {key for row in series for key in row}
            assert categories == set(decoded_signs)
            for key in categories:
                values = [row.get(key, 0) for row in series]
                assert all(value == 0 for value in third_differences(values[1:]))
            assert all(
                sum(value * decoded_signs[key] for key, value in row.items()) == 0
                for row in series
            )
    print("VERIFIED 192 affine balance patterns without CP-SAT")


if __name__ == "__main__":
    main()
