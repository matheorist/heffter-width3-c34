"""Generate and independently verify all admissible c=4 cases in a range."""

import argparse

from generate_c4_all import build


def verify(data):
    m = data["m"]
    arrays = data["arrays"]
    assert len(arrays) == 4
    assert all(len(array) == m for array in arrays)
    assert all(sum(row) == 0 for array in arrays for row in array)
    assert all(
        [sum(row[column] for row in array) for column in range(3)] == [0, 0, 0]
        for array in arrays
    )
    support = [abs(entry) for array in arrays for row in array for entry in row]
    assert sorted(support) == list(range(1, 12 * m + 1))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--maximum", type=int, default=1001)
    args = parser.parse_args()
    count = 0
    for m in range(5, args.maximum + 1, 2):
        if m % 3 == 0:
            continue
        verify(build(m))
        count += 1
        if count % 50 == 0:
            print(f"verified {count} cases; last m={m}", flush=True)
    print(f"VERIFIED all {count} admissible cases through m={args.maximum}")


if __name__ == "__main__":
    main()
