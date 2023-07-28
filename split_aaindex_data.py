import argparse

from pathlib import Path

import pickle

import re


MULTIPLE_SPACES = re.compile(r"\s+")


def process_args() -> argparse.Namespace:
    processor = argparse.ArgumentParser()
    processor.add_argument(
        "aaindex_raw_filepath",
        type = Path,
        help = "Path to the aaindex1.txt file",
    )
    return processor.parse_args()

def main(args: argparse.Namespace) -> int:
    filepath = Path("./data/aaindex1.txt")
    data = {}
    last_id = ""
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("H"):
                last_id = line.split()[1]
            if line.startswith("I"):
                data_line_1 = next(f)
                data_line_2 = next(f)
                header_line = list(line[1 :].strip().replace(" ", "").replace("/", ""))
                data_line_1 = MULTIPLE_SPACES.sub(r" ", data_line_1.strip()).split()
                data_line_2 = MULTIPLE_SPACES.sub(r" ", data_line_2.strip()).split()
                data_line = data_line_1 + data_line_2
                data[last_id] = dict(zip(header_line, data_line))
                assert len(header_line) == len(data_line)
    with open("aaindex_bin", "wb") as f:
        pickle.dump(data, f)

    return 0


if __name__ == "__main__":
    args = process_args()
    raise SystemExit(main(args))