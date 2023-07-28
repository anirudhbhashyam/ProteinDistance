import argparse

import collections

import itertools

import numpy as np

import os

import pandas as pd

import pickle

from pathlib import Path

import typing


CPD = Path(__file__).parents[0]
DATA_PATH = CPD / "data"
PathType = str | os.PathLike[str]


def process_args() -> argparse.Namespace:
    processor = argparse.ArgumentParser()
    processor.add_argument(
        "fasta_filepath",
        type = str,
        help = "The distance correlation between the two protein sequences.",
    )
    return processor.parse_args()


def get_aa_data(
    protein_sequence: str, 
    features: typing.Iterable[str], 
    index_data: dict[str, dict[str, float | int]],
) -> dict[str, list[float | int]]:
    aa_data = collections.defaultdict(list)
    for feat in features:
        for aa in protein_sequence:
            aa_data[feat].append(index_data[feat][aa])
    return aa_data


def seq_data_to_mat(aa_data: dict[str, list[float | int]]) -> np.ndarray:
    return np.array(list(aa_data.values()), dtype = np.float32)


def get_distance_correlation(seq_data_vec_1: np.ndarray, seq_data_vec_2: np.ndarray) -> np.ndarray:
    assert seq_data_vec_1.shape == seq_data_vec_2.shape, "Cannot compute distance matrix for sequences of different lengths."
    centered_vec_1 = seq_data_vec_1 - seq_data_vec_1.mean()
    centered_vec_2 = seq_data_vec_2 - seq_data_vec_2.mean()
    dxy = 1 - (centered_vec_1 @ centered_vec_2) / (np.linalg.norm(centered_vec_1) * np.linalg.norm(centered_vec_2))
    return dxy


def run(seq_1: str, seq_2: str, features: typing.Sequence[str], index_data: dict[str, dict[str, float | int]]) -> float:
    seq_data_1 = get_aa_data(seq_1, features, index_data)
    seq_data_2 = get_aa_data(seq_2, features, index_data)
    assert len(features) == len(seq_data_1.keys())
    assert len(features) == len(seq_data_2.keys())
    for v in seq_data_1.values():
        assert len(v) == len(seq_1)
    seq_data_mat_1 = seq_data_to_mat(seq_data_1)
    seq_data_mat_2 = seq_data_to_mat(seq_data_2)
    seq_data_mat_1_fft = np.apply_along_axis(np.fft.fft, 1, seq_data_mat_1)
    seq_data_mat_2_fft = np.apply_along_axis(np.fft.fft, 1, seq_data_mat_2)
    assert np.allclose(seq_data_mat_1_fft[0], np.fft.fft(seq_data_mat_1[0]))
    assert seq_data_mat_1_fft.shape[0] == 24
    assert seq_data_mat_1_fft.shape[1] == len(seq_1)
    assert np.allclose(seq_data_mat_2_fft[0], np.fft.fft(seq_data_mat_2[0]))
    assert seq_data_mat_2_fft.shape[0] == 24
    assert seq_data_mat_2_fft.shape[1] == len(seq_2)
    seq_data_vec_1 = np.absolute(seq_data_mat_1_fft.flatten()) ** 2
    seq_data_vec_2 = np.absolute(seq_data_mat_2_fft.flatten()) ** 2
    dxy = get_distance_correlation(seq_data_vec_1, seq_data_vec_2)
    return float(dxy)


def read_fasta(filepath: PathType) -> typing.Iterator[str]:
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith(">"):
                continue
            yield line.strip()


def main(args: argparse.Namespace) -> int:
    aa_features = pd.read_csv(DATA_PATH / "aa_features.csv")
    if not Path("aaindex_bin").exists():
        raise FileNotFoundError("aaindex_bin not found. Please run `python3 split_aaindex_data.py` to generate it.")
    with open("aaindex_bin", "rb") as f:
        aa_data = pickle.load(f)
    long_range_contact_data = pd.read_csv(DATA_PATH / "long_range_contacts.csv")
    long_range_contact_data = dict(zip(long_range_contact_data.columns, long_range_contact_data.iloc[0]))
    relative_connectivity_data = pd.read_csv(DATA_PATH / "relative_connectivity.csv")
    relative_connectivity_data = dict(zip(relative_connectivity_data.columns, relative_connectivity_data.iloc[0]))
    aa_data["Nl"] = long_range_contact_data
    aa_data["Rk"] = relative_connectivity_data
    # WARNING: Proscale_4 not present.
    features = aa_features["ID"].tolist()
    seq_1, seq_2 = itertools.islice(read_fasta(args.fasta_filepath), 2)
    dxy = run(seq_1, seq_2, features, aa_data)
    print(dxy)
    return 0


if __name__ == "__main__":
    args = process_args()
    raise SystemExit(main(args))
