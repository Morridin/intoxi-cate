"""
Here be all the miscellaneous utilities that fit nowhere else.
"""
import os
from functools import cache
from pathlib import Path
from typing import Generator

import pandas as pd
from Bio import SeqIO

from lib import config

__all__ = ["fasta_to_dataframe", "global_output", "get_cys_pattern", "get_threads"]


def _generate_fasta_records(fasta_path: Path) -> Generator[dict[str, str], None, None]:
    """
    Yields records from a fasta file.
    :param fasta_path: The path to the fasta file.
    :return: A generator yielding dictionaries consisting of ID and Sequence key.
    """
    for record in SeqIO.parse(fasta_path, 'fasta'):
        yield {'ID': record.id, 'Sequence': str(record.seq)}


def fasta_to_dataframe(fasta_path):
    """
    Transforms a fasta file to a pandas dataframe..
    :param fasta_path: The path to the fasta file.
    :return: A pandas dataframe containing the sequence data from the fasta file.
    """
    return pd.DataFrame(_generate_fasta_records(fasta_path))


def global_output(path: str | Path) -> Path:
    """
    Constructs the global output path.
    """
    if type(path) == str:
        path = path.strip()

    output_dir = Path(config.get("output_dir", "").strip())
    if not output_dir.exists():
        output_dir.mkdir(parents=True, exist_ok=True)

    return output_dir / path


def get_cys_pattern(seq):
    """
    ???
    :param seq:
    :return:
    """
    pattern = ""
    status = False
    if pd.isna(seq).empty and seq.count('C') >= 4:
        for char in seq:
            if char == "C":
                pattern += "C"
                status = True
            else:
                if status:
                    pattern += "-"
                    status = False
        if pattern[-1] == "-":
            pattern = pattern[0:-1]
    if pattern == "":
        pattern = None
    return pattern


@cache
def get_threads() -> int:
    """
    Retrieves the number of threads from config file, or if not available, tries to get the number by other means.
    :return: The number of threads available to the program.
    """

    threads = config.get("threads")
    if threads is None:
        threads = len(os.sched_getaffinity(0))

    return threads
