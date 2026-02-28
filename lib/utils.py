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


def _generate_fasta_records(fasta_path: Path, as_sequence: bool) -> Generator[dict[str, str], None, None] | Generator[
    dict[str, SeqIO.SeqRecord], None, None]:
    """
    Yields records from a fasta file.
    :param fasta_path: The path to the fasta file.
    :param as_sequence: Determines if the BioPython Sequence objects are kept as such or transformed to str.
    :return: A generator yielding dictionaries consisting of ID and Sequence key.
    """
    for record in SeqIO.parse(fasta_path, 'fasta'):
        if as_sequence:
            seq = record
        else:
            seq = str(record.seq)
        yield {'ID': record.id, 'Sequence': seq}


def fasta_to_dataframe(fasta_path, as_sequence=False):
    """
    Transforms a fasta file to a pandas dataframe..
    :param fasta_path: The path to the fasta file.
    :param as_sequence: Keep the BioPython sequences as objects if true, else convert to str.
    :return: A pandas dataframe containing the sequence data from the fasta file.
    """
    return pd.DataFrame(_generate_fasta_records(fasta_path, as_sequence))


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


def get_cys_pattern(seq: str) -> str | None:
    """
    Retrieves cysteine patterns within an aa sequence.
    :param seq: A string encoding a sequence of amino acids.
    :return: The cysteine pattern, if any or None else.
    """
    pattern = ""
    status = False
    if isinstance(seq, str) and seq.strip() and seq.count('C') >= 4:
        for char in seq:
            if char == "C":
                pattern += "C"
                status = True
            elif status:
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

def get_sequence_id(line: str) -> str:
    """
    Retrieves the Sequence ID from a FASTA/TMbed prediction file line.
    :param line: A line starting with a caret (>) and the sequence ID following.
    :return: The sequence ID with nothing else included.
    """
    if not line.startswith(">"):
        raise ValueError("Malformatted prediction file, expected '>' at the beginning of ID line.")
    seq_id, _, _ = line.lstrip(">").partition(" ")
    seq_id = seq_id.strip()
    if not seq_id:
        raise ValueError("Malformatted prediction file, expected '>' at the end of ID line.")
    return seq_id