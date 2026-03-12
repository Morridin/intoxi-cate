# Copyright 2026 Paul Zanner
# Some functions are adopted from DeToX, Ringeval et al., 2024, https://doi.org/10.1093/bib/bbae094
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
Here be all the miscellaneous utilities that fit nowhere else.
"""
import os
import subprocess
from functools import cache
from pathlib import Path
from typing import Generator

import pandas as pd
from Bio import SeqIO

from . import config

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

    output_dir = config.get_path("output_dir")
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

@cache
def ensure_mmseqs2(default: Path = None) -> str:
    """
    Checks if either the MMSeqs path given in config or as default leads to an existing and working MMSeqs installation
    or MMSeqs2 is available on PATH.
    :param default: An alternative path to look for an MMSeqs installation.
    :return: A string that can be used to call MMSeqs2. If MMseqs2 could not be detected, a FileNotFoundError is raised.
    """
    def _check_path(mmseqs: Path | None) -> tuple[bool, Path]:
        if not mmseqs:
            return False, Path("")

        if mmseqs.name != "mmseqs":
            mmseqs /= "mmseqs"

        try:
            subprocess.run([mmseqs, "version"], capture_output=True)
        except FileNotFoundError:
            return False, Path("")

        return True, mmseqs

    # First, check if something's within config.
    found, path = _check_path(config.get_path("mmseqs_path"))
    if found:
        return str(path.resolve())

    # Then, check the default value
    found, path = _check_path(default)
    if found:
        return str(path.resolve())

    # Finally, check if MMSeqs is on Path
    try:
        subprocess.run(["mmseqs", "version"], capture_output=True)
    except FileNotFoundError as e:
        raise e
    return "mmseqs"