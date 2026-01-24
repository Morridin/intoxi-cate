"""
Here be all the miscellaneous utilities that fit nowhere else.
"""
from os.path import join
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from lib import snakemake_checkpoints, config


def _get_signalp_splits(wildcards):
    checkpoint_output = snakemake_checkpoints.split_fasta.get(**wildcards).output[0]
    return expand(global_output("") / 'split_files/{i}.fasta',
           i=glob_wildcards(join(checkpoint_output, '{i}.fasta')).i)

def fasta_to_dataframe(fasta_ath):
    sequences = SeqIO.parse(open(fasta_ath), 'fasta')
    data = []
    for record in sequences:
        data.append({'ID': record.id, 'Sequence': str(record.seq)})
    return pd.DataFrame(data)


def global_output(path: str | Path):
    """
    Constructs the global output path.
    """
    output_dir = config.get("output_dir", "")
    if output_dir:
        return Path(f"{output_dir}/{path}")
    else:
        return Path(path)

def expand(path: Path | str, **kwargs) -> Path | tuple[Path, ...]:
    """
    Temp function stub to reduce error markings
    TODO: Remove
    :param path: -
    :param kwargs: -
    :return:
    """
    _ = path, kwargs
    raise NotImplementedError()

def glob_wildcards(path: str) -> str:
    """
    Temp function stub to reduce error markings
    TODO: Remove
    :param path: -
    :return:
    """
    _ = path
    raise NotImplementedError()

def aggregate_splits(wildcards):
    """
    ?
    :param wildcards:
    :return:
    """
    checkpoint_output = snakemake_checkpoints.split_fasta.get(**wildcards).output[0]
    return expand(global_output("")+"split_files/{i}_summary.signalp5",
        i=glob_wildcards(join(checkpoint_output, "{i}.faa")).i)

def get_cys_pattern(seq):
    """
    ???
    :param seq:
    :return:
    """
    pattern = ""
    status = False
    if not pd.isna(seq) and seq.count('C') >= 4:
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
