"""
This module bundles all the former rules around SignalP:

trim_peptides <- split_fasta (checkpoint) <- run_signalp <- filter_signalp_outputs

The public API consists of filter_signalp_outputs.
"""
import itertools
import subprocess
from functools import cache
from pathlib import Path
from typing import Iterable

import pandas as pd
from Bio import SeqIO

from lib import utils, config

__all__ = ["signalp"]


# ============================= Public functions ============================= #
@cache
def signalp(clustered_peptides: Path) -> pd.DataFrame:
    """
    Prepares the data for SignalP, runs SignalP and returns a filtered result.
    :return: A DataFrame containing the SignalP result filtered for such peptides that can pass through cell membranes.
    """
    chunk_size = 5000
    trimmed_peptides_dir = _trim_peptides(clustered_peptides, 50, chunk_size)

    output_dir = trimmed_peptides_dir.with_name("signalp_outputs")
    output_dir.mkdir(parents=True, exist_ok=True)

    signalp_outputs = [
        _run_signalp(file, output_dir / file.stem, chunk_size, config.get_path("signalp_path")) for file in
        trimmed_peptides_dir.iterdir()
    ]

    signalp_threshold = config.get("signalp_dvalue")
    if signalp_threshold is None:
        signalp_threshold = 0.7

    return _filter_signalp_outputs(signalp_outputs, signalp_threshold)


# ============================ Private functions ============================= #
def _trim_peptides(aa_sequences: Path, cut_off: int, chunk_size: int) -> Path:
    """
    Trims all peptides to only the first 50 aminoacids, as they are the only useful part for SignalP.
    Also, this function splits its output into multiple smaller file with up to `chunk_size` sequences each.
    :param aa_sequences: The path to a FASTA file holding clustered nucleotide sequences (see cluster_peptides).
        Given the expected contents, the file should wear the FAA suffix.
    :param cut_off: The maximum length of the trimmed sequences.
    :param chunk_size: Determines after how many sequences go into one output file.
    :return: The path to a directory containing FASTA files holding each up to `chunk_size` sequences, cut off after the 50th amino acid.
        The ending will be FAA corresponding to the file contents.
    """
    output_dir = utils.global_output("split_peptides_files")
    output_dir.mkdir(parents=True, exist_ok=True)

    trimmed_sequences = (record if len(record.seq) < cut_off else record[:cut_off] for record in
                         SeqIO.parse(aa_sequences, "fasta"))

    for index, batch in enumerate(itertools.batched(trimmed_sequences, chunk_size)):
        SeqIO.write(batch, output_dir / f"{index}.faa", "fasta")

    return output_dir


def _run_signalp(input_sequences: Path, output_prefix: Path, chunk_size: int, signalp_path: Path | None) -> Path:
    """
    Runs SignalP on the file `fasta_file` points to.
    Each file is expected to contain at most `chunk_size` sequences of eukaryotic origin.
    :param input_sequences: The path to the FASTA file containing the sequences on which SignalP shall be run.
    :param output_prefix: A prefix to be prepended to the output file name. This may also include directories.
    :param chunk_size: Describes how many sequences are expected in the input file and shall be processed per batch.
    :param signalp_path: The path to the SignalP binary, if SignalP cannot be called directly.
    :return: A TSV file containing the predicted secretory status for each input sequence.
    """
    if signalp_path is None or not (signalp_path / "signalp").exists():
        prefix = ""
    else:
        prefix = f"./{signalp_path}"

    subprocess.run(
        f"./signalp -batch {chunk_size} -fasta {"../" * 3}{input_sequences} -org euk -format short -verbose -prefix {"../" * 3}{output_prefix}",
        shell=True,
        cwd=prefix
    )

    return Path(f"{output_prefix}_summary.signalp5")


def _filter_signalp_outputs(files: Iterable[Path], threshold: float) -> pd.DataFrame:
    """
    Filters the files given by `files` for those peptides whose probability of signal peptides is greater than `threshold`.

    :return: A DataFrame containing only those lines from the input files that passed the filter.
    """
    data = [pd.read_csv(file, sep="\t", index_col=0, header=None, comment="#",
                        names=["ID", "signalp_prediction", "prob_signal", "prob_nosignal", "cutsite"]) for file in
            files]

    data = pd.concat(data)
    data = data[
        data[2] > threshold &
        ~data.astype(str).apply(lambda r: r.str.contains("?", regex=False))
        ]

    return data
