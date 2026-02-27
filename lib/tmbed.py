"""
This module bundles the functionality that previously culminated in several runs of SignalP: Detecting amino acid
sequences that contain signal peptides and are thus potentially toxins.
"""
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

from lib import config, utils

__all__ = ["tmbed"]


# ============================= Public functions ============================= #
def tmbed(clustered_peptides: Path, use_gpu: bool = True, cpu_fallback: bool = False) -> pd.DataFrame:
    """
    This function runs tmbed, collects the results and returns them as versatile DataFrame.
    :param clustered_peptides: The path to a FASTA file containing the amino acid sequences on which TMbed shall run.
    :param use_gpu: If you don't want to utilise the GPU to speed up the prediction process, set to False.
    :param cpu_fallback: If you want to fallback on the CPU instead of having the pipeline fail, set to True.
    :return: A DataFrame containing the following information: an ID column holding the sequence ID,
    the sequence itself in a second column, and the mature peptide (the peptide without the signal peptide) in a
    third column. The DataFrame only contains such sequences that TMbed could identify as featuring a signal peptide.
    """
    model_dir = config.get_path("tmbed_model_path")
    threads = utils.get_threads()

    with tempfile.NamedTemporaryFile(suffix=".tmbed", delete_on_close=False) as output_file:
        _run_tmbed(clustered_peptides, output_file.name, use_gpu, cpu_fallback, threads, model_dir)

        return utils.parse_tmbed_predictions(output_file.name, "signal")


# ============================ Private functions ============================= #
def _run_tmbed(clustered_peptides: Path, output_file_name: str, use_gpu: bool, cpu_fallback: bool, threads: int,
               model_dir: Path | None) -> None:
    """
    The actual runner for TMbed. Assembles the command from the given parameters and runs TMbed.
    :param clustered_peptides: A Path pointing to a FASTA file.
    :param output_file_name: A string representing the path of the output file, including the file itself.
    :param use_gpu: Whether TMbed shall utilise the GPU or not.
    :param cpu_fallback: Whether TMbed shall fall back on CPU if the GPU is not available.
    :param threads: The number of CPU threads to use if TMbed does not run on the GPU.
    :param model_dir: The directory where TMbed shall store its models for later reuse.
    """
    command = ["tmbed", "predict", "-f", clustered_peptides, "-p", output_file_name, "-t", f"{threads}"]
    if use_gpu:
        command.append("--use-gpu")
    else:
        command.append("--no-use-gpu")

    if cpu_fallback:
        command.append("--cpu-fallback")
    else:
        command.append("--no-cpu-fallback")

    if model_dir is not None:
        command += ["-m", model_dir]

    subprocess.run(
        command,
    )
