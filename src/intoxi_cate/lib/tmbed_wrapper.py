# Copyright 2026 Paul Zanner
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
This module bundles the functionality that previously culminated in several runs of SignalP: Detecting amino acid
sequences that contain signal peptides and are thus potentially toxins.
"""
import subprocess
import tempfile
from pathlib import Path
from typing import Generator, Callable

import pandas as pd

from . import config, utils

__all__ = ["detect_by_structure"]


# ============================= Public functions ============================= #
def detect_by_structure(clustered_peptides: Path) -> pd.DataFrame:
    """
    This function runs tmbed, collects the results and returns them as versatile DataFrame.
    :param clustered_peptides: The path to a FASTA file containing the amino acid sequences on which TMbed shall run.
    :return: A DataFrame containing the following information: an ID column holding the sequence ID,
    the sequence itself in a second column ("Sequence"), and the mature peptide (the peptide without the signal peptide) in a
    third column ("Mature Peptide"). Next, it contains a boolean column "Signal Peptide Predicted" that is `True` if TMbed
    predicted a signal peptide for the sequence and the predicted signal peptide is longer than the configured cut-off value.
    Also, the DataFrame contains a string column "Raw Prediction" containing TMbed's prediction string for the sequence,
    and finally, the numerical, 0-based position of the cleavage site along the sequence.
    The DataFrame contains all input sequences. You can filter for the sequences with signal peptides by selecting
    via the "Signal Peptide Predicted" column.
    """
    model_dir = config.get_path("tmbed_model_path")
    use_cpu = config.get("tmbed_use_cpu", True)
    use_gpu = config.get("tmbed_use_gpu", True)
    threads = utils.get_threads()

    with tempfile.NamedTemporaryFile(suffix=".tmbed", delete_on_close=False) as output_file:
        _run_tmbed(clustered_peptides, output_file.name, use_gpu, use_cpu, threads, model_dir)

        return pd.DataFrame(_generate_tmbed_pred_df_rows_signal_only(output_file.name))


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
    command = [
        "tmbed", "predict",
        "-f", clustered_peptides,
        "-p", output_file_name,
        "-t", f"{threads}",
        "--out-format", "4"
    ]

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

    subprocess.run(command, check=True)


def _generate_tmbed_pred_df_rows_signal_only(file: Path | str) -> Generator[dict[str, str], None, None]:
    """
    This Generator function yields DataFrame rows while iterating over a FASTA-like file that contains TMbed prediction output.
    :param file: The path to the file out of which to generate the DataFrame rows.
    :returns: A Generator yielding dictionaries with the keys "ID", "Sequence", "Signal Peptide Predicted", "Raw Prediction",
        "Mature Peptide" and "cutsite", corresponding to the DataFrame columns desscribed in `detect_by_structure`.
    """
    sp_threshold = config.get("signalpeptide_minlen")
    with open(file) as f:
        seq_id = ""
        sequence = ""

        for index, line in enumerate(f):
            index %= 3
            if not line:
                break
            if index == 0:
                if not line.startswith(">"):
                    raise ValueError("Malformatted prediction file, expected '>' at the beginning of ID line.")
                seq_id, _, _ = line.lstrip(">").partition(" ")
                seq_id = seq_id.strip()
                if not seq_id:
                    raise ValueError("Malformatted prediction file, expected '>' at the end of ID line.")
            if index == 1:
                sequence = line.strip()
            if index == 2:
                diff = len(line) - len(line.lstrip("S"))
                if diff >= sp_threshold:
                    yield {
                        "ID": seq_id,
                        "Sequence": sequence,
                        "Signal Peptide Predicted": True,
                        "Raw Prediction": line.strip(),
                        "Mature Peptide": sequence[diff:],
                        "cutsite": diff,
                    }
                else:
                    yield {
                        "ID": seq_id,
                        "Sequence": sequence,
                        "Signal Peptide Predicted": False,
                        "Raw Prediction": line.strip(),
                        "Mature Peptide": None,
                        "cutsite": None,
                    }
