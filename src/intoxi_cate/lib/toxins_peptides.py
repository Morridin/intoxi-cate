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
This module contains the functionality for DeToX's last filtering step where all the transmembrane peptides are removed
from the remaining list of toxin candidate peptides.
"""
import tempfile
from pathlib import Path
from typing import Generator

import pandas as pd
from Bio import SeqIO

from . import config, utils, tmbed_wrapper

__all__ = ["retrieve_candidate_toxins"]


# ============================= Public functions ============================= #
def retrieve_candidate_toxins(clustered_peptides: Path, toxins_blast_result: pd.DataFrame,
                              signal_peptides: pd.DataFrame) -> Path:
    """
    Builds a list of toxin candidates out of the data provided with the parameters.
    In this list, only those peptides are included that were either detected by structural features or sequence
    similarity to known toxins and that do not feature any transmembrane structures.
    :param clustered_peptides: A FASTA file containing amino acid sequences.
    :param toxins_blast_result: A DataFrame containing peptides that were detected by sequence similarity to other, known toxins.
        It is expected to contain a column 'ID' listing the sequence IDs of the contained peptides.
    :param signal_peptides: A DataFrame containing peptides with clear signal sequences.
        The only relevant property is that its index consists of the sequence IDs of the contained peptides.
    :return: A DataFrame holding those sequences that are candidate toxins.
    """
    secreted_peptides, non_secreted_peptides = _extract_secreted_peptides(signal_peptides, clustered_peptides)

    threads = utils.get_threads()
    use_cpu = config.get("tmbed_use_cpu", True)
    use_gpu = config.get("tmbed_use_gpu", True)
    model_dir = config.get_path("tmbed_model_path")

    with tempfile.NamedTemporaryFile(suffix=".tmbed", delete_on_close=False) as output_file:
        tmbed_wrapper.run_tmbed(secreted_peptides, output_file.name, use_gpu, use_cpu, threads, model_dir)
        tmbed_result = tmbed_wrapper.parse_tmbed_predictions(output_file.name, _generate_non_transmembrane_rows).set_index("ID")

    output_file = utils.global_output(config.get("basename") + "_candidate_toxins.fasta")

    with open(output_file, "w") as out_file:
        for seq in pd.concat((_filter_fasta_file(secreted_peptides, tmbed_result),
                              _filter_fasta_file(non_secreted_peptides, toxins_blast_result))).iterrows():
            SeqIO.write(seq[1].iloc[0], out_file, "fasta")

    return output_file


# ============================ Private functions ============================= #
def _extract_secreted_peptides(signal_peptides: pd.DataFrame, clustered_peptides: Path) -> tuple[Path, Path]:
    """
    This function filters a FASTA file with amino acid sequences against a list of secreted peptides and
    sorts them into 'secreted' and 'non-secreted' peptides.
    :param signal_peptides: A DataFrame containing only definitively secreted peptides. Its index is expected to contain the relevant Sequence IDs.
    :param clustered_peptides: The path to a FASTA file containing amino acid sequences (e.g. the output of `cluster_peptides`).
    :return: A tuple of two file paths with the first one pointing to the secreted peptides and the second one pointing to the non-secreted peptides.
    """
    secreted_peptides = utils.global_output(config.get("basename") + "_secreted_peptides.fasta")
    non_secreted_peptides = utils.global_output(config.get("basename") + "_non_secreted_peptides.fasta")

    with open(non_secreted_peptides, "w") as n_outfile, open(secreted_peptides, "w") as out_file:
        for seq in SeqIO.parse(clustered_peptides, "fasta"):
            if seq.id in signal_peptides.index:
                SeqIO.write(seq, out_file, "fasta")
            else:
                SeqIO.write(seq, n_outfile, "fasta")

    return secreted_peptides, non_secreted_peptides


def _filter_fasta_file(fasta_file: Path, filter_map: pd.DataFrame) -> pd.DataFrame:
    """
    Filters a FASTA file for the peptides present in the file given by the second parameter.
    All peptides in the first file, that are also present in the second file,
    will be copied into a pandas DataFrame and returned.

    :param fasta_file: The path to the FASTA file that shall be filtered.
    :param filter_map: A DataFrame containing the sequence IDs of the peptides to keep in the column 'ID'.
    :return: A DataFrame containing the peptides present in both files, consisting of the columns 'ID' and 'Sequence',
        with the latter holding original the Bio.Seq object from parsing the input FASTA file.
    """

    records = utils.fasta_to_dataframe(fasta_file, True).set_index("ID")
    return records.drop(records.join(filter_map, how="left_anti", rsuffix="_filter").index)


def _generate_non_transmembrane_rows(file: Path | str) -> Generator[dict[str, str], None, None]:
    with open(file) as f:
        seq_id = ""
        markers = set("bBhH")

        for index, line in enumerate(f):
            index %= 3
            if not line:
                break
            if index == 0:
                seq_id = utils.get_sequence_id(line)
            if index == 2:
                if markers.isdisjoint(line):
                    yield {"ID": seq_id}
