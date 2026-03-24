# Copyright 2026 Paul Zanner
# Adopted and redesigned from DeToX pipeline, Ringeval et al., 2024, https://doi.org/10.1093/bib/bbae094
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
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from . import config, utils

__all__ = ["retrieve_candidate_toxins"]


# ============================= Public functions ============================= #
def retrieve_candidate_toxins(peptides: pd.DataFrame, toxins_blast_result: pd.DataFrame) -> Path:
    """
    Builds a list of toxin candidates out of the data provided with the parameters.
    In this list, only those peptides are included that were either detected by structural features or sequence
    similarity to known toxins and that do not feature any transmembrane structures.
    :param peptides: A DataFrame containing at least the following columns: "ID" - str - the sequence ID, "Sequence" -
        str - the amino acid sequence itself, "Signal Peptide Predicted" - bool - whether there is a signal peptide
        present in the sequence or not, and "Raw Prediction" - str - a string with one letter per amino acid in the sequence,
        in which B/b mark detected trans-membrane beta-barrels and H/h trans-membrane alpha-helices.
    :param toxins_blast_result: A DataFrame containing peptides that were detected by sequence similarity to other, known toxins.
        It is expected to contain a column 'ID' listing the sequence IDs of the contained peptides.
    :return: A FASTA file containing those sequences that are considered candidate toxins.
    """
    is_secreted = peptides["Signal Peptide Predicted"]
    secreted_peptides = peptides[is_secreted]
    non_secreted_peptides = peptides[~is_secreted]

    # tm - trans-membrane
    is_non_tm = peptides["Raw Prediction"].apply(_has_no_transmembrane_region)
    non_tm_peptides = peptides[is_non_tm]

    secreted_candidates = _filter_fasta_file(secreted_peptides, non_tm_peptides)
    non_secreted_candidates = _filter_fasta_file(non_secreted_peptides, toxins_blast_result)

    output_file = utils.global_output(config.get("basename") + "_candidate_toxins.fasta")

    records = []
    for row in pd.concat((secreted_candidates, non_secreted_candidates)).itertuples():
        records.append(
            SeqRecord(
                Seq(row.Sequence),
                id=row.ID,
                description=""
            )
        )

    if records:
        with open(output_file, "w") as out_file:
            SeqIO.write(records, out_file, "fasta")

    return output_file


# ============================ Private functions ============================= #
def _filter_fasta_file(target: pd.DataFrame, filter_map: pd.DataFrame) -> pd.DataFrame:
    """
    Filters a DataFrame for any members of another and keeps only those that are present in the `filter_map` DataFrame.
    The filtering is not carried out in-place.

    :param target: A DataFrame with at least the column "ID" (str), which is to be filtered. All columns within this DataFrame
        are retained in the result.
    :param filter_map: A DataFrame containing a column "ID" (str) with the sequence IDs of the peptides to keep.
    :return: A DataFrame containing the all rows from the `target` DataFrame whose "ID" column has at least one match
        with the "ID" column of the `filter_map` DataFrame.
    """
    return target.merge(filter_map[["ID"]].drop_duplicates(), on="ID", how="inner")


def _has_no_transmembrane_region(aa_sequence: str) -> bool:
    """
    Decides whether a sequence has trans-membrane regions in it, based on its TMbed prediction string.
    :param aa_sequence: A string containing b/B for trans-membrane beta-barrels and h/H for trans-membrane alpha-helices.
    :return: True if a trans-membrane region was found in the prediction and False otherwise.
    """
    markers = set("bBhH")
    return markers.isdisjoint(aa_sequence)
