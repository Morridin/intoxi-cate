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
This module contains all functionality around running BLASTp and BLASTn.

Originally, this only included the task groups of running BLASTp on a toxins database and Swiss-Prot.
However, after switching to MMSeqs2 for this task, the command to replace BLASTn was equal to that of the former
BLASTp runs, so it was moved here.
"""
import math
import subprocess
import sys
import tempfile
from enum import IntEnum
from functools import cache
from pathlib import Path

import pandas as pd

from . import config, utils

__all__ = ["on_contaminants", "on_toxins", "on_uniprot"]


# ============================= Public functions ============================= #
@cache
def on_contaminants(contigs: Path) -> pd.DataFrame:
    """
    Runs MMSeqs search against contaminants database, replacing BLASTn for the sake of reducing dependencies.
    Goal of this function is to obtain a list of contigs that are similar in sequence to contigs known to be contaminants.
    :param contigs: A path leading to a FASTA file holding transcriptome data.
    :return:
    """
    db_file = config.get_path("contaminants")
    if db_file is None:
        print("Missing config value 'contaminants' pointing at the contaminants (blast) database file!",
              file=sys.stderr)
        exit(1)

    e_value = config.get("contamination_evalue", 1e-5)
    columns = ["ID", "c1", "c2", "c3"]

    return _run(contigs, db_file, e_value, SearchType.NUCLEOTIDE, columns)


@cache
def on_toxins(filtered_clustered_aa_sequences: Path) -> pd.DataFrame:
    """
    Runs MMSeqs2 search against the toxin database, replacing diamond blastp.
    The query are the peptides without any signal sequence.

    The function runs MMSeqs and extracts the fasta at the same time.

    :param filtered_clustered_aa_sequences: A filtered and clustered FASTA file containing amino acids.
    :return: A DataFrame containing the 'BLAST' results.
    """
    db_file = config.get_path("toxin_db")
    if db_file is None:
        print("Missing config value 'toxin_db' pointing at the toxins (blast) database file!", file=sys.stderr)
        exit(1)

    e_value = config.get("toxins_evalue", 1e-10)
    columns = ["qseqid", "toxinDB_sseqid", "toxinDB_pident", "toxinDB_evalue"]

    return _run(filtered_clustered_aa_sequences, db_file, e_value, SearchType.AMINO_ACID, columns)


@cache
def on_uniprot(toxin_candidates: Path) -> pd.DataFrame:
    """
    Runs MMSeqs search against UniProt/Swiss-Prot. This replaces the formerly used diamond blastp command.
    :param toxin_candidates: The proteins that are deemed potential toxins.
    :return: A DataFrame containing all those proteins that have a similarity match to a protein in Swiss-Prot
    """

    db_file = config.get_path("swissprot_db_path") or _download_uniprot()
    e_value = config.get("swissprot_evalue", 1e-10)
    columns = ["qseqid", "uniprot_sseqid", "uniprot_pident", "uniprot_evalue"]

    return _run(toxin_candidates, db_file, e_value, SearchType.AMINO_ACID, columns)


# ============================ Private functions ============================= #
class SearchType(IntEnum):
    NUCLEOTIDE = 1,
    AMINO_ACID = 3


def _run(aa_sequences: Path, db: Path, e_value: float, search_type: SearchType,
         columns: list[str]) -> pd.DataFrame:
    """
    Runs BLASTp against the FASTA file given in `aa_sequences` on the database provided in `db`.
    :param aa_sequences: The path to a FASTA file containing the query sequences for the BLASTp run.
    :param db: The database on which BLASTp is to operate
    :param e_value: The e_value for BLASTp
    :param columns: The column names for the returned DataFrame, as list of strings.
    :return: A pandas DataFrame containing the BLASTp results.
    """
    mmseqs_path = utils.ensure_mmseqs2()
    threads = utils.get_threads()
    memory = config.get("memory")

    blast_result: pd.DataFrame

    with tempfile.NamedTemporaryFile(suffix=".m8", delete_on_close=False) as result_file:
        command = [
            mmseqs_path, "easy-search",
            aa_sequences,
            db,
            result_file.name,
            utils.global_output("mmseqs"),
            "-s", "5.7",  # In preparation for later adjustments
            "-e", f"{e_value}",  # Replaces --evalue
            "--max-accept", "1",  # Replaces --max-target-seqs
            "--format-output", "query,target,pident,evalue",
            # replaces the --outfmt param (output is already in tabular format by default)
            "--threads", f"{threads}",
            "--split-memory-limit", f"{math.floor(memory * 0.8)}G",
            "--search-type", f"{search_type}",
        ]
        subprocess.run(command)

        blast_result = pd.read_csv(result_file.name,
                                   sep="\t",
                                   header=None,
                                   names=columns,
                                   index_col=0)
    return blast_result


@cache
def _download_uniprot() -> Path:
    """
    Downloads the latest release of the UniProt database
    """
    mmseqs_path = utils.ensure_mmseqs2()
    mmseqs_dir = utils.global_output("mmseqs")
    mmseqs_dir.mkdir(parents=True, exist_ok=True)

    db_file = mmseqs_dir / "uniprot.db"

    command = [
        mmseqs_path, "databases",
        "UniProtKB/Swiss-Prot",
        db_file,
        mmseqs_dir,
        "--threads", "16"
    ]
    subprocess.run(command)

    return db_file
