"""
This module contains all functionality around running BLASTn.

This contains the former rule blocks
build_toxin_blast_db <- blast_on_toxins
and
download_uniprot <- make_uniprot_blast_database <- blast_on_uniprot

It features a public API that yields the results of blast_on_toxins and blast_on_uniprot.
"""
import subprocess
import sys
import tempfile
from functools import cache
from pathlib import Path

import pandas as pd

from lib import config, utils

__all__ = ["blast_on_toxins", "blast_on_uniprot"]


# ============================= Public functions ============================= #
@cache
def blast_on_toxins(filtered_clustered_aa_sequences: Path) -> pd.DataFrame:
    """
    Runs MMSeqs2 search against the toxin database, replacing diamond blastp.
    The query are the peptides without any signal sequence.

    The function runs MMSeqs and extracts the fasta at the same time.

    :param filtered_clustered_aa_sequences: A filtered and clustered FASTA file containing amino acids.
    :return: A DataFrame containing the 'BLAST' results.
    """
    db_file = config.get_path("toxin_db")
    if db_file is None:
        print("Missing config value 'toxin_db' pointing at the Toxins blast database file!", file=sys.stderr)
        exit(1)

    mmseqs_path = utils.ensure_mmseqs2(Path("software/mmseqs/bin"))
    threads = utils.get_threads()
    e_value = config.get("toxins_evalue", 1e-10)

    columns = ["qseqid", "toxinDB_sseqid", "toxinDB_pident", "toxinDB_evalue"]

    return _run_blast(filtered_clustered_aa_sequences, db_file, e_value, threads, columns, mmseqs_path=mmseqs_path)


@cache
def blast_on_uniprot(toxin_candidates: Path) -> pd.DataFrame:
    """
    Runs MMSeqs search against UniProt/Swiss-Prot. This replaces the formerly used diamond blastp command.
    :param toxin_candidates: The proteins that are deemed potential toxins.
    :return: A DataFrame containing all those proteins that have a similarity match to a protein in Swiss-Prot
    """
    mmseqs_path = utils.ensure_mmseqs2(Path("software/mmseqs/bin"))

    db_file = config.get_path("swissprot_db_path") or _download_uniprot(mmseqs_path=mmseqs_path)

    threads = utils.get_threads()
    e_value = config.get("swissprot_evalue", 1e-10)

    columns = ["qseqid", "uniprot_sseqid", "uniprot_pident", "uniprot_evalue"]

    return _run_blast(toxin_candidates, db_file, e_value, threads, columns, mmseqs_path=mmseqs_path)


# ============================ Private functions ============================= #
def _run_blast(aa_sequences: Path, db: Path, e_value: float, threads: int, columns: list[str], *,
               mmseqs_path: Path) -> pd.DataFrame:
    """
    Runs BLASTp against the FASTA file given in `aa_sequences` on the database provided in `db`.
    :param aa_sequences: The path to a FASTA file containing the query sequences for the BLASTp run.
    :param db: The database on which BLASTp is to operate
    :param e_value: The e_value for BLASTp
    :param threads: The number of threads BLASTp shall use
    :param columns: The column names for the returned DataFrame, as list of strings.
    :return: A pandas DataFrame containing the BLASTp results.
    """
    blast_result: pd.DataFrame

    with tempfile.NamedTemporaryFile(suffix=".m8", delete_on_close=False) as result_file:
        command = [
            f"{mmseqs_path}/mmseqs", "easy-search",
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
        ]
        subprocess.run(command)

        blast_result = pd.read_csv(result_file.name,
                                   sep="\t",
                                   header=None,
                                   names=columns,
                                   index_col=0)
    return blast_result


@cache
def _download_uniprot(*, mmseqs_path: Path) -> Path:
    """
    Downloads the latest release of the UniProt database
    """
    mmseqs_dir = utils.global_output("mmseqs")
    mmseqs_dir.mkdir(parents=True, exist_ok=True)

    db_file = mmseqs_dir / "uniprot.db"

    command = [
        f"{mmseqs_path}/mmseqs", "databases",
        "UniProtKB/Swiss-Prot",
        db_file,
        mmseqs_dir,
        "--threads", "16"
    ]
    subprocess.run(command)

    return db_file
