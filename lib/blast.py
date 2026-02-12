"""
This module contains all functionality around running BLASTn.

This contains the former rule blocks
build_toxin_blast_db <- blast_on_toxins
and
download_uniprot <- make_uniprot_blast_database <- blast_on_uniprot

It features a public API that yields the results of blast_on_toxins and blast_on_uniprot.
"""
import subprocess
import tempfile
from functools import cache
from pathlib import Path

import httpx
import pandas as pd

from lib import config, utils

__all__ = ["blast_on_toxins", "blast_on_uniprot"]


# ======================== Block 0: Common Functions ========================= #
def _build_blast_db(source_file: Path, target_file: Path):
    """
    Builds a BLAST database using diamond.
    There is no return value as the target_path parameter would else be.
    :param source_file: The path to the (usually) FASTA file from which diamond shall build the database.
    :param target_file: The path to the file where diamond shall store the newly created DB.
    """
    subprocess.run(
        f"diamond makedb --in {source_file} --db {target_file}"
    )


def _run_blast(aa_sequences: Path, db: Path, e_value: float, threads: int, columns: str) -> pd.DataFrame:
    """
    Runs BLASTp against the FASTA file given in `aa_sequences` on the database provided in `db`.
    :param aa_sequences: The path to a FASTA file containing the query sequences for the BLASTp run.
    :param db: The database on which BLASTp is to operate
    :param e_value: The e_value for BLASTp
    :param threads: The number of threads BLASTp shall use
    :param columns: The column names for the returned DataFrame, as one, tab-separated string.
    :return: A pandas DataFrame containing the BLASTp results.
    """
    blast_result: pd.DataFrame

    with tempfile.NamedTemporaryFile(suffix=".tsv", delete_on_close=False) as result_file:
        result_file.write(columns)

        subprocess.run(
            f"diamond blastp -d {db} -q {aa_sequences} --evalue {e_value} --max-target-seqs 1 --threads {threads} --outfmt 6 qseqid sseqid pident evalue >> {result_file.name}",
            shell=True
        )

        blast_result = pd.read_csv(result_file, sep="\t", index_col=0)

    return blast_result


# ============================= Block 1: Toxins ============================== #
@cache
def blast_on_toxins(filtered_clustered_aa_sequences: Path) -> pd.DataFrame:
    """
    Runs BLASTp against the toxin database.
    The query are the peptides without any signal sequence.

    The rule runs blast and extracts the fasta at the same time.

    :param filtered_clustered_aa_sequences: A filtered and clustered FASTA file containing amino acids.
    :return: A DataFrame containing the BLAST results.
    """
    db_source: Path = config.get("toxin_db")
    db_file = utils.global_output(db_source.name + ".dmnd")
    _build_blast_db(db_source, db_file)

    threads = utils.get_threads()
    e_value = config.get("toxins_evalue")
    if e_value is None:
        e_value = 1e-10

    columns = "qseqid\ttoxinDB_sseqid\ttoxinDB_pident\ttoxinDB_evalue\n"

    return _run_blast(filtered_clustered_aa_sequences, db_file, e_value, threads, columns)


# ============================= Block 2: UniProt ============================= #
def _download_uniprot(url: str) -> Path:
    """
    Downloads the latest release of the UniProt database
    """
    db_dir = utils.global_output("databases/uniprot")
    db_dir.mkdir(parents=True, exist_ok=True)

    database = db_dir / "uniprot_sprot.fasta.gz"

    with httpx.stream("GET", url) as response:
        with open(database, "wb") as f:
            for data in response.iter_bytes():
                f.write(data)

    return database


@cache
def blast_on_uniprot(toxin_candidates: Path) -> pd.DataFrame:
    """

    :param toxin_candidates:
    :return:
    """
    url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"

    db_source = _download_uniprot(url)
    db_file = utils.global_output("uniprot_blast_db.dmnd")
    _build_blast_db(db_file, db_source)

    threads = utils.get_threads()
    e_value = config.get("swissprot_evalue")
    if e_value is None:
        e_value = 1e-10

    columns = "qseqid\tuniprot_sseqid\tuniprot_pident\tuniprot_evalue"

    return _run_blast(toxin_candidates, db_file, e_value, threads, columns)
