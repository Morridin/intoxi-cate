"""
This module bundles the hmmer related rules:
download_pfam <- run_hmmer <- parse_hmmsearch_output

The public API yields the results of parse_hmmsearch_output, but as pandas DataFrame.
"""
import gzip
import subprocess
from functools import cache
from pathlib import Path

import httpx
import pandas as pd

from lib import utils, config

__all__ = ["hmmer"]


def _download_pfam():
    """
    Loads an HMM for later use with hmmsearch from PFAM.
    :return: The path to the uncompressed HMM
    """
    db_dir = utils.global_output("databases/pfam")
    db_dir.mkdir(parents=True, exist_ok=True)

    pfam_db = db_dir / "Pfam-A.hmm"

    url = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"

    response = httpx.get(url)

    data = gzip.decompress(response.content)

    with open(pfam_db, "wb") as outfile:
        outfile.write(data)

    return pfam_db


def _run_hmmer(fasta_file: Path, pfam_db: Path) -> Path:
    """
    runs hmmer against the pfam database.
    :param fasta_file: A path to a file containing sequence data to run hmmer on.
    :param pfam_db: The HMM to run hmmer on.
    :return: The path to the file holding tabular domain output in TSV format.
    """
    log_out = utils.global_output(config.get("basename") + ".log")
    domain_tbl_out = utils.global_output(config.get("basename") + ".domtblout")

    threads = config.get("threads")

    subprocess.run(
        f"hmmsearch --cut_ga --cpu {threads} --domtblout {domain_tbl_out} -o {log_out} {pfam_db} {fasta_file}",
        shell=True,
        capture_output=True
    )

    return domain_tbl_out


def _parse_hmmsearch_output(domain_table: Path) -> pd.DataFrame:
    """
    parses and aggregates the hmmer output, uses the domtblout file
    :return: A DataFrame holding the aggregated output.
    """
    hmmer_domain_output = pd.read_csv(
        domain_table,
        comment="#",
        sep="\\s+",
        usecols=[0, 1, 2, 3, 4],
        names=["target name", "accession_t", "tlen", "query name", "accession_Q"]
    )

    aggregated_domains = hmmer_domain_output.groupby('target name')['query name'].apply(lambda x: set(x)).reset_index()
    aggregated_domains['pfam domains'] = aggregated_domains['query name'].apply(lambda x: "; ".join(x))

    return aggregated_domains.drop(["query name"], axis=1)

@cache
def hmmer(toxin_candidates: Path) -> pd.DataFrame:
    """
    Downloads the required HMM, runs HMMer on it against the sequences in toxin_candidates and returns the aggregated domains in a Pandas DataFrame.
    :param toxin_candidates: A Path pointing to a fasta file with toxin candidates to run HMMer against.
    :return: A DataFrame holding the aggregated domain results.
    """
    pfam_db = _download_pfam()
    domain_table = _run_hmmer(toxin_candidates, pfam_db)
    return _parse_hmmsearch_output(domain_table)
