"""
This module contains the functionality around contaminants and generating a first, clustered set of peptides.

The following former rules are included here:
build_contaminants_db <- blast_on_contaminants <- filter_contaminants <- detect_orfs <- drop_X <- cluster_peptides

The module features a public API that yields the results of cluster_peptides.
"""
import subprocess
import tempfile
from functools import cache
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from lib import config, utils

__all__ = ["cluster_peptides"]

from lib.assemble_transcriptome import get_transcriptome_db


def _build_contaminants_database(fasta_db: str) -> Path:
    """
    Builds blast database for the removal of contaminants
    :param fasta_db: A string representing the path to the fasta database.
    :return: The path to the blast database file.
    """
    blast_db = utils.global_output(fasta_db.split("/")[-1] + ".out")

    blast_db.touch()
    subprocess.run(
        f"makeblastdb -dbtype nucl -in {fasta_db} -out {blast_db}",
        shell=True
    )

    return blast_db


def _blast_on_contaminants(blast_db: Path, contigs: Path, e_value: float, threads: int) -> Path:
    """
    performs the actual blast of the contigs against the contaminants database
    :param blast_db: The path to the blast database file - the return value from build_contaminants_database
    :param contigs: The path leading to a fasta file holding transcriptome data.
    :param e_value: Contamination e-value
    :param threads: Number of threads available to this function
    :return: The path to the result of running blast.
    """
    blast_result = utils.global_output(config.get("basename") + ".blastsnuc.out")

    subprocess.run(
        f"blastn -db {blast_db} -query {contigs} -out {blast_result} -outfmt 6 -evalue {e_value} -max_target_seqs 1 -num_threads {threads}",
        shell=True
    )

    return blast_result


def _filter_contaminants(blast_result: Path, contigs: Path) -> Path:
    """
    performs the actual filtering
    :param blast_result: The return value from blast_on_contaminants - a path leading to the blastn results file.
    :param contigs: The path leading to a fasta file holding transcriptome data.
    :return: The path to the file containing the blastn result filtered from the contaminants.
    """
    filtered_contigs = utils.global_output(config.get("basename") + ".filtered.fasta")

    records = []
    with open(blast_result) as infile:
        for line in infile:
            line = line.rstrip()
            if line[0] != "#":
                records.append(line.split()[0])  # we recover the ID of the significan hits

    with open(filtered_contigs, "w") as outfile:
        for rec in SeqIO.parse(contigs, "fasta"):
            if rec.id not in records:
                SeqIO.write(rec, outfile, "fasta")

    return filtered_contigs


def _detect_orfs(nucleotide_sequences: Path, *, min_len=99, max_len=30_000_000, threads, mmseqs_path: Path) -> Path:
    """
    finds complete orfs within the input nucleotide sequences.
    DeTox were testing this with orfipy instead of orffinder to leverage multithreading
    :param nucleotide_sequences: The path leading to a fasta file holding transcriptome data/blastn results.
    :param min_len: minimum length in nucleotides an orf needs to have to be considered
    :param max_len: maximum length in nucleotides an orf needs to have to be considered
    :param threads: Number of threads available to this function
    :param mmseqs_path: The path to the folder containing the MMSeqs binary
    :return: The path leading to the faa file holding the nucleotide sequences that are considered to be complete orfs.
    """
    orfs_db = utils.global_output("mmseqs/" + config.get("basename") + "-orfs.db")
    orfs_db.resolve().parent.mkdir(parents=True, exist_ok=True)
    orfs_fasta = utils.global_output(config.get("basename") + "-orfs.faa")

    command = [
        f"{mmseqs_path}/mmseqs", "extractorfs",
        nucleotide_sequences,
        orfs_db,
        "--min-length", f"{min_len // 3}",
        "--max-length", f"{max_len // 3}",
        "--translate", "1",
        "--use-all-table-starts", "false",
        "--create-lookup", "1",
        "--threads", f"{threads}"
    ]
    subprocess.run(command)

    with tempfile.NamedTemporaryFile(suffix=".faa", delete_on_close=False) as temp:
        command = [
            f"{mmseqs_path}/mmseqs", "convert2fasta",
            orfs_db,
            temp.name
        ]
        subprocess.run(command)

        temp.seek(0)
        counter = 0
        previous_id = ""
        target = pd.read_csv(
            nucleotide_sequences.with_suffix(".db.lookup"),
            sep="\t",
            header=None,
            usecols=[0, 1],
            names=["id", "seq_id"],
            dtype=str
        )

        with open(orfs_fasta, "w") as result:
            for seq in SeqIO.parse(temp.name, "fasta"):
                seq_id = target[target["id"] == seq.id]["seq_id"].iloc[0]
                if seq.id == previous_id:
                    counter += 1
                else:
                    counter = 0
                    previous_id = seq.id

                seq.id = f"{seq_id}_ORF.{counter}"

                SeqIO.write(seq, result, "fasta")

    return orfs_fasta


def _drop_x(orfs_db: Path) -> Path:
    """
    Removes all entries from the ORFS DB file that contain at least one `X`
    :param orfs_db: The path to an faa file holding nucleotide sequences.
    :return: The path to the faa file that holds the filtered nucleotide sequences.
    """
    drop_sequence = utils.global_output(config.get('basename') + "_noX.faa")

    with open(drop_sequence, "w") as outfile:
        for seq in SeqIO.parse(orfs_db, "fasta"):
            if "X" not in seq.seq:
                SeqIO.write(seq, outfile, "fasta")

    return drop_sequence


def _cluster_peptides(aa_sequences: Path, clustering_threshold: float, max_memory: int, threads: int, *, mmseqs_path: Path) -> Path:
    """
    runs cd-hit on predicted peptide to remove excess redundancy
    :param aa_sequences: The path to an faa file holding nucleotide sequences.
    :param clustering_threshold: The sequence identity threshold that guides cluster generation.
        Sequences that are sufficiently similar to pass this threshold are clustered together.
    :param max_memory: The maximum amount of RAM that this function may use, in Gigabyte (GB)
    :param threads: The number of threads available to this function
    :param mmseqs_path: The path to the folder containing the MMSeqs binary
    :return: The path to the faa file that holds the clustered nucleotide sequences/the sequences representing the clusters.
    """
    clustered_aa_sequences = utils.global_output(config.get("basename") + ".clustered.faa")
    tmp_dir = utils.global_output("mmseqs")

    command = [
        f"{mmseqs_path}/mmseqs", "easy-linclust",
        aa_sequences,  # Input FASTA
        clustered_aa_sequences, # "Cluster Prefix"
        tmp_dir,
        "--min-seq-id", f"{clustering_threshold}",
        "--threads", f"{threads}",
        "--cluster-mode", "2", # Simulates CD-Hit's approach to clustering.
        "--dbtype", "1" # As we only expect AA sequences to arrive in this function, we can guide this a little bit.
        "--createdb-mode", "0", # As the input may contain multi-line sequences, this is safer (also, using more disk space is cheaper than rerunning the pipeline)
        "--remove-tmp-files", "false", # Helps debugging. Or just with understanding what the pipeline does.
    ] # The -d 40 flag is just for display purposes and thus not relevant for us.

    subprocess.run(command)

    return clustered_aa_sequences


@cache
def cluster_peptides(transcriptome: Path):
    """
    Runs the complete section that checks for contaminants and finally clusters peptides.
    :param transcriptome: The path to a FASTA file containing transcriptome data, e.g. the result of `assemble_transcriptome`
    :return: A FASTA file containing amino acid sequences of the representatives of the different peptide clusters.
    """
    threads = utils.get_threads()
    mmseqs_path = config.get_path("mmseqs_path") or Path("software/mmseqs/bin")
    try:
        subprocess.run([f"{mmseqs_path}/mmseqs", "version"])
    except FileNotFoundError:
        raise FileNotFoundError("Could not find MMSeqs binary!")


    aa_sequences = config.get("proteome_fasta")
    if aa_sequences is None:

        # If there is a config entry for the contaminants, work with that value, else just go with the transcriptome.
        contaminants = config.get("contaminants")
        if contaminants is not None:
            blast_db = _build_contaminants_database(contaminants)

            e_value = config.get("contamination_evalue")
            blast_result = _blast_on_contaminants(blast_db, transcriptome, e_value, threads)

            nucleotide_sequences = get_transcriptome_db(_filter_contaminants(blast_result, transcriptome), mmseqs_path=mmseqs_path)
        else:
            nucleotide_sequences = get_transcriptome_db(transcriptome, mmseqs_path=mmseqs_path)

        # Read length limits for open reading frames from config, but keep only if available.
        frame_size: dict[str, int] = {
            k: v for k, v in {
                "min_len": config.get("minlen"),
                "max_len": config.get("maxlen")
            }.items() if v is not None
        }

        aa_sequences = _detect_orfs(nucleotide_sequences, threads=threads, mmseqs_path=mmseqs_path, **frame_size)
    else:
        aa_sequences = get_transcriptome_db(aa_sequences)

    return _cluster_peptides(_drop_x(aa_sequences), config.get("clustering_threshold"),
                             config.get("memory"), threads, mmseqs_path=mmseqs_path)
