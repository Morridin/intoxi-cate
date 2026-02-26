"""
This module contains the functionality around contaminants and generating a first, clustered set of peptides.

The following former rules are included here:
build_contaminants_db <- blast_on_contaminants <- filter_contaminants <- detect_orfs <- drop_X <- cluster_peptides

The module features a public API that yields the results of cluster_peptides.
"""
import subprocess
from functools import cache
from pathlib import Path

from Bio import SeqIO

from lib import config, utils

__all__ = ["cluster_peptides"]


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


def _detect_orfs(nucleotide_sequences: Path, *, min_len=99, max_len=30_000_000, threads) -> Path:
    """
    finds complete orfs within the input nucleotide sequences.
    DeTox were testing this with orfipy instead of orffinder to leverage multithreading
    :param nucleotide_sequences: The path leading to a fasta file holding transcriptome data/blastn results.
    :param min_len: minimum length an orf needs to have to be considered
    :param max_len: maximum length an orf needs to have to be considered
    :param threads: Number of threads available to this function
    :return: The path leading to the faa file holding the nucleotide sequences that are considered to be complete orfs.
    """
    out_dir = utils.global_output(".")
    aa_sequences = Path(config.get("basename") + ".faa")

    subprocess.run(
        f"orfipy --procs {threads} --start ATG --partial-3 --partial-5 --pep {aa_sequences} --min {min_len} --max {max_len} {nucleotide_sequences} --outdir {out_dir}",
        shell=True
    )

    return out_dir / aa_sequences


def _drop_x(aa_sequences: Path) -> Path:
    """
    removes all entries from a nucleotide sequence file that contain at least one `X`
    :param aa_sequences: The path to an faa file holding nucleotide sequences.
    :return: The path to the faa file that holds the filtered nucleotide sequences.
    """
    drop_sequence = utils.global_output(config.get('basename') + "_noX.faa")

    with open(drop_sequence, "w") as outfile:
        for seq in SeqIO.parse(aa_sequences, "fasta"):
            if "X" not in seq.seq:
                SeqIO.write(seq, outfile, "fasta")

    return drop_sequence


def _cluster_peptides(aa_sequences: Path, clustering_threshold: float, max_memory: int, threads: int) -> Path:
    """
    runs cd-hit on predicted peptide to remove excess redundancy
    :param aa_sequences: The path to an faa file holding nucleotide sequences.
    :param clustering_threshold: The sequence identity threshold that guides cluster generation.
        Sequences that are sufficiently similar to pass this threshold are clustered together.
    :param max_memory: The maximum amount of RAM that this function may use, in Gigabyte (GB)
    :param threads: The number of threads available to this function
    :return: The path to the faa file that holds the clustered nucleotide sequences/the sequences representing the clusters.
    """
    filtered_aa_sequences = utils.global_output(config.get("basename") + ".clustered.fasta")

    subprocess.run(
        f"cd-hit -i {aa_sequences} -o {filtered_aa_sequences} -c {clustering_threshold} -M {max_memory * 1000} -T {threads} -d 40",
        shell=True
    )

    return filtered_aa_sequences

@cache
def cluster_peptides(transcriptome: Path):
    """
    Runs the complete section that checks for contaminants and finally clusters peptides.
    :param transcriptome: The path to a FASTA file containing transcriptome data, e.g. the result of `assemble_transcriptome`
    :return: A FASTA file containing amino acid sequences of the representatives of the different peptide clusters.
    """
    threads = utils.get_threads()

    aa_sequences = config.get("proteome_fasta")
    if aa_sequences is None:

        # If there is a config entry for the contaminants, work with that value, else just go with the transcriptome.
        contaminants = config.get("contaminants")
        if contaminants is not None:
            blast_db = _build_contaminants_database(contaminants)

            e_value = config.get("contamination_evalue")
            blast_result = _blast_on_contaminants(blast_db, transcriptome, e_value, threads)

            nucleotide_sequences = _filter_contaminants(blast_result, transcriptome)
        else:
            nucleotide_sequences = transcriptome

        # Read length limits for open reading frames from config, but keep only if available.
        frame_size: dict[str, int] = {
            k: v for k, v in {
                "min_len": config.get("minlen"),
                "max_len": config.get("maxlen")
            }.items() if v is not None
        }

        aa_sequences = _detect_orfs(nucleotide_sequences, threads=threads, **frame_size)

    return _cluster_peptides(_drop_x(aa_sequences), config.get("clustering_threshold"), config.get("memory"), threads)
