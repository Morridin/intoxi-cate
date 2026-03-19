# Copyright 2026 Paul Zanner
# Adopted and redesigned from DeToX, Ringeval et al., 2024, https://doi.org/10.1093/bib/bbae094
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
This module contains the functionality around contaminants and generating a first, clustered set of peptides.

The following former rules are included here:
build_contaminants_db <- blast_on_contaminants <- filter_contaminants <- detect_orfs <- drop_X <- cluster_peptides

The module features a public API that yields the results of cluster_peptides.
"""
import math
import subprocess
import tempfile
from functools import cache
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from . import config, utils, blast, transcriptome

__all__ = ["cluster_peptides"]


# ============================= Public functions ============================= #
@cache
def cluster_peptides(transcriptome_file: Path):
    """
    Runs the complete section that checks for contaminants and finally clusters peptides.
    :param transcriptome_file: The path to a FASTA file containing transcriptome data, e.g. the result of `assemble_transcriptome`
    :return: A FASTA file containing amino acid sequences of the representatives of the different peptide clusters.
    """
    threads = utils.get_threads()
    mmseqs_path = utils.ensure_mmseqs2()

    aa_sequences = config.get("proteome_fasta")
    if aa_sequences is None:

        # If there is a config entry for the contaminants, work with that value, else just go with the transcriptome.
        contaminants = config.get("contaminants")
        if contaminants is not None:
            blast_result = blast.on_contaminants(transcriptome_file)

            filtered_contigs = _filter_contaminants(blast_result, transcriptome_file)

            nucleotide_sequences = transcriptome.get_db(filtered_contigs, mmseqs_path=mmseqs_path)
        else:
            nucleotide_sequences = transcriptome.get_db(transcriptome_file, mmseqs_path=mmseqs_path)

        frame_size: dict[str, int] = {
            "min_len": config.get("minlen", 99),
            "max_len": config.get("maxlen", 45_000)
        }

        aa_sequences = _detect_orfs(nucleotide_sequences, threads=threads, mmseqs_path=mmseqs_path, **frame_size)
    else:
        aa_sequences = _drop_x(aa_sequences)

    return _cluster_peptides(aa_sequences, config.get("clustering_threshold", 0.99),
                             config.get("memory"), threads, mmseqs_path=mmseqs_path)


# ============================ Private functions ============================= #
def _filter_contaminants(blast_result: pd.DataFrame, contigs: Path) -> Path:
    """
    Filters the input contigs against the BLAST result and keeps only those reads/sequences that are not within the BLAST result.
    :param blast_result: The return value from blast_on_contaminants - a DataFrame containing the MMSeqs search result.
    :param contigs: The path leading to a fasta file holding transcriptome data.
    :return: The path to the file containing the blastn result filtered from the contaminants.
    """
    filtered_contigs = utils.global_output(config.get("basename") + ".filtered.fasta")

    with open(filtered_contigs, "w") as outfile:
        for rec in SeqIO.parse(contigs, "fasta"):
            if rec.id not in blast_result.index:
                SeqIO.write(rec, outfile, "fasta")

    return filtered_contigs


def _detect_orfs(nucleotide_sequences: Path, *, min_len=99, max_len=30_000_000, threads, mmseqs_path: str) -> Path:
    """
    Finds complete orfs within the input nucleotide sequences.
    Finds also such orfs with open 3' or 5' ends.
    Uses MMSeqs to do this.
    The result is also directly filtered from any sequence containing an unknown amino acid marked by X.
    :param nucleotide_sequences: The path leading to a fasta file holding transcriptome data/blastn results.
    :param min_len: minimum length in nucleotides an orf needs to have to be considered
    :param max_len: maximum length in nucleotides an orf needs to have to be considered
    :param threads: Number of threads available to this function
    :param mmseqs_path: A string containing a path to the MMSeqs executable such that copying it into terminal runs MMSeqs.
    :return: The path leading to the faa file holding the nucleotide sequences that are considered to be complete orfs.
    """
    orfs_db = utils.global_output("mmseqs/" + config.get("basename") + "-orfs.db")
    orfs_db.resolve().parent.mkdir(parents=True, exist_ok=True)
    orfs_fasta = utils.global_output(config.get("basename") + "-orfs.faa")

    command = [
        mmseqs_path, "extractorfs",
        nucleotide_sequences,
        orfs_db,
        "--min-length", f"{min_len // 3}",
        "--max-length", f"{max_len // 3}",
        "--translate", "1",
        "--use-all-table-starts", "false",
        "--create-lookup", "1",
        "--threads", f"{threads}"
    ]
    subprocess.run(command, check=True)

    with tempfile.NamedTemporaryFile(suffix=".faa", delete_on_close=False) as temp:
        command = [
            mmseqs_path, "convert2fasta",
            orfs_db,
            temp.name
        ]
        subprocess.run(command, check=True)

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

        id_map = dict(zip(target["id"], target["seq_id"]))

        with open(orfs_fasta, "w") as result:
            for seq in SeqIO.parse(temp.name, "fasta"):
                if "X" in seq.seq:
                    continue
                # Every sequence in the MMSeqs2 output FASTA either has a corresponding original (read) ID or the MMSeqs database is broken.
                # Hence, the program shall fail in this case.
                seq_id = id_map[seq.id]
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


def _cluster_peptides(aa_sequences: Path, min_sequence_identity: float, max_memory: int, threads: int, *,
                      mmseqs_path: str) -> Path:
    """
    Runs MMSeqs LinClust on predicted peptides to remove excess redundancy
    :param aa_sequences: The path to an faa file holding nucleotide sequences.
    :param min_sequence_identity: The sequence identity threshold that guides cluster generation.
        Sequences that are sufficiently similar to pass this threshold are clustered together.
    :param max_memory: The maximum amount of RAM that this function may use, in Gigabyte (GB)
    :param threads: The number of threads available to this function
    :param mmseqs_path: A string containing a path to the MMSeqs executable such that copying it into terminal runs MMSeqs.
    :return: The path to the faa file that holds the clustered nucleotide sequences/the sequences representing the clusters.
    """
    output_prefix = utils.global_output(config.get('basename'))  # MMSeqs adds file endings on their own.
    tmp_dir = utils.global_output("mmseqs")

    clustering_threshold = config.get('mmseqs_max_sequence_coverage', 0.8)

    command = [
        mmseqs_path, "easy-linclust",
        aa_sequences,  # Input FASTA
        output_prefix,  # Prefix that is prepended to all the output files (which are 3:
        # 1 FASTA with all seq ids (..._all_seqs.fasta),
        # 1 TSV mapping each seq ID to the ID of the cluster representative (..._cluster.tsv),
        # 1 FASTA containing only the cluster representatives (..._rep_seq.fasta))
        tmp_dir,
        "--min-seq-id", f"{min_sequence_identity}",
        "--threads", f"{threads}",
        "--cluster-mode", "2",  # Simulates CD-Hit's approach to clustering.
        "--dbtype", "1",  # As we only expect AA sequences to arrive in this function, we can guide this a little bit.
        "--createdb-mode", "0",
        "--cov-mode", "1",
        "-c", f"{clustering_threshold}",
        # As the input may contain multi-line sequences, this is safer (also, using more disk space is cheaper than rerunning the pipeline)
        "--remove-tmp-files", "false",  # Helps debugging. Or just with understanding what the pipeline does.
        "--split-memory-limit", f"{math.floor(max_memory * 0.8)}G",
    ]  # The -d 40 flag is just for display purposes and thus not relevant for us.

    subprocess.run(command, check=True)

    return output_prefix.with_name(output_prefix.name + "_rep_seq.fasta")
