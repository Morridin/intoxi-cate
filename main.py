import gzip
import itertools
import subprocess
import urllib.request
from os import makedirs
from os.path import commonpath
from pathlib import Path
from sys import exc_value

import pandas as pd
from typing import Iterable

from snakemake.io import expand, glob_wildcards, directory

from lib import config, utils
from lib.utils import global_output

from Bio import SeqIO


# configfile: "config.yaml"

def trim_reads(r1: Path, r2: Path | None, adapters: Path) -> dict[str, Path]:
    """
    This function trims the raw reads provided in paths r1 and r2 using what is provided in adapters.
    :param r1: The path to the forward paired-end or single-end reads file in FASTQ format.
    :param r2: The path to the reverse paired-end reads file in FASTQ format.
    :param adapters: The path to the adapters file in FASTA format.
    :return: The paths to the forward and reverse paired-end reads in FASTA format, as well as their unpaired
    counterparts, all sorted in a dict with the following keys: r1, (r1_unpaired, r2, r2_unpaired (only if r2 is not None))
    """
    assert r1 != "" and r2 != "", "Params r1 and r2 must not be empty!"

    base_path = Path("trimmed_reads")
    base_path.mkdir(parents=True, exist_ok=True)

    output = {"r1": global_output(base_path / config.get("R1").split("/")[-1])}
    num_threads = config.get("threads")

    if r2 is None:
        subprocess.run(
            f"trimmomatic SE -threads {num_threads} {r1} {output['r1']} "
            f"ILLUMINACLIP:{input.adapters}:2:40:15 LEADING:15 TRAILING:15 MINLEN:25 SLIDINGWINDOW:4:15".split(" ")
        )
        return output
    else:
        output |= {
            "r2": global_output(base_path / config.get("R2").split("/")[-1]),
            "r1_unpaired": global_output(base_path / ("unpaired." + config.get("R1").split("/")[-1])),
            "r2_unpaired": global_output(base_path / ("unpaired." + config.get("R2").split("/")[-1]))
        }
        subprocess.run(
            f"trimmomatic PE -threads {num_threads} {r1} {r2} "
            f"{output['r1']} {output['r1_unpaired']} {output['r2']} {output['r2_unpaired']} "
            f"ILLUMINACLIP:{adapters}:2:40:15 LEADING:15 TRAILING:15 MINLEN:25 SLIDINGWINDOW:4:15".split(" ")
        )
        return output


def assemble_transcriptome(r1: Path, r2: Path | None) -> Path:
    """
    Assembles a transcriptome if it is not provided. Uses Trinity
    In this case sequencing reads MUST be provided in the config.
    :param r1: the r1 return of trim_reads
    :param r2: the r2 return of trim_reads, if available.
    :return: The path to the file containing the assembled transcriptome
    """
    if r2 is None:
        reads = f"--single {r1}"
    else:
        reads = f"--left {r1} --right {r2}"

    memory = f"{config.get("memory")}G"
    seq_type = "fq"
    threads: int = config.get("threads")

    assembly = global_output("trinity_out_dir/Trinity.fasta")

    subprocess.run(
        f"Trinity --seqType {seq_type} {reads} --CPU {threads} --max_memory {memory} --output {assembly}"
    )

    return assembly


if 'contaminants' in config and config['contaminants'] not in [None, ""]:
    def build_contaminants_database(fasta_db: str = config.get("contaminants")) -> Path:
        """
        builds blast database for the removal of contaminants
        Todo: make this optional like in the original code
        :param fasta_db: A string representing the path to the fasta database.
        :return The path to the blast database file.
        """
        blast_db = global_output(fasta_db.split("/")[-1] + ".out")

        blast_db.touch()
        subprocess.run(
            f"makeblastdb -dbtype nucl -in {fasta_db} -out {blast_db}"
        )

        return blast_db


    def blast_on_contaminants(blast_db: Path, contigs: Path) -> Path:
        """
        performs the actual blast of the contigs against the contaminants database
        :param blast_db: The path to the blast database file - the return value from build_contaminants_database
        :param contigs: The path leading to a fasta file holding transcriptome data.
        :return: The path to the result of running blast.
        """
        blast_result = global_output(config.get("basename") + ".blastsnuc.out")

        e_value: float = config.get("contamination_evalue")
        threads: int = config.get("threads")

        subprocess.run(
            f"blastn -db {blast_db} -query {contigs} -out {blast_result} -outfmt 6 -evalue {e_value} -max_target_seqs 1 -num_threads {threads}"
        )

        return blast_result


    def filter_contaminants(contigs: Path, blast_result: Path) -> Path:
        """
        performs the actual filtering
        :param contigs: The path leading to a fasta file holding transcriptome data.
        :param blast_result: The return value from blast_on_contaminants - a path leading to the blastn results file.
        :return: The path to the file containing the blastn result filtered from the contaminants.
        """
        filtered_contigs = global_output(config.get("basename") + ".filtered.fasta")

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


def detect_orfs(nucleotide_sequences: Path, minlen=99, maxlen=30000000) -> Path:
    """
    finds complete orfs within the input nucleotide sequences.
    DeTox were testing this with orfipy instead of orffinder to leverage multithreading
    :param nucleotide_sequences: The path leading to a fasta file holding transcriptome data/blastn results.
    :param minlen: minimum length an orf needs to have to be considered
    :param maxlen: maximum length an orf needs to have to be considered
    :return: The path leading to the faa file holding the nucleotide sequences that are considered to be complete orfs.
    """
    aa_sequences = global_output(config.get("basename") + ".faa")

    threads = config.get("threads")

    subprocess.run(
        f"orfipy --procs {threads} --start ATG --partial-3 --partial-5 --pep {aa_sequences} --min {minlen} --max {maxlen} {nucleotide_sequences} --outdir ."
    )

    return nucleotide_sequences


def drop_X(aa_sequences: Path) -> Path:
    """
    removes all entries from a nucleotide sequence file that contain at least one `X`
    :param aa_sequences: The path to an faa file holding nucleotide sequences.
    :return: The path to the faa file that holds the filtered nucleotide sequences.
    """
    drop_sequence = global_output(config.get('basename') + "_noX.faa")

    with open(drop_sequence, "w") as outfile:
        for seq in SeqIO.parse(aa_sequences, "fasta"):
            if "X" not in seq.seq:
                SeqIO.write(seq, outfile, "fasta")

    return drop_sequence


def cluster_peptides(aa_sequences: Path, clustering_threshold: float, max_memory: int) -> Path:
    """
    runs cd-hit on predicted peptide to remove excess redundancy
    :param aa_sequences: The path to an faa file holding nucleotide sequences.
    :param clustering_threshold: The sequence identity threshold that guides cluster generation.
        Sequences that are sufficiently similar to pass this threshold are clustered together.
    :param max_memory: The maximum amount of RAM that this function may use, in Gigabyte (GB)
    :return: The path to the faa file that holds the clustered nucleotide sequences/the sequences representing the clusters.
    """
    filtered_aa_sequences = global_output(config.get("basename") + ".clustered.fasta")

    threads = config.get("threads")

    subprocess.run(
        f"cd-hit -i {aa_sequences} -o {filtered_aa_sequences} -c {clustering_threshold} -M {max_memory} -T {threads} -d 40"
    )

    return filtered_aa_sequences


def trim_peptides(aa_sequences: Path) -> Path:
    """
    trims all peptides to only the first 50 aminoacids, as they are the only useful part for signalp. This step improves load time.
    :param aa_sequences: The path to an faa file holding clustered nucleotide sequences (see cluster_peptides).
    :return: The path to an faa file holding the same nucleotide sequences, but cut off after the fiftieth aminoacid.
    """
    trimmed_sequences = global_output(config.get("basename") + ".trimmed.faa")

    with open(trimmed_sequences, "w") as outfile:
        for seq in SeqIO.parse(aa_sequences, "fasta"):
            outfile.write(f">{seq.id}\n{seq.seq[:50]}\n")

    return trimmed_sequences


def checkpoint__split_fasta(fasta_file: Path, chunk_size: int = 5000) -> Path:
    """
    ?
    :param chunk_size:
    :param fasta_file:
    :return:
    """
    split_dir = global_output("split_files")

    split_dir.mkdir(parents=True, exist_ok=True)

    threads = config.get("threads")

    subprocess.run(
        f"seqkit split2 -s {chunk_size} -O {split_dir} --by-size-prefix \"\"  -j {threads} {fasta_file}"
    )

    return split_dir


def run_signalp(fasta_file: Path, prefix: Path, signalp_path: Path) -> Path:
    """
    ?
    :param fasta_file: Default: global_output("")+"split_files/{i}.faa"
    :param prefix: Default: global_output("split_files/{i}"),
    :param signalp_path: Default: config.get("signalp_path")
    :return: ?
    """
    if not (signalp_path / "signalp").exists():
        subprocess.run(
            f"signalp -batch 5000 -fasta {fasta_file} -org euk -format short -verbose -prefix {prefix}"
        )
    else:
        subprocess.run(
            f"./{signalp_path}signalp -batch 5000 -fasta {fasta_file} -org euk -format short -verbose -prefix {prefix}"
        )
    return fasta_file.with_name(fasta_file.stem + "_summary.signalp5")


def filter_signalp_outputs(files: Iterable[Path], threshold: float = 0.7) -> Path:
    """
    filters the output of the multiple signalp runs and extracts only those with a probability of signal peptide greater than a threshold.
    Only one file should be produced from the multiple signalp files.
    Two outputs are expected: a filtered (or not?) table with the signalp results and a filtered fasta of only those peptides with a signal
    :return:
    """
    out_file = global_output(config.get("basename") + "_filtered_sigp.tsv")

    subprocess.run(
        f"awk -v b={threshold} -F'\t' '!/^#/ && !/\?/  && $3 > b' {" ".join({str(file) for file in files})} > {out_file}")

    return out_file


def extract_secreted_peptides(signalp_result: Path, fasta_file: Path) -> tuple[Path, Path]:
    """
    ?
    :param signalp_result: The output of running signalp and filtering its output
    :param fasta_file: The path to a file containing amino acid sequences as per the filtered output of cluster_peptides.
    :return:
    """
    secreted_peptides = global_output(config.get("basename") + "_secreted_peptides.fasta")
    non_secreted_peptides = global_output(config.get("basename") + "_non_secreted_peptides.fasta")

    with open(str(signalp_result)) as infile:
        records = []
        for line in infile:
            records.append(line.rstrip().split("\t")[0])
    with open(non_secreted_peptides, "w") as n_outfile:
        with open(secreted_peptides, "w") as outfile:
            for seq in SeqIO.parse(fasta_file, "fasta"):
                if seq.id in records:
                    SeqIO.write(seq, outfile, "fasta")
                else:
                    SeqIO.write(seq, n_outfile, "fasta")

    return secreted_peptides, non_secreted_peptides


def run_phobius(secreted_peptides: Path) -> Path:
    """
    Runs Phobius.
    #todo: remember to inform the user about the installation procedure. I added a dependency in the conda env with a convenient installation script
    :param secreted_peptides: The path to the file containing secreted peptides as generated by extract_secreted_peptides function.
    :return: The path to the predictions from Phobius.
    """
    table = global_output(config.get("basename") + "_phobius_predictions.tsv")
    subprocess.run(
        f"phobius.pl -short {secreted_peptides} | sed 's/\s\+/\t/g' | awk '$2 == 0' > {table}"
    )
    return table


def extract_non_tm_peptides(phobius_result: Path, fasta_file: Path) -> Path:
    """
    extracts non-TM peptides from the phobius output

    :param phobius_result: The path to the output of Phobius.
    :param fasta_file: The path to the file containing the secreted peptides as generated by extract_secreted_peptides function.
    :return: The path to the fasta file holding only non-TM peptides.
    """
    non_tm_peptides = global_output(config.get("basename") + "_non_TM_pep.fasta")

    records: list[str] = []
    with open(phobius_result) as infile:
        records = [line.rstrip().split("\t")[0] for line in infile]
    with open(non_tm_peptides, "w") as outfile:
        for seq in SeqIO.parse(fasta_file, "fasta"):
            if seq.id in records:
                SeqIO.write(seq, outfile, "fasta")

    return non_tm_peptides


def build_toxin_blast_db() -> Path:
    """
    builds a blast database for toxin prediction
    :return: The path to the generated database file
    """
    db = config.get("toxin_db")
    outfile = global_output(config.get("toxin_db").split("/")[-1] + ".dmnd")

    subprocess.run(
        f"diamond makedb --db {outfile} --in {db}"
    )

    return outfile


def blast_on_toxins(orf_fasta_clustered_file: Path, db_file: Path, e_value: float = 1e-10) -> Path:
    """
    Runs blastp against the toxin database.
    The query are the peptides without any signal sequence.

    The rule runs blast and extracts the fasta at the same time.
    Might be split in two rules for easier management.

    :param orf_fasta_clustered_file: A filtered and clustered fasta file containing amino acids.
    :param db_file: A blast database file.
    :param e_value: The evalue to use for blastp
    :return: The Path to a file containing the blast results.
    """
    threads = config.get("threads")

    blast_result = global_output(config.get("basename") + "_toxin_blast_results.tsv")

    build_header = f"echo \"qseqid\ttoxinDB_sseqid\ttoxinDB_pident\ttoxinDB_evalue\" > {blast_result}"
    command_line = f"{build_header} && diamond blastp -q {orf_fasta_clustered_file} --evalue {e_value} --max-target-seqs 1 --threads {threads} -d {db_file} --outfmt 6 qseqid sseqid pident evalue >> {blast_result}"

    subprocess.run(command_line, shell=True)

    return blast_result


def retrieve_orfs_with_blast_without_signalp(nonsec_fasta_file: Path, blast_toxins_result: Path) -> Path:
    """
    Filters `nonsec_fasta_file` for the peptides in `blast_toxins_result`.
    :param nonsec_fasta_file:
    :param blast_toxins_result:
    :return:
    """
    hits_fasta = global_output(config.get("basename") + "_toxins_by_similarity.fasta")

    records: list[str] = []
    with open(blast_toxins_result) as infile:
        records = [line.rstrip().split("\t")[0] for line in infile]
    with open(hits_fasta, "w") as outfile:
        for seq in SeqIO.parse(nonsec_fasta_file, "fasta"):
            if seq.id in records:
                SeqIO.write(seq, outfile, "fasta")

    return hits_fasta


def retrieve_candidate_toxins(non_tm_peptides: Path, hits_fasta: Path) -> Path:
    """
    this rule just creates a fasta from the positive hits in the toxin similarity and structure searches.
    :return:
    """
    output = global_output(config.get("basename") + "_candidate_toxins.fasta")

    subprocess.run(
        f"cat {non_tm_peptides} {hits_fasta} > {output}",
    )

    return output


def download_pfam():
    """
    ?
    :return:
    """
    db_dir = global_output("databases/pfam")
    db_dir.mkdir(parents=True, exist_ok=True)

    pfam_db = global_output("databases/pfam/Pfam-A.hmm")

    url = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"

    response = urllib.request.urlopen(url)
    zipped_data = response.read()

    data = gzip.decompress(zipped_data)

    with open(pfam_db, "wb") as outfile:
        outfile.write(data)

    return pfam_db


def run_hmmer(fasta_file: Path, pfam_db: Path) -> tuple[Path, Path]:
    """
    runs hmmer against the pfam database.
    :param fasta_file:
    :param pfam_db:
    :return:
    """
    tblout = global_output(config.get("basename") + ".tblout")
    domtblout = global_output(config.get("basename") + ".domtblout")

    threads = config.get("threads")

    subprocess.run(
        f"hmmsearch --cut_ga --cpu {threads} --domtblout {domtblout} --tblout {tblout} {pfam_db} {fasta_file}"
    )

    return tblout, domtblout


def parse_hmmsearch_output(domtblout: Path) -> Path:
    """
    parses and aggregates the hmmer output, uses the domtblout file
    :return:
    """

    filtered_table = global_output(config.get("basename") + ".domtblout.tsv")

    df_domtblout = pd.read_csv(
        domtblout,
        comment="#",
        delim_whitespace=True,
        usecols=[0, 1, 2, 3, 4],
        names=["target name", "accession_t", "tlen", "query name", "accession_Q"]
    )

    aggregated_domains = df_domtblout.groupby('target name')['query name'].apply(list).reset_index()
    aggregated_domains['pfam domains'] = aggregated_domains['query name'].apply(lambda x: "; ".join(set(x)))
    aggregated_domains.to_csv(f"{filtered_table}", sep="\t", index=False)

    return filtered_table


def run_wolfpsort(candidate_toxins: Path) -> Path:
    """
    runs wolfpsort on secreted peptides inferred by signalp
    :param candidate_toxins: Fasta file
    :return:
    """
    output = global_output(config.get("basename") + "_secreted_wolfpsort_prediction.tsv")

    wps_path = config.get("wolfPsort_path")
    awk = "awk '{print $1\"\t\"$2}'"

    subprocess.run(
        f"{wps_path} animal < {candidate_toxins} | {awk} > {output}",
    )

    return output


def detect_repeated_aa(fasta_file: Path, threshold: int = 5) -> Path:
    """
    this rule looks at the fasta aminoacid sequences in input and produces a table. The table reports whether some kind of repeated pattern is found in the sequences (up to 3AA long). The default threshold for repetition is 5. The input is processed with biopython
    :param fasta_file:
    :param threshold:
    :return:
    """
    repeated_aa = global_output(config.get("basename") + "_repeated_aa.tsv")

    secreted = utils.fasta_to_dataframe(f"{fasta_file}")
    secreted["Repeats1"] = secreted.apply(lambda x: _find_repetition(1, x["Sequence"], threshold), axis=1)
    secreted["Repeats2"] = secreted.apply(lambda x: _find_repetition(2, x["Sequence"], threshold), axis=1)
    secreted["Repeats3"] = secreted.apply(lambda x: _find_repetition(3, x["Sequence"], threshold), axis=1)
    secreted["Repeats"] = secreted["Repeats1"] + secreted["Repeats2"] + secreted["Repeats3"]
    secreted['RepeatsTypes'] = secreted['Repeats'].apply(lambda t: [n for (n, _) in t])
    secreted['RepeatsLengths'] = secreted['Repeats'].apply(lambda t: [n for (_, n) in t])
    secreted['RepeatsLengths'] = [','.join(map(str, l)) for l in secreted['RepeatsLengths']]
    secreted['RepeatsTypes'] = [','.join(map(str, l)) for l in secreted['RepeatsTypes']]
    secreted = secreted.drop(columns=["Repeats", "Repeats1", "Repeats2", "Repeats3"])
    secreted.to_csv(f"{repeated_aa}", index=False, sep='\t')

    return repeated_aa


def _find_repetition(size: int, seq: pd.Series, threshold: int) -> list:
    repetition = []
    for cdl in range(0, size):
        sub = [seq[i:i + size] for i in range(cdl, len(seq), size)]
        groups = itertools.groupby(sub)
        result = [(label, sum(1 for _ in group)) for label, group in groups]
        for elem, nbRep in result:
            if int(nbRep) >= threshold:
                repetition.append((elem, nbRep))
    return repetition


### optional rules


def run_salmon(transcriptome: Path, r1: Path, r2: Path | None = None):
    """
    ?
    :return:
    """

    quant_dir = global_output(config.get("basename") + "_quant")
    quant_dir.mkdir(parents=True, exist_ok=True)
    index = global_output("salmon.idx")
    index.mkdir(parents=True, exist_ok=True)
    quantification = global_output(config.get('basename') + "_quant/quant.sf")

    threads = config.get("threads")

    subprocess.run(
        f"salmon index -t {transcriptome} -i {index} -p {threads}"
    )
    if r2 is not None:
        command = f"salmon quant -i {index} -l A -1 {r1} -2 {r2} --validateMappings -o {quant_dir}"
    else:
        command = f"salmon quant -i {index} -l A -r {r1} --validateMappings -o {quant_dir}"

    subprocess.run(command)

    return quantification


def download_uniprot():
    db_dir = global_output("databases/uniprot")
    db_dir.mkdir(parents=True, exist_ok=True)

    database = global_output("databases/uniprot/uniprot_sprot.fasta.gz")

    url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"

    urllib.request.urlretrieve(url, database)


def make_uniprot_blast_database(fasta_file: Path):
    """
    builds a blast database from the uniprot fasta
    :return: 
    """
    db_file = global_output("uniprot_blast_db.dmnd")
    subprocess.run(
        f"diamond makedb --in {fasta_file} --db {db_file}"
    )


def blast_on_uniprot(fasta_file: Path, db_file: Path, e_value: float = 1e-10) -> Path:
    """
    run blast against the uniprot database, return only the best hit
    :type e_value: 
    :type db_file: 
    :type fasta_file: 
    :return: 
    """
    blast_result = global_output(config.get("basename") + "_uniprot_blast_results.tsv")

    threads = config.get('threads')

    with open(blast_result, "w") as f:
        f.write("qseqid\tuniprot_sseqid\tuniprot_pident\tuniprot_evalue")

    subprocess.run(
        f"diamond blastp -d {db_file} -q {fasta_file} --evalue {e_value} --outfmt 6 qseqid sseqid pident evalue --max-target-seqs 1 --threads {threads} >> {blast_result}"
    )

    return blast_result


# TODO: follow this comment for the rule that will wraps everything up and create the final table.
#  -> Also, in my opinion these peptides should be marked with a warning flag in the output,
#  specifying which issue affects them (e.g. “this peptide lacks a signal peptide”,
#  “this peptide contains a transmembrane domain”, etc.)


# TODO: try to run signalp during the split rule to avoid problems.
#  issue: if the process is interrupted abnormally during the run the rule is almost certain to
#  misbehave and rerun the whole thing

# this is the list with all the expected output to be put in the final table,
# will be filled depending on the configuration file.
'''
outputs = [
    (rules.run_wolfpsort.output if config['wolfpsort'] else []),
    rules.parse_hmmsearch_output.output,
    rules.blast_on_toxins.output.blast_result,
    (rules.blast_on_uniprot.output.blast_result if config['swissprot'] else []),
    rules.detect_repeated_aa.output.repeated_aa,
]
'''


def build_output_table(wolfpsort: Path, uniprot_blast: Path, tpm_threshold: int = 1000,
                       cys_pattern: bool = False, quant: Path = None):
    """
    this rule merges the tabular output of the other rules and merges it in a single table.
    It uses the outputs list defined above.
    :param cys_pattern:
    :param quant: Quantifications file from running salmon
    :param wolfpsort: Output from running WolfPSort
    :param uniprot_blast: Output from running `blast_on_uniprot`
    :param tpm_threshold: Threshold for tpm
    :return:
    """
    fasta_file = retrieve_candidate_toxins()
    signalp_result = filter_signalp_outputs()
    hmm_search = parse_hmmsearch_output()
    toxins_blast = blast_on_toxins()
    repeated_aa = detect_repeated_aa()

    extra = [wolfpsort, hmm_search, toxins_blast, uniprot_blast, repeated_aa]

    output = global_output(config.get("basename") + "_toxins.tsv")

    seq_df = utils.fasta_to_dataframe(fasta_file)
    signalp_df = pd.read_csv(
        signalp_result,
        sep="\t",
        names=["ID", "signalp_prediction", "prob_signal", "prob_nosignal", "cutsite"]
    )
    df = seq_df.merge(signalp_df, on="ID", how="left")

    for file in extra:
        dfi = pd.read_csv(file, sep="\t")
        if "Sequence" in dfi.columns:
            dfi = dfi.drop(columns=["Sequence"])
        new_cols = [col for col in dfi.columns]
        new_cols[0] = "ID"
        dfi.columns = new_cols
        df = df.merge(dfi, how="left", on="ID")

    if cys_pattern:
        df["cutsite"] = df["cutsite"].fillna("")
        df['cut_site_position'] = df['cutsite'].apply(
            lambda x: int(x.split(" ")[2].split("-")[-1][:-1]) if "pos:" in x else -1)
        df['mature_peptide'] = df.apply(
            lambda x: x['Sequence'][x['cut_site_position']:] if x['cut_site_position'] > 0 else None, axis=1)
        df['Cys_pattern'] = df['mature_peptide'].apply(lambda x: utils.get_cys_pattern(x) if pd.notna(x) else None)

    if quant is not None:
        df["contig"] = df['ID'].apply(lambda x: x.split("_ORF")[0])
        q = pd.read_csv(f"{quant}", sep="\t")
        new_cols = [col for col in q.columns]
        new_cols[0] = "contig"
        q.columns = new_cols
        df = df.merge(q, how="left", on="contig")
        df = df.drop(['EffectiveLength', 'NumReads'], axis=1)

    df = df.assign(Rating="")
    df['Rating'] = df.apply(
        lambda row: str(row['Rating'] + 'S') if pd.notna(row['signalp_prediction']) else str(row['Rating'] + '*'),
        axis=1
    )
    df['Rating'] = df.apply(
        lambda row: str(row['Rating'] + 'B') if pd.notna(row['toxinDB_sseqid']) else row['Rating'],
        axis=1
    )
    if 'Cys_pattern' in df.columns:
        df['Rating'] = df.apply(
            lambda row: str(row['Rating'] + 'C') if pd.notna(row['Cys_pattern']) else row['Rating'],
            axis=1
        )
    if 'TPM' in df.columns:
        df['Rating'] = df.apply(
            lambda row: str(row['Rating'] + 'T') if (float(row['TPM']) >= float(f"{tpm_threshold}")) else row['Rating'],
            axis=1
        )
    df['Rating'] = df.apply(
        lambda row: str(row['Rating'] + 'D') if pd.notna(row['pfam domains']) else row['Rating'],
        axis=1
    )
    if 'uniprot_sseqid' in df.columns:
        df['Rating'] = df.apply(
            lambda row: str(row['Rating'] + '!') if pd.notna(row['uniprot_sseqid']) and pd.isna(
                row['toxinDB_sseqid']) else row['Rating'],
            axis=1
        )
    df = df.drop(['cut_site_position', 'query name'], axis=1)
    df.rename(columns={'k': 'wolfpsort_prediction'}, inplace=True)
    df.drop_duplicates().to_csv(f"{output}", sep='\t', index=False)


if __name__ == "__main__":
    build_output_table();
