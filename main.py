import itertools
import subprocess
from pathlib import Path

import pandas as pd

from lib.utils import global_output
from lib import config, utils, assemble_transcriptome, cluster_peptides, hmmer, blast_on_toxins, blast_on_uniprot, \
    signalp, retrieve_candidate_toxins


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
    signalp_result = signalp()
    hmm_search = hmmer()
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
    # --- Validate Config ----------------------------------------------------------
    has_reads = config.get("R1") is not None
    has_transcripts = config.get("transcriptome") is not None
    has_proteome = config.get("proteome_fasta") is not None
    quant = config.get("quant", True)

    # --- 1. Proteome OR 2. Reads (+/- transcriptome)  -----------------------------
    if has_proteome and (has_reads or has_transcripts):
        raise ValueError("Choose between providing only a proteome or a transcriptome/read")
    if not (has_reads or has_transcripts or has_proteome):
        raise ValueError(
            "No entry point provided in config.yaml. Please provide either a proteome, a transcriptome or reads.")

    if has_proteome:
        if quant:
            raise ValueError("quant must be false in config.yaml if you provide proteome_fasta.")
        if config.get("R1") is not None or config.get("R2") is not None:
            raise ValueError("You must not fill in R1/R2 in config.yaml when proteome_fasta is used.")

    transcriptome = assemble_transcriptome()
    clustered_peptides = cluster_peptides(transcriptome)

    toxins_blast_result = blast_on_toxins(clustered_peptides)

    signalp_result = signalp(clustered_peptides)

    toxin_candidates = retrieve_candidate_toxins(clustered_peptides, toxins_blast_result, signalp_result)

    hmmer_result = hmmer(toxin_candidates)
    uniprot_blast_result = blast_on_uniprot(toxin_candidates)

    build_output_table()
