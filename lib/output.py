"""
This module bundles all former rules that are only for used for direct pipeline output generation.

Those include:

detect_repeated_aa, run_wolfpsort

and

build_output_table

that relies on the aforementioned two and six other former rules and forms the public API of this module.
"""
import itertools
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

from lib import config, utils

__all__ = ["build_output_table"]


# ============================= Public functions ============================= #
def build_output_table(candidate_toxins: Path, hmmer: pd.DataFrame, toxins_blast_result: pd.DataFrame,
                       signal_p_result: pd.DataFrame, uniprot_blast_result: pd.DataFrame = None,
                       salmon_result: Path = None) -> Path:
    """
    Builds the output table
    :param candidate_toxins: A FASTA file containing toxin candidates.
    :param hmmer: A DataFrame with the results from HMMer.
    :param toxins_blast_result: A DataFrame with the results of BLASTn on the toxins file(s).
    :param uniprot_blast_result: A DataFrame with the results of BLASTn on the uniprot database.
    :param signal_p_result: A DataFrame with the results of running SignalP.
    :param salmon_result: A DataFrame with the results of running Salmon.
    :return: The Path to the TSV file with the combined output of the entire pipeline.
    """
    if config.get("wolfpsort", False):
        wps_path = config.get("wolfPsort_path")
        wolf_p_result = _run_wolfpsort(candidate_toxins, wps_path)
    else:
        wolf_p_result = None

    repetition_threshold = config.get("repeated_aa_threshold", 5)

    repeated_aa = _detect_repeated_aa(candidate_toxins, repetition_threshold)

    final_output = utils.global_output(config.get("basename") + "_toxins.tsv")

    tpm_threshold = config.get("TPMthreshold", 1000)

    cys_pattern = config.get("cys_pattern", False)

    _build_output_table(final_output, hmmer, toxins_blast_result, repeated_aa, candidate_toxins, signal_p_result,
                        uniprot_blast_result, salmon_result, wolf_p_result, tpm_threshold=tpm_threshold,
                        cys_pattern=cys_pattern)

    return final_output


# ============================ Private functions ============================= #
def _run_wolfpsort(candidate_toxins: Path, wolf_p_sort_path: Path) -> pd.DataFrame:
    """
    Runs wolfpsort on secreted peptides inferred by SignalP
    :param candidate_toxins: FASTA file
    :return: DataFrame with results.
    """
    awk = "awk '{print $1\"\t\"$2}'"

    with tempfile.NamedTemporaryFile(suffix=".tsv", delete_on_close=False) as output:
        subprocess.run(
            f"{wolf_p_sort_path} animal < {candidate_toxins} | {awk} > {output}",
        )
        return pd.read_csv(output, sep="\t", index_col=0)


def _detect_repeated_aa(candidate_toxins: Path, threshold: int) -> pd.DataFrame:
    """
    This rule looks at the fasta aminoacid sequences in input and produces a table.
    The table reports whether some kind of repeated pattern is found in the sequences (up to 3AA long).
    The default threshold for repetition is 5.
    The input is processed with biopython
    :param candidate_toxins:
    :param threshold:
    :return:
    """
    secreted = utils.fasta_to_dataframe(f"{candidate_toxins}")
    secreted["Repeats1"] = secreted.apply(lambda x: _find_repetition(1, x["Sequence"], threshold), axis=1)
    secreted["Repeats2"] = secreted.apply(lambda x: _find_repetition(2, x["Sequence"], threshold), axis=1)
    secreted["Repeats3"] = secreted.apply(lambda x: _find_repetition(3, x["Sequence"], threshold), axis=1)
    secreted["Repeats"] = secreted["Repeats1"] + secreted["Repeats2"] + secreted["Repeats3"]
    secreted['RepeatsTypes'] = secreted['Repeats'].apply(lambda t: [n for (n, _) in t])
    secreted['RepeatsLengths'] = secreted['Repeats'].apply(lambda t: [n for (_, n) in t])
    secreted['RepeatsLengths'] = [','.join(map(str, l)) for l in secreted['RepeatsLengths']]
    secreted['RepeatsTypes'] = [','.join(map(str, l)) for l in secreted['RepeatsTypes']]
    secreted = secreted.drop(columns=["Repeats", "Repeats1", "Repeats2", "Repeats3"])
    return secreted


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


def _build_output_table(output_file: Path, hmmer: pd.DataFrame, toxins_blast_result: pd.DataFrame,
                        repeated_aa: pd.DataFrame, candidate_toxins: Path, signal_p_result: pd.DataFrame,
                        uniprot_blast_result: pd.DataFrame = None, salmon_result: Path = None,
                        wolf_p_sort: pd.DataFrame = None, *, tpm_threshold: int,
                        cys_pattern: bool):
    """
    this rule merges the tabular output of the other rules and merges it in a single table.
    It uses the outputs list defined above.
    :param cys_pattern:
    :param tpm_threshold: Threshold for tpm
    :param output_file: A TSV file summing up the entire pipeline result.
    """

    extra = [x for x in (wolf_p_sort, hmmer, toxins_blast_result, uniprot_blast_result, repeated_aa) if x is not None]

    seq_df = utils.fasta_to_dataframe(candidate_toxins)
    df = seq_df.merge(signal_p_result, on="ID", how="left")

    for dfi in extra:
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

    if salmon_result is not None:
        df["contig"] = df['ID'].apply(lambda x: x.split("_ORF")[0])
        q = pd.read_csv(salmon_result, sep="\t")
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
    df.drop_duplicates().to_csv(f"{output_file}", sep='\t', index=False)
