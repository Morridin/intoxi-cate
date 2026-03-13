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
This module is the CLI runner for Intoxi-Cate.
"""
import argparse
import io
import shutil
import stat
import tempfile
import zipfile
from pathlib import Path
import sys

import requests

from .lib import assemble_transcriptome, cluster_peptides, blast, detect_by_structure, \
    retrieve_candidate_toxins, hmmer, run_salmon, build_output_table
from .lib._config import Config, _WOLF_PSORT_GITHUB_URL, _WOLF_PSORT_INSTALL_DIR


def accept_bool(arg: str) -> bool:
    if arg.lower() == "true":
        return True
    elif arg.lower() == "false":
        return False
    else:
        raise ValueError(f"Expected a boolean value, found: {arg}")


def validate_file(file_path: Path | None, message_detail: str) -> None:
    """
    Subtask of `validate_args`. Checks if a given path is None or exists and terminates the program otherwise.
    :param file_path: The path to validate.
    :param message_detail: The name of the file to be displayed in the possible error message when exiting.
    :return: Nothing.
    """
    if file_path is not None and not file_path.is_file():
        print(f"Could not find {message_detail} file. Please ensure you provided a correct path.", file=sys.stderr)
        exit(2)


def validate_args(config: Config) -> None:
    """
    Checks that the arguments handed over are valid for a pipeline run.
    If they are invalid, this function prints an error message and terminates the program with exit code 2.
    :param config: A pre-assembled config to validate.
    :return: Nothing.
    """
    has_r1 = config.get("R1") is not None
    has_r2 = config.get("R2") is not None
    has_reads = has_r1 or has_r2
    has_transcripts = config.get("transcriptome") is not None
    has_proteome = config.get("proteome_fasta") is not None

    # Sanitise input files
    if has_proteome and (has_reads or has_transcripts):
        print("Choose between providing only a proteome or a transcriptome/reads", file=sys.stderr)
        exit(2)
    if not (has_reads or has_transcripts or has_proteome):
        print("No input data specified. Please provide either a proteome, a transcriptome or reads.", file=sys.stderr)
        exit(2)
    if has_reads and not has_r1:
        print("Your reads data is incomplete. If you have only single-end reads, provide `R1` only.", file=sys.stderr)
        exit(2)

    toxins_db = config.get_path("toxin_db")
    if toxins_db is None:
        print("The toxins database file is required. Please provide one.", file=sys.stderr)
        exit(2)
    validate_file(toxins_db, f"toxins database")

    validate_file(config.get_path("contaminants"), "contaminants database")

    mmseqs_path = config.get_path("mmseqs_path")
    if mmseqs_path is not None and mmseqs_path.name != "mmseqs":
        mmseqs_path /= "mmseqs"
    validate_file(mmseqs_path, "MMSeqs binary")

    if config.get_path("output_dir") == Path("/dev/null"):
        print("The output directory is specified as /dev/null. This run is pointless. Anyways ...", file=sys.stderr)

    if config.get("quant") and not has_reads:
        print(
            "`quant` requires reads to be available, so `R1` and, depending on your data, `R2` must be provided as well.",
            file=sys.stderr)
        exit(2)

    if config.get("swissprot"):
        validate_file(config.get_path("swissprot_db_path"), "Swiss-Prot database")

    if not config.get("tmbed_use_cpu") and not config.get("tmbed_use_gpu"):
        print("You do not allow TMbed to use either CPU or GPU. Hence, this pipeline cannot run successfully.",
              file=sys.stderr)
        exit(2)

    tpm_threshold = config.get("TPMthreshold")
    if tpm_threshold is not None:
        try:
            config.config["TPMthreshold"] = float(tpm_threshold)
        except ValueError:
            print("TPMthreshold must be a real number!", file=sys.stderr)
            exit(2)


def handle_wolf_psort(wolf_psort_path: Path = None) -> Path:
    """
    Validates the WoLF PSORT binary path and installs the program if the path is the default value.
    If installation or validation fail, the function will terminate the program.
    :param wolf_psort_path: The expected path to the WoLF PSORT binary.
    :return: The path to the actual binary.
    """
    # Take care of WoLF PSORT installation if not yet installed
    if wolf_psort_path is None:
        wolf_psort_path = Config._DEFAULT["wolfPsort_path"]

    if wolf_psort_path == Config._DEFAULT["wolfPsort_path"]:
        if not wolf_psort_path.is_file():
            print("Downloading and installing WoLF PSORT ...")
            wolf_psort_path.unlink(True) # Delete anything that was there before, especially if it was a directory.

            response = requests.get(_WOLF_PSORT_GITHUB_URL)
            with tempfile.TemporaryDirectory() as t:
                with zipfile.ZipFile(io.BytesIO(response.content)) as zip_ref:
                    zip_ref.extractall(t)
                source = Path(t) / "WoLFPSort-master"
                shutil.copytree(source, _WOLF_PSORT_INSTALL_DIR, dirs_exist_ok=True)
            print(f"Done!\nYou can now find WoLF PSort in {_WOLF_PSORT_INSTALL_DIR}.\n")
            # Make it executable...
            previous_mode = wolf_psort_path.stat().st_mode
            wolf_psort_path.chmod(previous_mode | stat.S_IXUSR)

    if not wolf_psort_path.is_file():
        print("Path does not lead to a valid WoLF PSORT binary. Please check the path you entered.", sys.stderr)
        exit(2)

    return wolf_psort_path


def parse_args() -> Config:
    """
    Parses the command line arguments and validates them.
    """
    parser = argparse.ArgumentParser(
        description="A pythonic pipeline for annotation of de novo toxins.\n \n"
                    "For more information, especially on default values, please refer to the ReadMe file. Setting the "
                    "--config option will overwrite the expected location of the config file. Setting other options will "
                    "overwrite the corresponding config keys.\n \n"
                    "In case of conflicting values between values in a config file and options set here, the program will exit with code 1.",
        usage="%(prog)s [options]",
        add_help=False
    )
    general_group = parser.add_argument_group(
        "General options"
    )
    general_group.add_argument(
        "-h", "--help",
        action="help",
        help="Show this help message and exit",
    )
    general_group.add_argument(
        "-V", "--version",
        action="version",
        version="%(prog)s 0.1.2"
    )
    general_group.add_argument(
        "-c", "--config",
        default="config.yaml",
        help="Path to configuration file in YAML format. Default value is `config.yaml`. If no config file can be found, the default config is used as base.",
        metavar="PATH",
        type=Path,
    )

    config_group = parser.add_argument_group(
        "Configuration settings",
        "The following options overwrite values set in any config file set by the --config option (or its default value).",
    )
    config_group.add_argument(
        "--proteome_fasta",
        help="Path to a FASTA file containing a ready-to-use proteome.",
        metavar="PATH",
        type=Path
    )
    config_group.add_argument(
        "--r1",
        dest="R1",
        help="Path to a file containing raw reads from an Illumina device, forward read file.",
        metavar="PATH",
        type=Path,
    )
    config_group.add_argument(
        "--r2",
        dest="R2",
        help="Path to a file containing raw reads from an Illumina device, reverse read file.",
        metavar="PATH",
        type=Path,
    )
    config_group.add_argument(
        "--transcriptome",
        help="Path to a FASTA file containing an assembled transcriptome.",
        metavar="PATH",
        type=Path
    )
    config_group.add_argument(
        "--toxin_db",
        help="Path to a FASTA file containing a toxin database.",
        metavar="PATH",
        type=Path
    )
    config_group.add_argument(
        "--basename",
        help="The prefix that is prepended to all persistent output files.",
        metavar="STR",
        type=str
    )
    config_group.add_argument(
        "--contaminants",
        help="Path to a FASTA file containing sequences of known contaminants, to be used as contaminants database.",
        metavar="PATH",
        type=Path
    )
    config_group.add_argument(
        "--contamination_evalue",
        help="The e-value threshold for contaminants removal using MMSeqs2.",
        metavar="FLOAT",
        type=float
    )
    config_group.add_argument(
        "--clustering_threshold",
        help="The peptide clustering threshold for MMSeqs2.",
        metavar="FLOAT",
        type=float
    )
    config_group.add_argument(
        "--cys_pattern",
        help="Analyse output sequences for cysteine patters.",
        metavar="BOOL",
        type=accept_bool
    )
    config_group.add_argument(
        "--maxlen",
        help="Maximum length an ORF may have to be detected.",
        metavar="INT",
        type=int
    )
    config_group.add_argument(
        "--memory",
        help="The maximum amount of memory for Intoxi-Cate to use, in GB. Note: This limit will only be applied to certain steps within the pipeline.",
        metavar="INT",
        type=int
    )
    config_group.add_argument(
        "--minlen",
        help="Minimum length an ORF may have to be detected.",
        metavar="INT",
        type=int
    )
    config_group.add_argument(
        "--mmseqs_path",
        help="The path to an alternative MMSeqs installation",
        metavar="PATH",
        type=Path
    )
    config_group.add_argument(
        "--output_dir",
        help="The path to the directory where all persistent outputs will be stored.",
        metavar="PATH",
        type=Path
    )
    config_group.add_argument(
        "--pfam_db_path",
        help="The path to an already downloaded Pfam database file.",
        metavar="PATH",
        type=Path
    )
    config_group.add_argument(
        "--quant",
        help="Whether to quantise the raw reads data.",
        metavar="BOOL",
        type=accept_bool
    )
    config_group.add_argument(
        "--repeated_aa_threshold",
        help="The number of amino acid repetitions required for a sequence to be flagged.",
        metavar="INT",
        type=int
    )
    config_group.add_argument(
        "--swissprot",
        help="Whether to perform an alignment search against Swiss-Prot.",
        metavar="BOOL",
        type=accept_bool
    )
    config_group.add_argument(
        "--swissprot_evalue",
        help="The e-value for the alignment search against Swiss-Prot.",
        metavar="FLOAT",
        type=float
    )
    config_group.add_argument(
        "--swissprot_db_path",
        help="The path to an already downloaded Swiss-Prot database file (FASTA or MMSeqs format).",
        metavar="PATH",
        type=Path
    )
    config_group.add_argument(
        "--threads",
        help="Maximum number of threads to use. This value is not sanity-checked.",
        metavar="INT",
        type=int
    )
    config_group.add_argument(
        "--tmbed_model_path",
        help="The path to a directory where TMbed shall store its embeddings persistently.",
        metavar="PATH",
        type=Path
    )
    config_group.add_argument(
        "--tmbed_use_cpu",
        help="Whether TMbed shall fallback on CPU if no GPU is available.",
        metavar="BOOL",
        type=accept_bool
    )
    config_group.add_argument(
        "--tmbed_use_gpu",
        help="Whether TMbed shall use the GPU if available.",
        metavar="BOOL",
        type=accept_bool
    )
    config_group.add_argument(
        "--toxins_evalue",
        help="The e-value threshold for toxin detection using MMSeqs2.",
        metavar="FLOAT",
        type=float
    )
    config_group.add_argument(
        "--TPMthreshold",
        help="The TPM threshold for a sequence to be flagged.",
        metavar="FLOAT",
        type=float
    )
    config_group.add_argument(
        "--wolfpsort",
        help="Whether to perform sub-cellular localisation of the output sequences with WoLF PSORT.",
        metavar="BOOL",
        type=accept_bool
    )

    args = parser.parse_args()

    config_path = args.config or "config.yaml"

    if config_path.exists():
        config = Config()
        config.set_config_file(config_path)
    else:
        config = Config()

    del args.config

    for key, value in vars(args).items():
        if value is not None:
            config.config[key] = value

    validate_args(config)
    
    if config.get("wolfpsort"):
        config.config["wolfPsort_path"] = handle_wolf_psort(config.get_path("wolfPsort_path"))

    return config

def app(config: Config):
    """
    Runs Intoxi-Cate.
    :param config: The configuration for the current run of Intoxi-Cate.
    """

    transcriptome = assemble_transcriptome()
    clustered_peptides = cluster_peptides(transcriptome)

    toxins_blast_result = blast.on_toxins(clustered_peptides)

    signal_peptides = detect_by_structure(clustered_peptides).set_index("ID")

    toxin_candidates = retrieve_candidate_toxins(clustered_peptides, toxins_blast_result, signal_peptides)

    hmmer_result = hmmer(toxin_candidates)

    if config.get("swissprot", False):
        uniprot_blast_result = blast.on_uniprot(toxin_candidates).reset_index()
    else:
        uniprot_blast_result = None

    if config.get("quant"):
        salmon_result = run_salmon()
    else:
        salmon_result = None

    print(
        f"\n\n"
        f"# ============================================================================ #\n"
        f"#                              Pipeline complete                               #\n"
        f"# ============================================================================ #\n"
        f"The final pipeline output can be found under {
        build_output_table(toxin_candidates, hmmer_result, toxins_blast_result.reset_index(), signal_peptides, uniprot_blast_result, salmon_result)
        }"
    )

def run():
    """
    The CLI entry point for Intoxi-Cate.
    """
    config = parse_args()
    app(config)