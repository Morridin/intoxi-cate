"""
This module encapsulates the two rules trim_reads <- assemble_transcriptome

The public API exposes the results of assemble_transcriptome in such a way that it is executed only once per Pipeline run.

Also, the API exposes run_salmon as it needs the outputs of trim_reads which is not any longer publicly available to code.

Currently, four methods are known to rely on this output.
"""
import subprocess
from functools import cache
from pathlib import Path

from lib import utils, config

__all__ = ["assemble_transcriptome", "run_salmon"]

@cache
def _trim_reads(r1: Path, r2: Path | None, threads: int) -> dict[str, Path]:
    """
    This function trims the raw reads provided in paths r1 and r2 using what is provided in adapters.
    :param threads: Maximum number of threads to use
    :param r1: The path to the forward paired-end or single-end reads file in FASTQ format.
    :param r2: The path to the reverse paired-end reads file in FASTQ format.
    :return: The paths to the forward and reverse paired-end reads in FASTA format, as well as their unpaired
    counterparts, all sorted in a dict with the following keys: r1, (r1_unpaired, r2, r2_unpaired (only if r2 is not None))
    """
    assert r1 != "" and r2 != "", "Params r1 and r2 must not be empty!"

    base_path = utils.global_output("trimmed_reads")
    base_path.mkdir(parents=True, exist_ok=True)

    output = {"r1": base_path / r1.name}

    if r2 is not None:
        output |= {
            "r2": base_path / r2.name,
            "r1_unpaired": base_path / f"unpaired.{r1.name}", # Currently not used
            "r2_unpaired": base_path / f"unpaired.{r2.name}"  # Currently not used
        }
        command_extension = f" -I {r2} -O {output["r2"]} --unpaired1 {output["r1_unpaired"]} --unpaired2 {output["r2_unpaired"]} --detect_adapter_for_pe "
    else:
        command_extension = " "

    subprocess.run(
        f"fastp -i {r1} -o {output["r1"]}{command_extension}"
        f"-w {threads} --cut_front --cut_front_window_size 1 --cut_tail --cut_tail_window_size 1 "
        f"--cut_right --cut_right_window_size 4 --cut_mean_quality 15 --length_required 25"
    )

    return output

def _assemble_transcriptome(r1: Path, r2: Path | None) -> Path:
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

    assembly = utils.global_output("trinity_out_dir/Trinity.fasta")

    subprocess.run(
        f"Trinity --seqType {seq_type} {reads} --CPU {threads} --max_memory {memory} --output {assembly}"
    )

    return assembly

@cache
def assemble_transcriptome() -> Path:
    """
    This function assembles the transcriptome if it is not provided via config.
    Else, it will instead return the file path from the config: If available, the value at "transcriptome" is returned, else "proteome_fasta".
    :param r1: The path to the forward paired-end or single-end reads file in FASTQ format.
    :param r2: The path to the reverse paired-end reads file in FASTQ format.
    :param threads: Maximum number of threads to use, optional. If left as None, the value from the config file is used.
    If that is not possible, the program tries to estimate this on its own.
    :return: The path to a FASTA file containing an assembled transcriptome.
    """
    threads = utils.get_threads()

    transcriptome = config.get("transcriptome")
    if transcriptome is None:
        transcriptome = config.get("proteome_fasta")
    if transcriptome is None:
        reads = _trim_reads(config.get("R1"), config.get("R2"), threads)
        transcriptome = _assemble_transcriptome(reads["r1"], reads.get("r2"))
    return transcriptome


@cache
def run_salmon(r1: Path, r2: Path | None, threads: int = None) -> Path:
    """
    Runs salmon on the reads and assembled transcriptome to quantify <?>
    :return:
    """

    if threads is None:
        threads = utils.get_threads()
    transcriptome = assemble_transcriptome(r1, r2, threads)

    base_dir = utils.global_output(config.get('basename'))

    quant_dir = base_dir.with_name(base_dir.name + "_quant")
    quantification = quant_dir / "quant.sf"

    index = utils.global_output("salmon.idx")

    quant_dir.mkdir(parents=True, exist_ok=True)
    index.mkdir(parents=True, exist_ok=True)

    subprocess.run(
        f"salmon index -t {transcriptome} -i {index} -p {threads}"
    )

    if r2 is not None:
        reads = f"-1 {r1} -2 {r2}"
    else:
        reads = f"-r {r1}"

    subprocess.run(
        f"salmon quant -i {index} -l A {reads} --validateMappings -o {quant_dir}"
    )

    return quantification