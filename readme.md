# Intoxi-Cate
[![Conda Install](https://github.com/Morridin/intoxi-cate/actions/workflows/conda-install.yml/badge.svg)](https://github.com/Morridin/intoxi-cate/actions/workflows/conda-install.yml)

Intoxi-Cate is an improved version of the DeToX pipeline [[1](#cite1)] by Ringeval et al.

The software provides an end-to-end pipeline for detecting (de novo) toxins in protein transcripts.
It also features some preliminary steps to build up protein transcripts from assembled transcriptomes or even raw
sequence reads.

As DeTox [[1](#cite1)], it includes detection by sequence similarity as well as by structural features.

While keeping the general structure of the pipeline as is, we replaced several tools the pipeline uses with freely
available open-source software. This not only improves its usability as Intoxi-Cate can be installed with just a little
more than a `conda install` command, but also its overall performance, as the replacement software generates equal or
better results at the same or less runtime.

Also, thanks to the translation into actual Python and modularisation into task-specific blocks it is now easier to
maintain or further improve the software by replacing additional components.

Intoxi-Cate is currently only tested on Linux, thus we cannot give any guarantees that it will run correctly on other
platforms.

Intoxi-Cate uses MMSeqs2 [[2](#cite2)] to replace both BLASTp and BLASTn used in this pipeline.
To further reduce the number of dependencies, it also replaces orfipy in the pipeline section responsible for
translating an assembled transcriptome into peptides.

To replace SignalP and Phobius which are both access-restricted to scientific use, we use TMbed [[3](#cite3)], which is
available under Apache 2.0 license.

## Installation

The installation is rather simple in the end:
All you need to do is downloading the software, create a new environment based on the environment file included in the
repository and maybe check that everything is working as expected.

### 1. Clone Repository

Clone the repository to your local machine:

  ```bash
  git clone https://github.com/Morridin/intoxi-cate.git
  ```

### 2. Create and Activate Conda Environment

First, you need conda in order to get all the dependencies Intoxi-Cate needs to work installed along with the pipeline.
In most cases, Miniforge or Miniconda should be sufficient.

#### 2.1 Check that conda is available

To check whether conda is already installed on your system, enter

```bash
conda -V
```

into your terminal.
If you receive an error message, conda is likely to be not installed on your system.
Then, proceed with step [2.2.](#22-install-conda)

If you see something like `conda 26.1.0`, conda is installed on your system and ready for use.
Please, proceed with step [2.3](#23-create-conda-environment)

#### 2.2 Install conda

To install either of the aforementioned options, please visit the conda website and follow along the instructions in
the [installation guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html).
The link points the Linux version of the guide as Linux is the only supported platform.

Some notes on that guide:

- In case of miniforge, you can't verify the installer checksum as miniforge doesn't offer a fingerprint to compare to.
- Instead of executing step 5 as in the guide, you can run the following command, if your shell is bash:
  ```bash
  source ~/.bashrc
  ```
- If you don't like to always have a conda environment activated, you can run
  ```bash
  conda config --set auto_activate_base false
  ```
  to deactivate this behaviour.

#### 2.3 Create conda environment

As final step installing Intoxi-Cate, you need to create a conda-environment.

In order to create the environment, please first switch to the directory into which Intoxi-Cate was downloaded.
Git will most likely download Intoxi-Cate into a directory `intoxi-cate` in the current directory, and report this as very first line on the screen.

Once you are in the correct directory, please run
```bash
conda env create -f intoxi-cate.yaml
```

This will likely take some time (expect 5 - 15 minutes).
If you want to speed this up, check if you have mamba available on your machine (included in miniforge) by running
```bash
mamba --version
```
As previously when checking for the conda availability, you have mamba available if you see a version number.

If mamba is available, you can replace `conda` in the environment creation command above with `mamba`.

Now, once the command has successfully finished, you're ready to use Intoxi-Cate!

## Usage

### Activate Environment and Prepare Configuration

1. Activate the `intoxi-cate` environment.
   ```bash
   conda activate intoxi-cate
   ```
2. Copy `config.yaml` to your working directory and populate it with the required information.

### Run Pipeline

The pipeline can be run with

```
intoxi-cate
```

You can specify a `--config` option pointing to an alternative location of a config file other than `config.yaml` in the
current directory.
Also, you can specify all the options from within the config file as parameters.
Just add two dashes in front of the config key.

For a complete overview of all available options, run

```
intoxi-cate --help
```

or

```
intoxi-cate -h
```

### Configuration of settings and options

The pipeline expects the config file to be of the format `key: value` and to contain the following keys, as listed in the tables below.

#### Required keys

The following table lists all configuration keys that necessary for the pipeline to run.
Out of the group `R1`/`R2`, `transcriptome` and `proteome_fasta`, only one key is required and `proteome_fasta` excludes the other two.

| Config Key       | Description                                                                                                                                                                                                                                                                                   | 
|:-----------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `proteome_fasta` | If set, neither `transcriptome` or `R1` may be set. If you have your data already as proteome available, set this value to the path to your proteome file in FASTA format.<br/>As Salmon needs raw reads, you can't have this set to a value and `quant` to true at the same time.            | 
| `R1`             | If you do not provide a value for either `transcriptome` or `proteome`, this value **must** be set with the path to a file containing raw reads in FASTQ format. <br/>If only `R1` is set, the reads are assumed to be single-end, if `R2` is set as well, they are assumed to be paired-end. |       
| `R2`             | If you have paired-end raw read data, set this to the path to the reverse paired-end file in FASTQ format.                                                                                                                                                                                    |        
| `transcriptome`  | If you have the transcriptome already assembled in FASTA format, provide this parameter with the path to that file. <br/>If this parameter and `proteome` are not provided, the assembly will be performed using the `R1` and `R2` parameters.                                                |     
| `toxin_db`       | The path to your toxins database file in FASTA format.                                                                                                                                                                                                                                        |  

The following table lists all values that can optionally be set to adjust the pipeline behaviour.
Some options are only effective if others are set, some exclude other options.

| Config Key                 | Description                                                                                                                                                                                                                                                             | Default                                                        | 
|:---------------------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:---------------------------------------------------------------|
| `basename`                 | With this string, the names of all files generated by the pipeline will start. <br/>If this value is not set, the file names will be left empty, except for their extensions.                                                                                           | `""` (empty string)                                            |
| `contaminants`             | Setting this value to the path of a contaminants database file in FASTA format activates contaminant filtering within the pipeline.                                                                                                                                     | (none)                                                         |
| `contamination_evalue`     | The e-value threshold for contamination removal using MMSeqs2.                                                                                                                                                                                                          | 1E-5                                                           |
| `clustering_threshold`     | The clustering threshold for MMSeqs2 in the peptide clustering step when Intoxi-Cate is run without `proteome` parameter set.                                                                                                                                           | 0.99                                                           |                                                                                                                                    
| `cys_pattern`              | Set to `True` if you want Intoxi-Cate to perform cysteine pattern analysis using DeTox's in-house script for that purpose.                                                                                                                                              | False                                                          |
| `maxlen`                   | The maximum length a section within an NA sequence may have to be considered an ORF, in nucleotides.                                                                                                                                                                    | 45 000                                                         |
| `memory`                   | Add a value to limit the maximum amount of RAM in GigaByte available to the pipeline. <br/>If not set, Intoxi-Cate will use as much as it needs up to machine limits.<br/>This value is not sanitised!                                                                  | (none)                                                         |
| `minlen`                   | The minimum length a section within an NA sequence needs to have to be considered an ORF, in nucleotides.                                                                                                                                                               | 99                                                             |
| `mmseqs_sequence_coverage` | The threshold of sequence coverage above which two sequences are clustered, given sufficient sequence identity. <br/>A sequence contained to at least this fraction in another is clusterd to that one.<br/>To get the closest results to DeTox, use the default value. | 0.0                                                            |
| `mmseqs_path`              | If you prefer to use an already existing MMSeqs2 installation on your system <br/>instead of the one that's installed by conda together with Intoxi-Cate, fill this setting with the path to the MMSeqs binary.                                                         | (none)                                                         | 
| `output_dir`               | Add a value to specify a dedicated directory for the output files. The directory does not need to exist in beforehand.<br/>If not set, all files the pipeline produces will be put into the current directory.                                                          | `.`                                                            |
| `pfam_db_path`             | The path to your Pfam database file in .hmm format. If not provided, Intoxi-Cate downloads the database.                                                                                                                                                                | (none)                                                         |
| `quant`                    | If set (to `True`), Intoxi-Cate performs transcript quantification using Salmon. <br/>This option requires the `R1` parameter to have a value and, depending on the reads data, the `R2` parameter to have a value as well.                                             | False                                                          |
| `repeated_aa_threshold`    | Set this value to change the lower threshold for how many repetitions of 1, 2 or 3 amino acid long frames are flagged as problematic repetitions.                                                                                                                       | 5                                                              |
| `signalpeptide_minlen`     | Sequences with a signal peptide below this length are considered false positives in terms of having a signal peptide and thus discarded.                                                                                                                                |  10                                                            |
| `swissprot`                | Set to `True` to perform functional annotation using SwissProt.                                                                                                                                                                                                         | False                                                          |
| `swissprot_evalue`         | The e-value threshold for SwissProt annotation using MMSeqs2. Only needed at all if `swissprot` is set to `True`.                                                                                                                                                       | 1E-10                                                          |
| `swissprot_db_path`        | The path to your SwissProt database file in either (gzip-compressed) FASTA format or in the MMSeqs database format. <br/>If no path is provided, MMSeqs will download the database when running a search against SwissProt.                                             | (none)                                                         |
| `threads`                  | Add a value to limit how many CPU threads the pipeline may use. Should not be more than the maximum number of threads available to your machine. <br/>If not set, the pipeline will use all available threads.                                                          | (none)                                                         |
| `tmbed_model_path`         | Add a value to specify a directory where TMbed shall store its ProtT5 model persistently. <br/>If no value is set, TMbed will download the model each time it is executed.                                                                                              | (none)                                                         |
| `tmbed_use_cpu`            | If you want Intoxi-Cate to fail instead of falling back on CPU when running TMbed, set this to `False`.                                                                                                                                                                 | True                                                           | 
| `tmbed_use_gpu`            | Set to `False` if your system is not equipped with a GPU capable of running TMbed (it should have at least 12 GB VRAM [[3](#cite3)]).                                                                                                                                   | True                                                           |
| `toxins_evalue`            | The e-value threshold for toxin annotation using MMSeqs2.                                                                                                                                                                                                               | 1E-10                                                          |
| `TPMthreshold`             | The TPM threshold for a sequence to be flagged "T".                                                                                                                                                                                                                     | 1000                                                           |
| `wolfpsort`                | Set to `True` if you want Intoxi-Cate to perform subcellular localization prediction as per WoLF PSORT with TMbed.                                                                                                                                                      | False                                                          |
| `wolfPsort_path`           | If you want to use an already existing installation of WoLF PSORT in your program (or an entirely different program with similar output), plase add the path to and including the binary.                                                                               | `$CONDA_PREFIX/share/WoLFPSort-master/bin/runWolfPsortSummary` |

## A note on WoLF PSORT
WoLF PSORT is licensed under GPL v3 which conflicts with Intoxi-Cate's license. 
Hence, WoLF PSORT is not distributed as part of Intoxi-Cate but only downloaded on the user's behalf when the user sets the `wolfpsort` option to `True` while not providing a path to an existing WoLF PSORT binary.
This automatic download is only performed for the user's convenience.

I am more than happy to replace this rather skewed mechanic with a conda dependency if WoLF PSORT is added to conda in the future.
Alternatively, I'm also happy to replace it with another program that performs the same job at least equally well and is either available on PyPI, conda or under a license compatible to Apache 2.0.

## References

<p id="cite1">[1] RINGEVAL, A., S. FARHAT, A. FEDOSOV, M. GERDOL, S. GRECO, L. MARY, M. V. MODICA and N.
PUILLANDRE, 2024. DeTox: a pipeline for the detection of toxins in venomous organisms. In: Briefings in Bioinformatics.
2024, vol. 25, no. 2, bbae094. <a href="https://doi.org/10.1093/bib/bbae094">DOI 10.1093/bib/bbae094</a>. PMID: 38493344.</p>
<p id="cite2">[2] STEINEGGER, Martin and Johannes SÖDING, 2017. MMseqs2 enables sensitive protein sequence
searching for the analysis of massive data sets. In: Nature Biotechnology. 2017, vol. 35, no. 11, pp. 1026–1028.
<a href="https://doi.org/10.1038/nbt.3988">DOI 10.1038/nbt.3988</a>. ISSN 1546-1696.</p>
<p id="cite3">[3] BERNHOFER, Michael and Burkhard ROST, 2022. TMbed: transmembrane proteins predicted through language 
model embeddings. In: BMC Bioinformatics. 2022, vol. 23, no. 1, p. 326. DOI 
<a href="https://doi.org/10.1186/s12859-022-04873-x">DOI 10.1186/s12859-022-04873-x</a>.</p>




