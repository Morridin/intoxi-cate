# Intoxi-Cate

Intoxi-Cate is an improved version of the DeToX pipeline [[1](#cite1)] by Ringeval et al.

The software provides an end-to-end pipeline for detecting (de novo) toxins in protein transcripts.
It also features some preliminary steps to build up protein transcripts from assembled transcriptomes or even raw
sequence reads.

As DeTox [[1](#cite1)], it includes detection by sequence similarity as well as by structural features.

While keeping the general structure of the pipeline as is, we replaced several tools the pipeline uses with freely
available open-source software. This not only improves its usability as Intoxi-Cate can be installed with just a little
more than a `conda install` command, but also its overall performance, as the replacement software generates equal or
better results at the same or less runtime.

Also, thanks to the translation into actual Python and modularisation into task-specific blocks it is now
easier to maintain or further improve the software by replacing additional components.

Intoxi-Cate is currently only tested on Linux, thus we cannot give any guarantees that it will run correctly on other
platforms.

Intoxi-Cate uses MMSeqs2 [[2](#cite2)] to replace both BLASTp and BLASTn used in this pipeline.
To further reduce the number of dependencies, it also replaces orfipy [[3](#cite3)] in the pipeline section responsible
for translating an assembled transcriptome into peptides.

To replace SignalP and Phobius which are both access-restricted to scientific use [[4](#cite4), [5](#cite5)], we use
TMbed [[6](#cite6)], which is available under Apache 2.0 license.

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
- Instead of executing step 5 as in the guide, you can instead run the following command, if your shell is bash:
  ```bash
  source ~/.bashrc
  ```
- If you don't like to always have a conda environment activated, we suggest you run
  ```bash
  conda config --set auto_activate_base false
  ```

#### 2.3 Create conda environment

As final step installing Intoxi-Cate, you need to create a conda-environment.

In order to create the environment, please run

```bash
conda env create -f intoxi-cate.yaml
```

This will likely take some time (expect 5 - 15 minutes).
If you want to speed this up, if you have mamba available on your machine (included in miniforge) by running

```bash
mamba --version
```

As previously when checking for the conda availability, you have mamba available if you see a version number.

If mamba is available, please replace `conda` in the command above with `mamba` to speed up the proccess.

Now, you're ready to use the pipeline!

## Usage

### Activate Environment and Prepare Configuration

1. Activate the `intoxi-cate` environment.
   ```bash
   conda activate intoxi-cate
   ```
2. Copy `config.yaml` to your working directory and populate it with required information.

### Run Pipeline

The pipeline can be run with

```
intoxi-cate
```

You can specify a `--config` option pointing to an alternative location of a config file other than `config.yaml` in the
current directory.
Also, you can specify more or less any of the options from within the config file as parameters.

For a complete overview of all available options, run

```
intoxi-cate --help
```

or

```
intoxi-cate -h
```

### Configuration of settings and options

The pipeline expects the config file to be of the format `key: value` and to contain the following keys:

#### Required keys

The following table lists all configuration keys that necessary for the pipeline to run.
Out of the group `Rx`, `transcriptome` and `proteome_fasta`, only one key is required.

| Config Key       | Description                                                                                                                                                                                                                                                                                   | 
|:-----------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `proteome_fasta` | If set, neither `transcriptome` or `R1` may be set. If you have your data already as proteome available, set this value to the path to your proteome file in FASTA format.<br/>As Salmon needs raw reads, you can't have this set to a value and `quant` to true at the same time.            | 
| `R1`             | If you do not provide a value for either `transcriptome` or `proteome`, this value **must** be set with the path to a file containing raw reads in FASTQ format. <br/>If only `R1` is set, the reads are assumed to be single-end, if `R2` is set as well, they are assuemd to be paired-end. |       
| `R2`             | If you have paired-end raw read data, set this to the path to the reverse paired-end file in FASTQ format.                                                                                                                                                                                    |        
| `transcriptome`  | If you have the transcriptome already assembled in FASTA format, provide this parameter with the path to that file. <br/>If this parameter and `proteome` are not provided, the assembly will be performed using the `R1` and `R2` parameters.                                                |     
| `toxin_db`       | The path to your toxins database file in FASTA format.                                                                                                                                                                                                                                        |  

The following table lists all values that can optionally be set to adjust the pipeline behaviour.
Some options are only effective if others are set, some exclude other options.

| Config Key             | Description                                                                                                                                                                                                                           |    Default | 
|:-----------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------:|
| `basename`             | Optional. With this string, the names of all files generated by the pipeline will start. <br/>If this value is not set, the file names will be left empty, except for their extensions.                                               |     (none) |
| `contaminants`         | Optional. Setting this value to the path of a contaminants database file in FASTA format activates contaminant filtering within the pipeline.                                                                                         |     (none) |
| `contamination_evalue` | Optional. The e-value threshold for contamination removal using MMSeqs2.                                                                                                                                                              |       1E-5 |
| `clustering_threshold` | Optional. The clustering threshold for MMSeqs2 in the peptide clustering step when Intoxi-Cate is run without `proteome` parameter set.                                                                                               |       0.99 |                                                                                                                                    
| `cys_pattern`          | Optional. Set to `True` if you want Intoxi-Cate to perform cysteine pattern analysis using DeTox's in-house script for that purpose.                                                                                                  |      False |
| `maxlen`               | Optional. The maximum length a section within an NA sequence may have to be considered an ORF, in nucleotides.                                                                                                                        | 30 000 000 |
| `memory`               | Optional. Add a value to limit the maximum amount of RAM in GigaByte available to the pipeline. <br/>If not set, Intoxi-Cate will use as much as it needs up to machine limits.                                                       |     (none) |
| `minlen`               | Optional. The minimum length a section within an NA sequence needs to have to be considered an ORF, in nucleotides.                                                                                                                   |         99 |
| `output_dir`           | Optional. Add a value to specify a dedicated directory for the output files. The directory does not need to exist in beforehand.<br/>If not set, all files the pipeline produces will be put into the current directory.              |        `.` |
| `pfam_db_path`         | Optional. The path to your Pfam database file in .hmm format. If not provided, Intoxi-Cate downloads the database.                                                                                                                    |     (none) |
| `quant`                | Optional. If set (to `True`), Intoxi-Cate performs transcript quantification using Salmon. <br/>This option requires the `R1` parameter to have a value and, depending on the reads data, the `R2` parameter to have a value as well. |      False |
| `swissprot`            | Optional. Set to `True` to perform functional annotation using SwissProt.                                                                                                                                                             |      False |
| `swissprot_evalue`     | Optional. The e-value threshold for SwissProt annotation using MMSeqs2. Only needed at all if `swissprot` is set to `True`.                                                                                                           |      1E-10 |
| `swissprot_db_path`    | Optional. The path to your SwissProt database file in either (gzip-compressed) FASTA format or in the MMSeqs database format. <br/>If no path is provided, MMSeqs will download the database when running a search against SwissProt. |     (none) |
| `threads`              | Optional. Add a value to limit how many CPU threads the pipeline may use. Should not be more than the maximum number of threads avalailable to your machine. <br/>If not set, the pipeline will use all available threads.            |     (none) |
| `toxins_evalue`        | Optional. The e-value threshold for toxin annotation using MMSeqs2.                                                                                                                                                                   |      1E-10 |
| `TPMthreshold`         | Optional. The TPM threshold for a sequence to be flagged "T".                                                                                                                                                                         |       1000 |
| `wolfpsort`            | Optional. Set to `True` if you want Intoxi-Cate to perform subcellular localization prediction as per WoLF PSORT with TMbed.                                                                                                          |      False |

## References

<p id="cite1">[1] RINGEVAL, A., S. FARHAT, A. FEDOSOV, M. GERDOL, S. GRECO, L. MARY, M. V. MODICA und N.
PUILLANDRE, 2024. DeTox: a pipeline for the detection of toxins in venomous organisms. In: Briefings in Bioinformatics.
2024, vol. 25, no. 2, bbae094. <a href="https://doi.org/10.1093/bib/bbae094">DOI 10.1093/bib/bbae094</a>. PMID: 38493344.</p>
<p id="cite2">[2] STEINEGGER, Martin and Johannes SÖDING, 2017. MMseqs2 enables sensitive protein sequence
searching for the analysis of massive data sets. In: Nature Biotechnology. 2017, vol. 35, no. 11, pp. 1026–1028.
<a href="https://doi.org/10.1038/nbt.3988">DOI 10.1038/nbt.3988</a>. ISSN 1546-1696.</p>
<p id="cite3">[3] SINGH, Urminder and Eve Syrkin WURTELE, 2021. orfipy: a fast and flexible tool for extracting
ORFs. In: Bioinformatics. 2021, vol. 37, no. 18, pp. 3019–3020.
<a href="https://doi.org/10.1093/bioinformatics/btab090">DOI 10.1093/bioinformatics/btab090</a>.</p>
<p id="cite4">[4] ALMAGRO ARMENTEROS, José Juan, Konstantinos D. TSIRIGOS, Casper Kaae SØNDERBY, Thomas Nordahl PETERSEN, 
Ole WINTHER, Søren BRUNAK, Gunnar von HEIJNE and Henrik NIELSEN, 2019. SignalP 5.0 improves signal peptide predictions 
using deep neural networks. In: <em>Nature Biotechnology</em>. 2019, vol. 37, pp. 420–423. 
<a href="https://doi.org/10.1038/s41587-019-0036-z">DOI 10.1038/s41587-019-0036-z</a>.</p>
<p id="cite5">[5] KÄLL, Lukas, Anders KROGH und Erik L.L. SONNHAMMER, 2004. A Combined Transmembrane Topology and Signal 
Peptide Prediction Method. In: <em>Journal of Molecular Biology</em>. 2004, vol. 338, no. 5, pp. 1027–1036. 
<a href="https://doi.org/10.1016/j.jmb.2004.03.016">DOI 10.1016/j.jmb.2004.03.016</a>.</p>
<p id="cite6">[6] BERNHOFER, Michael and Burkhard ROST, 2022. TMbed: transmembrane proteins predicted through language 
model embeddings. In: BMC Bioinformatics. 2022, vol. 23, no. 1, p. 326. DOI 
<a href="https://doi.org/10.1186/s12859-022-04873-x">DOI 10.1186/s12859-022-04873-x</a>.</p>
