"""
This folder shall later contain everything that is not part of the following:
- User Interface (whichever will be implemented)
- Main program logic (what you execute when you run Intoxi-Cate).
- Data

This leaves us with the following:
- Miscellaneous utils and stuff
- Algorithms the pipeline needs
- All the data handling routines
"""
from config import config
from assemble_transcriptome import *
from generate_peptides import *
global snakemake_checkpoints

# Set parts of lib
__all__ = ["config", "utils"]

# Add the different modules
__all__ += ["assemble_transcriptome", "run_salmon", "cluster_peptides"]

# Temporary stuff I later want to get rid of
__all__ += ["snakemake_checkpoints"]