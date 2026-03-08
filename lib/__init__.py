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
from .config import config
from .transcriptome import assemble_transcriptome, run_salmon
from .generate_peptides import cluster_peptides
from .toxins_peptides import retrieve_candidate_toxins
from .tmbed_wrapper import detect_by_structure
from .hmmer import hmmer
from .output import build_output_table
