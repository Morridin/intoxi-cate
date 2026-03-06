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
from .blast import on_toxins, on_uniprot
from .tmbed import detect_by_structure
from .toxins_peptides import retrieve_candidate_toxins
from .hmmer import hmmer
from .output import build_output_table


# Set parts of lib
__all__ = ["config", "utils", "transcriptome", "generate_peptides", "blast", "tmbed", "toxins_peptides", "output"]

# Add the different modules
__all__ += ["assemble_transcriptome", "cluster_peptides", "blast_on_toxins", "detect_by_structure",
            "retrieve_candidate_toxins", "hmmer", "blast_on_uniprot", "run_salmon", "build_output_table"]
