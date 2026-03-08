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
