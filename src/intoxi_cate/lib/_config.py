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
This file contains the logic required for configuration reading and controlling.
"""
import os
from pathlib import Path
from typing import TypeVar, Optional, Any

import yaml

T = TypeVar("T")

__all__ = ["config"]

_WOLF_PSORT_GITHUB_URL = "https://github.com/fmaguire/WoLFPSort/archive/refs/heads/master.zip"
_WOLF_PSORT_INSTALL_DIR = Path(os.environ.get("CONDA_PREFIX", "")) / "bin" / "WoLFPSort"

class Config:
    _instance = None
    _DEFAULT = {
        "proteome_fasta": None,
        "R1": None,
        "R2": None,
        "transcriptome": None,
        "toxin_db": None,
        "basename": "",
        "contaminants": None,
        "contamination_evalue": 1E-5,
        "clustering_threshold": 0.99,
        "cys_pattern": False,
        "maxlen": 45_000,
        "memory": None,
        "minlen": 99,
        "mmseqs_sequence_coverage": 0.0,
        "mmseqs_path": None,
        "output_dir": Path(""),
        "pfam_db_path": None,
        "quant": False,
        "repeated_aa_threshold": 5,
        "signalpeptide_minlen": 10,
        "swissprot": False,
        "swissprot_evalue": 1E-10,
        "swissprot_db_path": None,
        "threads": None,
        "tmbed_model_path": None,
        "tmbed_use_cpu": True,
        "tmbed_use_gpu": True,
        "toxins_evalue": 1E-10,
        "TPMthreshold": 1000.0,
        "wolfpsort": False,
        "wolfPsort_path": _WOLF_PSORT_INSTALL_DIR / "bin" / "runWolfPsortSummary"
    }

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance.config = cls._DEFAULT.copy()
        return cls._instance

    def set_config_file(self, config_path: Path) -> None:
        """
        This function sets a new config file and loads its data to the config handler.
        The config file is expected to be formatted as YAML, with a dictionary structure (aka key-value paris) at root level.
        If the file given by `config_path` does not exist, or can't be read or parsed, errors will be raised.
        :param config_path: The path to the config file holding the configuration data valid from this point.
        """
        with open(config_path, "r") as f:
            self.config = self._DEFAULT | dict(yaml.safe_load(f))

    def get(self, key: str, default: Optional[T] = None) -> Optional[T]:
        """
        This function searches the config file given by config_path and returns the value associated with key in this file.
        If key is not found in config_path, or is an empty string the default value is returned.
        If no default is set, None is returned.
        Exception to this is an empty string as default value which is returned if else None had been returned.

        :param default: A value to return if `key` is not found in config_path
        :param key: The configuration key to look for.
        :return: The value associated with `key` in the file denoted by config_path.
        """
        output = self.config.get(key, default)

        if output is None or output == "":
            return default

        return output

    def get_path(self, key: str) -> Optional[Path]:
        """
        This method does a config lookup as `get` does, but instead of returning a found value directly,
        this method attempts to transform said value into a Path object.
        If no value can be found or the transformation to a Path object fails, None is returned,
        else a Path object pointing to whatever value was found, interpreted as file system path.
        :param key: The configuration key to look for.
        :return: A path to a transformable value found in the config file, if any, else None.
        """
        output = self.config.get(key)
        if output is None:
            return None
        if isinstance(output, Path):
            return output
        return Path(str(output).strip())

config = Config()
