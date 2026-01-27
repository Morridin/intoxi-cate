"""
This file contains the logic required for configuration reading and controlling.
"""
from pathlib import Path
from typing import TypeVar

import yaml

T = TypeVar("T")

__all__ = ["config"]

class Config:
    def __init__(self, config_path: Path = Path("config.yaml")):
        self.config_path = config_path
        with open(config_path, "r") as f:
            self.config = dict(yaml.safe_load(f))

    def set_config_file(self, config_path: Path) -> None:
        """
        This function sets a new config file and loads its data to the config handler.
        The config file is expected to be formatted as YAML, with a dictionary structure (aka key-value paris) at root level.
        If the file given by `config_path` does not exist, or can't be read or parsed, errors will be raised.
        :param config_path: The path to the config file holding the configuration data valid from this point.
        """
        self.config_path = config_path
        with open(config_path, "r") as f:
            self.config = dict(yaml.safe_load(f))

    def get(self, key: str, default: T = None) -> T | None:
        """
        This function searches the config file given by config_path and returns the value associated with key in this file.
        If key is not found in config_path, or is an empty string, None is returned.

        :param default: A value to return if `key` is not found in config_path
        :param key: The configuration key to look for.
        :return: The value associated with `key` in the file denoted by config_path.
        """
        output = self.config.get(key, default)

        if output == "":
            return None

        return output


config = Config()
