"""
This file contains the logic required for configuration reading and controlling.
"""
from pathlib import Path
from typing import TypeVar, Optional

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

    def get(self, key: str, default: Optional[T]) -> Optional[T]:
        """
        This function searches the config file given by config_path and returns the value associated with key in this file.
        If key is not found in config_path, or is an empty string, None is returned.
        Exception to this is an empty string as default value which is returned if else None had been returned.

        :param default: A value to return if `key` is not found in config_path
        :param key: The configuration key to look for.
        :return: The value associated with `key` in the file denoted by config_path.
        """
        output = self.config.get(key, default)

        if output == "" != default:
            return None

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
        if type(output) != str:
            return None
        return Path(output.strip())

config = Config()
