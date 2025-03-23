"""
This file contains the logic required for configuration reading and controlling.
"""
from pathlib import Path
from typing import Any


def get(key: str, config_path: str | Path = None) -> Any:
    """
    This function searches the config file given by config_path and returns the value associated with key in this file.
    If key is not found in config_path, a key error will be raised.
    If key exists multiple times in config_path, the behaviour is undefined.
    If config_path is not given, the file is assumed to reside in the current working directory and to be named
    'config.yaml'.
    The file specified by config_path is assumed to be formatted as YAML.
    TODO: Determine, if YAML is really a good choice

    :param key: The configuration key to look for.
    :param config_path: The path to the config file in which this function shall search.
    :return: The value associated with key in the file denoted by config_path.
    :raises KeyError: If key does not exist in config_path.
    """
    _ = key, config_path
    raise NotImplementedError()
