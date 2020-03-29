# -*- coding: utf-8 -*-

"""Utilities for Bio2BEL HGNC."""

__all__ = [
    'VERSION',
    'get_version',
]

VERSION = '0.3.1-dev'


def get_version() -> str:
    """Return the software version of Bio2BEL HGNC."""
    return VERSION
