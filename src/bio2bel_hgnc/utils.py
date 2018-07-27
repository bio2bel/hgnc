# -*- coding: utf-8 -*-

"""Utilities for Bio2BEL HGNC."""

from .constants import VERSION

__all__ = [
    'get_version',
]


def get_version() -> str:
    """Return the software version of Bio2BEL HGNC."""
    return VERSION
