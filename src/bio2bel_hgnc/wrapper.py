# -*- coding: utf-8 -*-

"""Wrap PyHGNC functionality."""

from pyhgnc.manager.database import DbManager
from pyhgnc.manager.query import QueryManager


class BaseManager(DbManager, QueryManager):
    """Wraps the PyHGNC managers."""
