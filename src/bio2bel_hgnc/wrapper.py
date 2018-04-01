# -*- coding: utf-8 -*-

from pyhgnc.manager.database import DbManager
from pyhgnc.manager.query import QueryManager


class BaseManager(DbManager, QueryManager):
    """Wraps the relevant PyHGNC functions"""
