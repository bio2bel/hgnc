# -*- coding: utf-8 -*-

"""Enrichment functions for BEL graphs"""

from pybel.constants import NAMESPACE


def get_node(graph, node, manager=None):
    """Gets a node from the PyHGNC database, whether it has a HGNC, RGD, MGI, or EG identifier.

    :param pybel.BELGraph graph:
    :param tuple node: A PyBEL node tuple
    :rtype: pyhgnc.models.HGNC
    """
    raise NotImplementedError


def add_metadata(graph, node, manager=None):
    """Add to the node based on PyHGNC's database

    :param pybel.BELGraph graph:
    :param tuple node: A PyBEL node tuple
    :param optional[pyhgnc.QueryManager] manager: A PyHGNC database manager
    """
    raise NotImplementedError
