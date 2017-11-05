# -*- coding: utf-8 -*-

"""Enrichment functions for BEL graphs"""

from pybel.constants import NAMESPACE


def get_node(graph, node, manager=None):
    """Gets a node from the PyHGNC database, whether it has a HGNC, RGD, MGI, or EG identifier.

    :param pybel.BELGraph graph: A BEL graph
    :param tuple node: A PyBEL node tuple
    :param optional[pyhgnc.QueryManager] manager: A PyHGNC database manager
    :rtype: pyhgnc.manager.models.HGNC
    """
    raise NotImplementedError


def add_metadata(graph, node, manager=None):
    """Add to the node based on PyHGNC's database

    :param pybel.BELGraph graph: A BEL Graph
    :param tuple node: A PyBEL node tuple
    :param optional[pyhgnc.QueryManager] manager: A PyHGNC database manager
    """
    raise NotImplementedError


def add_node_equivalencies(grpah, node, manager=None, add_leaves=False):
    """Given an HGNC node, add all equivalencies to EG and UniProt

    :param pybel.BELGraph graph: A BEL graph
    :param tuple node: A PyBEL node tuple
    :param optional[pyhgnc.QueryManager] manager: A PyHGNC database manager
    :param bool add_leaves: Should equivalencies that are not already in the graph be added?
    """
    raise NotImplementedError


def add_node_orthologies(graph, node, manager=None, add_leaves=False):
    """Given a node that's HGNC, add orthology relationships

    :param pybel.BELGraph graph: A BEL graph
    :param tuple node: A PyBEL node tuple
    :param optional[pyhgnc.QueryManager] manager: A PyHGNC database manager
    :param bool add_leaves: Should orthologs that are not already in the graph be added?
    """
    raise NotImplementedError


def add_orthologies(graph, manager=None, add_leaves=False):
    """Adds orthology relationships for all HGNC nodes to all of the available data for nodes that are in the
    graph

    :param pybel.BELGraph graph: A BEL graph
    :param optional[pyhgnc.QueryManager] manager: A PyHGNC database manager
    :param bool add_leaves: Should orthologs that are not already in the graph be added?
    """
    for node, data in graph.nodes(data=True):
        if NAMESPACE in data and data[NAMESPACE] == 'HGNC':
            add_node_orthologies(graph, node, manager=manager, add_leaves=add_leaves)
