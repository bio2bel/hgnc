# -*- coding: utf-8 -*-

"""Enrichment functions for BEL graphs"""

from pybel.constants import FUNCTION, GENE, NAMESPACE

from .manager import Manager

__all__ = [
    'get_node',
    'add_node_orthologies',
    'add_orthologies',
    'add_node_central_dogma',
    'add_central_dogma',
]


def get_node(graph, node, manager=None):
    """Gets a node from the database, whether it has a HGNC, RGD, MGI, or EG identifier.

    :param pybel.BELGraph graph: A BEL graph
    :param tuple node: A PyBEL node tuple
    :param Optional[str or bio2bel.Manager] manager: A database manager
    :rtype: pyhgnc.manager.models.HGNC
    """
    if manager is None:
        manager = Manager()
    return manager.get_node(graph, node)


def add_node_orthologies(graph, node, manager=None, add_leaves=False):
    """Given a node that's HGNC, add orthology relationships

    :param pybel.BELGraph graph: A BEL graph
    :param tuple node: A PyBEL node tuple
    :param Optional[str or bio2bel_hgnc.Manager] manager: A database manager
    :param bool add_leaves: Should orthologs that are not already in the graph be added?
    """
    raise NotImplementedError


def add_orthologies(graph, manager=None, add_leaves=False):
    """Adds orthology relationships for all HGNC nodes to all of the available data for nodes that are in the
    graph

    :param pybel.BELGraph graph: A BEL graph
    :param Optional[str or bio2bel_hgnc.Manager] manager: A database manager
    :param bool add_leaves: Should orthologs that are not already in the graph be added?
    """
    if manager is None:
        manager = Manager()

    for node, data in graph.nodes(data=True):
        if NAMESPACE in data and data[NAMESPACE] == 'HGNC':
            add_node_orthologies(graph, node, manager=manager, add_leaves=add_leaves)


def add_node_equivalencies(graph, node, manager=None):
    """Given an HGNC node, add all equivalencies to EG and UniProt

    :param pybel.BELGraph graph: A BEL graph
    :param tuple node: A PyBEL node tuple
    :param Optional[bio2bel_hgnc.Manager] manager: A database manager
    """
    if manager is None:
        manager = Manager()
    manager.add_node_equivalencies(graph, node)


def add_node_central_dogma(graph, node, manager=None):
    """Uses the field :py:attr:`pyhgnc.manager.models.HGNC.locus_type` to determine if a gene node should be extended
    with nothing, an miRNA, an RNA, and possibly then a protein.

    :param pybel.BELGraph graph: A BEL graph
    :param tuple node: A PyBEL node tuple
    :param Optional[str or bio2bel_hgnc.Manager] manager: A database manager
    """
    if manager is None:
        manager = Manager()
    return manager.add_central_dogma(graph, node)


def add_central_dogma(graph, manager=None):
    """Add central dogma information for all gene nodes when possible

    :param pybel.BELGraph graph: A BEL graph
    :param Optional[str or bio2bel_hgnc.Manager] manager: A database manager
    """
    if manager is None:
        manager = Manager()

    for node, data in graph.nodes(data=True):
        if data[FUNCTION] == GENE:
            add_node_central_dogma(graph, node, manager=manager)
