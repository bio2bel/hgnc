# -*- coding: utf-8 -*-

"""Enrichment functions for BEL graphs"""

from pybel.constants import FUNCTION, GENE, NAMESPACE

from .manager import Manager

__all__ = [
    'get_node',
    'add_metadata',
    'add_node_orthologies',
    'add_orthologies',
    'add_node_central_dogma',
    'add_central_dogma',
]


def get_node(graph, node, connection=None):
    """Gets a node from the database, whether it has a HGNC, RGD, MGI, or EG identifier.

    :param pybel.BELGraph graph: A BEL graph
    :param tuple node: A PyBEL node tuple
    :param Optional[str or bio2bel.Manager] connection: A database manager
    :rtype: pyhgnc.manager.models.HGNC
    """
    manager = Manager.ensure(connection=connection)
    return manager.get_node(graph, node)


def add_metadata(graph, node, manager=None):
    """Add to the node based on PyHGNC's database

    :param pybel.BELGraph graph: A BEL Graph
    :param tuple node: A PyBEL node tuple
    :param Optional[str or bio2bel_hgnc.Manager] manager: A database manager
    """
    raise NotImplementedError


def add_node_orthologies(graph, node, manager=None, add_leaves=False):
    """Given a node that's HGNC, add orthology relationships

    :param pybel.BELGraph graph: A BEL graph
    :param tuple node: A PyBEL node tuple
    :param Optional[str or bio2bel_hgnc.Manager] manager: A database manager
    :param bool add_leaves: Should orthologs that are not already in the graph be added?
    """
    raise NotImplementedError


def add_orthologies(graph, connection=None, add_leaves=False):
    """Adds orthology relationships for all HGNC nodes to all of the available data for nodes that are in the
    graph

    :param pybel.BELGraph graph: A BEL graph
    :param Optional[str or bio2bel_hgnc.Manager] connection: A database manager
    :param bool add_leaves: Should orthologs that are not already in the graph be added?
    """
    manager = Manager.ensure(connection)

    for node, data in graph.nodes(data=True):
        if NAMESPACE in data and data[NAMESPACE] == 'HGNC':
            add_node_orthologies(graph, node, manager=manager, add_leaves=add_leaves)


def add_node_equivalencies(graph, node, connection=None, add_leaves=False):
    """Given an HGNC node, add all equivalencies to EG and UniProt

    :param pybel.BELGraph graph: A BEL graph
    :param tuple node: A PyBEL node tuple
    :param Optional[str or bio2bel_hgnc.Manager] connection: A database manager
    """
    manager = Manager.ensure(connection)
    manager.add_node_equivalencies(graph, node, add_leaves=add_leaves)


def add_node_central_dogma(graph, node, connection=None):
    """Uses the field :py:attr:`pyhgnc.manager.models.HGNC.locus_type` to determine if a gene node should be extended
    with nothing, an miRNA, an RNA, and possibly then a protein.

    :param pybel.BELGraph graph: A BEL graph
    :param tuple node: A PyBEL node tuple
    :param Optional[str or bio2bel_hgnc.Manager] connection: A database manager
    """
    manager = Manager.ensure(connection)
    return manager.add_central_dogma(graph, node)


def add_central_dogma(graph, connection=None):
    """Add central dogma information for all gene nodes when possible

    :param pybel.BELGraph graph: A BEL graph
    :param Optional[str or bio2bel_hgnc.Manager] connection: A database manager
    """
    connection = Manager.ensure(connection)

    for node, data in graph.nodes(data=True):
        if data[FUNCTION] == GENE:
            add_node_central_dogma(graph, node, connection=connection)
