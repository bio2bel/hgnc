# -*- coding: utf-8 -*-

from pybel.constants import FUNCTION, GENE, IDENTIFIER, IS_A, NAME, NAMESPACE
from pybel.dsl import gene
from pyhgnc.manager.database import DbManager
from pyhgnc.manager.query import QueryManager
from .constants import GENE_FAMILY_KEYWORD
from .models import HGNC

__all__ = [
    'Manager',
]


def _deal_with_nonsense(results):
    if not results:
        return

    if 1 < len(results):
        raise ValueError

    return results[0]


def family_to_gene(family):
    """Converts a PyHGNC Gene Family model to a PyBEL gene

    :param pyhgnc.manager.models.GeneFamily family: An HGNC Gene Family model
    :rtype: pybel.dsl.gene
    """
    return gene(
        namespace=GENE_FAMILY_KEYWORD,
        identifier=family.family_identifier,
        name=family.family_name
    )


class Manager(DbManager, QueryManager):
    """An extended version of the PyHGNC manager to have useful functions"""

    def populate(self, *args, **kwargs):
        self.db_import(*args, **kwargs)

    def drop_all(self):
        self._drop_tables()

    def get_gene_by_hgnc_symbol(self, hgnc_symbol):
        """Gets a gene by the symbol

        :param str hgnc_symbol: The HGNC gene symbol
        :rtype: Optional[pyhgnc.manager.models.HGNC]
        """
        results = self.hgnc(symbol=hgnc_symbol)
        return _deal_with_nonsense(results)

    def get_gene_by_hgnc_id(self, hgnc_id):
        """Gets a gene by the identifier

        :param str hgnc_id: The HGNC gene identifier
        :rtype: Optional[pyhgnc.manager.models.HGNC]
        """
        results = self.hgnc(identifier=hgnc_id)
        return _deal_with_nonsense(results)

    def get_gene_by_entrez_id(self, entrez_id):
        """Gets a HGNC gene by its Entrez gene identifier

        :param str entrez_id: The Entrez gene identifier
        :rtype: Optional[pyhgnc.manager.models.HGNC]
        """
        results = self.hgnc(entrez=entrez_id)
        return _deal_with_nonsense(results)

    def get_gene_by_mgi_id(self, mgi_id):
        """Gets a HGNC gene by an orthologous MGI identifier

        :param str mgi_id: MGI identifier
        :rtype: Optional[pyhgnc.manager.models.HGNC]
        """
        raise NotImplementedError

    def get_gene_by_mgi_symbol(self, mgi_symbol):
        """Gets a HGNC gene by an orthologous MGI gene symbol

        :param str mgi_symbol: MGI gene symbol
        :rtype: Optional[pyhgnc.manager.models.HGNC]
        """
        # TODO how to deal with getting the MGI name to MGI identifiers mapping? Should this be part of PyHGNC, or something different?

        raise NotImplementedError

    def get_gene_by_rgd_id(self, rgd_id):
        """Gets a HGNC gene by an orthologous RGD identifier

        :param str rgd_id: RGD identifier
        :rtype: Optional[pyhgnc.manager.models.HGNC]
        """
        raise NotImplementedError

    def get_gene_by_rgd_symbol(self, rgd_symbol):
        """Gets a HGNC gene by an orthologous RGD identifier

        :param str rgd_symbol: RGD gene symbol
        :rtype: Optional[pyhgnc.manager.models.HGNC]
        """
        # TODO how to deal with getting the RGD name to RGD identifiers mapping? Should this be part of PyHGNC, or something different?

        raise NotImplementedError

    def get_node(self, graph, node):
        """Gets a node from the PyHGNC database, whether it has a HGNC, RGD, MGI, or EG identifier.

        :param pybel.BELGraph graph: A BEL graph
        :param tuple node: A PyBEL node tuple
        :param Optional[pyhgnc.manager.query.QueryManager] manager: A PyHGNC database manager
        :rtype: pyhgnc.manager.models.HGNC
        :raises: KeyError
        """
        data = graph.node[node]

        if NAMESPACE not in data:
            raise KeyError

        namespace = data[NAMESPACE]

        if namespace == 'HGNC':
            if IDENTIFIER in data:
                return self.get_gene_by_hgnc_id(data[IDENTIFIER])
            elif NAME in data:
                return self.get_gene_by_hgnc_symbol(data[NAME])
            raise KeyError

        if namespace in {'ENTREZ', 'EGID', 'EG'}:
            if IDENTIFIER in data:
                return self.get_gene_by_entrez_id(data[IDENTIFIER])
            elif NAME in data:
                return self.get_gene_by_entrez_id(data[NAME])
            raise KeyError

        if namespace in {'MGI'}:
            if IDENTIFIER in data:
                return self.get_gene_by_mgi_id(data[IDENTIFIER])
            elif NAME in data:
                return self.get_gene_by_mgi_symbol(data[NAME])
            raise KeyError

        if namespace in {'RGD'}:
            if IDENTIFIER in data:
                return self.get_gene_by_rgd_id(data[IDENTIFIER])
            elif NAME in data:
                return self.get_gene_by_rgd_symbol(data[NAME])
            raise KeyError

    def enrich_genes_with_families(self, graph):
        """Enrich genes in the BEL graph with their families

        :param pybel.BELGraph graph: The BEL graph to enrich
        """
        for n, data in graph.nodes(data=True):
            if data[FUNCTION] != GENE:
                continue

            if data.get(NAMESPACE) != 'HGNC':
                continue

            if IDENTIFIER in data:
                m = self.get_gene_by_hgnc_id(data[IDENTIFIER])
            elif NAME in data:
                m = self.get_gene_by_hgnc_symbol(data[NAME])
            else:
                raise ValueError

            for family in m.gene_families:
                graph.add_unqualified_edge(n, family_to_gene(family), IS_A)

    def get_family_by_id(self, family_identifier):
        """Gets a gene family by its identifier

        :param str family_identifier: The identifier of a HGNC Gene Family
        :rtype: Optional[GeneFamily]
        """
        results = self.gene_family(family_identifier=family_identifier)
        return _deal_with_nonsense(results)

    def get_family_by_name(self, family_name):
        """Gets a gene family by its name

        :param str family_name: The name of a HGNC Gene Family
        :rtype: Optional[GeneFamily]
        """
        results = self.gene_family(family_name=family_name)
        return _deal_with_nonsense(results)

    def enrich_families_with_genes(self, graph):
        """Enrich gene families in the BEL graph with their member genes

        :param pybel.BELGraph graph: The BEL graph to enrich
        """
        for n, data in graph.nodes(data=True):
            if data[FUNCTION] != GENE:
                continue

            if data.get(NAMESPACE) != GENE_FAMILY_KEYWORD:
                continue

            if IDENTIFIER in data:
                m = self.get_family_by_id(data[IDENTIFIER])
            elif NAME in data:
                m = self.get_family_by_name(data[NAME])
            else:
                raise ValueError

            for g in m.hgncs:
                graph.add_unqualified_edge(gene(namespace='HGNC', name=g.symbol, identifier=g.identifier), n, IS_A)

    def build_hgnc_id_symbol_mapping(self):
        """Builds a mapping from HGNC identifier to HGNC symbol

        :rtype: dict[str,str]
        """
        return {
            str(identifier): symbol
            for identifier, symbol in self.session.query(HGNC.identifier, HGNC.symbol).all()
        }

    def build_hgnc_symbol_id_mapping(self):
        """Builds a mapping from HGNC symbol to HGNC identifier

        :rtype: dict[str,str]
        """
        return {
            symbol: str(identifier)
            for symbol, identifier in self.session.query(HGNC.symbol, HGNC.identifier).all()
        }

    @staticmethod
    def ensure(connection=None):
        """
        :param Optional[str or Manager] connection:
        :rtype: Manager
        """
        if connection is None or isinstance(connection, str):
            return Manager(connection=connection)
        return connection

    def __repr__(self):
        return '<{} connection={}>'.format(self.__class__.__name__, self.engine.url)
