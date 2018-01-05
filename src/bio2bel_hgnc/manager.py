# -*- coding: utf-8 -*-

from pybel.constants import IDENTIFIER, NAME, NAMESPACE
from pyhgnc.manager.database import DbManager
from pyhgnc.manager.query import QueryManager

__all__ = [
    'Manager',
]


def _deal_with_nonsense(results):
    if not results:
        return

    if 1 < len(results):
        raise ValueError

    return results[0]


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

    @staticmethod
    def ensure(connection=None):
        if connection is None or isinstance(connection, str):
            return Manager(connection=connection)
        return connection

    def __repr__(self):
        return '<{} connection={}>'.format(self.__class__.__name__, self.engine.url)
