# -*- coding: utf-8 -*-

from pyhgnc import QueryManager


def _deal_with_nonsense(results):
    if not results:
        return

    if 1 < len(results):
        raise ValueError

    return results[0]


class Manager(QueryManager):
    """An extended version of the PyHGNC manager to have useful functions"""

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
