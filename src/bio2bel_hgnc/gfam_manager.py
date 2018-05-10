# -*- coding: utf-8 -*-

"""Bio2BEL HGNC Gene Family Manager."""

import logging

from bio2bel.namespace_manager import NamespaceManagerMixin
from pybel.manager.models import NamespaceEntry
from .models import Base, GeneFamily
from .wrapper import BaseManager

log = logging.getLogger(__name__)

__all__ = [
    'Manager',
    'main',
]


class Manager(NamespaceManagerMixin, BaseManager):
    """Bio2BEL HGNC Manager"""

    module_name = 'gfam'
    namespace_model = GeneFamily

    @property
    def _base(self):
        return Base

    def populate(self, *args, **kwargs):
        raise NotImplementedError

    def is_populated(self):
        return 0 < self._count_model(self.namespace_model)

    def _create_namespace_entry_from_model(self, gene_family, namespace):
        return NamespaceEntry(
            encoding='G',
            identifier=gene_family.family_identifier,
            name=gene_family.family_name,
            namespace=namespace
        )

    def _get_identifier(self, gene_family):
        return gene_family.family_identifier


main = Manager.get_cli()

if __name__ == '__main__':
    main()
