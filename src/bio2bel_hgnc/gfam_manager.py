# -*- coding: utf-8 -*-

"""Bio2BEL HGNC Gene Family Manager."""

import logging

from bio2bel import AbstractManager
from bio2bel.manager.namespace_manager import BELNamespaceManagerMixin
from pybel.manager.models import NamespaceEntry
from .models import Base, GeneFamily
from .wrapper import BaseManager

log = logging.getLogger(__name__)

__all__ = [
    'Manager',
    'main',
]


class Manager(AbstractManager, BELNamespaceManagerMixin, BaseManager):
    """Bio2BEL HGNC Manager"""

    module_name = 'gfam'
    namespace_model = GeneFamily

    identifiers_reccommended = 'HGNC gene family'
    identifiers_pattern = '^\d+$'
    identifiers_miriam = 'MIR:00000573'
    identifiers_namespace = 'hgnc.genefamily'
    identifiers_url = 'http://identifiers.org/hgnc.genefamily/'

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
