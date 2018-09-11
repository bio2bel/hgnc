# -*- coding: utf-8 -*-

"""Bio2BEL HGNC Gene Family Manager."""

import logging
from typing import Mapping

from bio2bel import AbstractManager
from bio2bel.manager.namespace_manager import BELNamespaceManagerMixin
from pybel import BELGraph
from pybel.manager.models import Namespace, NamespaceEntry
from .models import Base, GeneFamily
from .wrapper import BaseManager

log = logging.getLogger(__name__)

__all__ = [
    'Manager',
    'main',
]


class Manager(AbstractManager, BELNamespaceManagerMixin, BaseManager):
    """Bio2BEL HGNC Manager."""

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

    def populate(self, *args, **kwargs):  # noqa: D102
        raise NotImplementedError

    def is_populated(self) -> bool:
        """Check if the gene families have been populated."""
        return 0 < self._count_model(self.namespace_model)

    def _create_namespace_entry_from_model(self, gene_family: GeneFamily, namespace: Namespace) -> NamespaceEntry:
        return NamespaceEntry(
            encoding='G',
            identifier=str(gene_family.family_identifier),
            name=gene_family.family_name,
            namespace=namespace,
        )

    def _get_identifier(self, gene_family: GeneFamily) -> str:
        return str(gene_family.family_identifier)

    def add_namespace_to_graph(self, graph: BELGraph):
        """Add this manager's namespace to the graph."""
        namespace = self.upload_bel_namespace()
        graph.namespace_url[namespace.keyword] = namespace.url

    def summarize(self) -> Mapping[str, int]:
        """Summarize the database."""
        return dict(families=self._count_model(GeneFamily))


main = Manager.get_cli()

if __name__ == '__main__':
    main()
