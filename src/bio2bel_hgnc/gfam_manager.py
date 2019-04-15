# -*- coding: utf-8 -*-

"""Bio2BEL HGNC Gene Family Manager."""

import logging
from typing import List, Mapping

from tqdm import tqdm

from bio2bel import AbstractManager
from bio2bel.manager.bel_manager import BELManagerMixin
from bio2bel.manager.namespace_manager import BELNamespaceManagerMixin
from pybel import BELGraph
from pybel.manager.models import Namespace, NamespaceEntry
from .model_utils import family_to_bel, gene_to_bel
from .models import (
    Base, GeneFamily, gene_gene_family,
)
from .wrapper import BaseManager

__all__ = [
    'Manager',
    'main',
]

log = logging.getLogger(__name__)


class Manager(AbstractManager, BELNamespaceManagerMixin, BELManagerMixin, BaseManager):
    """Human gene-gene family memberships."""

    _base = Base
    module_name = 'hgnc.genefamily'

    namespace_model = GeneFamily
    edge_model = gene_gene_family
    identifiers_recommended = 'HGNC gene family'
    identifiers_pattern = r'^\d+$'
    identifiers_miriam = 'MIR:00000573'
    identifiers_namespace = 'hgnc.genefamily'
    identifiers_url = 'http://identifiers.org/hgnc.genefamily/'

    def populate(self, *args, **kwargs):  # noqa: D102
        raise NotImplementedError

    def is_populated(self) -> bool:
        """Check if the gene families have been populated."""
        return 0 < self._count_model(self.namespace_model)

    def _create_namespace_entry_from_model(self, gene_family: GeneFamily, namespace: Namespace) -> NamespaceEntry:
        return NamespaceEntry(
            encoding=self._get_encoding(gene_family),
            identifier=str(gene_family.family_identifier),
            name=gene_family.family_name,
            namespace=namespace,
        )

    def _get_identifier(self, gene_family: GeneFamily) -> str:
        return str(gene_family.family_identifier)

    def _get_name(self, gene_family: GeneFamily) -> str:
        return gene_family.family_name

    def _get_encoding(self, gene_family: GeneFamily) -> str:
        return 'GRP'

    def _get_namespace_name_to_encoding(self, **kwargs) -> Mapping[str, str]:
        return {
            model.family_name: 'GRP'
            for model in self._iterate_namespace_models(**kwargs)
        }

    def _get_namespace_identifier_to_encoding(self, **kwargs) -> Mapping[str, str]:
        return {
            str(model.family_identifier): 'GRP'
            for model in self._iterate_namespace_models(**kwargs)
        }

    def add_namespace_to_graph(self, graph: BELGraph):
        """Add this manager's namespace to the graph."""
        namespace = self.upload_bel_namespace()
        graph.namespace_url[namespace.keyword] = namespace.url

    def summarize(self) -> Mapping[str, int]:
        """Summarize the database."""
        return dict(
            families=self.count_families(),
        )

    def normalize_families(self, graph: BELGraph) -> None:
        """Normalize the families in the graph."""
        raise NotImplementedError

    def count_families(self) -> int:
        """Count all gene families."""
        return self._count_model(GeneFamily)

    def list_families(self) -> List[GeneFamily]:
        """List all gene families."""
        return self._list_model(GeneFamily)

    def to_bel(self) -> BELGraph:
        """Export gene family definitions as a BEL graph."""
        graph = BELGraph(
            name='HGNC Gene Family Memberships',
            version='1.0.0'
        )

        for family in tqdm(self.list_families(), total=self.count_families(),
                           desc='Mapping gene family definitions to BEL'):
            for human_gene in family.hgncs:
                graph.add_is_a(gene_to_bel(human_gene), family_to_bel(family))

        return graph


main = Manager.get_cli()

if __name__ == '__main__':
    main()
