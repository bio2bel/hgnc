# -*- coding: utf-8 -*-

"""Bio2BEL HGNC Manager."""

import logging
from collections import Counter
from typing import Dict, Iterable, Optional, Tuple

import click
import networkx as nx
from bio2bel import AbstractManager
from bio2bel.manager.bel_manager import BELManagerMixin
from bio2bel.manager.flask_manager import FlaskMixin
from bio2bel.manager.namespace_manager import BELNamespaceManagerMixin
from pybel import BELGraph
from pybel.constants import FUNCTION, GENE, IDENTIFIER, MIRNA, NAME, NAMESPACE, PROTEIN, RNA, VARIANTS
from pybel.dsl import gene as gene_dsl, mirna as mirna_dsl, protein as protein_dsl, rna as rna_dsl
from pybel.dsl.nodes import BaseEntity
from pybel.manager.models import NamespaceEntry
from tqdm import tqdm

from .constants import MODULE_NAME, encodings
from .gfam_manager import Manager as GfamManager
from .model_utils import (
    add_central_dogma, family_to_bel, gene_to_bel, uniprot_to_bel,
)
from .models import Base, GeneFamily, HumanGene, MouseGene, RatGene, UniProt
from .wrapper import BaseManager

log = logging.getLogger(__name__)

__all__ = [
    'Manager',
]

_func_to_dsl = {
    GENE: gene_dsl,
    RNA: rna_dsl,
    PROTEIN: protein_dsl,
    MIRNA: mirna_dsl,
}


def _deal_with_nonsense(results):
    """
    :param list[X] results:
    :rtype: Optional[X]
    """
    if not results:
        return

    if 1 < len(results):
        raise ValueError('{}'.format(results))

    return results[0]


class Manager(AbstractManager, FlaskMixin, BELManagerMixin, BELNamespaceManagerMixin, BaseManager):
    """Bio2BEL HGNC Manager."""

    module_name = MODULE_NAME

    namespace_model = HumanGene
    identifiers_recommended = 'HGNC'
    identifiers_pattern = '^((HGNC|hgnc):)?\d{1,5}$'
    identifiers_miriam = 'MIR:00000080'
    identifiers_namespace = 'hgnc'
    identifiers_url = 'http://identifiers.org/hgnc/'

    flask_admin_models = [HumanGene, GeneFamily, UniProt, MouseGene, RatGene]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._hgnc_symbol_entrez_id_mapping = {}

    @property
    def _base(self):
        return Base

    def is_populated(self):
        """Check if the database is already populated.

        :rtype: bool
        """
        return 0 < self.count_human_genes()

    def populate(self, silent=False, hgnc_file_path=None, use_hcop=True, hcop_file_path=None, low_memory=False):
        """Populate the database."""
        json_data = self.load_hgnc_json(hgnc_file_path=hgnc_file_path)
        self.insert_hgnc(hgnc_dict=json_data, silent=silent, low_memory=low_memory)

        if use_hcop:
            self.insert_hcop(silent=silent, hcop_file_path=hcop_file_path)

    #: Clobber this PyHGNC function so it doesn't accidentally get called
    def db_import(self, silent=False, hgnc_file_path=None, hcop_file_path=None, low_memory=False):
        raise NotImplemented('call manager.populate instead')

    #: Clobber this PyHGNC function so it doesn't accidentally get called
    def _drop_tables(self):
        raise NotImplemented('call manager.drop_all instead')

    #####################
    # Summary Functions #
    #####################

    def count_human_genes(self):
        """Count the number of human genes in the database."""
        return self._count_model(HumanGene)

    def count_families(self):
        """Count the number of human gene families in the database."""
        return self._count_model(GeneFamily)

    def count_mouse_genes(self):
        """Count the number of mouse genes in the database."""
        return self._count_model(MouseGene)

    def count_rat_genes(self):
        """Count the number of rat genes in the database."""
        return self._count_model(RatGene)

    def count_uniprots(self):
        """Count the number of UniProt proteins in the database."""
        return self._count_model(UniProt)

    def summarize(self):
        """Summarize the database.

        :rtype: dict[str,int]
        """
        return dict(
            human_genes=self.count_human_genes(),
            rat_genes=self.count_rat_genes(),
            mouse_genes=self.count_mouse_genes(),
            families=self.count_families(),
            uniprots=self.count_uniprots()
        )

    def get_gene_by_hgnc_symbol(self, hgnc_symbol: str) -> Optional[HumanGene]:
        """Get a human gene by HGNC symbol.

        :param str hgnc_symbol: The HGNC gene symbol
        """
        results = self.hgnc(symbol=hgnc_symbol)
        return _deal_with_nonsense(results)

    def get_gene_by_hgnc_id(self, hgnc_id: str) -> Optional[HumanGene]:
        """Get a human gene by HGNC identifier.

        :param str hgnc_id: The HGNC gene identifier
        :rtype: Optional[bio2bel_hgnc.models.HumanGene]
        """
        results = self.hgnc(identifier=int(hgnc_id))  # it's actually an integer in the database
        return _deal_with_nonsense(results)

    def get_gene_by_entrez_id(self, entrez_id):
        """Get a human gene by its Entrez gene identifier.

        :param str entrez_id: The Entrez gene identifier
        :rtype: Optional[bio2bel_hgnc.models.HumanGene]
        """
        results = self.hgnc(entrez=entrez_id)
        return _deal_with_nonsense(results)

    def get_gene_by_ensembl_id(self, ensembl_id):
        """Get a human gene by its ENSEMBL gene identifier.

        :param str ensembl_id: The ENSEMBL gene identifier
        :rtype: Optional[bio2bel_hgnc.models.HumanGene]
        """
        results = self.hgnc(ensembl_gene=ensembl_id)
        return _deal_with_nonsense(results)

    def get_gene_by_uniprot_id(self, uniprot_id):
        """Get a human gene by its UniProt gene identifier.

        :param str uniprot_id: The UniProt gene identifier
        :rtype: Optional[bio2bel_hgnc.models.HumanGene]
        """
        return self.hgnc(uniprotid=uniprot_id)

    def get_gene_by_mgi_id(self, mgi_id):
        """Get a human gene by an orthologous MGI identifier.

        :param str mgi_id: MGI identifier
        :rtype: Optional[bio2bel_hgnc.models.HumanGene]
        """
        results = self.mgd(mgdid=mgi_id)
        mouse_gene = _deal_with_nonsense(results)

        if mouse_gene is None:
            return

        human_genes = mouse_gene.hgncs

        if len(human_genes) > 1:
            log.warning('multiple human genes mapped to mgi_id:%s: %s', mgi_id, human_genes)
            return

        return human_genes[0]

    def get_gene_by_rgd_id(self, rgd_id):
        """Get a human gene by an orthologous RGD identifier.

        :param str rgd_id: RGD identifier
        :rtype: Optional[bio2bel_hgnc.models.HumanGene]
        """
        results = self.rgd(rgdid=rgd_id)
        rat_gene = _deal_with_nonsense(results)

        if rat_gene is None:
            return

        human_genes = rat_gene.hgncs

        if len(human_genes) > 1:
            log.warning('multiple human genes mapped to rgd_id:%s: %s', rgd_id, human_genes)
            return

        return human_genes[0]

    def get_node(self, graph, node) -> Optional[HumanGene]:
        """Get a node from the database, whether it has a HGNC, RGD, MGI, or EG identifier.

        :param pybel.BELGraph graph: A BEL graph
        :param node: A PyBEL node tuple
        :type node: tuple or BaseEntity
        :raises: KeyError
        """
        if isinstance(node, BaseEntity):
            data = node
        else:
            data = graph.node[node]

        if NAMESPACE not in data:
            return

        namespace = data[NAMESPACE]
        identifer = data.get(IDENTIFIER)
        name = data.get(NAME)

        if namespace.lower() == 'hgnc':
            if identifer is not None:
                return self.get_gene_by_hgnc_id(identifer)
            elif name is not None:
                return self.get_gene_by_hgnc_symbol(name)
            raise KeyError

        if namespace.lower() in {'entrez', 'egid', 'eg', 'ncbigene'}:
            if identifer is not None:
                return self.get_gene_by_entrez_id(identifer)
            elif name is not None:
                return self.get_gene_by_entrez_id(name)
            raise KeyError

        if namespace.lower() in {'mgi'}:
            if identifer is not None:
                return self.get_gene_by_mgi_id(identifer)
            elif name is not None:
                return
            raise KeyError

        if namespace == 'MGIID':
            if name is None:
                raise KeyError
            return self.get_gene_by_mgi_id(name)

        if namespace.lower() in {'rgd'}:
            if identifer is not None:
                return self.get_gene_by_rgd_id(identifer)
            elif name is not None:
                return
            raise KeyError

        if namespace == 'RGDID':
            if name is None:
                raise KeyError
            return self.get_gene_by_rgd_id(name)

    def add_namespace_to_graph(self, graph):
        """Add this manager's namespace to the graph.

        :param pybel.BELGraph graph:
        """
        namespace = self.upload_bel_namespace()
        graph.namespace_url[namespace.keyword] = namespace.url

        gfam_manager = GfamManager(engine=self.engine, session=self.session)
        gfam_manager.add_namespace_to_graph(graph)

    def _iter_genes(self, graph) -> Iterable[Tuple[tuple, dict, HumanGene]]:
        for node_tuple, node_data in graph.nodes(data=True):
            human_gene = self.get_node(graph, node_tuple)
            if human_gene is not None:
                yield node_tuple, node_data, human_gene

    def normalize_genes(self, graph) -> None:
        """Add identifiers to all HGNC genes.

        :param pybel.BELGraph graph: The BEL graph to enrich
        """
        mapping = {}

        for node_tuple, node_data, human_gene in self._iter_genes(graph):
            dsl = gene_to_bel(human_gene, func=node_data[FUNCTION], variants=node_data.get(VARIANTS))
            graph.node[node_tuple] = dsl
            mapping[node_tuple] = dsl.as_tuple()

        nx.relabel_nodes(graph, mapping, copy=False)

    def enrich_genes_with_equivalences(self, graph: BELGraph) -> None:
        """Enrich genes with their corresponding UniProt."""
        self.add_namespace_to_graph(graph)

        if 'uniprot' not in graph.namespace_url:
            graph.namespace_pattern[
                'uniprot'] = '^([A-N,R-Z][0-9]([A-Z][A-Z, 0-9][A-Z, 0-9][0-9]){1,2})|([O,P,Q][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9])(\.\d+)?$'

        for node_tuple, node_data, human_gene in self._iter_genes(graph):
            func = node_data[FUNCTION]

            if human_gene.entrez:
                graph.add_equivalence(node_tuple, _func_to_dsl[func](
                    namespace='ncbigene',
                    name=human_gene.symbol,
                    identifier=str(human_gene.entrez)
                ))

            if func in {PROTEIN, RNA, GENE}:
                for uniprot in human_gene.uniprots:
                    graph.add_equivalence(node_tuple, uniprot_to_bel(uniprot))

            if func in {RNA, GENE}:
                if human_gene.mirbase:
                    mirbase_rna_node = rna_dsl(
                        namespace='mirbase',
                        identifier=str(human_gene.mirbase)
                    )
                    graph.add_equivalence(node_tuple, mirbase_rna_node)

    def enrich_genes_with_families(self, graph: BELGraph) -> None:
        """Enrich genes in the BEL graph with their families."""
        self.add_namespace_to_graph(graph)
        for node_tuple, node_data, human_gene in self._iter_genes(graph):
            for family in human_gene.gene_families:
                graph.add_is_a(node_tuple, family_to_bel(family, node_data[FUNCTION]))

    def get_family_by_id(self, family_identifier: str) -> Optional[GeneFamily]:
        """Get a gene family by its hgnc.genefamily identifier.

        :param family_identifier: The identifier of a HGNC Gene Family
        """
        results = self.gene_family(family_identifier=family_identifier)
        return _deal_with_nonsense(results)

    def get_family_by_name(self, family_name: str) -> Optional[GeneFamily]:
        """Get a gene family by its name.

        :param family_name: The name of a HGNC Gene Family
        """
        results = self.gene_family(family_name=family_name)
        return _deal_with_nonsense(results)

    def _enrich_hgnc_with_entrez_equivalences(self, graph, node):
        """

        :param pybel.BELGraph graph:
        :param node:
        :return: the hash of the edge added
        :rtype: str
        """
        data = graph.node[node]

        namespace = data.get(NAMESPACE)

        if namespace != 'HGNC':
            return

        func = data[FUNCTION]
        name = data[NAME]
        entrez = self.hgnc_symbol_entrez_id_mapping[name]

        return graph.add_equivalence(node, _func_to_dsl[func](
            namespace='ncbigene',
            name=name,
            identifier=str(entrez)
        ))

    def enrich_hgnc_with_entrez_equivalences(self, graph):
        """Add equivalent Entrez nodes for all HGNC genes.

        :param pybel.BELGraph graph: The BEL graph to enrich
        """
        self.add_namespace_to_graph(graph)

        for node in graph.nodes():
            self._enrich_hgnc_with_entrez_equivalences(graph, node)

    def enrich_families_with_genes(self, graph):
        """Enrich gene families in the BEL graph with their member genes.

        :param pybel.BELGraph graph: The BEL graph to enrich
        """
        self.add_namespace_to_graph(graph)

        for gene_family_node, data in graph.nodes(data=True):
            if data[FUNCTION] != GENE:
                continue

            if data.get(NAMESPACE).lower() not in {'gfam', 'hgnc.family', 'hgnc.genefamily'}:
                continue

            if IDENTIFIER in data:
                gene_family_model = self.get_family_by_id(data[IDENTIFIER])
            elif NAME in data:
                gene_family_model = self.get_family_by_name(data[NAME])
            else:
                raise ValueError

            if gene_family_model is None:
                log.info('family not found: %s', data)
                continue

            for human_gene in gene_family_model.hgncs:
                graph.add_is_a(gene_to_bel(human_gene), gene_family_node)

    """ Mapping dictionaries"""

    def _get_identifier(self, human_gene: HumanGene) -> str:
        """Get the identifier from a human gene SQLAlchemy model.

        :rtype: str
        """
        return str(human_gene.identifier)

    def build_entrez_id_symbol_mapping(self) -> Dict[str, str]:
        """Build a mapping from Entrez gene identifier to HGNC gene symbols."""
        return dict(self.session.query(HumanGene.entrez, HumanGene.symbol).all())

    @property
    def hgnc_symbol_entrez_id_mapping(self) -> Dict[str, str]:
        """Get a mapping from Entrez gene identifiers to HGNC gene symbols."""
        if not self._hgnc_symbol_entrez_id_mapping:
            self._hgnc_symbol_entrez_id_mapping = self.build_hgnc_symbol_entrez_id_mapping()

        return self._hgnc_symbol_entrez_id_mapping

    def build_hgnc_symbol_entrez_id_mapping(self) -> Dict[str, str]:
        """Build a mapping from HGNC symbol to ENTREZ identifier.

        :rtype: dict[str,str]
        """
        return {
            symbol: identifier
            for symbol, identifier in self.session.query(HumanGene.symbol, HumanGene.entrez).all()
        }

    def build_hgnc_id_symbol_mapping(self):
        """Build a mapping from HGNC identifiers to HGNC gene symbols.

        :rtype: dict[str,str]
        """
        return {
            str(identifier): symbol
            for identifier, symbol in self.session.query(HumanGene.identifier, HumanGene.symbol).all()
        }

    def build_hgnc_symbol_id_mapping(self):
        """Build a mapping from HGNC gene symbols to HGNC identifiers.

        :rtype: dict[str,str]
        """
        return {
            symbol: str(identifier)
            for symbol, identifier in self.session.query(HumanGene.symbol, HumanGene.identifier).all()
        }

    def build_hgnc_symbol_uniprot_ids_mapping(self):
        """Build a mapping from HGNC gene symbols to UniProt identifiers.

        :rtype: dict[str,set[str]]
        """
        return {
            symbol: uniprot_id
            for symbol, uniprot_id in self.session.query(HumanGene.symbol, UniProt.uniprotid).all()
        }

    def build_hgnc_id_uniprot_ids_mapping(self):
        """Build a mapping from HGNC identifier to UniProt identifiers.

        :rtype: dict[str,set[str]]
        """
        return {
            hgnc_id: uniprot_id
            for hgnc_id, uniprot_id in self.session.query(HumanGene.identifier, UniProt.uniprotid).all()
        }

    def build_uniprot_id_hgnc_id_mapping(self):
        """Build a mapping from UniProt identifiers to HGNC identifiers.

        :rtype: dict[str,str]
        """
        return {
            uniprot_id: hgnc_id
            for hgnc_id, uniprot_id in self.session.query(HumanGene.identifier, UniProt.uniprotid).all()
        }

    def build_uniprot_id_hgnc_symbol_mapping(self):
        """Build a mapping from UniProt identifiers to HGNC gene symbols.

        :rtype: dict[str,str]
        """
        return {
            uniprot_id: symbol
            for symbol, uniprot_id in self.session.query(HumanGene.symbol, UniProt.uniprotid).all()
        }

    def get_all_hgnc_symbols(self):
        """Return the set of HGNC gene symbols in the database.

        :rtype: set[str]
        """
        return {
            _deal_with_nonsense(symbol)
            for symbol in self.session.query(HumanGene.symbol).all()
        }

    def _get_gene_encodings(self):
        """Get the name to encoding dictionary for HGNC gene names/

        :rtype: dict[str,str]
        """
        return {
            symbol: encodings.get(locus_type, 'GRP')
            for symbol, locus_type in self.session.query(HumanGene.symbol, HumanGene.locus_type).all()
        }

    def list_families(self):
        """List families in the database.

        :rtype: list[GeneFamily]
        """
        return self._list_model(GeneFamily)

    def list_human_genes(self):
        """List human genes in the database.

        :rtype: list[HumanGene]
        """
        return self._list_model(HumanGene)

    def to_bel(self):
        """Export gene family definitions as a BEL graph.

        :rtype: pybel.BELGraph
        """
        graph = BELGraph(
            name='HGNC Central Dogma and Gene Family Definitions',
            version='1.0.0',  # FIXME use database version
            authors='HUGO Gene Nomenclature Consortium',
            description='Gene transcription, translation, and memberships to Gene Families.',
            contact='charles.hoyt@scai.fraunhofer.de',
        )

        hgnc_namespace = self.upload_bel_namespace()
        log.info('using default namespace: %s at %s', hgnc_namespace, hgnc_namespace.url)
        graph.namespace_url[hgnc_namespace.keyword] = hgnc_namespace.url

        gfam_manager = GfamManager(connection=self.connection)
        gfam_namespace = gfam_manager.upload_bel_namespace()
        log.info('using default namespace: %s at %s', gfam_namespace, gfam_namespace.url)
        graph.namespace_url[gfam_namespace.keyword] = gfam_namespace.url

        for human_gene in tqdm(self.list_human_genes(), total=self.count_human_genes(),
                               desc='Mapping central dogma to BEL'):
            graph.add_node_from_data(gene_to_bel(human_gene))
            add_central_dogma(graph, human_gene)

        for family in tqdm(self.list_families(), total=self.count_families(),
                           desc='Mapping gene family definitions to BEL'):
            for human_gene in family.hgncs:
                graph.add_is_a(gene_to_bel(human_gene), family_to_bel(family))

        return graph

    def get_pathway_size_distribution(self):
        """

        :rtype: dict[str,int]
        """
        return Counter({
            family.family_name: len(family.hgncs)
            for family in self.session.query(GeneFamily)
        })

    def get_all_hgnc_symbols_family(self):
        """Get all Gene symbols that appear in gene families.

        :rtype: set[str]
        """
        return {
            human_gene.symbol
            for family in self.session.query(GeneFamily)
            for human_gene in family.hgncs
        }

    def add_central_dogma(self, graph, node):
        """Add the central dogma of biology.

        :param graph:
        :param node:
        """
        data = graph.node[node]
        if VARIANTS in data:
            return

        human_gene = self.get_node(graph, node)
        encoding = encodings.get(human_gene.locus_type, 'GRP')

        if 'M' in encoding:
            mirna = gene_to_bel(human_gene, func=MIRNA)
            graph.add_transcription(gene_to_bel(human_gene), mirna)
            return mirna

        if 'R' not in encoding:
            return

        rna = gene_to_bel(human_gene, func=RNA)
        graph.add_transcription(gene_to_bel(human_gene), rna)

        if 'P' not in encoding:
            return rna
        else:
            protein = gene_to_bel(human_gene, func=PROTEIN)
            graph.add_translation(rna, protein)
            return protein

    def _create_namespace_entry_from_model(self, human_gene, namespace):
        return NamespaceEntry(
            encoding=encodings.get(human_gene.locus_type, 'GRP'),
            identifier=human_gene.identifier,
            name=human_gene.symbol,
            namespace=namespace
        )

    @staticmethod
    def _cli_add_populate(main):
        """Override default method to make it possible to add more flags.

        :type main: click.Group
        :rtype: click.Group
        """

        @main.command()
        @click.option('--reset', is_flag=True)
        @click.option('--skip-hcop', is_flag=True)
        @click.pass_obj
        def populate(manager, reset, skip_hcop):
            """Populate the database."""

            if reset:
                log.info('Deleting the previous instance of the database')
                manager.drop_all()
                log.info('Creating new models')
                manager.create_all()

            manager.populate(use_hcop=(not skip_hcop))

        return main
