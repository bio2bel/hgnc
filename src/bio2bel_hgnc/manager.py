# -*- coding: utf-8 -*-

"""Bio2BEL HGNC Manager."""

import logging
from collections import Counter
from typing import Iterable, List, Mapping, Optional, Set, Tuple

import click
from networkx import relabel_nodes
from tqdm import tqdm

import pybel.dsl
from bio2bel import AbstractManager
from bio2bel.manager.bel_manager import BELManagerMixin
from bio2bel.manager.flask_manager import FlaskMixin
from bio2bel.manager.namespace_manager import BELNamespaceManagerMixin
from pybel import BELGraph
from pybel.constants import FUNCTION, GENE, IDENTIFIER, MIRNA, NAME, NAMESPACE, PROTEIN, RNA, VARIANTS
from pybel.dsl import BaseEntity, CentralDogma, FUNC_TO_DSL, rna as rna_dsl
from pybel.manager.models import Namespace, NamespaceEntry
from .constants import ENCODINGS, ENTREZ, MODULE_NAME
from .gfam_manager import Manager as GfamManager
from .model_utils import add_central_dogma, family_to_bel, gene_to_bel, uniprot_to_bel
from .models import (
    AliasName, AliasSymbol, Base, Enzyme, GeneFamily, HumanGene, MouseGene, RatGene, UniProt, gene_enzyme,
    gene_mouse_gene, gene_rat_gene,
)
from .wrapper import BaseManager

log = logging.getLogger(__name__)

__all__ = [
    'Manager',
]

UNIPROT_RE = r'^([A-N,R-Z][0-9]([A-Z][A-Z, 0-9][A-Z, 0-9][0-9]){1,2})|([O,P,Q][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9])(\.\d+)?$'

GENE_FAMILY_NAMESPACES = {'gfam', 'hgnc.family', 'hgnc.genefamily'}


def _deal_with_nonsense(results):
    """Fix lists that only have one thing inside them.

    :param list[X] results:
    :rtype: Optional[X]
    """
    if not results:
        return

    if 1 < len(results):
        raise ValueError('{}'.format(results))

    return results[0]


class Manager(AbstractManager, FlaskMixin, BELManagerMixin, BELNamespaceManagerMixin, BaseManager):
    """Human gene nomenclature and orthologies to mouse and rat."""

    module_name = MODULE_NAME
    _base = Base

    namespace_model = HumanGene
    edge_model = [gene_rat_gene, gene_enzyme, gene_mouse_gene]
    identifiers_recommended = 'HGNC'
    identifiers_pattern = r'^((HGNC|hgnc):)?\d{1,5}$'
    identifiers_miriam = 'MIR:00000080'
    identifiers_namespace = 'hgnc'
    identifiers_url = 'http://identifiers.org/hgnc/'

    flask_admin_models = [HumanGene, GeneFamily, UniProt, MouseGene, RatGene, AliasName, AliasSymbol, Enzyme]

    def __init__(self, *args, **kwargs):  # noqa: D107
        super().__init__(*args, **kwargs)
        self._hgnc_symbol_entrez_id_mapping = {}

    def is_populated(self) -> bool:
        """Check if the database is already populated."""
        return 0 < self.count_human_genes()

    def populate(self, silent=False, hgnc_file_path=None, use_hcop=False, hcop_file_path=None, low_memory=False):
        """Populate the database."""
        json_data = self.load_hgnc_json(hgnc_file_path=hgnc_file_path)
        self.insert_hgnc(hgnc_dict=json_data, silent=silent, low_memory=low_memory)

        if use_hcop:
            self.insert_hcop(silent=silent, hcop_file_path=hcop_file_path)

    #: Clobber this PyHGNC function so it doesn't accidentally get called
    def db_import(self, silent=False, hgnc_file_path=None, hcop_file_path=None, low_memory=False):  # noqa: D102
        raise NotImplementedError('call manager.populate instead')

    #: Clobber this PyHGNC function so it doesn't accidentally get called
    def _drop_tables(self):  # noqa: D102
        raise NotImplementedError('call manager.drop_all instead')

    #####################
    # Summary Functions #
    #####################

    def count_human_genes(self) -> int:
        """Count the number of human genes in the database."""
        return self._count_model(HumanGene)

    def count_families(self) -> int:
        """Count the number of human gene families in the database."""
        return self._count_model(GeneFamily)

    def count_mouse_genes(self) -> int:
        """Count the number of mouse genes in the database."""
        return self._count_model(MouseGene)

    def count_rat_genes(self) -> int:
        """Count the number of rat genes in the database."""
        return self._count_model(RatGene)

    def count_proteins(self) -> int:
        """Count the number of UniProt proteins in the database."""
        return self._count_model(UniProt)

    def summarize(self) -> Mapping[str, int]:
        """Summarize the database."""
        return dict(
            human_genes=self.count_human_genes(),
            rat_genes=self.count_rat_genes(),
            mouse_genes=self.count_mouse_genes(),
            families=self.count_families(),
            proteins=self.count_proteins(),
        )

    def get_gene_by_hgnc_symbol(self, hgnc_symbol: str) -> Optional[HumanGene]:
        """Get a human gene by HGNC symbol.

        :param hgnc_symbol: The HGNC gene symbol
        """
        results = self.hgnc(symbol=hgnc_symbol)
        return _deal_with_nonsense(results)

    def get_gene_by_hgnc_id(self, hgnc_id: str) -> Optional[HumanGene]:
        """Get a human gene by HGNC identifier.

        :param hgnc_id: The HGNC gene identifier
        """
        results = self.hgnc(identifier=int(hgnc_id))  # it's actually an integer in the database
        return _deal_with_nonsense(results)

    def get_gene_by_entrez_id(self, entrez_id: str) -> Optional[HumanGene]:
        """Get a human gene by its Entrez gene identifier.

        :param entrez_id: The Entrez gene identifier
        """
        results = self.hgnc(entrez=entrez_id)
        return _deal_with_nonsense(results)

    def get_gene_by_ensembl_id(self, ensembl_id: str) -> Optional[HumanGene]:
        """Get a human gene by its ENSEMBL gene identifier.

        :param ensembl_id: The ENSEMBL gene identifier
        """
        results = self.hgnc(ensembl_gene=ensembl_id)
        return _deal_with_nonsense(results)

    def get_gene_by_uniprot_id(self, uniprot_id: str) -> Optional[HumanGene]:
        """Get a human gene by its UniProt gene identifier.

        :param uniprot_id: The UniProt gene identifier
        """
        return self.hgnc(uniprotid=uniprot_id)

    def get_gene_by_mgi_id(self, mgi_id: str) -> Optional[HumanGene]:
        """Get a human gene by an orthologous MGI identifier.

        :param mgi_id: MGI identifier
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

    def get_gene_by_rgd_id(self, rgd_id: str) -> Optional[HumanGene]:
        """Get a human gene by an orthologous RGD identifier.

        :param rgd_id: RGD identifier
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

    def get_enzyme_by_ec_number(self, ec_number: str) -> Optional[HumanGene]:
        """Get a enzyme by its associated EC number.

        :param ec_number: EC number
        """
        return self.session.query(Enzyme).filter(Enzyme.ec_number == ec_number).one_or_none()

    def get_hgnc_from_alias_symbol(self, alias_symbol: str) -> Optional[HumanGene]:
        """Get HGNC from alias symbol.

        :param alias_symbol: alias symbol
        """
        query_result = self.session.query(AliasSymbol).filter(AliasSymbol.alias_symbol == alias_symbol).all()

        if not query_result:
            return None

        return query_result[0].hgnc

    def get_node(self, node: BaseEntity) -> Optional[HumanGene]:
        """Get a node from the database, whether it has a HGNC, RGD, MGI, or EG identifier.

        :param node: The node to look for
        :raises: KeyError
        """
        if NAMESPACE not in node:
            return

        namespace = node[NAMESPACE]
        identifier = node.get(IDENTIFIER)
        name = node.get(NAME)

        if namespace.lower() in {'hgnc'}:
            return self._get_node_handle_hgnc(identifier, name)

        if namespace.lower() in {'entrez', 'egid', 'eg', 'ncbigene'}:
            return self._get_node_handle_entrez(identifier, name)

        # if namespace.lower() in {'mgi'}:
        #     return self._get_node_handle_mgi(identifier, name)
        #
        # if namespace.lower() == 'mgiid':
        #     return self._get_node_handle_mgiid(identifier, name)
        #
        # if namespace.lower() in {'rgd'}:
        #     return self._get_node_handle_rgd(identifier, name)
        #
        # if namespace.lower() in {'rgdid'}:
        #     return self._get_node_handle_rgdid(identifier, name)

    def _get_node_handle_hgnc(self, identifier, name) -> Optional[HumanGene]:
        if identifier is not None:
            return self.get_gene_by_hgnc_id(identifier)
        elif name is not None:
            return self.get_gene_by_hgnc_symbol(name)
        raise KeyError

    def _get_node_handle_entrez(self, identifier, name) -> Optional[HumanGene]:
        if identifier is not None:
            return self.get_gene_by_entrez_id(identifier)
        elif name is not None:
            return self.get_gene_by_entrez_id(name)
        raise KeyError

    def _get_node_handle_mgi(self, identifier, name) -> Optional[HumanGene]:
        if identifier is not None:
            return self.get_gene_by_mgi_id(identifier)
        elif name is not None:
            return
        raise KeyError

    def _get_node_handle_rgd(self, identifier, name) -> Optional[HumanGene]:
        if identifier is not None:
            return self.get_gene_by_rgd_id(identifier)
        elif name is not None:
            return
        raise KeyError

    def _get_node_handle_mgiid(self, _, name):
        if name is None:
            raise KeyError
        return self.get_gene_by_mgi_id(name)

    def _get_node_handle_rgdid(self, _, name):
        if name is None:
            raise KeyError
        return self.get_gene_by_rgd_id(name)

    def add_namespace_to_graph(self, graph: BELGraph):
        """Add this manager's namespace to the graph."""
        namespace = self.upload_bel_namespace()
        graph.namespace_url[namespace.keyword] = namespace.url

        gfam_manager = GfamManager(engine=self.engine, session=self.session)
        gfam_manager.add_namespace_to_graph(graph)

    def iter_genes(self, graph: BELGraph, use_tqdm: bool = False) -> Iterable[Tuple[BaseEntity, HumanGene]]:
        """Iterate over pairs of BEL nodes and HGNC genes."""
        if use_tqdm:
            it = tqdm(graph, desc='HGNC Genes')
        else:
            it = graph

        for node in it:
            human_gene = self.get_node(node)
            if human_gene is not None:
                yield node, human_gene

    def normalize_genes(self, graph: BELGraph, use_tqdm: bool = False) -> None:
        """Add identifiers to all HGNC genes."""
        mapping = {}
        for node_data, human_gene in list(self.iter_genes(graph, use_tqdm=use_tqdm)):
            mapping[node_data] = gene_to_bel(human_gene, func=node_data.function, variants=node_data.get(VARIANTS))

        # FIXME what about when an HGNC node appears in a fusion, complex, or composite?

        relabel_nodes(graph, mapping, copy=False)

    def enrich_genes_with_equivalences(self, graph: BELGraph) -> None:
        """Enrich genes with their corresponding UniProt."""
        self.add_namespace_to_graph(graph)

        if 'uniprot' not in graph.namespace_url:
            graph.namespace_pattern['uniprot'] = UNIPROT_RE

        for node, human_gene in list(self.iter_genes(graph)):
            func = node.function

            if human_gene.entrez:
                graph.add_equivalence(node, FUNC_TO_DSL[func](
                    namespace=ENTREZ,
                    name=human_gene.symbol,
                    identifier=str(human_gene.entrez)
                ))

            if func in {PROTEIN, RNA, GENE}:
                for uniprot in human_gene.uniprots:
                    graph.add_equivalence(node, uniprot_to_bel(uniprot))

            if func in {RNA, GENE}:
                if human_gene.mirbase:
                    mirbase_rna_node = rna_dsl(
                        namespace='mirbase',
                        identifier=str(human_gene.mirbase)
                    )
                    graph.add_equivalence(node, mirbase_rna_node)

    def enrich_genes_with_families(self, graph: BELGraph) -> None:
        """Enrich genes in the BEL graph with their families."""
        self.add_namespace_to_graph(graph)
        for node_data, human_gene in list(self.iter_genes(graph)):
            for family in human_gene.gene_families:
                graph.add_is_a(node_data, family_to_bel(family, node_data[FUNCTION]))

    def get_family_by_id(self, family_identifier: str) -> Optional[GeneFamily]:
        """Get a gene family by its hgnc.genefamily identifier, if it exists."""
        results = self.gene_family(family_identifier=family_identifier)
        return _deal_with_nonsense(results)

    def get_family_by_name(self, family_name: str) -> Optional[GeneFamily]:
        """Get a gene family by its name, if it exists."""
        results = self.gene_family(family_name=family_name)
        return _deal_with_nonsense(results)

    def _enrich_hgnc_with_entrez_equivalences(self, graph: BELGraph, data: BaseEntity) -> Optional[str]:
        namespace = data.get(NAMESPACE)

        if namespace.lower() not in {'hgnc'}:
            return

        func = data[FUNCTION]
        name = data[NAME]
        entrez = self.hgnc_symbol_entrez_id_mapping[name]

        return graph.add_equivalence(data, FUNC_TO_DSL[func](
            namespace=ENTREZ,
            name=name,
            identifier=str(entrez),
        ))

    def enrich_hgnc_with_entrez_equivalences(self, graph: BELGraph):
        """Add equivalent Entrez nodes for all HGNC genes."""
        self.add_namespace_to_graph(graph)

        for _, data in graph.nodes(data=True):
            self._enrich_hgnc_with_entrez_equivalences(graph, data)

    def enrich_families_with_genes(self, graph: BELGraph):
        """Enrich gene families in the BEL graph with their member genes."""
        self.add_namespace_to_graph(graph)

        for gene_family_node in list(graph):
            if gene_family_node[FUNCTION] != GENE:
                continue

            if gene_family_node.get(NAMESPACE).lower() not in GENE_FAMILY_NAMESPACES:
                continue

            if IDENTIFIER in gene_family_node:
                gene_family_model = self.get_family_by_id(gene_family_node[IDENTIFIER])
            elif NAME in gene_family_node:
                gene_family_model = self.get_family_by_name(gene_family_node[NAME])
            else:
                raise ValueError

            if gene_family_model is None:
                log.info('family not found: %s', gene_family_node)
                continue

            for human_gene in gene_family_model.hgncs:
                graph.add_is_a(gene_to_bel(human_gene), gene_family_node)

    """Mapping dictionaries"""

    @staticmethod
    def _get_name(human_gene: HumanGene) -> str:
        return str(human_gene.symbol)

    @staticmethod
    def _get_identifier(human_gene: HumanGene) -> str:
        """Get the identifier from a human gene SQLAlchemy model."""
        return str(human_gene.identifier)

    @staticmethod
    def _get_encoding(human_gene: HumanGene) -> str:
        """Get the BEL encoding for a Human gene."""
        return ENCODINGS.get(human_gene.locus_type, 'GRP')

    def build_entrez_id_to_hgnc_symbol_mapping(self) -> Mapping[str, str]:
        """Build a mapping from Entrez gene identifier to HGNC gene symbols."""
        return {
            str(entrez_id): hgnc_symbol
            for entrez_id, hgnc_symbol in self.session.query(HumanGene.entrez, HumanGene.symbol).all()
        }

    def build_entrez_id_to_hgnc_id_mapping(self) -> Mapping[str, str]:
        """Build a mapping from Entrez gene identifier to HGNC identifier."""
        return {
            str(entrez_id): str(hgnc_id)
            for entrez_id, hgnc_id in self.session.query(HumanGene.entrez, HumanGene.identifier).all()
        }

    @property
    def hgnc_symbol_entrez_id_mapping(self) -> Mapping[str, str]:
        """Get a mapping from Entrez gene identifiers to HGNC gene symbols."""
        if not self._hgnc_symbol_entrez_id_mapping:
            self._hgnc_symbol_entrez_id_mapping = self.build_hgnc_symbol_entrez_id_mapping()

        return self._hgnc_symbol_entrez_id_mapping

    def build_hgnc_symbol_entrez_id_mapping(self) -> Mapping[str, str]:
        """Build a mapping from HGNC symbol to ENTREZ identifier."""
        return {
            hgnc_symbol: str(entrez_id)
            for hgnc_symbol, entrez_id in self.session.query(HumanGene.symbol, HumanGene.entrez).all()
        }

    def build_hgnc_id_symbol_mapping(self) -> Mapping[str, str]:
        """Build a mapping from HGNC identifiers to HGNC gene symbols."""
        return {
            str(hgnc_id): hgnc_symbol
            for hgnc_id, hgnc_symbol in self.session.query(HumanGene.identifier, HumanGene.symbol).all()
        }

    def build_hgnc_symbol_id_mapping(self) -> Mapping[str, str]:
        """Build a mapping from HGNC gene symbols to HGNC identifiers."""
        return {
            hgnc_symbol: str(hgnc_id)
            for hgnc_symbol, hgnc_id in self.session.query(HumanGene.symbol, HumanGene.identifier).all()
        }

    def build_uniprot_id_hgnc_id_mapping(self) -> Mapping[str, str]:
        """Build a mapping from UniProt identifiers to HGNC identifiers."""
        return {
            uniprot_id: str(hgnc_id)
            for hgnc_id, uniprot_id in self.session.query(HumanGene.identifier, UniProt.uniprotid).all()
        }

    def build_uniprot_id_hgnc_symbol_mapping(self) -> Mapping[str, str]:
        """Build a mapping from UniProt identifiers to HGNC gene symbols."""
        return {
            uniprot_id: symbol
            for symbol, uniprot_id in self.session.query(HumanGene.symbol, UniProt.uniprotid).all()
        }

    def get_all_hgnc_symbols(self) -> Set[str]:
        """Return the set of HGNC gene symbols in the database."""
        return {
            gene.symbol
            for gene in self.session.query(HumanGene.symbol).all()
        }

    def _get_gene_encodings(self) -> Mapping[str, str]:
        """Get the name to encoding dictionary for HGNC gene names."""
        return {
            symbol: ENCODINGS.get(locus_type, 'GRP')
            for symbol, locus_type in self.session.query(HumanGene.symbol, HumanGene.locus_type).all()
        }

    def list_families(self) -> List[GeneFamily]:
        """List families in the database."""
        return self._list_model(GeneFamily)

    def list_human_genes(self) -> List[HumanGene]:
        """List human genes in the database."""
        return self._list_model(HumanGene)

    def to_bel(self) -> BELGraph:
        """Export gene family definitions as a BEL graph."""
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
            gene_dsl = gene_to_bel(human_gene)
            graph.add_node_from_data(gene_dsl)
            add_central_dogma(graph, human_gene)

            for mouse_gene in human_gene.mgds:
                graph.add_orthology(gene_dsl, pybel.dsl.Gene('mgi', identifier=str(mouse_gene.mgdid)))
            for rat_gene in human_gene.rgds:
                graph.add_orthology(gene_dsl, pybel.dsl.Gene('rgd', identifier=str(rat_gene.rgdid)))

            protein_dsl = gene_to_bel(human_gene, PROTEIN)
            for enzyme in human_gene.enzymes:
                graph.add_is_a(protein_dsl, pybel.dsl.Protein('ec-code', enzyme.ec_number))

        for family in tqdm(self.list_families(), total=self.count_families(),
                           desc='Mapping gene family definitions to BEL'):
            for human_gene in family.hgncs:
                graph.add_is_a(gene_to_bel(human_gene), family_to_bel(family))

        return graph

    def get_pathway_size_distribution(self) -> Counter:
        """Get the pathway size distribution."""
        return Counter({
            family.family_name: len(family.hgncs)
            for family in self.session.query(GeneFamily)
        })

    def get_all_hgnc_symbols_family(self) -> Set[str]:
        """Get all Gene symbols that appear in gene families."""
        return {
            human_gene.symbol
            for family in self.session.query(GeneFamily)
            for human_gene in family.hgncs
        }

    def add_central_dogma(self, graph: BELGraph, node: BaseEntity) -> Optional[CentralDogma]:
        """Add the central dogma of biology."""
        if VARIANTS in node:
            return

        human_gene = self.get_node(node)
        encoding = ENCODINGS.get(human_gene.locus_type, 'GRP')

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

    def _create_namespace_entry_from_model(self, human_gene: HumanGene, namespace: Namespace) -> NamespaceEntry:
        return NamespaceEntry(
            encoding=self._get_encoding(human_gene),
            identifier=human_gene.identifier,
            name=human_gene.symbol,
            namespace=namespace
        )

    @staticmethod
    def _cli_add_populate(main: click.Group) -> click.Group:  # noqa: D202
        """Override default method to make it possible to add more flags."""

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
