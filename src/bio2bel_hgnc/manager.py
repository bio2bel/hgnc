# -*- coding: utf-8 -*-

"""Bio2BEL HGNC Manager."""

import logging
from typing import Iterable, List, Mapping, Optional, Set, Tuple, TypeVar

from networkx import relabel_nodes
from tqdm import tqdm

import pybel.dsl
import pyobo
from bio2bel.compath import CompathManager
from pybel import BELGraph
from pybel.constants import FUNCTION, MIRNA, NAME, NAMESPACE, PROTEIN, RNA, VARIANTS
from pybel.dsl import BaseEntity, CentralDogma, FUNC_TO_DSL
from pybel.manager.models import Namespace, NamespaceEntry
from .constants import ENCODINGS, ENTREZ, MODULE_NAME
from .gfam_manager import Manager as GfamManager
from .model_utils import add_central_dogma, family_to_bel, gene_to_bel
from .models import Base, GeneFamily, HumanGene, MouseGene, RatGene, human_mouse, human_rat

__all__ = [
    'Manager',
]

logger = logging.getLogger(__name__)

UNIPROT_RE = r'^([A-N,R-Z][0-9]([A-Z][A-Z, 0-9][A-Z, 0-9][0-9]){1,2})|([O,P,Q][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9])(\.\d+)?$'
GENE_FAMILY_NAMESPACES = {'gfam', 'hgnc.family', 'hgnc.genefamily'}

X = TypeVar('X')


class Manager(CompathManager):
    """Human gene nomenclature and orthologies to mouse and rat."""

    module_name = MODULE_NAME
    _base = Base

    protein_model = namespace_model = HumanGene
    pathway_model = GeneFamily
    edge_model = [human_mouse, human_rat]
    identifiers_recommended = 'HGNC'
    identifiers_pattern = r'^((HGNC|hgnc):)?\d{1,5}$'
    identifiers_miriam = 'MIR:00000080'
    identifiers_namespace = 'hgnc'
    identifiers_url = 'http://identifiers.org/hgnc/'

    flask_admin_models = [HumanGene, GeneFamily]

    def __init__(self, *args, **kwargs):  # noqa: D107
        super().__init__(*args, **kwargs)
        self._hgnc_symbol_entrez_id_mapping = {}

    def is_populated(self) -> bool:
        """Check if the database is already populated."""
        return 0 < self.count_human_genes()

    def populate(self):  # noqa: C901
        """Populate the database."""
        gfam = pyobo.get('hgnc.genefamily')
        gfam_id_to_model = {
            term.identifier: GeneFamily(identifier=term.identifier, name=term.name)
            for term in gfam
        }
        mgi = pyobo.get('mgi')
        mgi_id_to_model = {
            term.identifier: MouseGene(mgi_id=term.identifier, mgi_symbol=term.name)
            for term in mgi
        }
        rgd = pyobo.get('rgd')
        rgd_id_to_model = {
            term.identifier: RatGene(rgd_id=term.identifier, rgd_symbol=term.name)
            for term in rgd
        }
        hgnc = pyobo.get('hgnc')
        hgnc_id_to_entrez_id = hgnc.get_filtered_xrefs_mapping('ncbigene')
        for term in hgnc:
            try:
                entrez_id = hgnc_id_to_entrez_id[term.identifier]
            except KeyError:
                logger.warning('could not map hgnc:%s to entrez', term.identifier)
                continue

            rat_genes = []
            mouse_genes = []
            for xref in term.xrefs:
                if xref.prefix == 'rgd':
                    try:
                        rat_gene = rgd_id_to_model[xref.identifier]
                    except KeyError:
                        logger.warning('could not map hgnc:%s to rgd:%s', term.identifier, xref.identifier)
                    else:
                        rat_genes.append(rat_gene)
                elif xref.prefix == 'mgi':
                    try:
                        mouse_gene = mgi_id_to_model[xref.identifier]
                    except KeyError:
                        logger.warning('could not map hgnc:%s to mgi:%s', term.identifier, xref.identifier)
                    else:
                        mouse_genes.append(mouse_gene)

            model = HumanGene(
                hgnc_id=term.identifier,
                hgnc_symbol=term.name,
                entrez_id=entrez_id,
                gene_families=[
                    gfam_id_to_model[parent.identifier]
                    for parent in term.parents
                ],
                rat_genes=rat_genes,
                mouse_genes=mouse_genes,
            )
            self.session.add(model)
        self.session.commit()

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

    def summarize(self) -> Mapping[str, int]:
        """Summarize the database."""
        return dict(
            human_genes=self.count_human_genes(),
            rat_genes=self.count_rat_genes(),
            mouse_genes=self.count_mouse_genes(),
            families=self.count_families(),
        )

    def get_gene_by_hgnc_symbol(self, hgnc_symbol: str) -> Optional[HumanGene]:
        """Get a human gene by HGNC symbol.

        :param hgnc_symbol: The HGNC gene symbol
        """
        return self.session.query(HumanGene).filter(HumanGene.hgnc_symbol == hgnc_symbol).one_or_none()

    def get_gene_by_hgnc_id(self, hgnc_id: str) -> Optional[HumanGene]:
        """Get a human gene by HGNC identifier.

        :param hgnc_id: The HGNC gene identifier
        """
        return self.session.query(HumanGene).filter(HumanGene.hgnc_id == hgnc_id).one_or_none()

    def get_gene_by_entrez_id(self, entrez_id: str) -> Optional[HumanGene]:
        """Get a human gene by its Entrez gene identifier.

        :param entrez_id: The Entrez gene identifier
        """
        return self.session.query(HumanGene).filter(HumanGene.entrez_id == entrez_id).one_or_none()

    def get_gene_by_mgi_id(self, mgi_id: str) -> Optional[HumanGene]:
        """Get a human gene by an orthologous MGI identifier.

        :param mgi_id: MGI identifier
        """
        mouse_gene = self.session.query(MouseGene).filter(MouseGene.mgi_id == mgi_id).one_or_none()
        if mouse_gene is None:
            return

        human_genes = mouse_gene.human_genes

        if len(human_genes) > 1:
            logger.warning('multiple human genes mapped to mgi_id:%s: %s', mgi_id, human_genes)
            return

        return human_genes[0]

    def get_gene_by_rgd_id(self, rgd_id: str) -> Optional[HumanGene]:
        """Get a human gene by an orthologous RGD identifier.

        :param rgd_id: RGD identifier
        """
        mouse_gene = self.session.query(RatGene).filter(RatGene.rgd_id == rgd_id).one_or_none()
        if mouse_gene is None:
            return

        human_genes = mouse_gene.human_genes

        if len(human_genes) > 1:
            logger.warning('multiple human genes mapped to rgd_id:%s: %s', rgd_id, human_genes)
            return

        return human_genes[0]

    def get_node(self, node: BaseEntity) -> Optional[HumanGene]:
        """Get a node from the database.

        :param node: The node to look for
        :raises: KeyError
        """
        if not isinstance(node, CentralDogma):
            return

        namespace = node.namespace
        if namespace is None:
            return

        identifier = node.identifier
        name = node.name

        if namespace.lower() in {'hgnc'}:
            return self._get_node_handle_hgnc(identifier, name)

    def _get_node_handle_hgnc(self, identifier, name) -> Optional[HumanGene]:
        if identifier is not None:
            return self.get_gene_by_hgnc_id(identifier)
        elif name is not None:
            return self.get_gene_by_hgnc_symbol(name)
        raise KeyError

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
        for node, human_gene in list(self.iter_genes(graph)):
            func = node.function

            if human_gene.entrez_id:
                entrez_node = FUNC_TO_DSL[func](
                    namespace=ENTREZ,
                    name=human_gene.hgnc_symbol,
                    identifier=str(human_gene.entrez_id)
                )
                graph.add_equivalence(node, entrez_node)

            # if func == PROTEIN:
            #     for uniprot in human_gene.uniprots:
            #         graph.add_equivalence(node, uniprot_to_bel(uniprot))

            # if func == RNA:
            #     if human_gene.mirbase:
            #         mirbase_rna_node = pybel.dsl.Rna(
            #             namespace='mirbase',
            #             identifier=str(human_gene.mirbase),
            #         )
            #         graph.add_equivalence(node, mirbase_rna_node)

    def enrich_genes_with_families(self, graph: BELGraph) -> None:
        """Enrich genes in the BEL graph with their families."""
        self.add_namespace_to_graph(graph)
        for node_data, human_gene in list(self.iter_genes(graph)):
            for family in human_gene.gene_families:
                graph.add_is_a(node_data, family_to_bel(family, node_data[FUNCTION]))

    def get_family_by_id(self, hgnc_genefamily_id: str) -> Optional[GeneFamily]:
        """Get a gene family by its hgnc.genefamily identifier, if it exists."""
        return self.session.query(GeneFamily).filter(GeneFamily.identifier == hgnc_genefamily_id).one_or_none()

    def get_family_by_name(self, family_name: str) -> Optional[GeneFamily]:
        """Get a gene family by its name, if it exists."""
        return self.session.query(GeneFamily).filter(GeneFamily.name == family_name).one_or_none()

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
            if not isinstance(gene_family_node, pybel.dsl.Gene):
                continue

            namespace = gene_family_node.namespace
            if namespace is None or namespace.lower() not in GENE_FAMILY_NAMESPACES:
                continue

            identifier = gene_family_node.identifier
            name = gene_family_node.name

            if identifier:
                gene_family_model = self.get_family_by_id(identifier)
            elif name:
                gene_family_model = self.get_family_by_name(name)
            else:
                raise ValueError

            if gene_family_model is None:
                logger.info('family not found: %s', gene_family_node)
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
            contact='cthoyt@gmail.com',
        )

        hgnc_namespace = self.upload_bel_namespace()
        logger.info('using default namespace: %s at %s', hgnc_namespace, hgnc_namespace.url)
        graph.namespace_url[hgnc_namespace.keyword] = hgnc_namespace.url

        gfam_manager = GfamManager(connection=self.connection)
        gfam_namespace = gfam_manager.upload_bel_namespace()
        logger.info('using default namespace: %s at %s', gfam_namespace, gfam_namespace.url)
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
