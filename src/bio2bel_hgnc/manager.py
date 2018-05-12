# -*- coding: utf-8 -*-

"""Bio2BEL HGNC Manager."""

import logging
from collections import Counter

import click
from tqdm import tqdm

from bio2bel.namespace_manager import NamespaceManagerMixin
from pybel import BELGraph
from pybel.constants import FUNCTION, GENE, IDENTIFIER, NAME, NAMESPACE, NAMESPACE_DOMAIN_GENE, PROTEIN, RNA
from pybel.dsl import gene as gene_dsl, protein as protein_dsl, rna as rna_dsl
from pybel.manager.models import NamespaceEntry
from .constants import GENE_FAMILY_KEYWORD, MODULE_NAME, encodings
from .gfam_manager import Manager as GfamManager
from .model_utils import *
from .models import Base, GeneFamily, HumanGene, MouseGene, RatGene, UniProt
from .wrapper import BaseManager

log = logging.getLogger(__name__)

__all__ = [
    'Manager',
]

_func_to_dsl = {
    GENE: gene_dsl,
    RNA: rna_dsl,
    PROTEIN: protein_dsl
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


def _write_gene_bel_namespace_helper(values, file):
    """
    :param dict[str,str] values:
    :param file:
    """
    from pybel.resources import write_namespace
    write_namespace(
        namespace_name='HGNC Human Gene Names',
        namespace_keyword='HGNC',
        namespace_domain=NAMESPACE_DOMAIN_GENE,
        author_name='Charles Tapley Hoyt',
        citation_name='HGNC?',
        values=values,
        functions='G',
        file=file
    )


def _write_gene_families_bel_namespace_helper(values, file):
    """
    :param list[str] values:
    :param file:
    """
    from pybel.resources import write_namespace
    write_namespace(
        namespace_name='HGNC Gene Families',
        namespace_keyword='GFAM',
        namespace_domain=NAMESPACE_DOMAIN_GENE,
        author_name='Charles Tapley Hoyt',
        citation_name='HGNC?',
        values=values,
        functions='G',
        file=file
    )


class Manager(NamespaceManagerMixin, BaseManager):
    """Bio2BEL HGNC Manager"""

    module_name = MODULE_NAME
    flask_admin_models = [HumanGene, GeneFamily, UniProt, MouseGene, RatGene]
    namespace_model = HumanGene

    def __init__(self, connection=None):
        super().__init__(connection=connection)

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

    def count_human_genes(self):
        return self._count_model(HumanGene)

    def count_families(self):
        return self._count_model(GeneFamily)

    def count_mouse_genes(self):
        return self._count_model(MouseGene)

    def count_rat_genes(self):
        return self._count_model(RatGene)

    def count_uniprots(self):
        return self._count_model(UniProt)

    def summarize(self):
        """Returns a summary dictionary over the content of the database

        :rtype: dict[str,int]
        """
        return dict(
            human_genes=self.count_human_genes(),
            rat_genes=self.count_rat_genes(),
            mouse_genes=self.count_mouse_genes(),
            families=self.count_families(),
            uniprots=self.count_uniprots()
        )

    def get_gene_by_hgnc_symbol(self, hgnc_symbol):
        """Gets a gene by the symbol

        :param str hgnc_symbol: The HGNC gene symbol
        :rtype: Optional[bio2bel_hgnc.models.HumanGene]
        """
        results = self.hgnc(symbol=hgnc_symbol)
        return _deal_with_nonsense(results)

    def get_gene_by_hgnc_id(self, hgnc_id):
        """Gets a gene by the identifier

        :param str hgnc_id: The HGNC gene identifier
        :rtype: Optional[bio2bel_hgnc.models.HumanGene]
        """
        results = self.hgnc(identifier=hgnc_id)
        return _deal_with_nonsense(results)

    def get_gene_by_entrez_id(self, entrez_id):
        """Gets a HGNC gene by its Entrez gene identifier

        :param str entrez_id: The Entrez gene identifier
        :rtype: Optional[bio2bel_hgnc.models.HumanGene]
        """
        results = self.hgnc(entrez=entrez_id)
        return _deal_with_nonsense(results)

    def get_gene_by_mgi_id(self, mgi_id):
        """Gets a HGNC gene by an orthologous MGI identifier

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

    def get_gene_by_mgi_symbol(self, mgi_symbol):
        """Gets a HGNC gene by an orthologous MGI gene symbol

        :param str mgi_symbol: MGI gene symbol
        :rtype: Optional[bio2bel_hgnc.models.HumanGene]
        """
        # TODO this data is not in HGNC...
        raise NotImplementedError

    def get_gene_by_rgd_id(self, rgd_id):
        """Gets a HGNC gene by an orthologous RGD identifier

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

    def get_gene_by_rgd_symbol(self, rgd_symbol):
        """Gets a HGNC gene by an orthologous RGD identifier

        :param str rgd_symbol: RGD gene symbol
        :rtype: Optional[bio2bel_hgnc.models.HumanGene]
        """
        # TODO this data is not in HGNC...
        raise NotImplementedError

    def get_node(self, graph, node):
        """Gets a node from the database, whether it has a HGNC, RGD, MGI, or EG identifier.

        :param pybel.BELGraph graph: A BEL graph
        :param tuple node: A PyBEL node tuple
        :rtype: bio2bel.models.HGNC
        :raises: KeyError
        """
        data = graph.node[node]

        if NAMESPACE not in data:
            raise KeyError

        namespace = data[NAMESPACE]
        identifer = data.get(IDENTIFIER)
        name = data.get(NAME)

        if namespace == 'HGNC':
            if identifer is not None:
                return self.get_gene_by_hgnc_id(identifer)
            elif name is not None:
                return self.get_gene_by_hgnc_symbol(name)
            raise KeyError

        if namespace in {'ENTREZ', 'EGID', 'EG'}:
            if identifer is not None:
                return self.get_gene_by_entrez_id(identifer)
            elif name is not None:
                return self.get_gene_by_entrez_id(name)
            raise KeyError

        if namespace in {'MGI'}:
            if identifer is not None:
                return self.get_gene_by_mgi_id(identifer)
            elif name is not None:
                return self.get_gene_by_mgi_symbol(name)
            raise KeyError

        if namespace == 'MGIID':
            if name is None:
                raise KeyError
            return self.get_gene_by_mgi_id(name)

        if namespace in {'RGD'}:
            if identifer is not None:
                return self.get_gene_by_rgd_id(identifer)
            elif name is not None:
                return self.get_gene_by_rgd_symbol(name)
            raise KeyError

        if namespace == 'RGDID':
            if name is None:
                raise KeyError
            return self.get_gene_by_rgd_id(name)

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

            if m is None:
                log.info('gene not found: %s', data)
                continue

            for family in m.gene_families:
                graph.add_is_a(n, family_to_bel(family))

    def get_family_by_id(self, family_identifier):
        """Gets a gene family by its identifier

        :param str family_identifier: The identifier of a HGNC Gene Family
        :rtype: Optional[bio2bel_hgnc.models.GeneFamily]
        """
        results = self.gene_family(family_identifier=family_identifier)
        return _deal_with_nonsense(results)

    def get_family_by_name(self, family_name):
        """Gets a gene family by its name

        :param str family_name: The name of a HGNC Gene Family
        :rtype: Optional[bio2bel_hgnc.models.GeneFamily]
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
            namespace='ENTREZ',
            name=name,
            identifier=str(entrez)
        ))

    def enrich_hgnc_with_entrez_equivalences(self, graph):
        """For all HGNC genes, adds their Entrez equivalent nodes

        :param pybel.BELGraph graph: The BEL graph to enrich
        """
        for node in graph.nodes():
            self._enrich_hgnc_with_entrez_equivalences(graph, node)

    def enrich_families_with_genes(self, graph):
        """Enrich gene families in the BEL graph with their member genes

        :param pybel.BELGraph graph: The BEL graph to enrich
        """
        for gene_family_node, data in graph.nodes(data=True):
            if data[FUNCTION] != GENE:
                continue

            if data.get(NAMESPACE) != GENE_FAMILY_KEYWORD:
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

            for gene in gene_family_model.hgncs:
                graph.add_is_a(gene_to_bel(gene), gene_family_node)

    """ Mapping dictionaries"""

    def _get_identifier(self, human_gene):
        return human_gene.identifier

    def build_entrez_id_symbol_mapping(self):
        """Builds a mapping from ENTREZ identifier to HGNC symbol

        :rtype: dict[str,str]
        """
        return {
            identifier: symbol
            for identifier, symbol in self.session.query(HumanGene.entrez, HumanGene.symbol).all()
        }

    @property
    def hgnc_symbol_entrez_id_mapping(self):
        if not self._hgnc_symbol_entrez_id_mapping:
            self._hgnc_symbol_entrez_id_mapping = self.build_hgnc_symbol_entrez_id_mapping()

        return self._hgnc_symbol_entrez_id_mapping

    def build_hgnc_symbol_entrez_id_mapping(self):
        """Builds a mapping from HGNC symbol to ENTREZ identifier

        :rtype: dict[str,str]
        """
        return {
            symbol: identifier
            for symbol, identifier in self.session.query(HumanGene.symbol, HumanGene.entrez).all()
        }

    def build_hgnc_id_symbol_mapping(self):
        """Builds a mapping from HGNC identifier to HGNC symbol

        :rtype: dict[str,str]
        """
        return {
            str(identifier): symbol
            for identifier, symbol in self.session.query(HumanGene.identifier, HumanGene.symbol).all()
        }

    def build_hgnc_symbol_id_mapping(self):
        """Builds a mapping from HGNC symbol to HGNC identifier

        :rtype: dict[str,str]
        """
        return {
            symbol: str(identifier)
            for symbol, identifier in self.session.query(HumanGene.symbol, HumanGene.identifier).all()
        }

    def build_hgnc_symbol_uniprot_ids_mapping(self):
        """Builds mapping from HGNC symbol to UniProt identifiers

        :rtype: dict[str,set[str]]
        """
        return {
            symbol: uniprot_id
            for symbol, uniprot_id in self.session.query(HumanGene.symbol, UniProt.uniprotid).all()
        }

    def build_hgnc_id_uniprot_ids_mapping(self):
        """Builds mapping from HGNC identifier to UniProt identifiers

        :rtype: dict[str,set[str]]
        """
        return {
            hgnc_id: uniprot_id
            for hgnc_id, uniprot_id in self.session.query(HumanGene.identifier, UniProt.uniprotid).all()
        }

    def build_uniprot_id_hgnc_id_mapping(self):
        """Builds mapping from UniProt identifiers to HGNC identifier

        :rtype: dict[str,str]
        """
        return {
            uniprot_id: hgnc_id
            for hgnc_id, uniprot_id in self.session.query(HumanGene.identifier, UniProt.uniprotid).all()
        }

    def build_uniprot_id_hgnc_symbol_mapping(self):
        """Builds mapping from UniProt identifier to HGNC symbol

        :rtype: dict[str,str]
        """
        return {
            uniprot_id: symbol
            for symbol, uniprot_id in self.session.query(HumanGene.symbol, UniProt.uniprotid).all()
        }

    def get_all_hgnc_symbols(self):
        """Returns the set of hgnc symbols in the database

        :rtype: set
        """
        return {
            _deal_with_nonsense(symbol)
            for symbol in self.session.query(HumanGene.symbol).all()
        }

    def _get_gene_encodings(self):
        """Gets the name to encoding dictionary for HGNC gene names

        :rtype: dict[str,str]
        """
        return {
            symbol: encodings.get(locus_type, 'GRP')
            for symbol, locus_type in self.session.query(HumanGene.symbol, HumanGene.locus_type).all()
        }

    def write_gene_bel_namespace(self, file):
        values = self._get_gene_encodings()
        _write_gene_bel_namespace_helper(values, file)

    def deploy_gene_bel_namespace(self):
        """Creates and deploys the Gene Names Namespace

        :rtype: Optional[str]
        """
        from pybel.resources.deploy import deploy_namespace
        from pybel.resources.arty import get_today_arty_namespace

        file_name = get_today_arty_namespace('hgnc')

        with open(file_name, 'w') as file:
            self.write_gene_bel_namespace(file)

        return deploy_namespace(file_name, module_name='hgnc')

    def write_gene_family_bel_namespace(self, file):
        values = [name for name, in self.session.query(GeneFamily.family_name).all()]
        _write_gene_families_bel_namespace_helper(values=values, file=file)

    def deploy_gene_family_bel_namespace(self):
        """Creates and deploys the Gene Families Namespace

        :rtype: Optional[str]
        """
        from pybel.resources.deploy import deploy_namespace
        from pybel.resources.arty import get_today_arty_namespace

        file_name = get_today_arty_namespace('gfam')

        with open(file_name, 'w') as file:
            self.write_gene_family_bel_namespace(file)

        return deploy_namespace(file_name, 'gfam')

    def list_families(self):
        """
        :rtype: list[GeneFamily]
        """
        return self._list_model(GeneFamily)

    def list_human_genes(self):
        """
        :rtype: list[HumanGene]
        """
        return self._list_model(HumanGene)

    def to_bel(self):
        """Export gene family definitions as a BEL graph

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

        for family in tqdm(self.list_families(), total=self.count_families(),
                           desc='Mapping gene family definitions to BEL'):
            for gene in family.hgncs:
                graph.add_is_a(
                    gene_to_bel(gene),
                    family_to_bel(family)
                )

        for human_gene in tqdm(self.list_human_genes(), total=self.count_human_genes(),
                               desc='Mapping central dogma to BEL'):
            encoding = encodings.get(human_gene.locus_type, 'GRP')

            if 'R' in encoding:
                rna = gene_to_rna_to_bel(human_gene)
                graph.add_transcription(gene_to_bel(human_gene), rna)

                if 'P' in encoding:
                    graph.add_translation(rna, gene_to_protein_to_bel(human_gene))

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
        """Gets all Gene symbols in gene families

        :rtype: set[str]
        """
        res = (
            gene.symbol
            for family in self.session.query(GeneFamily)
            for gene in family.hgncs
        )

        return set(res)

    def add_central_dogma(self, graph, node):
        """

        :param graph:
        :param node:
        :return:
        """
        if node not in graph:
            raise ValueError

        human_gene = self.get_node(graph, node)
        encoding = encodings.get(human_gene.locus_type, 'GRP')

        if 'R' in encoding:
            graph.add_unqualified_edge(node, )

    def add_node_equivalencies(self, graph, node, add_leaves=True):
        """Given an HGNC node, add equivalencies found in the database.

         - Entrez
         - UniProt
         - miRBase

        :param pybel.BELGraph graph: A BEL graph
        :param tuple node: A PyBEL node tuple
        :param bool add_leaves: Should equivalencies that are not already in the graph be added?
        """
        gene = self.get_node(graph, node)
        if gene is None:
            return

        for uniprot in gene.uniprots:
            graph.add_translation(
                gene_to_bel(gene),
                gene_to_rna_to_bel(gene)
            )
            graph.add_translation(
                gene_to_rna_to_bel(gene),
                gene_to_protein_to_bel(gene),
            )
            graph.add_equivalence(
                gene_to_protein_to_bel(gene),
                uniprot_to_pybel(uniprot)
            )

        if gene.entrez:
            graph.add_equivalence(
                node,
                gene_dsl(namespace='ENTREZ', identifier=str(gene.entrez))
            )

        if gene.mirbase:
            graph.add_translation(
                gene_to_bel(gene),
                gene_to_rna_to_bel(gene)
            )
            graph.add_equivalence(
                gene_to_rna_to_bel(gene),
                rna_dsl(namespace='MIRBASE', identifier=str(gene.entrez))
            )

    def _create_namespace_entry_from_model(self, human_gene, namespace):
        return NamespaceEntry(
            encoding=encodings.get(human_gene.locus_type, 'GRP'),
            identifier=human_gene.identifier,
            name=human_gene.symbol,
            namespace=namespace
        )

    @staticmethod
    def _cli_add_populate(main):
        """Overrides default method to make it possible to add more flags"""

        @main.command()
        @click.option('--reset', is_flag=True)
        @click.option('--skip-hcop', is_flag=True)
        @click.pass_obj
        def populate(manager, reset, skip_hcop):
            """Populates the database"""

            if reset:
                log.info('Deleting the previous instance of the database')
                manager.drop_all()
                log.info('Creating new models')
                manager.create_all()

            manager.populate(use_hcop=(not skip_hcop))

        return main
