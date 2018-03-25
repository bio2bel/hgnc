# -*- coding: utf-8 -*-

import logging
from collections import Counter

from pybel import BELGraph, to_bel
from pybel.constants import FUNCTION, GENE, IDENTIFIER, IS_A, NAME, NAMESPACE, NAMESPACE_DOMAIN_GENE, PROTEIN, RNA
from pybel.dsl import gene as gene_dsl, protein as protein_dsl, rna as rna_dsl
from pybel.resources import get_latest_arty_namespace, write_namespace
from pybel.resources.arty import get_today_arty_knowledge, get_today_arty_namespace
from pybel.resources.deploy import deploy_knowledge, deploy_namespace

from bio2bel.abstractmanager import AbstractManager
from pyhgnc.manager.database import DbManager
from pyhgnc.manager.models import Base
from pyhgnc.manager.query import QueryManager
from .constants import GENE_FAMILY_KEYWORD, MODULE_NAME, encodings
from .models import GeneFamily, HGNC, UniProt

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


def gene_to_pybel(gene):
    """Converts a PyHGNC Gene to a PyBEL gene

    :param pyhgnc.manager.models.HGNC gene:  A PyHGNC Gene model
    :rtype: pybel.dsl.gene
    """
    return gene_dsl(
        namespace='HGNC',
        name=str(gene.symbol),
        identifier=str(gene.identifier)
    )


def family_to_pybel(family):
    """Converts a PyHGNC Gene Family model to a PyBEL gene

    :param pyhgnc.manager.models.GeneFamily family: A PyHGNC Gene Family model
    :rtype: pybel.dsl.gene
    """
    return gene_dsl(
        namespace='GFAM',
        identifier=str(family.family_identifier),
        name=str(family.family_name)
    )


def _write_gene_bel_namespace_helper(values, file):
    """
    :param dict[str,str] values:
    :param file:
    """
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


class _PyHGNCManager(DbManager, QueryManager):
    pass


class Manager(AbstractManager, _PyHGNCManager):
    """An extended version of the PyHGNC manager to have useful functions"""

    module_name = MODULE_NAME

    @property
    def base(self):
        return Base

    def populate(self, silent=False, hgnc_file_path=None, hcop_file_path=None, low_memory=False):
        json_data = self.load_hgnc_json(hgnc_file_path=hgnc_file_path)
        self.insert_hgnc(hgnc_dict=json_data, silent=silent, low_memory=low_memory)
        self.insert_hcop(silent=silent, hcop_file_path=hcop_file_path)

    #: Clobber this PyHGNC function so it doesn't accidentally get called
    def db_import(self, silent=False, hgnc_file_path=None, hcop_file_path=None, low_memory=False):
        raise NotImplemented('call manager.populate instead')

    #: Clobber this PyHGNC function so it doesn't accidentally get called
    def _drop_tables(self):
        raise NotImplemented('call manager.drop_all instead')

    def count_genes(self):
        return self.session.query(HGNC).count()

    def count_families(self):
        return self.session.query(GeneFamily).count()

    def count_uniprots(self):
        return self.session.query(UniProt).count()

    def summarize(self):
        """Returns a summary dictionary over the content of the database

        :rtype: dict[str,int]
        """
        return dict(
            genes=self.count_genes(),
            families=self.count_families(),
            uniprots=self.count_uniprots()
        )

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
        """Gets a HGNC gene by its Entrez gene identifier

        :param str entrez_id: The Entrez gene identifier
        :rtype: Optional[pyhgnc.manager.models.HGNC]
        """
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
                graph.add_unqualified_edge(n, family_to_pybel(family), IS_A)

    def get_family_by_id(self, family_identifier):
        """Gets a gene family by its identifier

        :param str family_identifier: The identifier of a HGNC Gene Family
        :rtype: Optional[pyhgnc.manager.models.GeneFamily]
        """
        results = self.gene_family(family_identifier=family_identifier)
        return _deal_with_nonsense(results)

    def get_family_by_name(self, family_name):
        """Gets a gene family by its name

        :param str family_name: The name of a HGNC Gene Family
        :rtype: Optional[pyhgnc.manager.models.GeneFamily]
        """
        results = self.gene_family(family_name=family_name)
        return _deal_with_nonsense(results)

    def enrich_hgnc_with_entrez_equivalences(self, graph):
        """For all HGNC genes, adds their Entrez equivalent nodes

        :param pybel.BELGraph graph: The BEL graph to enrich
        """
        hgnc_symbol_to_entrez = self.build_hgnc_symbol_entrez_id_mapping()

        for gene_node, data in graph.nodes(data=True):
            namespace = data.get(NAMESPACE)

            if namespace != 'HGNC':
                continue

            func = data[FUNCTION]
            name = data[NAME]
            entrez = hgnc_symbol_to_entrez[name]

            graph.add_equivalence(gene_node, _func_to_dsl[func](
                namespace='ENTREZ',
                name=name,
                identifier=str(entrez)
            ))

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
                graph.add_is_a(gene_to_pybel(gene), gene_family_node)

    """ Mapping dictionaries"""

    def build_entrez_id_symbol_mapping(self):
        """Builds a mapping from ENTREZ identifier to HGNC symbol

        :rtype: dict[str,str]
        """
        return {
            identifier: symbol
            for identifier, symbol in self.session.query(HGNC.entrez, HGNC.symbol).all()
        }

    def build_hgnc_symbol_entrez_id_mapping(self):
        """Builds a mapping from HGNC symbol to ENTREZ identifier

        :rtype: dict[str,str]
        """
        return {
            symbol: identifier
            for symbol, identifier in self.session.query(HGNC.symbol, HGNC.entrez).all()
        }

    def build_hgnc_id_symbol_mapping(self):
        """Builds a mapping from HGNC identifier to HGNC symbol

        :rtype: dict[str,str]
        """
        return {
            str(identifier): symbol
            for identifier, symbol in self.session.query(HGNC.identifier, HGNC.symbol).all()
        }

    def build_hgnc_symbol_id_mapping(self):
        """Builds a mapping from HGNC symbol to HGNC identifier

        :rtype: dict[str,str]
        """
        return {
            symbol: str(identifier)
            for symbol, identifier in self.session.query(HGNC.symbol, HGNC.identifier).all()
        }

    def build_hgnc_symbol_uniprot_ids_mapping(self):
        """Builds mapping from HGNC symbol to UniProt identifiers

        :rtype: dict[str,set[str]]
        """
        return {
            symbol: uniprot_id
            for symbol, uniprot_id in self.session.query(HGNC.symbol, UniProt.uniprotid).all()
        }

    def build_hgnc_id_uniprot_ids_mapping(self):
        """Builds mapping from HGNC identifier to UniProt identifiers

        :rtype: dict[str,set[str]]
        """
        return {
            hgnc_id: uniprot_id
            for hgnc_id, uniprot_id in self.session.query(HGNC.identifier, UniProt.uniprotid).all()
        }

    def build_uniprot_id_hgnc_id_mapping(self):
        """Builds mapping from UniProt identifiers to HGNC identifier

        :rtype: dict[str,str]
        """
        return {
            uniprot_id: hgnc_id
            for hgnc_id, uniprot_id in self.session.query(HGNC.identifier, UniProt.uniprotid).all()
        }

    def build_uniprot_id_hgnc_symbol_mapping(self):
        """Builds mapping from UniProt identifier to HGNC symbol

        :rtype: dict[str,str]
        """
        return {
            uniprot_id: symbol
            for symbol, uniprot_id in self.session.query(HGNC.symbol, UniProt.uniprotid).all()
        }

    def get_all_hgnc_symbols(self):
        """Returns the set of hgnc symbols in PyHGNC

        :rtype: set
        """
        return {
            _deal_with_nonsense(symbol)
            for symbol in self.session.query(HGNC.symbol).all()
        }

    def _get_gene_encodings(self):
        """Gets the name to encoding dictionary for HGNC gene names

        :rtype: dict[str,str]
        """
        return {
            symbol: encodings.get(locus_type, 'GRP')
            for symbol, locus_type in self.session.query(HGNC.symbol, HGNC.locus_type).all()
        }

    def write_gene_bel_namespace(self, file):
        values = self._get_gene_encodings()
        _write_gene_bel_namespace_helper(values, file)

    def deploy_gene_bel_namespace(self):
        """Creates and deploys the Gene Names Namespace

        :rtype: Optional[str]
        """
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
        file_name = get_today_arty_namespace('gfam')

        with open(file_name, 'w') as file:
            self.write_gene_family_bel_namespace(file)

        return deploy_namespace(file_name, 'gfam')

    def export_gene_families_to_bel_graph(self):
        """Export gene family definitions as a BEL graph

        :rtype: pybel.BELGraph
        """
        graph = BELGraph(
            name='HGNC Gene Family Definitions',
            version='1.0.0',  # FIXME use database version
            authors='HGNC Consortium',
            description='Gene memberships to gene families',
            contact='charles.hoyt@scai.fraunhofer.de',
        )

        graph.namespace_url.update({
            'HGNC': get_latest_arty_namespace('hgnc'),
            'GFAM': get_latest_arty_namespace('gfam')
        })

        for family in self.session.query(GeneFamily).all():
            for gene in family.hgncs:
                graph.add_is_a(
                    gene_to_pybel(gene),
                    family_to_pybel(family)
                )

        return graph

    def write_gene_family_bel(self, file):
        """Writes the Gene Family memberships as BEL

        :param file:
        """
        graph = self.export_gene_families_to_bel_graph()
        to_bel(graph=graph, file=file)

    def deploy_gene_family_bel(self):
        """Writes the Gene family memberships as BEL and deploys"""
        name = 'hgnc-gene-family-membership'
        file_name = get_today_arty_knowledge(name)

        with open(file_name, 'w') as file:
            self.write_gene_family_bel(file)

        return deploy_knowledge(file_name, name)

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
