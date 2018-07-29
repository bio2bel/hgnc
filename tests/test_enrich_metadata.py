# -*- coding: utf-8 -*-

"""Tests for enrichment methods."""

import logging
import unittest
from pybel import BELGraph
from pybel.constants import EQUIVALENT_TO, ORTHOLOGOUS, RELATION, TRANSCRIBED_TO, TRANSLATED_TO, unqualified_edge_code
from pybel.dsl import gene, mirna, protein, rna

from bio2bel_hgnc.constants import GENE_FAMILY_KEYWORD
from bio2bel_hgnc.enrich import (
    add_node_orthologies, get_node,
)
from tests.cases import TemporaryCacheMixin

log = logging.getLogger(__name__)

translate_code = unqualified_edge_code[TRANSLATED_TO]
transcribe_code = unqualified_edge_code[TRANSCRIBED_TO]
equivalence_code = unqualified_edge_code[EQUIVALENT_TO]


class TestEnrich(TemporaryCacheMixin):
    """Test enrichment methods."""

    def help_check_cd33_model(self, model):
        """Check if the given model is CD33.

        :param pyhgnc.manager.models.HGNC model: The result from a search of the PyHGNC database
        """
        self.assertIsNotNone(model)
        self.assertEqual('1659', str(model.identifier))
        self.assertEqual('CD33', model.symbol)
        self.assertEqual('CD33 molecule', model.name)

    def test_pyhgnc_loaded(self):
        """"Test that data has been loaded properly."""
        cd33_results = self.manager.hgnc(symbol='CD33')
        self.assertIsNotNone(cd33_results)

        self.assertEqual(1, len(cd33_results))

        cd33_model = cd33_results[0]
        self.help_check_cd33_model(cd33_model)

    def test_get_hgnc_node(self):
        """Test getting a node by name from the database."""
        graph = BELGraph()

        cd33_tuple = graph.add_node_from_data(protein(name='CD33', namespace='HGNC'))

        cd33_model = get_node(graph, cd33_tuple, manager=self.manager)
        self.help_check_cd33_model(cd33_model)

    def test_get_hgnc_id_node(self):
        """Test getting a node by identifier from the database."""
        graph = BELGraph()

        cd33_tuple = graph.add_node_from_data(protein(identifier='1659', namespace='HGNC'))

        cd33_model = get_node(graph, cd33_tuple, manager=self.manager)
        self.help_check_cd33_model(cd33_model)

    @unittest.skip('HGNC does not have RGD symbols')
    def test_get_rgd_node(self):
        """Test getting a node by RGD name."""
        graph = BELGraph()

        # CD33's MGI counterpart's identifier
        cd33_tuple = graph.add_node_from_data(protein(name='Cd33', namespace='RGD'))

        cd33_model = get_node(graph, cd33_tuple, manager=self.manager)
        self.help_check_cd33_model(cd33_model)

    def test_get_rgd_id_node(self):
        """Test getting a node by RGD identifier."""
        graph = BELGraph()

        # CD33's RGD counterpart's identifier
        cd33_tuple = graph.add_node_from_data(protein(identifier='1596020', namespace='RGD'))

        cd33_model = get_node(graph, cd33_tuple, manager=self.manager)
        self.help_check_cd33_model(cd33_model)

    @unittest.skip('HGNC does not have MGI symbol information')
    def test_get_mgi_node(self):
        """Test getting a node by MGI name."""
        graph = BELGraph()

        # CD33's MGI counterpart's identifier
        cd33_tuple = graph.add_node_from_data(protein(name='Cd33', namespace='MGI'))

        cd33_model = get_node(graph, cd33_tuple, manager=self.manager)
        self.help_check_cd33_model(cd33_model)

    def test_get_mgi_id_node(self):
        """Test getting a node by MGI identifier."""
        graph = BELGraph()

        # CD33's MGI counterpart's identifier
        cd33_tuple = graph.add_node_from_data(protein(identifier='99440', namespace='MGI'))

        cd33_model = get_node(graph, cd33_tuple, manager=self.manager)
        self.help_check_cd33_model(cd33_model)

    def test_get_entrez_node(self):
        """Test getting a node by Entrez Gene identifier."""
        graph = BELGraph()

        # CD33's MGI counterpart's identifier
        cd33_tuple = graph.add_node_from_data(protein(identifier='945', namespace='ENTREZ'))

        cd33_model = get_node(graph, cd33_tuple, manager=self.manager)
        self.help_check_cd33_model(cd33_model)

    def test_add_equivalency(self):
        """Test that CD33 identified by HGNC and entrez can be equivalenced."""
        graph = BELGraph()

        cd33_hgnc = gene(name='CD33', namespace='HGNC')
        cd33_eg = gene(identifier='945', namespace='ENTREZ')

        cd33_hgnc_tuple = graph.add_node_from_data(cd33_hgnc)
        cd33_eg_tuple = graph.add_node_from_data(cd33_eg)

        # self.assertEqual(2, graph.number_of_nodes(), msg='wrong initial number of nodes')
        # self.assertEqual(0, graph.number_of_edges(), msg='wrong initial number of edges')

        self.manager.add_node_equivalencies(graph, cd33_hgnc)

        # self.assertEqual(2, graph.number_of_nodes(), msg='nodes: {}'.format(list(graph)))
        # self.assertEqual(2, graph.number_of_edges(), msg='adding equivalence added wrong number of edges')

        self.assertIn(cd33_eg_tuple, graph[cd33_hgnc_tuple])
        v = list(graph[cd33_hgnc_tuple][cd33_eg_tuple].values())[0]
        self.assertIn(RELATION, v)
        self.assertEqual(EQUIVALENT_TO, v[RELATION])

        self.assertIn(cd33_hgnc_tuple, graph.edge[cd33_eg_tuple])
        v = list(graph[cd33_eg_tuple][cd33_hgnc_tuple].values())[0]
        self.assertIn(RELATION, v)
        self.assertEqual(EQUIVALENT_TO, v[RELATION])

    @unittest.skip(reason='Not implemented yet!')
    def test_add_orthology(self):
        """Test adding orthologies when both nodes are identified by name."""
        graph = BELGraph()

        cd33_hgnc = protein(name='CD33', namespace='HGNC')
        cd33_mgi = protein(name='Cd33', namespace='MGI')

        graph.add_node_from_data(cd33_hgnc)
        graph.add_node_from_data(cd33_mgi)

        self.assertEqual(2, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())

        add_node_orthologies(graph, cd33_hgnc, manager=self.manager, add_leaves=False)

        self.assertEqual(2, graph.number_of_nodes())
        self.assertEqual(2, graph.number_of_edges())

        self.assertIn(cd33_mgi.as_tuple(), graph.edge[cd33_hgnc.as_tuple()])
        v = list(graph.edge[cd33_hgnc.as_tuple()][cd33_mgi.as_tuple()].values())[0]
        self.assertIn(RELATION, v)
        self.assertEqual(ORTHOLOGOUS, v[RELATION])

        self.assertIn(cd33_hgnc.as_tuple(), graph.edge[cd33_mgi.as_tuple()])
        v = list(graph.edge[cd33_mgi.as_tuple()][cd33_hgnc.as_tuple()].values())[0]
        self.assertIn(RELATION, v)
        self.assertEqual(ORTHOLOGOUS, v[RELATION])

    @unittest.skip(reason='Not implemented yet!')
    def test_add_orthology_target_has_id(self):
        """Test adding orthologies when the target is identified by its database identifier."""
        graph = BELGraph()

        cd33_hgnc = protein(name='CD33', namespace='HGNC')
        cd33_mgi_id = protein(identifier='99440', namespace='MGI')

        graph.add_node_from_data(cd33_hgnc)
        graph.add_node_from_data(cd33_mgi_id)

        self.assertEqual(2, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())

        add_node_orthologies(graph, cd33_hgnc, manager=self.manager, add_leaves=False)

        self.assertEqual(2, graph.number_of_nodes())
        self.assertEqual(2, graph.number_of_edges())

        self.assertIn(cd33_hgnc.as_tuple(), graph.edge[cd33_mgi_id.as_tuple()])
        v = list(graph.edge[cd33_mgi_id.as_tuple()][cd33_hgnc.as_tuple()].values())[0]
        self.assertIn(RELATION, v)
        self.assertEqual(ORTHOLOGOUS, v[RELATION])

        self.assertIn(cd33_mgi_id.as_tuple(), graph.edge[cd33_hgnc.as_tuple()])
        v = list(graph.edge[cd33_hgnc.as_tuple()][cd33_mgi_id.as_tuple()].values())[0]
        self.assertIn(RELATION, v)
        self.assertEqual(ORTHOLOGOUS, v[RELATION])

    @unittest.skip('needs serious reinvestigation')
    def test_add_mirna(self):
        """Test inferring the central dogma for an miRNA."""
        graph = BELGraph()
        mir489_gene = gene(namespace='HGNC', name='MIR489', identifier='32074')
        graph.add_node_from_data(mir489_gene)

        self.assertEqual(1, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())

        result = self.manager.add_central_dogma(graph, mir489_gene)
        self.assertIsNotNone(result)

        self.assertEqual(2, graph.number_of_nodes())
        self.assertEqual(1, graph.number_of_edges())

        mir489_mirna = mirna(namespace='HGNC', name='MIR489', identifier='32074')
        self.assertTrue(graph.has_node_with_data(mir489_mirna))

        # Check doesn't add protein
        mir489_protein = mirna(namespace='HGNC', name='MIR489', identifier='32074')
        self.assertFalse(graph.has_node_with_data(mir489_protein))

    @unittest.skip('need to reinvestigate assignment of RNA')
    def test_add_rna(self):
        """Test inferring the central dogma for an RNA."""
        graph = BELGraph()
        mir503hg_gene = gene(namespace='HGNC', name='MIR503HG', identifier='28258')
        graph.add_node_from_data(mir503hg_gene)

        self.assertEqual(1, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())

        result = self.manager.add_central_dogma(graph, mir503hg_gene)
        self.assertIsNotNone(result)
        self.assertIsInstance(result, mirna)

        self.assertEqual(2, graph.number_of_nodes())
        self.assertEqual(1, graph.number_of_edges())

        mir503hg_rna = rna(namespace='HGNC', name='MIR503HG', identifier='28258')
        self.assertTrue(graph.has_node_with_data(mir503hg_rna))

        # Check doesn't add protein
        mir503hg_protein = protein(namespace='HGNC', name='MIR503HG', identifier='28258')
        self.assertFalse(graph.has_node_with_data(mir503hg_protein))

    def test_add_protein(self):
        """Test inferring the central dogma for a protein."""
        graph = BELGraph()
        cd33_gene = gene(name='CD33', namespace='HGNC')
        graph.add_node_from_data(cd33_gene)

        self.assertEqual(1, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())

        result = self.manager.add_central_dogma(graph, cd33_gene)
        self.assertIsNotNone(result)
        self.assertIsInstance(result, protein)

        self.assertEqual(3, graph.number_of_nodes())
        self.assertEqual(2, graph.number_of_edges())

        cd33_rna = rna(name='CD33', namespace='HGNC')
        self.assertTrue(graph.has_node_with_data(cd33_rna), msg='Nodes: {}'.format(list(graph)))

        self.assertIn(cd33_rna.as_tuple(), graph.edge[cd33_gene.as_tuple()])
        self.assertIn(transcribe_code, graph.edge[cd33_gene.as_tuple()][cd33_rna.as_tuple()])

        cd33_protein = protein(name='CD33', namespace='HGNC')
        self.assertTrue(graph.has_node_with_data(cd33_protein))

        self.assertIn(cd33_protein.as_tuple(), graph.edge[cd33_rna.as_tuple()])
        self.assertIn(translate_code, graph.edge[cd33_rna.as_tuple()][cd33_protein.as_tuple()])

    def test_enrich_genes_with_families(self):
        """Tests enriching genes with their families."""
        graph = BELGraph()
        cd33_gene = gene(name='CD33', namespace='HGNC')
        graph.add_node_from_data(cd33_gene)

        self.assertEqual(1, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())

        self.manager.enrich_genes_with_families(graph)

        self.assertEqual(4, graph.number_of_nodes())
        self.assertEqual(3, graph.number_of_edges())

        for x in ["CD molecules", "V-set domain containing", "Sialic acid binding Ig like lectins"]:
            g = gene(namespace=GENE_FAMILY_KEYWORD, name=x)
            self.assertTrue(graph.has_node_with_data(g))
            self.assertIn(g.as_tuple(), graph.edge[cd33_gene.as_tuple()])

    def test_enrich_family_with_genes(self):
        """Test enriching families with their genes."""
        graph = BELGraph()
        f = gene(name="CD molecules", namespace=GENE_FAMILY_KEYWORD)
        graph.add_node_from_data(f)

        self.assertEqual(1, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())

        self.manager.enrich_families_with_genes(graph)

        self.assertEqual(3, graph.number_of_nodes())
        self.assertEqual(2, graph.number_of_edges())

        for x in ["CD33", 'CD34']:
            g = gene(namespace='HGNC', name=x)
            self.assertTrue(graph.has_node_with_data(g))
            self.assertIn(f.as_tuple(), graph.edge[g.as_tuple()])


if __name__ == '__main__':
    unittest.main()
