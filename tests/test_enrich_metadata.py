# -*- coding: utf-8 -*-

"""Tests for enrichment methods."""

import logging
import unittest

from bio2bel_hgnc import Manager
from bio2bel_hgnc.constants import ENTREZ, HGNC, HGNC_GENE_FAMILY
from bio2bel_hgnc.models import HumanGene
from pybel import BELGraph
from pybel.constants import EQUIVALENT_TO, ORTHOLOGOUS, RELATION
from pybel.dsl import gene, mirna, protein, rna
from tests.cases import TemporaryCacheMixin

log = logging.getLogger(__name__)

protein_hgnc_cd33 = protein(name='CD33', namespace=HGNC)
gene_hgnc_cd33 = protein_hgnc_cd33.get_rna().get_gene()
protein_hgnc_1659 = protein(identifier='1659', namespace=HGNC)
cd33_hgnc_name_as_identifier = protein(name='1659', namespace=HGNC)
cd33_hgnc_name_and_identifier = protein(name='CD33', identifier='1659', namespace=HGNC)

cd34_hgnc_name_and_identifier = protein(name='CD34', identifier='1662', namespace=HGNC)

cd33_rgd_name = protein(name='Cd33', namespace='RGD')
cd33_rgd_id = protein(identifier='1596020', namespace='RGD')

cd33_mgi_name = protein(name='Cd33', namespace='MGI')
cd33_mgi_identifier = protein(identifier='99440', namespace='MGI')

cd33_entrez = protein(identifier='945', name='CD33', namespace=ENTREZ)

cd_family = gene(name="CD molecules", namespace=HGNC_GENE_FAMILY)

mir489_gene = gene(namespace=HGNC, name='MIR489', identifier='32074')


class TestEnrich(TemporaryCacheMixin):
    """Test enrichment methods."""

    manager: Manager

    def help_check_cd33_model(self, model: HumanGene):
        """Check if the given model is CD33."""
        self.assertIsNotNone(model)
        self.assertEqual('1659', str(model.identifier))
        self.assertEqual('CD33', model.symbol)
        self.assertEqual('CD33 molecule', model.name)

    def test_pyhgnc_loaded(self):
        """Test that data has been loaded properly."""
        cd33_results = self.manager.hgnc(symbol='CD33')
        self.assertIsNotNone(cd33_results)

        self.assertEqual(1, len(cd33_results))

        cd33_model = cd33_results[0]
        self.help_check_cd33_model(cd33_model)

    def test_get_hgnc_node(self):
        """Test getting a node by name from the database."""
        cd33_model = self.manager.get_node(protein_hgnc_cd33)
        self.help_check_cd33_model(cd33_model)

    def test_get_hgnc_id_node(self):
        """Test getting a node by identifier from the database."""
        cd33_model = self.manager.get_node(protein_hgnc_1659)
        self.help_check_cd33_model(cd33_model)

    @unittest.skip('HGNC does not have RGD symbols')
    def test_get_rgd_node(self):
        """Test getting a node by RGD name."""
        cd33_model = self.manager.get_node(cd33_rgd_name)
        self.help_check_cd33_model(cd33_model)

    def test_get_rgd_id_node(self):
        """Test getting a node by RGD identifier."""
        cd33_model = self.manager.get_node(cd33_rgd_id)
        self.help_check_cd33_model(cd33_model)

    @unittest.skip('HGNC does not have MGI symbol information')
    def test_get_mgi_node(self):
        """Test getting a node by MGI name."""
        cd33_model = self.manager.get_node(cd33_mgi_name)
        self.help_check_cd33_model(cd33_model)

    def test_get_mgi_id_node(self):
        """Test getting a node by MGI identifier."""
        cd33_model = self.manager.get_node(cd33_mgi_identifier)
        self.help_check_cd33_model(cd33_model)

    def test_get_entrez_node(self):
        """Test getting a node by Entrez Gene identifier."""
        cd33_model = self.manager.get_node(cd33_entrez)
        self.help_check_cd33_model(cd33_model)

    def test_add_equivalency(self):
        """Test that CD33 identified by HGNC and entrez can be equivalenced."""
        graph = BELGraph()
        graph.add_node_from_data(protein_hgnc_cd33)

        self.assertEqual(1, graph.number_of_nodes(), msg='wrong initial number of nodes')
        self.assertEqual(0, graph.number_of_edges(), msg='wrong initial number of edges')

        self.manager.enrich_genes_with_equivalences(graph)

        # self.assertEqual(2, graph.number_of_nodes(), msg='nodes: {}'.format(list(graph)))
        # self.assertEqual(2, graph.number_of_edges(), msg='adding equivalence added wrong number of edges')

        self.assertIn(cd33_entrez, set(graph[protein_hgnc_cd33]))
        v = list(graph[protein_hgnc_cd33][cd33_entrez].values())[0]
        self.assertIn(RELATION, v)
        self.assertEqual(EQUIVALENT_TO, v[RELATION])

        self.assertIn(protein_hgnc_cd33, set(graph[cd33_entrez]))
        v = list(graph[cd33_entrez][protein_hgnc_cd33].values())[0]
        self.assertIn(RELATION, v)
        self.assertEqual(EQUIVALENT_TO, v[RELATION])

    @unittest.skip(reason='Not implemented yet!')
    def test_add_orthology(self):
        """Test adding orthologies when both nodes are identified by name."""
        graph = BELGraph()
        graph.add_node_from_data(protein_hgnc_cd33)
        graph.add_node_from_data(cd33_mgi_name)

        self.assertEqual(2, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())

        self.manager.add_node_orthologies(graph, protein_hgnc_cd33, add_leaves=False)

        self.assertEqual(2, graph.number_of_nodes())
        self.assertEqual(2, graph.number_of_edges())

        self.assertIn(cd33_mgi_name, graph[protein_hgnc_cd33])
        v = list(graph[protein_hgnc_cd33][cd33_mgi_name].values())[0]
        self.assertIn(RELATION, v)
        self.assertEqual(ORTHOLOGOUS, v[RELATION])

        self.assertIn(protein_hgnc_cd33, graph[cd33_mgi_name])
        v = list(graph[cd33_mgi_name][protein_hgnc_cd33].values())[0]
        self.assertIn(RELATION, v)
        self.assertEqual(ORTHOLOGOUS, v[RELATION])

    @unittest.skip(reason='Not implemented yet!')
    def test_add_orthology_target_has_id(self):
        """Test adding orthologies when the target is identified by its database identifier."""
        graph = BELGraph()
        graph.add_node_from_data(protein_hgnc_cd33)
        graph.add_node_from_data(cd33_mgi_identifier)

        self.assertEqual(2, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())

        self.manager.add_node_orthologies(graph, protein_hgnc_cd33, add_leaves=False)

        self.assertEqual(2, graph.number_of_nodes())
        self.assertEqual(2, graph.number_of_edges())

        self.assertIn(protein_hgnc_cd33, graph[cd33_mgi_identifier])
        v = list(graph[cd33_mgi_identifier][protein_hgnc_cd33].values())[0]
        self.assertIn(RELATION, v)
        self.assertEqual(ORTHOLOGOUS, v[RELATION])

        self.assertIn(cd33_mgi_identifier, graph[protein_hgnc_cd33])
        v = list(graph[protein_hgnc_cd33][cd33_mgi_identifier].values())[0]
        self.assertIn(RELATION, v)
        self.assertEqual(ORTHOLOGOUS, v[RELATION])

    @unittest.skip('needs serious reinvestigation')
    def test_add_mirna(self):
        """Test inferring the central dogma for an miRNA."""
        graph = BELGraph()
        graph.add_node_from_data(mir489_gene)

        self.assertEqual(1, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())

        result = self.manager.add_central_dogma(graph, mir489_gene)
        self.assertIsNotNone(result)

        self.assertEqual(2, graph.number_of_nodes())
        self.assertEqual(1, graph.number_of_edges())

        mir489_mirna = mirna(namespace=HGNC, name='MIR489', identifier='32074')
        self.assertIn(mir489_mirna, graph)

        # Check doesn't add protein
        mir489_protein = mirna(namespace=HGNC, name='MIR489', identifier='32074')
        self.assertNotIn(mir489_protein, graph)

    @unittest.skip('need to reinvestigate assignment of RNA')
    def test_add_rna(self):
        """Test inferring the central dogma for an RNA."""
        graph = BELGraph()
        mir503hg_gene = gene(namespace=HGNC, name='MIR503HG', identifier='28258')
        graph.add_node_from_data(mir503hg_gene)

        self.assertEqual(1, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())

        result = self.manager.add_central_dogma(graph, mir503hg_gene)
        self.assertIsNotNone(result)
        self.assertIsInstance(result, mirna)

        self.assertEqual(2, graph.number_of_nodes())
        self.assertEqual(1, graph.number_of_edges())

        mir503hg_rna = rna(namespace=HGNC, name='MIR503HG', identifier='28258')
        self.assertIn(mir503hg_rna, graph)

        # Check doesn't add protein
        mir503hg_protein = protein(namespace=HGNC, name='MIR503HG', identifier='28258')
        self.assertNotIn(mir503hg_protein, graph)

    def test_add_protein(self):
        """Test inferring the central dogma for a protein."""
        graph = BELGraph()
        cd33_gene = gene(name='CD33', namespace=HGNC)
        graph.add_node_from_data(cd33_gene)

        self.assertEqual(1, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())

        result = self.manager.add_central_dogma(graph, cd33_gene)
        self.assertIsNotNone(result)
        self.assertIsInstance(result, protein)

        self.assertEqual(3, graph.number_of_nodes(), msg='Nodes:\n{}'.format("\n".join(map(str, graph))))
        self.assertEqual(2, graph.number_of_edges())

        cd33_rna = rna(name='CD33', namespace=HGNC)
        self.assertIn(cd33_rna, graph, msg='Nodes: {}'.format(list(graph)))

        self.assertIn(cd33_rna, graph[cd33_gene])
        # self.assertIn(transcribe_code, graph[cd33_gene][cd33_rna])

        cd33_protein = protein(name='CD33', namespace=HGNC)
        self.assertIn(cd33_protein, graph)

        self.assertIn(cd33_protein, graph[cd33_rna])
        # self.assertIn(translate_code, graph[cd33_rna][cd33_protein])

    def test_enrich_genes_with_families(self):
        """Tests enriching genes with their families."""
        graph = BELGraph()
        graph.add_node_from_data(gene_hgnc_cd33)

        self.assertEqual(1, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())

        self.manager.enrich_genes_with_families(graph)

        self.assertEqual(4, graph.number_of_nodes())
        self.assertEqual(3, graph.number_of_edges())

        families = [
            ('CD molecules', '471'),
            ('V-set domain containing', '590'),
            ('Sialic acid binding Ig like lectins', '745'),
        ]
        for name, identifier in families:
            g = gene(namespace=HGNC_GENE_FAMILY, name=name, identifier=identifier)
            self.assertIn(g, graph, msg='Nodes: {}'.format(list(graph)))
            self.assertIn(g, graph[gene_hgnc_cd33])

    def test_enrich_family_with_genes(self):
        """Test enriching families with their genes."""
        graph = BELGraph()
        graph.add_node_from_data(cd_family)

        self.assertEqual(1, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())

        self.manager.enrich_families_with_genes(graph)

        self.assertEqual(3, graph.number_of_nodes())
        self.assertEqual(2, graph.number_of_edges())

        for g in [cd33_hgnc_name_and_identifier, cd34_hgnc_name_and_identifier]:
            g = g.get_rna().get_gene()
            self.assertIn(g, graph, msg=f'Nodes: {list(graph)}')
            self.assertIn(cd_family, graph[g])


if __name__ == '__main__':
    unittest.main()
