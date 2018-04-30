# -*- coding: utf-8 -*-

import logging
import os
import tempfile
import unittest

import pyhgnc
from bio2bel_hgnc import Manager
from bio2bel_hgnc.constants import GENE_FAMILY_KEYWORD
from bio2bel_hgnc.enrich import (
    add_metadata, add_node_central_dogma, add_node_equivalencies, add_node_orthologies, get_node,
)
from pybel import BELGraph
from pybel.constants import EQUIVALENT_TO, ORTHOLOGOUS, RELATION, TRANSCRIBED_TO, TRANSLATED_TO, unqualified_edge_code
from pybel.dsl import gene, mirna, protein, rna
from pybel.tokens import node_to_tuple
from tests.constants import hcop_test_path, hgnc_test_path

log = logging.getLogger(__name__)

translate_code = unqualified_edge_code[TRANSLATED_TO]
transcribe_code = unqualified_edge_code[TRANSCRIBED_TO]
equivalence_code = unqualified_edge_code[EQUIVALENT_TO]


class TemporaryCacheMixin(unittest.TestCase):
    """
    :type file_handle: int
    :type path: str
    :type manager: pyhgnc.manager.query.QueryManager
    """
    file_handle, path, manager = None, None, None

    @classmethod
    def setUpClass(cls):
        """Create and populate temporary PyHGNC cache"""
        cls.file_handle, cls.path = tempfile.mkstemp()
        cls.connection = 'sqlite:///' + cls.path
        log.info('Test generated connection string %s', cls.connection)

        pyhgnc.update(
            connection=cls.connection,
            hgnc_file_path=hgnc_test_path,
            hcop_file_path=hcop_test_path,
        )

        cls.manager = Manager(connection=cls.connection)

    @classmethod
    def tearDownClass(cls):
        """Deletes the temporary PyHGNC cache"""
        cls.manager.session.close()
        os.close(cls.file_handle)
        os.remove(cls.path)


class TestEnrich(TemporaryCacheMixin):
    def help_check_cd33_model(self, model):
        """Checks if the given model is CD33

        :param pyhgnc.manager.models.HGNC model: The result from a search of the PyHGNC database
        """
        self.assertIsNotNone(model)
        self.assertEqual('1659', str(model.identifier))
        self.assertEqual('CD33', model.symbol)
        self.assertEqual('CD33 molecule', model.name)

    def test_pyhgnc_loaded(self):
        cd33_results = self.manager.hgnc(symbol='CD33')
        self.assertIsNotNone(cd33_results)

        self.assertEqual(1, len(cd33_results))

        cd33_model = cd33_results[0]
        self.help_check_cd33_model(cd33_model)

    def test_get_hgnc_node(self):
        graph = BELGraph()

        cd33_tuple = graph.add_node_from_data(protein(name='CD33', namespace='HGNC'))

        cd33_model = get_node(graph, cd33_tuple, connection=self.manager)
        self.help_check_cd33_model(cd33_model)

    def test_get_hgnc_id_node(self):
        graph = BELGraph()

        cd33_tuple = graph.add_node_from_data(protein(identifier='1659', namespace='HGNC'))

        cd33_model = get_node(graph, cd33_tuple, connection=self.manager)
        self.help_check_cd33_model(cd33_model)

    @unittest.skip('HGNC does not have RGD symbols')
    def test_get_rgd_node(self):
        graph = BELGraph()

        # CD33's MGI counterpart's identifier
        cd33_tuple = graph.add_node_from_data(protein(name='Cd33', namespace='RGD'))

        cd33_model = get_node(graph, cd33_tuple, connection=self.manager)
        self.help_check_cd33_model(cd33_model)

    def test_get_rgd_id_node(self):
        graph = BELGraph()

        # CD33's RGD counterpart's identifier
        cd33_tuple = graph.add_node_from_data(protein(identifier='1596020', namespace='RGD'))

        cd33_model = get_node(graph, cd33_tuple, connection=self.manager)
        self.help_check_cd33_model(cd33_model)

    @unittest.skip('HGNC does not have MGI symbol information')
    def test_get_mgi_node(self):
        graph = BELGraph()

        # CD33's MGI counterpart's identifier
        cd33_tuple = graph.add_node_from_data(protein(name='Cd33', namespace='MGI'))

        cd33_model = get_node(graph, cd33_tuple, connection=self.manager)
        self.help_check_cd33_model(cd33_model)

    def test_get_mgi_id_node(self):
        graph = BELGraph()

        # CD33's MGI counterpart's identifier
        cd33_tuple = graph.add_node_from_data(protein(identifier='99440', namespace='MGI'))

        cd33_model = get_node(graph, cd33_tuple, connection=self.manager)
        self.help_check_cd33_model(cd33_model)

    def test_get_entrez_node(self):
        graph = BELGraph()

        # CD33's MGI counterpart's identifier
        cd33_tuple = graph.add_node_from_data(protein(identifier='945', namespace='ENTREZ'))

        cd33_model = get_node(graph, cd33_tuple, connection=self.manager)
        self.help_check_cd33_model(cd33_model)

    def test_add_metadata(self):
        graph = BELGraph()

        cd33_test = protein(namespace='HGNC', name='CD33')
        graph.add_node_from_data(cd33_test)

        self.assertIn(cd33_test.as_tuple(), graph, msg='Graph is missing CD33 protein node')
        self.assertIsNone(graph.get_node_label(cd33_test.as_tuple()), msg='CD33 should not have label information')
        self.assertIsNone(graph.get_node_identifier(cd33_test.as_tuple()), msg='CD33 should not have identifier information')

        add_metadata(graph, cd33_test.as_tuple(), manager=self.manager)

        self.assertIn(cd33_test.as_tuple(), graph, msg='Graph somehow lost CD33 protein node')

        self.assertEqual('1659', graph.get_node_identifier(cd33_test.as_tuple()), msg='Graph should be enriched with identifier')

    def test_add_equivalency(self):
        graph = BELGraph()
        cd33_hgnc_tuple = graph.add_node_from_data(protein(name='CD33', namespace='HGNC'))
        cd33_eg_tuple = graph.add_node_from_data(protein(identifier='945', namespace='ENTREZ'))

        self.assertEqual(2, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())

        add_node_equivalencies(graph, cd33_hgnc_tuple, connection=self.manager, add_leaves=False)

        self.assertEqual(2, graph.number_of_nodes())
        self.assertEqual(2, graph.number_of_edges())

        self.assertIn(cd33_eg_tuple, graph.edge[cd33_hgnc_tuple])
        v = list(graph.edge[cd33_hgnc_tuple][cd33_eg_tuple].values())[0]
        self.assertIn(RELATION, v)
        self.assertEqual(EQUIVALENT_TO, v[RELATION])

        self.assertIn(cd33_hgnc_tuple, graph.edge[cd33_eg_tuple])
        v = list(graph.edge[cd33_eg_tuple][cd33_hgnc_tuple].values())[0]
        self.assertIn(RELATION, v)
        self.assertEqual(EQUIVALENT_TO, v[RELATION])

    def test_add_equivalency_with_leaves(self):
        graph = BELGraph()
        cd33_hgnc_tuple = graph.add_node_from_data(protein(name='CD33', namespace='HGNC'))

        self.assertEqual(1, graph.number_of_nodes())

        add_node_equivalencies(graph, cd33_hgnc_tuple, connection=self.manager, add_leaves=True)

        cd33_entrez = protein(identifier='945', namespace='ENTREZ')

        self.assertTrue(graph.has_node_with_data(cd33_entrez))

        cd33_eg_node = node_to_tuple(cd33_entrez)

        self.assertIn(cd33_eg_node, graph.edge[cd33_hgnc_tuple])
        self.assertIn(equivalence_code, graph.edge[cd33_hgnc_tuple][cd33_eg_node])

    def test_add_orthology(self):
        """Tests adding orthologies when both nodes are identified by name"""
        graph = BELGraph()
        cd33_hgnc_tuple = graph.add_node_from_data(protein(name='CD33', namespace='HGNC'))
        cd33_mgi_tuple = graph.add_node_from_data(protein(name='Cd33', namespace='MGI'))

        self.assertEqual(2, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())

        add_node_orthologies(graph, cd33_hgnc_tuple, manager=self.manager, add_leaves=False)

        self.assertEqual(2, graph.number_of_nodes())
        self.assertEqual(2, graph.number_of_edges())

        self.assertIn(cd33_mgi_tuple, graph.edge[cd33_hgnc_tuple])
        v = list(graph.edge[cd33_hgnc_tuple][cd33_mgi_tuple].values())[0]
        self.assertIn(RELATION, v)
        self.assertEqual(ORTHOLOGOUS, v[RELATION])

        self.assertIn(cd33_hgnc_tuple, graph.edge[cd33_mgi_tuple])
        v = list(graph.edge[cd33_mgi_tuple][cd33_hgnc_tuple].values())[0]
        self.assertIn(RELATION, v)
        self.assertEqual(ORTHOLOGOUS, v[RELATION])

    def test_add_orthology_target_has_id(self):
        """Tests adding orthologies when the target is identified by its database identiier"""
        graph = BELGraph()
        cd33_hgnc_tuple = graph.add_node_from_data(protein(name='CD33', namespace='HGNC'))
        cd33_mgi_id_tuple = graph.add_node_from_data(protein(identifier='99440', namespace='MGI'))

        self.assertEqual(2, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())

        add_node_orthologies(graph, cd33_hgnc_tuple, manager=self.manager, add_leaves=False)

        self.assertEqual(2, graph.number_of_nodes())
        self.assertEqual(2, graph.number_of_edges())

        self.assertIn(cd33_hgnc_tuple, graph.edge[cd33_mgi_id_tuple])
        v = list(graph.edge[cd33_mgi_id_tuple][cd33_hgnc_tuple].values())[0]
        self.assertIn(RELATION, v)
        self.assertEqual(ORTHOLOGOUS, v[RELATION])

        self.assertIn(cd33_mgi_id_tuple, graph.edge[cd33_hgnc_tuple])
        v = list(graph.edge[cd33_hgnc_tuple][cd33_mgi_id_tuple].values())[0]
        self.assertIn(RELATION, v)
        self.assertEqual(ORTHOLOGOUS, v[RELATION])

    def test_add_mirna(self):
        graph = BELGraph()
        mir489_gene = gene(namespace='HGNC', name='MIR489', identifier='32074')
        mir489_gene_tuple = graph.add_node_from_data(mir489_gene)

        self.assertEqual(1, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())

        add_node_central_dogma(graph, mir489_gene_tuple, connection=self.manager)

        self.assertEqual(2, graph.number_of_nodes())
        self.assertEqual(1, graph.number_of_edges())

        mir489_mirna = mirna(namespace='HGNC', name='MIR489', identifier='32074')
        self.assertTrue(graph.has_node_with_data(mir489_mirna))

        # Check doesn't add protein
        mir489_protein = mirna(namespace='HGNC', name='MIR489', identifier='32074')
        self.assertFalse(graph.has_node_with_data(mir489_protein))

    def test_add_rna(self):
        graph = BELGraph()
        mir503hg_gene = gene(namespace='HGNC', name='MIR503HG', identifier='28258')
        graph.add_node_from_data(mir503hg_gene)

        self.assertEqual(1, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())

        add_node_central_dogma(graph, mir503hg_gene.as_tuple(), connection=self.manager)

        self.assertEqual(2, graph.number_of_nodes())
        self.assertEqual(1, graph.number_of_edges())

        mir503hg_rna = rna(namespace='HGNC', name='MIR503HG', identifier='28258')
        self.assertTrue(graph.has_node_with_data(mir503hg_rna))

        # Check doesn't add protein
        mir503hg_protein = protein(namespace='HGNC', name='MIR503HG', identifier='28258')
        self.assertFalse(graph.has_node_with_data(mir503hg_protein))

    def test_add_protein(self):
        graph = BELGraph()
        cd33_gene = gene(name='CD33', namespace='HGNC')
        cd33_gene_tuple = graph.add_node_from_data(cd33_gene)

        self.assertEqual(1, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())

        add_node_central_dogma(graph, cd33_gene_tuple, connection=self.manager)

        self.assertEqual(3, graph.number_of_nodes())
        self.assertEqual(2, graph.number_of_edges())

        cd33_rna = rna(name='CD33', namespace='HGNC')
        self.assertTrue(graph.has_node_with_data(cd33_rna))

        cd33_rna_tuple = node_to_tuple(cd33_rna)
        self.assertIn(cd33_rna_tuple, graph.edge[cd33_gene_tuple])
        self.assertIn(transcribe_code, graph.edge[cd33_gene_tuple][cd33_rna_tuple])

        cd33_protein = protein(name='CD33', namespace='HGNC')
        self.assertTrue(graph.has_node_with_data(cd33_protein))

        cd33_protein_tuple = node_to_tuple(cd33_protein)
        self.assertIn(cd33_protein_tuple, graph.edge[cd33_rna_tuple])
        self.assertIn(translate_code, graph.edge[cd33_rna_tuple][cd33_protein_tuple])

    def test_enrich_genes_with_families(self):
        """For genes, adds their corresponding families"""
        graph = BELGraph()
        cd33_gene = gene(name='CD33', namespace='HGNC')
        cd33_gene_tuple = graph.add_node_from_data(cd33_gene)

        self.assertEqual(1, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())

        self.manager.enrich_genes_with_families(graph)

        self.assertEqual(4, graph.number_of_nodes())
        self.assertEqual(3, graph.number_of_edges())

        for x in ["CD molecules", "V-set domain containing", "Sialic acid binding Ig like lectins"]:
            g = gene(namespace=GENE_FAMILY_KEYWORD, name=x)
            self.assertTrue(graph.has_node_with_data(g))
            self.assertIn(g.as_tuple(), graph.edge[cd33_gene_tuple])

    def test_enrich_family_with_genes(self):
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

    def test_enrich_families_with_genes(self):
        """For gene families, adds their member genes"""
        raise NotImplementedError
