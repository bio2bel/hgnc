# -*- coding: utf-8 -*-

import logging
import os
import tempfile
import unittest

import pyhgnc
from bio2bel_hgnc.enrich import add_metadata, add_node_equivalencies, add_node_orthologies, get_node
from pybel import BELGraph
from pybel.constants import EQUIVALENT_TO, ORTHOLOGOUS, RELATION, unqualified_edge_code
from pybel.dsl.nodes import protein
from pybel.examples import sialic_acid_graph
from pybel.examples.sialic_acid_example import cd33
from pybel.parser.canonicalize import po_to_tuple
from pyhgnc import QueryManager
from tests.constants import hcop_test_path, hgnc_test_path

log = logging.getLogger(__name__)

equivalence_code = unqualified_edge_code[EQUIVALENT_TO]


class TestEnrich(unittest.TestCase):
    """
    :type manager: QueryManager
    """

    @classmethod
    def setUpClass(cls):
        """Create and populate temporary PyHGNC cache"""
        cls.fd, cls.path = tempfile.mkstemp()
        cls.connection = 'sqlite:///' + cls.path
        log.info('Test generated connection string %s', cls.connection)

        pyhgnc.update(
            connection=cls.connection,
            hgnc_file_path=hgnc_test_path,
            hcop_file_path=hcop_test_path,
        )

        cls.manager = QueryManager(connection=cls.connection)

    @classmethod
    def tearDownClass(cls):
        """Deletes the temporary PyHGNC cache"""
        cls.manager.session.close()
        os.close(cls.fd)
        os.remove(cls.path)

    def help_check_cd33_model(self, model):
        """Checks if the given model is CD33

        :param pyhgnc.manager.models.HGNC model: The result from a search of the PyHGNC database
        """
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

        cd33_model = get_node(sialic_acid_graph, cd33_tuple, manager=self.manager)
        self.help_check_cd33_model(cd33_model)

    def test_get_hgnc_id_node(self):
        graph = BELGraph()

        cd33_tuple = graph.add_node_from_data(protein(identifier='1659', namespace='HGNC'))

        cd33_model = get_node(sialic_acid_graph, cd33_tuple, manager=self.manager)
        self.help_check_cd33_model(cd33_model)

    def test_get_rgd_node(self):
        graph = BELGraph()

        # CD33's MGI counterpart's identifier
        cd33_tuple = graph.add_node_from_data(protein(name='Cd33', namespace='RGD'))

        cd33_model = get_node(sialic_acid_graph, cd33_tuple, manager=self.manager)
        self.help_check_cd33_model(cd33_model)

    def test_get_rgd_id_node(self):
        graph = BELGraph()

        # CD33's RGD counterpart's identifier
        cd33_tuple = graph.add_node_from_data(protein(identifier='1596020', namespace='RGD'))

        cd33_model = get_node(sialic_acid_graph, cd33_tuple, manager=self.manager)
        self.help_check_cd33_model(cd33_model)

    def test_get_mgi_node(self):
        graph = BELGraph()

        # CD33's MGI counterpart's identifier
        cd33_tuple = graph.add_node_from_data(protein(name='Cd33', namespace='MGI'))

        cd33_model = get_node(sialic_acid_graph, cd33_tuple, manager=self.manager)
        self.help_check_cd33_model(cd33_model)

    def test_get_mgi_id_node(self):
        graph = BELGraph()

        # CD33's MGI counterpart's identifier
        cd33_tuple = graph.add_node_from_data(protein(identifier='99440', namespace='MGI'))

        cd33_model = get_node(sialic_acid_graph, cd33_tuple, manager=self.manager)
        self.help_check_cd33_model(cd33_model)

    def test_get_entrez_node(self):
        graph = BELGraph()

        # CD33's MGI counterpart's identifier
        cd33_tuple = graph.add_node_from_data(protein(identifier='945', namespace='EG'))

        cd33_model = get_node(sialic_acid_graph, cd33_tuple, manager=self.manager)
        self.help_check_cd33_model(cd33_model)

    def test_add_metadata(self):
        cd33_tuple = sialic_acid_graph.add_node_from_data(cd33)

        self.assertIn(cd33_tuple, sialic_acid_graph)
        self.assertIsNone(sialic_acid_graph.get_node_label(cd33_tuple))
        self.assertIsNone(sialic_acid_graph.get_node_identifier(cd33_tuple))

        add_metadata(sialic_acid_graph, cd33_tuple, manager=self.manager)

        self.assertIn(cd33_tuple, sialic_acid_graph)
        self.assertIsNotNone(sialic_acid_graph.get_node_label(cd33_tuple))

        self.assertEqual('CD33 molecule', sialic_acid_graph.get_node_label(cd33_tuple))
        self.assertEqual('1659', sialic_acid_graph.get_node_identifier(cd33_tuple))

    def test_add_equivalency(self):
        graph = BELGraph()
        cd33_hgnc_tuple = graph.add_node_from_data(protein(name='CD33', namespace='HGNC'))
        cd33_eg_tuple = graph.add_node_from_data(protein(identifier='945', namespace='EG'))

        self.assertEqual(2, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())

        add_node_equivalencies(graph, cd33_hgnc_tuple, manager=self.manager, add_leaves=False)

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

        add_node_equivalencies(graph, cd33_hgnc_tuple, manager=self.manager, add_leaves=False)

        cd33_ed = protein(identifier='945', namespace='EG')

        self.assertTrue(graph.has_node_with_data(cd33_ed))

        cd33_eg_node = po_to_tuple(cd33_ed)

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
