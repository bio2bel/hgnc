# -*- coding: utf-8 -*-

import logging
import os
import tempfile
import unittest

import pyhgnc
from bio2bel_hgnc.enrich import add_metadata, get_node
from pybel import BELGraph
from pybel.dsl.nodes import protein
from pybel.examples import sialic_acid_graph
from pybel.examples.sialic_acid_example import cd33
from pyhgnc import QueryManager
from tests.constants import hcop_test_path, hgnc_test_path

log = logging.getLogger(__name__)


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
