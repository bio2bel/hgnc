# -*- coding: utf-8 -*-

"""Test cases for Bio2BEL HGNC."""

from bio2bel.testing import AbstractTemporaryCacheClassMixin

from bio2bel_hgnc import Manager
from tests.constants import hcop_test_path, hgnc_test_path


class TemporaryCacheMixin(AbstractTemporaryCacheClassMixin):
    """A test case that has a pre-populated HGNC database."""

    Manager = Manager

    @classmethod
    def populate(cls):
        """Populate the database."""
        cls.manager.populate(
            hgnc_file_path=hgnc_test_path,
            hcop_file_path=hcop_test_path,
        )
