# -*- coding: utf-8 -*-

"""Test cases for Bio2BEL HGNC."""

from bio2bel.testing import AbstractTemporaryCacheClassMixin
from bio2bel_hgnc import Manager
from tests.constants import hcop_test_path, hgnc_test_path


class TemporaryCacheMixin(AbstractTemporaryCacheClassMixin):
    Manager = Manager

    @classmethod
    def populate(cls):
        cls.manager.populate(
            hgnc_file_path=hgnc_test_path,
            hcop_file_path=hcop_test_path,
        )
