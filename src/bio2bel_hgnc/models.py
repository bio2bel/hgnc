# -*- coding: utf-8 -*-

"""SQLAlchemy models for Bio2BEL HGNC.

Note: currently wraps models from :mod:`PyHGNC`.
"""

from pyhgnc.manager.models import AliasName, AliasSymbol, Base, Enzyme, GeneFamily, UniProt
from pyhgnc.manager.models import HGNC as HumanGene  # noqa: N811
from pyhgnc.manager.models import MGD as MouseGene  # noqa: N811
from pyhgnc.manager.models import RGD as RatGene  # noqa: N811

__all__ = [
    'AliasSymbol',
    'AliasName',
    'Base',
    'Enzyme',
    'GeneFamily',
    'HumanGene',
    'MouseGene',
    'RatGene',
    'UniProt',
]
