# -*- coding: utf-8 -*-

"""SQLAlchemy models for Bio2BEL HGNC.

Note: currently wraps models from :mod:`PyHGNC`.
"""

from pyhgnc.manager.models import AliasName, AliasSymbol, Base, Enzyme, GeneFamily, UniProt
from pyhgnc.manager.models import HGNC as HumanGene  # noqa: N811
from pyhgnc.manager.models import MGD as MouseGene  # noqa: N811
from pyhgnc.manager.models import RGD as RatGene  # noqa: N811
from pyhgnc.manager.models import (
    hgnc_enzyme as gene_enzyme,
    hgnc_gene_family as gene_gene_family,
    hgnc_mgd as gene_mouse_gene,
    hgnc_rgd as gene_rat_gene,
)


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
    'gene_gene_family',
    'gene_enzyme',
    'gene_mouse_gene',
    'gene_rat_gene',
]
