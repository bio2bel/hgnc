# -*- coding: utf-8 -*-

"""Exports HGNC Equivalences"""

from . import integrate
from .integrate import *

__all__ = (
    integrate.__all__
)

__version__ = '0.0.1-dev'

__title__ = 'bio2bel_hgnc'
__description__ = "A package for converting HGNC to BEL"
__url__ = 'https://github.com/bio2bel/hgnc'

__author__ = 'Charles Tapley Hoyt'
__email__ = 'charles.hoyt@scai.fraunhofer.de'

__license__ = 'Apache 2.0 License'
__copyright__ = 'Copyright (c) 2017 Charles Tapley Hoyt'