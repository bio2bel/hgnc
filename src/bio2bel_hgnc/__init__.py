# -*- coding: utf-8 -*-

"""A package for converting HGNC to BEL."""

from .gfam_manager import Manager as FamilyManager  # noqa: F401
from .manager import Manager  # noqa: F401
from .utils import get_version  # noqa: F401

__version__ = '0.2.1-dev'

__title__ = 'bio2bel_hgnc'
__description__ = "A package for converting HGNC to BEL"
__url__ = 'https://github.com/bio2bel/hgnc'

__author__ = 'Charles Tapley Hoyt'
__email__ = 'charles.hoyt@scai.fraunhofer.de'

__license__ = 'MIT License'
__copyright__ = 'Copyright (c) 2017-2018 Charles Tapley Hoyt'
