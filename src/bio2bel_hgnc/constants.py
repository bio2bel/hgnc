# -*- coding: utf-8 -*-

"""Enrichment functions for BEL graphs"""

from bio2bel.utils import get_data_dir

VERSION = '0.1.0-dev'

MODULE_NAME = 'hgnc'
#: The default directory where PyBEL files, including logs and the  default cache, are stored. Created if not exists.
DATA_DIR = get_data_dir(MODULE_NAME)

GENE_FAMILY_KEYWORD = 'GFAM'

#: Encodings from https://www.genenames.org/cgi-bin/statistics
encodings = {
    # protein-coding gene
    'gene with protein product': 'GRP',
    # non-coding RNA
    'RNA, Y': 'GR',
    'RNA, cluster': 'GR',
    'RNA, long non-coding': 'GR',
    'RNA, micro': 'GM',
    'RNA, misc': 'GR',
    'RNA, ribosomal': 'GR',
    'RNA, small cytoplasmic': 'GR',
    'RNA, small nuclear': 'GR',
    'RNA, small nucleolar': 'GR',
    'RNA, transfer': 'GR',
    'RNA, vault': 'GR',
    # phenotype
    'phenotype only': 'G',
    # pseudogene
    'T cell receptor pseudogene': 'G',
    'immunoglobulin pseudogene': 'G',
    'pseudogene': 'G',
    # other
    'T-cell receptor pseudogene': 'G',
    'T-cell receptor gene': 'G',
    'T cell receptor gene': 'G',
    'complex locus constituent': 'G',
    'endogenous retrovirus': 'G',
    'fragile site': 'G',
    'immunoglobulin gene': 'G',
    'protocadherin': 'G',
    'readthrough': 'G',
    'region': 'G',
    'transposable element': 'G',
    'unknown': 'GRP',
    'virus integration site': 'G',
}
