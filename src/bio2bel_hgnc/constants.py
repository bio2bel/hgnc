# -*- coding: utf-8 -*-

"""Enrichment functions for BEL graphs"""

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
