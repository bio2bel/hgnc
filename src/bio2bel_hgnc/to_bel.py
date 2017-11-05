# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import pandas as pd
from pybel.utils import ensure_quotes
from pybel_tools.document_utils import write_boilerplate
from pybel.resources.defaults import CONFIDENCE
from pybel.resources.arty import get_latest_arty_namespace

url = 'http://www.genenames.org/cgi-bin/download?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_prev_sym&col=gd_aliases&col=gd_pub_chrom_map&col=gd_date_name_change&col=gd_pub_acc_ids&col=gd_enz_ids&col=gd_pub_eg_id&col=gd_other_ids&col=gd_other_ids_list&col=md_prot_id&col=md_mgd_id&col=md_rgd_id&status=Approved&status=Entry+Withdrawn&status_opt=2&where=&order_by=gd_app_sym_sort&format=text&limit=&hgnc_dbtag=on&submit=submit'


def get_data():
    """Gets the gene information from HGNC

    :rtype: pandas.DataFrame
    """
    return pd.read_csv(url, sep='\t')


def write_hgnc_equivalences_boilerplate(file):
    """Writes the header of the HGNC equivalences file

    :param file file: A write-enabled file or file-like
    """
    write_boilerplate(
        name='HGNC Equivalences',
        description="""This document contains the equivalence information from HGNC to RGD, MGI, UniProt, EC, Entrez and more""",
        authors='Charles Tapley Hoyt',
        contact='charles.hoyt@scai.fraunhofer.de',
        licenses='Creative Commons by 4.0',
        copyright='Copyright (c) 2017 Charles Tapley Hoyt. All rights reserved',
        namespace_url={
            'HGNC': get_latest_arty_namespace('hgnc-human-genes'),
            'UniProt': get_latest_arty_namespace('uniprot'),
            'MGI': get_latest_arty_namespace('mgi-mouse-genes'),
            'RGD': get_latest_arty_namespace('rgd-rat-genes'),
            'EGID': get_latest_arty_namespace('entrez-gene-ids'),
        },
        namespace_patterns={
            'EC': ".+"
        },
        annotations_dict={'Confidence': CONFIDENCE},
        file=file
    )


columns = [
    'Approved Symbol',
    'Entrez Gene ID',
    'UniProt ID(supplied by UniProt)',
    'Mouse Genome Database ID(supplied by MGI)',
    'Rat Genome Database ID(supplied by RGD)',
    'Enzyme IDs',
]


def write_hgnc_equivalences_body(file=None):
    """Writes the body of the HGNC equivalences file

    :param file file: A write-enabled file or file-like
    """
    df = get_data()

    print('SET Citation = {"PubMed","HGNC","25361968"}', file=file)
    print('SET Evidence = "HGNC Definitions"', file=file)
    print('SET Confidence = "Axiomatic"', file=file)

    for _, hgnc, egid, upid, mgi, rgd, ecs in df[columns].itertuples():

        hgnc = ensure_quotes(hgnc)

        if pd.notnull(egid):
            print('g(HGNC:{}) eq g(EGID:{})'.format(hgnc, str(int(egid))), file=file)

        if pd.notnull(upid):
            print('p(HGNC:{}) eq p(UNIPROT:{})'.format(hgnc, upid), file=file)

        # print('g(HGNC:{}) orthologousTo g(MGI:{})'.format(hgnc, mgi), file=file)
        # print('g(HGNC:{}) orthologousTo g(RGD:{})'.format(hgnc, rgd), file=file)

        if pd.notnull(ecs):
            for ec in str(ecs).strip().split(','):
                print('p(EC:{}) hasMember p(HGNC:{})'.format(ensure_quotes(ec), hgnc), file=file)


def write_hgnc_equivalences(file=None):
    """Writes the HGNC equivalences to Entrez, UniProt, and Enzyme Commission.

    :param file file: A write-enabled file or file-like
    """
    write_hgnc_equivalences_boilerplate(file)
    write_hgnc_equivalences_body(file)


if __name__ == '__main__':
    desktop = os.path.join(os.path.expanduser('~'), 'Desktop', 'hgnc-equivalences.bel')

    with open(desktop, 'w') as f:
        write_hgnc_equivalences(f)
