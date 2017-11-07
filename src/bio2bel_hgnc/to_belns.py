# -*- coding: utf-8 -*-

import logging

from bio2bel_hgnc.constants import encodings
from bio2bel_hgnc.utils import ensure_pyhgnc_manager
from pybel.constants import NAMESPACE_DOMAIN_GENE
from pybel.resources.arty import get_today_arty_namespace
from pybel.resources.definitions.definitions_write import write_namespace
from pybel.resources.deploy import deploy_namespace

log = logging.getLogger(__name__)

MODULE_NAME = 'hgnc'


def write_belns(file=None, manager=None):
    """Writes an HGNC namespace

    :param file file: A write-enabled file or file-like. Defaults to standard out.
    :param pyhgnc.QueryManager manager: A PyHGNC database manager
    """
    manager = ensure_pyhgnc_manager(connection=manager)

    log.info('getting entries')

    all_entries = manager.hgnc()

    values = {
        entry.symbol: encodings.get(entry.locus_type, 'GRP')
        for entry in all_entries
    }

    log.info('writing namespace')

    write_namespace(
        namespace_name='HGNC Human Genes',
        namespace_keyword='HGNC',
        namespace_domain=NAMESPACE_DOMAIN_GENE,
        author_name='Charles Tapley Hoyt',
        citation_name='HGNC?',  # FIXME
        values=values,
        file=file
    )


def deploy_to_arty():
    """Gets the data, writes BEL namespace"""
    file_name = get_today_arty_namespace(MODULE_NAME)

    with open(file_name, 'w') as file:
        write_belns(file)

    namespace_deploy_success = deploy_namespace(file_name, MODULE_NAME)

    if not namespace_deploy_success:
        log.warning('did not redeploy')


if __name__ == '__main__':
    import os

    with open(os.path.expanduser('~/Desktop/hgnc.belns'), 'w') as f:
        write_belns(file=f)
