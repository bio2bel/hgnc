# -*- coding: utf-8 -*-

import logging
import sys

import click

from .constants import DEFAULT_CACHE_CONNECTION
from .manager import Manager

log = logging.getLogger(__name__)


def set_debug(level):
    logging.basicConfig(level=level, format="%(asctime)s - %(levelname)s - %(message)s")


def set_debug_param(debug):
    if debug == 1:
        set_debug(20)
    elif debug == 2:
        set_debug(10)


@click.group()
@click.option('-c', '--connection', help='Defaults to {}'.format(DEFAULT_CACHE_CONNECTION))
@click.pass_context
def main(ctx, connection):
    """Convert HGNC to BEL"""
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    ctx.obj = Manager(connection=connection)


@main.command()
@click.option('-v', '--debug', count=True, help="Turn on debugging.")
@click.pass_obj
def populate(manager, debug):
    """Populate the database"""
    set_debug_param(debug)
    manager.create_all()
    manager.populate()


@main.command()
@click.option('-v', '--debug', count=True, help="Turn on debugging.")
@click.option('-y', '--yes', is_flag=True, help="Automatically answer yes")
@click.pass_obj
def drop(manager, debug, yes):
    """Drops the database"""
    set_debug_param(debug)
    if yes or click.confirm('Do you really want to delete the database?'):
        manager.drop_all()


@main.command()
@click.pass_obj
def summarize(manager):
    """Summarize the contents of the database"""
    click.echo('Genes: {}'.format(manager.count_genes()))
    click.echo('Families: {}'.format(manager.count_families()))
    click.echo('UniProt Entries: {}'.format(manager.count_uniprots()))


@main.command()
@click.option('-v', '--debug', count=True, help='Turn on debugging')
@click.option('-p', '--port')
@click.option('-h', '--host')
@click.pass_obj
def web(manager, debug, port, host):
    """Run the web app"""
    set_debug_param(debug)

    from .web import get_app
    app = get_app(connection=manager, url='/')
    app.run(host=host, port=port, debug=debug)


@main.command()
@click.option('-v', '--debug', count=True, help='Turn on debugging')
@click.option('-o', '--output', type=click.File('w'), default=sys.stdout)
@click.pass_obj
def write_gene_belns(manager, debug, output):
    """Write gene namespace"""
    set_debug_param(debug)

    manager.write_gene_bel_namespace(output)


@main.command()
@click.option('-v', '--debug', count=True, help='Turn on debugging')
@click.pass_obj
def deploy_gene_belns(manager, debug):
    """Deploy gene namespace"""
    set_debug_param(debug)
    manager.deploy_gene_bel_namespace()


@main.command()
@click.option('-v', '--debug', count=True, help='Turn on debugging')
@click.option('-o', '--output', type=click.File('w'), default=sys.stdout)
@click.pass_obj
def write_gene_family_belns(manager, debug, output):
    """Write gene family namespace"""
    set_debug_param(debug)

    manager.write_gene_family_bel_namespace(output)


@main.command()
@click.option('-v', '--debug', count=True, help='Turn on debugging')
@click.pass_obj
def deploy_gene_family_belns(manager, debug):
    """Deploy gene family namespace"""
    set_debug_param(debug)

    manager.deploy_gene_family_bel_namespace()


@main.command()
@click.option('-v', '--debug', count=True, help='Turn on debugging')
@click.option('-o', '--output', type=click.File('w'), default=sys.stdout)
@click.pass_obj
def write_gene_family_bel(manager, debug, output):
    """Write gene family BEL script"""
    set_debug_param(debug)

    manager.write_gene_family_bel(output)


@main.command()
@click.option('-v', '--debug', count=True, help='Turn on debugging')
@click.pass_obj
def deploy_gene_family_bel(manager, debug):
    """Deploy gene family BEL script"""
    set_debug_param(debug)

    gene_namespace_url = manager.deploy_gene_bel_namespace()
    if gene_namespace_url:
        click.echo('Deployed gene namespace to: {}'.format(gene_namespace_url))

    family_namespace_url = manager.deploy_gene_family_bel_namespace()
    if family_namespace_url:
        click.echo('Deployed family namespace to {}'.format(family_namespace_url))

    family_bel_url = manager.deploy_gene_family_bel()
    if family_bel_url:
        click.echo('Deployed family BEL script to {}'.format(family_bel_url))


@main.command()
@click.option('-v', '--debug', count=True, help='Turn on debugging')
@click.pass_obj
def upload_gene_family_bel(manager, debug):
    """Uploads to PyBEL network store"""
    set_debug_param(debug)

    from pybel.manager import Manager as PyBELManager
    pybel_manager = PyBELManager(connection=manager.connection)

    log.info('exporting graph')
    graph = manager.export_gene_families_to_bel_graph()
    log.info('inserting graph')
    pybel_manager.insert_graph(graph)


if __name__ == '__main__':
    main()
