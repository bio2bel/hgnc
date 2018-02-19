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


@click.group(help='Convert HGNC to BEL. Default connection at {}'.format(DEFAULT_CACHE_CONNECTION))
def main():
    pass


@main.command()
@click.option('-c', '--connection', help="Defaults to {}".format(DEFAULT_CACHE_CONNECTION))
@click.option('-v', '--debug', count=True, help="Turn on debugging.")
def populate(connection, debug):
    """Populate the database"""
    set_debug_param(debug)

    manager = Manager(connection=connection)
    log.info('connected to %s', manager.engine.url)

    manager.create_all()
    manager.populate()


@main.command()
@click.option('-c', '--connection', help='Defaults to {}'.format(DEFAULT_CACHE_CONNECTION))
@click.option('-v', '--debug', count=True, help="Turn on debugging.")
@click.option('-y', '--yes', is_flag=True, help="Automatically answer yes")
def drop(connection, debug, yes):
    """Drops the database"""
    set_debug_param(debug)

    if yes or click.confirm('Do you really want to delete the database?'):
        m = Manager(connection=connection)
        m.drop_all()


@main.command()
@click.option('-c', '--connection', help='Defaults to {}'.format(DEFAULT_CACHE_CONNECTION))
@click.option('-v', '--debug', count=True, help='Turn on debugging')
@click.option('-p', '--port')
@click.option('-h', '--host')
def web(connection, debug, port, host):
    """Run the web app"""
    set_debug_param(debug)

    from .web import create_app
    app = create_app(connection=connection, url='/')
    app.run(host=host, port=port, debug=debug)


@main.command()
@click.option('-c', '--connection', help='Defaults to {}'.format(DEFAULT_CACHE_CONNECTION))
@click.option('-v', '--debug', count=True, help='Turn on debugging')
@click.option('-o', '--output', type=click.File('w'), default=sys.stdout)
def write_gene_belns(connection, debug, output):
    """Write gene namespace"""
    set_debug_param(debug)

    manager = Manager(connection=connection)
    manager.write_gene_bel_namespace(output)


@main.command()
@click.option('-c', '--connection', help='Defaults to {}'.format(DEFAULT_CACHE_CONNECTION))
@click.option('-v', '--debug', count=True, help='Turn on debugging')
def deploy_gene_belns(connection, debug):
    """Deploy gene namespace"""
    set_debug_param(debug)

    manager = Manager(connection=connection)
    manager.deploy_gene_bel_namespace()


@main.command()
@click.option('-c', '--connection', help='Defaults to {}'.format(DEFAULT_CACHE_CONNECTION))
@click.option('-v', '--debug', count=True, help='Turn on debugging')
@click.option('-o', '--output', type=click.File('w'), default=sys.stdout)
def write_gene_family_belns(connection, debug, output):
    """Write gene family namespace"""
    set_debug_param(debug)

    manager = Manager(connection=connection)
    manager.write_gene_family_bel_namespace(output)


@main.command()
@click.option('-c', '--connection', help='Defaults to {}'.format(DEFAULT_CACHE_CONNECTION))
@click.option('-v', '--debug', count=True, help='Turn on debugging')
def deploy_gene_family_belns(connection, debug):
    """Deploy gene family namespace"""
    set_debug_param(debug)

    manager = Manager(connection=connection)
    manager.deploy_gene_family_bel_namespace()


@main.command()
@click.option('-c', '--connection', help='Defaults to {}'.format(DEFAULT_CACHE_CONNECTION))
@click.option('-v', '--debug', count=True, help='Turn on debugging')
@click.option('-o', '--output', type=click.File('w'), default=sys.stdout)
def write_gene_family_bel(connection, debug, output):
    """Write gene family BEL script"""
    set_debug_param(debug)

    manager = Manager(connection=connection)
    manager.write_gene_family_bel(output)


@main.command()
@click.option('-c', '--connection', help='Defaults to {}'.format(DEFAULT_CACHE_CONNECTION))
@click.option('-v', '--debug', count=True, help='Turn on debugging')
def deploy_gene_family_bel(connection, debug):
    """Deploy gene family BEL script"""
    set_debug_param(debug)

    manager = Manager(connection=connection)

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
@click.option('-c', '--connection', help='Defaults to {}'.format(DEFAULT_CACHE_CONNECTION))
@click.option('-v', '--debug', count=True, help='Turn on debugging')
def upload_gene_family_bel(connection, debug):
    """Uploads to PyBEL network store"""
    set_debug_param(debug)

    from pybel.manager import Manager as PyBELManager

    hgnc_manager = Manager(connection=connection)
    pybel_manager = PyBELManager(connection=connection)

    log.info('exporting graph')
    graph = hgnc_manager.export_gene_families_to_bel_graph()
    log.info('inserting graph')
    pybel_manager.insert_graph(graph)


if __name__ == '__main__':
    main()
