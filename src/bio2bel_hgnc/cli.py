# -*- coding: utf-8 -*-

import logging

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
@click.pass_obj
def drop(manager, debug):
    """Drops the database"""
    set_debug_param(debug)
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


if __name__ == '__main__':
    main()
