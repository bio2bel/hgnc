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
def main():
    """HGNC to BEL"""


@main.command()
@click.option('-c', '--connection', help="Defaults to {}".format(DEFAULT_CACHE_CONNECTION))
@click.option('-v', '--debug', count=True, help="Turn on debugging.")
def populate(connection, debug):
    """Populate the database"""

    set_debug_param(debug)

    manager = Manager(connection=connection)
    manager.create_all()
    manager.populate()


@main.command()
@click.option('-c', '--connection', help='Defaults to {}'.format(DEFAULT_CACHE_CONNECTION))
@click.option('-v', '--debug', count=True, help="Turn on debugging.")
def drop(connection, debug):
    """Drops the database"""

    set_debug_param(debug)

    m = Manager(connection=connection)
    m.drop_all()


@main.command()
@click.option('-c', '--connection', help='Defaults to {}'.format(DEFAULT_CACHE_CONNECTION))
@click.option('-v', '--debug', is_flag=True)
@click.option('-p', '--port')
@click.option('-h', '--host')
def web(connection, debug, port, host):
    """Run the web app"""
    from .web import create_app
    app = create_app(connection=connection, url='/')
    app.run(host=host, port=port, debug=debug)


if __name__ == '__main__':
    main()
