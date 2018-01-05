# -*- coding: utf-8 -*-

import click

from .constants import DEFAULT_CACHE_CONNECTION
from .manager import Manager


@click.group()
def main():
    """HGNC to BEL"""


@main.command()
@click.option('-c', '--connection', help="Defaults to {}".format(DEFAULT_CACHE_CONNECTION))
def populate(connection):
    """Populate the database"""
    manager = Manager(connection=connection)
    manager.populate()


@main.command()
@click.option('-c', '--connection', help='Defaults to {}'.format(DEFAULT_CACHE_CONNECTION))
def drop(connection):
    """Drops the database"""
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
