# -*- coding: utf-8 -*-

import click

import pyhgnc
from .constants import DEFAULT_CACHE_CONNECTION


@click.group()
def main():
    """Bio2BEL HGNC Utilities"""


@main.command()
@click.option('-c', '--connection', help="Defaults to {}".format(DEFAULT_CACHE_CONNECTION))
def populate(connection):
    """Populate the database"""
    pyhgnc.update(connection=connection)


if __name__ == '__main__':
    main()
