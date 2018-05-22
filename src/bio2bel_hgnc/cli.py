# -*- coding: utf-8 -*-

import sys

import click
import logging

from .manager import Manager

log = logging.getLogger(__name__)

main = Manager.get_cli()


@main.command()
@click.option('-o', '--output', type=click.File('w'), default=sys.stdout)
@click.pass_obj
def write_gene_belns(manager, output):
    """Write gene namespace"""
    manager.write_gene_bel_namespace(output)


@main.command()
@click.option('-o', '--output', type=click.File('w'), default=sys.stdout)
@click.pass_obj
def write_gene_family_belns(manager, output):
    """Write gene family namespace"""
    manager.write_gene_family_bel_namespace(output)


if __name__ == '__main__':
    main()
