# -*- coding: utf-8 -*-

import logging
import sys

import click

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
@click.pass_obj
def deploy_gene_belns(manager):
    """Deploy gene namespace"""
    manager.deploy_gene_bel_namespace()


@main.command()
@click.option('-o', '--output', type=click.File('w'), default=sys.stdout)
@click.pass_obj
def write_gene_family_belns(manager, output):
    """Write gene family namespace"""
    manager.write_gene_family_bel_namespace(output)


@main.command()
@click.pass_obj
def deploy_gene_family_belns(manager):
    """Deploy gene family namespace"""
    manager.deploy_gene_family_bel_namespace()


@main.command()
@click.option('-o', '--output', type=click.File('w'), default=sys.stdout)
@click.pass_obj
def write_gene_family_bel(manager, output):
    """Write gene family BEL script"""
    manager.write_gene_family_bel(output)


@main.command()
@click.pass_obj
def deploy_gene_family_bel(manager):
    """Deploy gene family BEL script"""
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
@click.pass_obj
def upload_gene_family_bel(manager):
    """Uploads to PyBEL network store"""
    from pybel.manager import Manager as PyBELManager
    pybel_manager = PyBELManager(connection=manager.connection)

    log.info('exporting graph')
    graph = manager.export_gene_families_to_bel_graph()
    log.info('inserting graph')
    pybel_manager.insert_graph(graph)


if __name__ == '__main__':
    main()
