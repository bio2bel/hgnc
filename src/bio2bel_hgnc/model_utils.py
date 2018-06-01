# -*- coding: utf-8 -*-

"""Bio2BEL HGNC Model utilities."""

from pybel.dsl import gene as gene_dsl, mirna as mirna_dsl, protein as protein_dsl, rna as rna_dsl

__all__ = [
    'gene_to_bel',
    'gene_to_rna_to_bel',
    'gene_to_mirna_to_bel',
    'gene_to_protein_to_bel',
    'family_to_bel',
    'uniprot_to_pybel',
]


def gene_to_bel(gene):
    """Converts a Gene to a PyBEL gene

    :param bio2bel_hgnc.models.HumanGene gene:  A Gene model
    :rtype: pybel.dsl.gene
    """
    return gene_dsl(
        namespace='HGNC',
        name=str(gene.symbol),
        identifier=str(gene.identifier)
    )


def gene_to_rna_to_bel(gene):
    """Converts a Gene to a PyBEL gene

    :param bio2bel_hgnc.models.HumanGene gene:  A Gene model
    :rtype: pybel.dsl.gene
    """
    return rna_dsl(
        namespace='HGNC',
        name=str(gene.symbol),
        identifier=str(gene.identifier)
    )


def gene_to_mirna_to_bel(gene):
    """Converts a Gene to a PyBEL miRNA

    :param bio2bel_hgnc.models.HumanGene gene:  A Gene model
    :rtype: pybel.dsl.mirna
    """
    return mirna_dsl(
        namespace='HGNC',
        name=str(gene.symbol),
        identifier=str(gene.identifier)
    )


def gene_to_protein_to_bel(gene):
    """Converts a Gene to a PyBEL gene

    :param bio2bel_hgnc.models.HumanGene gene:  A Gene model
    :rtype: pybel.dsl.gene
    """
    return protein_dsl(
        namespace='HGNC',
        name=str(gene.symbol),
        identifier=str(gene.identifier)
    )


def family_to_bel(family):
    """Converts a Gene Family model to a PyBEL gene

    :param bio2bel_hgnc.models.GeneFamily family: A Gene Family model
    :rtype: pybel.dsl.gene
    """
    return gene_dsl(
        namespace='GFAM',
        identifier=str(family.family_identifier),
        name=str(family.family_name)
    )


def uniprot_to_pybel(uniprot):
    """

    :param bio2bel_hgnc.models.UniProt uniprot:
    :return:
    """
    return protein_dsl(
        namespace='UNIPROT',
        identifier=str(uniprot.uniprotid)
    )
