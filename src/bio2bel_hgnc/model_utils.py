# -*- coding: utf-8 -*-

"""Bio2BEL HGNC Model utilities."""

from pybel.constants import MIRNA, PROTEIN, RNA
from pybel.dsl import gene as gene_dsl, mirna as mirna_dsl, protein as protein_dsl, rna as rna_dsl
from pybel.dsl.nodes import CentralDogma

from .constants import encodings

__all__ = [
    'gene_to_bel',
    'gene_to_rna_to_bel',
    'gene_to_mirna_to_bel',
    'gene_to_protein_to_bel',
    'family_to_bel',
    'uniprot_to_bel',
]


def gene_to_bel(human_gene, func=None, variants=None) -> CentralDogma:
    """Convert a Gene to a PyBEL gene.

    :param bio2bel_hgnc.models.HumanGene human_gene:  A Gene model
    :rtype: pybel.dsl.gene
    """
    if func == PROTEIN:
        return gene_to_protein_to_bel(human_gene)
    elif func == RNA:
        return gene_to_rna_to_bel(human_gene)
    elif func == MIRNA:
        return gene_to_mirna_to_bel(human_gene)

    return gene_dsl(
        namespace='hgnc',
        name=str(human_gene.symbol),
        identifier=str(human_gene.identifier),
        variants=variants,
    )


def gene_to_rna_to_bel(gene) -> rna_dsl:
    """Converts a Gene to a PyBEL gene

    :param bio2bel_hgnc.models.HumanGene gene:  A Gene model
    :rtype: pybel.dsl.gene
    """
    return rna_dsl(
        namespace='hgnc',
        name=str(gene.symbol),
        identifier=str(gene.identifier)
    )


def gene_to_mirna_to_bel(gene):
    """Converts a Gene to a PyBEL miRNA

    :param bio2bel_hgnc.models.HumanGene gene:  A Gene model
    :rtype: pybel.dsl.mirna
    """
    return mirna_dsl(
        namespace='hgnc',
        name=str(gene.symbol),
        identifier=str(gene.identifier)
    )


def gene_to_protein_to_bel(gene):
    """Converts a Gene to a PyBEL gene

    :param bio2bel_hgnc.models.HumanGene gene:  A Gene model
    :rtype: pybel.dsl.gene
    """
    return protein_dsl(
        namespace='hgnc',
        name=str(gene.symbol),
        identifier=str(gene.identifier)
    )


def family_to_bel(family, func=None):
    """Converts a Gene Family model to a PyBEL gene

    :param bio2bel_hgnc.models.GeneFamily family: A Gene Family model
    :rtype: pybel.dsl.gene
    """
    if func == PROTEIN:
        dsl = protein_dsl
    elif func == RNA:
        dsl = rna_dsl
    elif func == MIRNA:
        dsl = mirna_dsl
    else:
        dsl = gene_dsl

    return dsl(
        namespace='hgnc.genefamily',
        identifier=str(family.family_identifier),
        name=str(family.family_name)
    )


def uniprot_to_bel(uniprot) -> protein_dsl:
    """

    :param bio2bel_hgnc.models.UniProt uniprot:
    :return:
    """
    return protein_dsl(
        namespace='uniprot',
        name=str(uniprot.uniprotid),
        identifier=str(uniprot.uniprotid),
    )


def add_central_dogma(graph, human_gene):
    """Add the corresponding protein and"""
    encoding = encodings.get(human_gene.locus_type, 'GRP')

    if 'M' in encoding:
        mirna = gene_to_mirna_to_bel(human_gene)
        graph.add_transcription(gene_to_bel(human_gene), mirna)

    elif 'R' in encoding:
        rna = gene_to_rna_to_bel(human_gene)
        graph.add_transcription(gene_to_bel(human_gene), rna)

        if 'P' in encoding:
            graph.add_translation(rna, gene_to_protein_to_bel(human_gene))
