# -*- coding: utf-8 -*-

"""Bio2BEL HGNC Model utilities."""

from typing import List, Optional

from pybel import BELGraph
from pybel.constants import MIRNA, PROTEIN, RNA
from pybel.dsl import CentralDogma, FUNC_TO_DSL, Variant, gene as gene_dsl, protein as protein_dsl
from .constants import encodings
from .models import GeneFamily, HumanGene, UniProt

__all__ = [
    'gene_to_bel',
    'family_to_bel',
    'uniprot_to_bel',
]


def gene_to_bel(human_gene: HumanGene, func: Optional[str] = None,
                variants: Optional[List[Variant]] = None) -> CentralDogma:
    """Convert a Gene to a PyBEL gene.

    :param human_gene:  A Gene model
    :rtype: pybel.dsl.gene
    """
    dsl = FUNC_TO_DSL[func] if func else gene_dsl

    rv = dsl(
        namespace='hgnc',
        name=str(human_gene.symbol),
        identifier=str(human_gene.identifier),
    )

    if variants is not None:
        return rv.with_variants(variants)

    return rv


def family_to_bel(family: GeneFamily, func: Optional[str] = None,
                  variants: Optional[List[Variant]] = None) -> CentralDogma:
    """Convert a Gene Family model to a PyBEL gene.

    :param family: A Gene Family model
    :rtype: pybel.dsl.gene
    """
    dsl = FUNC_TO_DSL[func] if func else gene_dsl

    rv = dsl(
        namespace='hgnc.genefamily',
        identifier=str(family.family_identifier),
        name=str(family.family_name)
    )

    if variants is not None:
        return rv.with_variants(variants)

    return rv


def uniprot_to_bel(uniprot: UniProt) -> protein_dsl:
    """Convert the uniprot model to BEL.

    :param bio2bel_hgnc.models.UniProt uniprot:
    """
    return protein_dsl(
        namespace='uniprot',
        name=str(uniprot.uniprotid),
        identifier=str(uniprot.uniprotid),
    )


def add_central_dogma(graph: BELGraph, human_gene: HumanGene):
    """Add the corresponding protein and/or RNA."""
    encoding = encodings.get(human_gene.locus_type, 'GRP')

    if 'M' in encoding:
        mirna = gene_to_bel(human_gene, func=MIRNA)
        graph.add_transcription(gene_to_bel(human_gene), mirna)

    elif 'R' in encoding:
        rna = gene_to_bel(human_gene, func=RNA)
        graph.add_transcription(gene_to_bel(human_gene), rna)

        if 'P' in encoding:
            graph.add_translation(rna, gene_to_bel(human_gene, func=PROTEIN))
