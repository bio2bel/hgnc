# -*- coding: utf-8 -*-

"""Bio2BEL HGNC Model utilities."""

from typing import List, Optional

from pybel import BELGraph
from pybel.constants import MIRNA, PROTEIN, RNA
from pybel.dsl import CentralDogma, FUNC_TO_DSL, Variant, gene as gene_dsl, protein as protein_dsl
from .constants import ENCODINGS, HGNC, HGNC_GENE_FAMILY, UNIPROT
from .models import GeneFamily, HumanGene, UniProt

__all__ = [
    'gene_to_bel',
    'family_to_bel',
    'uniprot_to_bel',
]


def gene_to_bel(human_gene: HumanGene, func: Optional[str] = None,
                variants: Optional[List[Variant]] = None) -> CentralDogma:
    """Convert a Gene to a PyBEL gene."""
    dsl = FUNC_TO_DSL[func] if func else gene_dsl

    rv = dsl(
        namespace=HGNC,
        name=human_gene.symbol,
        identifier=str(human_gene.identifier),
    )

    if variants is not None:
        return rv.with_variants(variants)

    return rv


def family_to_bel(family: GeneFamily, func: Optional[str] = None,
                  variants: Optional[List[Variant]] = None) -> CentralDogma:
    """Convert a Gene Family model to a PyBEL gene."""
    dsl = FUNC_TO_DSL[func] if func else gene_dsl

    rv = dsl(
        namespace=HGNC_GENE_FAMILY,
        identifier=str(family.family_identifier),
        name=str(family.family_name)
    )

    if variants is not None:
        return rv.with_variants(variants)

    return rv


def uniprot_to_bel(uniprot: UniProt) -> protein_dsl:
    """Convert a UniProt model to BEL."""
    return protein_dsl(
        namespace=UNIPROT,
        name=str(uniprot.uniprotid),
        identifier=str(uniprot.uniprotid),
    )


def add_central_dogma(graph: BELGraph, human_gene: HumanGene) -> None:
    """Add the corresponding protein and/or RNA."""
    encoding = ENCODINGS.get(human_gene.locus_type, 'GRP')

    if 'M' in encoding:
        mirna = gene_to_bel(human_gene, func=MIRNA)
        graph.add_transcription(gene_to_bel(human_gene), mirna)

    elif 'R' in encoding:
        rna = gene_to_bel(human_gene, func=RNA)
        graph.add_transcription(gene_to_bel(human_gene), rna)

        if 'P' in encoding:
            graph.add_translation(rna, gene_to_bel(human_gene, func=PROTEIN))
