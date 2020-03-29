# -*- coding: utf-8 -*-

"""SQLAlchemy models for Bio2BEL HGNC."""

from __future__ import annotations

from sqlalchemy import Column, ForeignKey, Integer, String, Table
from sqlalchemy.ext.declarative import DeclarativeMeta, declarative_base
from sqlalchemy.orm import relationship

from bio2bel.compath import CompathPathwayMixin, CompathProteinMixin
from .constants import MODULE_NAME

__all__ = [
    'Base',
    'GeneFamily',
    'HumanGene',
    'MouseGene',
    'RatGene',
    'human_mouse',
    'human_rat',
]

HUMAN_GENE_TABLE_NAME = f'{MODULE_NAME}_humanGene'
HUMAN_RAT_TABLE_NAME = f'{MODULE_NAME}_humanGene_ratGene'
RAT_GENE_TABLE_NAME = f'{MODULE_NAME}_ratGene'
HUMAN_MOUSE_TABLE_NAME = f'{MODULE_NAME}_humanGene_mouseGene'
MOUSE_GENE_TABLE_NAME = f'{MODULE_NAME}_mouseGene'
GENE_FAMILY_TABLE_NAME = f'{MODULE_NAME}_geneFamily'
GENE_TO_FAMILY_TABLE_NAME = f'{MODULE_NAME}_humanGene_geneFamily'

Base: DeclarativeMeta = declarative_base()

human_mouse = Table(
    HUMAN_MOUSE_TABLE_NAME,
    Base.metadata,
    Column('human_gene_id', Integer, ForeignKey(f'{HUMAN_GENE_TABLE_NAME}.id'), primary_key=True),
    Column('mouse_gene_id', Integer, ForeignKey(f'{MOUSE_GENE_TABLE_NAME}.id'), primary_key=True),
)

human_rat = Table(
    HUMAN_RAT_TABLE_NAME,
    Base.metadata,
    Column('human_gene_id', Integer, ForeignKey(f'{HUMAN_GENE_TABLE_NAME}.id'), primary_key=True),
    Column('rat_gene_id', Integer, ForeignKey(f'{RAT_GENE_TABLE_NAME}.id'), primary_key=True),
)

human_family = Table(
    GENE_TO_FAMILY_TABLE_NAME,
    Base.metadata,
    Column('human_gene_id', Integer, ForeignKey(f'{HUMAN_GENE_TABLE_NAME}.id'), primary_key=True),
    Column('gene_family_id', Integer, ForeignKey(f'{GENE_FAMILY_TABLE_NAME}.id'), primary_key=True),
)


class HumanGene(Base, CompathProteinMixin):
    """A SQLAlchemy model for a human gene."""

    __tablename__ = HUMAN_GENE_TABLE_NAME

    id = Column(Integer, primary_key=True)

    entrez_id = Column(String(255), doc='entrez id of the protein')
    hgnc_id = Column(String(255), doc='HGNC id of the protein')
    hgnc_symbol = Column(String(255), doc='HGN symbol of the protein')


class MouseGene(Base):
    """A SQLAlchemy model for a mouse gene."""

    __tablename__ = MOUSE_GENE_TABLE_NAME

    id = Column(Integer, primary_key=True)

    entrez_id = Column(String(255), doc='entrez id of the protein')
    mgi_id = Column(String(255), doc='MGI id of the protein')
    mgi_symbol = Column(String(255), doc='MGI symbol of the protein')

    human_genes = relationship(
        HumanGene,
        secondary=human_mouse,
        backref='mouse_genes',
    )


class RatGene(Base):
    """A SQLAlchemy model for an rat gene."""

    __tablename__ = RAT_GENE_TABLE_NAME

    id = Column(Integer, primary_key=True)

    entrez_id = Column(String(255), doc='entrez id of the protein')
    rgd_id = Column(String(255), doc='RGD id of the protein')
    rgd_symbol = Column(String(255), doc='RGD symbol of the protein')

    human_genes = relationship(
        HumanGene,
        secondary=human_rat,
        backref='rat_genes',
    )


class GeneFamily(CompathPathwayMixin, Base):
    """A SQLAlchemy model for an HGNC Gene family."""

    __tablename__ = GENE_FAMILY_TABLE_NAME

    id = Column(Integer, primary_key=True)

    identifier = Column(String(255), doc='HGNC gene family id of the protein')
    symbol = Column(String(255), doc='HGNC gene family symbol of the protein')
    name = Column(String(255), doc='HGNC gene family name of the protein')

    proteins = relationship(
        HumanGene,
        secondary=human_family,
        backref='gene_families',
    )
