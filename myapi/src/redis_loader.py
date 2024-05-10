#!/usr/bin/env python3
"""
Author : Ken Youens-Clark <kyclark@gmail.com>
Date   : 2024-05-09
Purpose: Rock the Casbah
"""

import argparse
import codecs
import csv
import gffutils
import gzip
import os
from collections import defaultdict
from urllib.request import urlopen, urlparse
from redisearch import RediSearchLoader
from configparser import ConfigParser
from typing import List, NamedTuple, TextIO


class Args(NamedTuple):
    """Command-line arguments"""

    load_type: str
    sequence_types: List[str]
    save: bool
    chunk_size: int
    genus: str
    species: str
    strain: str
    gene_gff: TextIO
    chr_gff: TextIO
    gfa: TextIO
    redis_host: str
    redis_port: int
    redis_db: str
    redis_password: str


# --------------------------------------------------
def get_args() -> Args:
    """Get command-line arguments"""

    # Read configuration for global settings
    config_file = "./config.ini"
    config = {}
    if os.path.isfile(config_file):
        conf = ConfigParser(interpolation=None)
        conf.read(config_file)
        config = conf["DEFAULT"]

    parser = argparse.ArgumentParser(
        description="Rock the Casbah",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--load-type",
        choices=["new", "reload", "append"],
        default="new",
        help="Load type",
    )

    parser.add_argument(
        "--sequence-types",
        choices=["chromosome", "supercontig", "chloroplast", "mitochondrion"],
        default="chromosome",
        help="Sequence types",
        nargs="+",
    )

    parser.add_argument(
        "--save",
        action="store_true",
        help="Save to disk",
    )

    parser.add_argument(
        "--chunk-size",
        type=int,
        default=100,
        help="Chunk size for Redis batch processing",
    )

    parser.add_argument(
        "--genus",
        required=True,
        help="Genus",
    )

    parser.add_argument(
        "--species",
        required=True,
        help="Species",
    )

    parser.add_argument("--strain", help="Strain")

    parser.add_argument(
        "--gene-gff",
        required=True,
        type=argparse.FileType("rt"),
        help="GFF(.gz) file containing gene records",
    )

    parser.add_argument(
        "--chromosome-gff",
        required=True,
        type=argparse.FileType("rt"),
        help="GFF(.gz) file containing chromosome/supercontig records",
    )

    parser.add_argument(
        "--gfa",
        required=True,
        type=argparse.FileType("rt"),
        help="GFA(.gz) file containing gene-gene family associations",
    )

    parser.add_argument(
        "--rhost", default=config.get("redis_host"), help="Redis host"
    )

    parser.add_argument(
        "--rport",
        type=int,
        default=config.get("redis_port"),
        help="Redis port",
    )

    parser.add_argument(
        "--rdb", default=config.get("redis_db"), help="Redis db"
    )

    parser.add_argument(
        "--rpassword",
        default=config.get("redis_password"),
        help="Redis password",
    )

    args = parser.parse_args()

    return Args(
        load_type=args.load_type,
        sequence_types=args.sequence_types,
        save=args.save,
        chunk_size=args.chunk_size,
        genus=args.genus,
        species=args.species,
        strain=args.strain,
        gene_gff=args.gene_gff,
        chr_gff=args.chromosome_gff,
        gfa=args.gfa,
        redis_host=args.rhost,
        redis_port=args.rport,
        redis_db=args.rdb,
        redis_password=args.rpassword,
    )


# --------------------------------------------------
def main() -> None:
    """Make a jazz noise here"""

    args = get_args()

    loader = RediSearchLoader(
        host=args.redis_host,
        port=args.redis_port,
        db=args.redis_db,
        password=args.redis_password,
        load_type=args.load_type,
        sequence_types=args.sequence_types,
        no_save=not args.save,
        chunk_size=args.chunk_size,
    )

    chromosome_names = transfer_chromosomes(
        loader, args.genus, args.species, args.chromosome_gff
    )

    transfer_genes(loader, args.gene_gff, args.gfa, chromosome_names)


# --------------------------------------------------
def transfer_chromosomes(loader, genus, species, chromosome_gff):
    """
    Loads chromosomes from a GFF file into a RediSearch database.

    Parameters:
      loader (RediSearchLoader): The loader to use to load data into
        RediSearch.
      genus (str): The genus of the chromosomes being loaded.
      species (str): The species of the chromosomes being loaded.
      chromosome_gff (str): The local path or URL to the GFF to load chromosomes from.

    Returns:
      set[str]: A set containing the names of all the chromosomes that were
        loaded.
    """

    # create chromosome SQLLite database from chromosomal GFF file
    gffchr_db = gffutils.create_db(
        chromosome_gff,
        ":memory:",
        force=True,
        keep_order=True,
    )

    # index the chromosomes
    chromosome_names = set()
    for chr in gffchr_db.features_of_type(
        ("chromosome", "supercontig"), order_by="attributes"
    ):
        name = chr.seqid
        length = chr.end
        chromosome_names.add(name)
        loader.indexChromosome(name, length, genus, species)

    return chromosome_names


def transfer_genes(loader, gene_gff, gfa, chromosome_names):
    """
    Loads genes from a GFF file into a RediSearch database.

    Parameters:
      loader (RediSearchLoader): The loader to use to load data into
        RediSearch.
      gene_gff (str): The local path or URL to the GFF to load genes from.
      gfa (str): The local path or URL to a GFA file containing gene family
        associations for the genes being loaded.
      chromosome_names (set[str]): A containing the names of all the chromosomes
        that have been loaded.
    """

    # create gene SQLLite database from gene GFF file
    gffgene_db = gffutils.create_db(
        gene_gff, ":memory:", force=True, keep_order=True
    )

    # index all the genes in the db
    gene_lookup = dict()
    chromosome_genes = defaultdict(list)
    for gffgene in gffgene_db.features_of_type("gene", order_by="attributes"):
        chr_name = gffgene.seqid
        if chr_name in chromosome_names:
            strand = 0
            if gffgene.strand == "+":
                strand = 1
            if gffgene.strand == "-":
                strand = -1
            gene = {
                "name": gffgene.id,
                "fmin": gffgene.start,
                "fmax": gffgene.end,
                "strand": strand,
                "family": "",
            }
            gene_lookup[gffgene.id] = gene
            chromosome_genes[chr_name].append(gene)

    # deal with family assignments (for non-orphans) from GFA
    with open(gfa, "rb") if urlparse(gfa).scheme == "" else urlopen(
        gfa
    ) as fileobj:
        tsv = gzip.GzipFile(fileobj=fileobj) if gfa.endswith("gz") else fileobj
        for line in csv.reader(
            codecs.iterdecode(tsv, "utf-8"), delimiter="\t"
        ):
            # skip comment and metadata lines
            if line[0].startswith("#") or line[0] == "ScoreMeaning":
                continue
            gene_id = line[0]
            if gene_id in gene_lookup:
                gene = gene_lookup[gene_id]
                genefamily_id = line[1]
                gene["family"] = genefamily_id

    # index the genes
    for chr_name, genes in chromosome_genes.items():
        loader.indexChromosomeGenes(chr_name, genes)


# --------------------------------------------------
if __name__ == "__main__":
    main()
