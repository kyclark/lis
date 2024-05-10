#!/usr/bin/env bash

OUTDIR=data

[[ ! -d "$OUTDIR" ]] && mkdir -p "$OUTDIR"

cd "$OUTDIR"

wget https://github.com/legumeinfo/microservices/raw/main/tests/data/genome.gff3.gz 
wget https://github.com/legumeinfo/microservices/raw/main/tests/data/genes.gff3.gz 
wget https://github.com/legumeinfo/microservices/raw/main/tests/data/gfa.tsv.gz

echo "Done, see \"$OUTDIR\""
