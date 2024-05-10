req:
	python3 -m pip install -r requirements.txt

data:
	./src/get_data.sh

load:
	./src/redis_loader.py --genus Glycine --species max --strain Wm82 \
        --chromosome-gff data/genome.gff3.gz \
        --gene-gff data/genes.gff3.gz \
        --gfa data/gfa.tsv.gz
