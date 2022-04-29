
all: indexlr \
	btllib_compile

clean:
	rm -f indexlr btllib/python/*

# Compile btllib wrappers
btllib_compile: btllib/python/btllib.py

btllib/python/btllib.py:
	cd btllib && ./compile-wrappers

indexlr: btllib/recipes/indexlr.cpp
	btllib/compile
	cp btllib/install/bin/indexlr .


check: check-index-fasta check-index-pos-fasta check-index-pos-strand-fasta check-index-pos-strand-seq-fasta check-index-fastq

check-index-fasta: all
	./indexlr --long -t16 -k100 -w5 --id data/mt.fa | diff -q - data/mt.physlr.tsv.good

check-index-pos-fasta: all
	./indexlr --long -t16 -k100 -w5 --id --pos data/mt.fa | diff -q - data/mt.pos.physlr.tsv.good

check-index-pos-strand-fasta: all
	./indexlr --long -t16 -k100 -w5 --id --pos --strand data/mt.fa | diff -q - data/mt.pos.strand.physlr.tsv.good

check-index-pos-strand-seq-fasta: all
	./indexlr --long -t16 -k100 -w5 --id --pos --strand --seq data/mt.fa | diff -q - data/mt.pos.strand.seq.physlr.tsv.good

check-index-fastq: all
	./indexlr --long -t16 -k100 -w5 --id --bx data/tiny.fq | diff -q - data/tiny.physlr.tsv.good
