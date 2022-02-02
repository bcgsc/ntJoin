CXXFLAGS=-O3 -Wall -Wextra -Werror -Wno-unknown-pragmas -std=c++11 -pthread -Ibtllib/include

all: indexlr \
	btllib_compile

clean:
	rm -f indexlr btllib/python/*

# Compile btllib wrappers
btllib_compile: btllib/python/btllib.py

btllib/python/btllib.py:
	cd btllib && ./compile-wrappers

# Check the C++ source code for errors.
lint: clang-format clang-tidy

# Check the C++ source code for errors with clang-tidy.
clang-tidy:
	clang-tidy -warnings-as-errors='*' *.cpp -- -std=c++11 -x c++ -Ibtllib/include

# Check the C++ source code for white-space errors with clang-format.
clang-format:
	for i in *.cpp; do clang-format -style=file $$i >$$i.fixed; done
	for i in *.cpp; do diff -su $$i $$i.fixed && rm -f $$i.fixed; done
	if ls *.fixed; then exit 1; fi

# Fix white-space errors with clang-format.
fix:
	for i in *.cpp; do clang-format -i -style=file $$i; done

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
