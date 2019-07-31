CXXFLAGS=-g -O3 -Wall -Werror -Wno-error=unknown-pragmas -std=c++11 -pthread -I./include -I./vendor -fopenmp
LDLIBS=-lz

all: physlr-indexlr

clean:
	rm -f physlr-indexlr

# Check the C++ source code for errors.
lint: clang-format clang-tidy

# Check the C++ source code for errors with clang-tidy.
clang-tidy:
	clang-tidy -warnings-as-errors='*' include/*.h *.cc -- -std=c++11 -x c++ -I./include -I./vendor

# Check the C++ source code for white-space errors with clang-format.
clang-format:
	for i in *.cc include/*.h; do clang-format -style=file $$i >$$i.fixed; done
	for i in *.cc include/*.h; do diff -su $$i $$i.fixed && rm -f $$i.fixed; done
	if ls *.fixed; then exit 1; fi

# Fix white-space errors with clang-format.
fix:
	for i in *.cc; do clang-format -i -style=file $$i; done

check: check-physlr-index-fasta check-physlr-index-pos-fasta check-physlr-index-pos-strand-fasta check-physlr-index-fastq

check-physlr-index-fasta: all
	./physlr-indexlr -t16 -k100 -w5 data/mt.fa | diff -q - data/mt.physlr.tsv.good

check-physlr-index-pos-fasta: all
	./physlr-indexlr -t16 -k100 -w5 --pos data/mt.fa | diff -q - data/mt.pos.physlr.tsv.good

check-physlr-index-pos-strand-fasta: all
	./physlr-indexlr -t16 -k100 -w5 --pos --strand data/mt.fa | diff -q - data/mt.pos.strand.physlr.tsv.good

check-physlr-index-fastq: all
	./physlr-indexlr -t16 -k100 -w5 data/tiny.fq | diff -q - data/tiny.physlr.tsv.good
