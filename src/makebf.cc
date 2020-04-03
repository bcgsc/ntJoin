#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#if _OPENMP
#include <omp.h>
#endif

#include "btl_bloomfilter/BloomFilter.hpp"
#include "btl_bloomfilter/vendor/ntHashIterator.hpp"

static void
printErrorMsg(const std::string& progname, const std::string& msg)
{
	std::cerr << progname << ": " << msg << "\nTry 'physlr-makebf --help' for more information.\n";
}

static void
printUsage(const std::string& progname)
{
	std::cout << "Usage:  " << progname
	          << "  -k K [-v] [-o FILE] FILE...\n\n"
	             "  -k K       use K as k-mer size\n"
	             "  -b B       Bloom filter size in Bytes\n"
	             "  -v         enable verbose output\n"
	             "  -o FILE    write Bloom filter to FILE [required]\n"
	             "  -t N       use N number of threads [1]\n"
	             "  --help     display this help and exit\n"
	             "  FILE       space separated list of ntHits tsv output files\n";
}

static void
printBloomStats(BloomFilter& bloom, ostream& os)
{
	os << "Bloom filter stats:"
	   << "\n\t#counters               = " << bloom.getFilterSize()
	   << "\n\t#size (B)               = " << bloom.sizeInBytes()
	   << "\n\tpopcount                = " << bloom.getPop()
	   << "\n\tFPR                     = " << setprecision(3) << 100.f * bloom.getFPR() << "%"
	   << "\n";
}

int
main(int argc, char* argv[])
{

	auto progname = "physlr-makebf";
	int c;
	int optindex = 0;
	static int help = 0;
	unsigned k = 0;
	uint64_t filterSize = 0;
	bool verbose = false;
	unsigned t = 1;
	unsigned hashNum = 1;
	bool failed = false;
	bool k_set = false;
	char* end = nullptr;
	std::string outfile;
	static const struct option longopts[] = { { "help", no_argument, &help, 1 },
		                                      { nullptr, 0, nullptr, 0 } };
	while ((c = getopt_long(argc, argv, "b:k:o:vt:", longopts, &optindex)) != -1) {
		switch (c) {
		case 0:
			break;
		case 'b':
			filterSize = strtoul(optarg, &end, 10) * 8;
			break;
		case 'k':
			k_set = true;
			k = strtoul(optarg, &end, 10);
			break;
		case 'o':
			outfile.assign(optarg);
			break;
		case 'v':
			verbose = true;
			break;
		case 't':
			t = strtoul(optarg, &end, 10);
			break;
		default:
			exit(EXIT_FAILURE);
		}
	}
	std::vector<std::string> infiles(&argv[optind], &argv[argc]);
	if (argc < 2) {
		printUsage(progname);
		exit(EXIT_FAILURE);
	}
	if (help != 0) {
		printUsage(progname);
		exit(EXIT_SUCCESS);
	} else if (!k_set) {
		printErrorMsg(progname, "missing option -- 'k'");
		failed = true;
	} else if (k == 0) {
		printErrorMsg(progname, "option has incorrect argument -- 'k'");
		failed = true;
	} else if (filterSize == 0) {
		printErrorMsg(progname, "option has incorrect argument -- 'b");
		failed = true;
	} else if (infiles.empty()) {
		printErrorMsg(progname, "missing file operand");
		failed = true;
	}
	if (failed) {
		exit(EXIT_FAILURE);
	}

	std::vector<std::string> kmerVect;
	kmerVect.reserve(10000);
	if (verbose) {
		std::cerr << "Collecting Kmers to insert into bloom filter" << std::endl;
	}
	for (auto& infile : infiles) {
		infile == "-" ? "/dev/stdin" : infile;
		std::ifstream infileStream(infile);
		std::string kmer;
		int count;
		while (infileStream >> kmer >> count) {
			kmerVect.push_back(kmer);
		}
	}

	BloomFilter bloomFilter(filterSize, hashNum, k);
	if (verbose) {
		std::cerr << "Made Bloom filter with:\n"
		          << "kmer size                 = " << k << "\n"
		          << "Bloom filter size         = " << filterSize << "\n"
		          << std::endl;
		std::cerr << "Inserting " << kmerVect.size() << " kmers into Bloom filter"
		          << "Using " << t << " threads." << std::endl;
	}

	uint64_t counter = 0;
#pragma omp parallel for num_threads(t)
	for (auto vectIt = kmerVect.begin(); vectIt < kmerVect.end(); ++vectIt) {
		ntHashIterator itr(*vectIt, hashNum, k);
		while (itr != ntHashIterator::end()) {
			bloomFilter.insert(*itr);
			++itr;
		}
		if (verbose) {
#pragma omp critical
			{
				counter++;
				if (counter % 100000 == 0) {
					std::cerr << "Processed " << counter << " kmers." << std::endl;
				}
			}
		}
	}

	printBloomStats(bloomFilter, std::cerr);
	bloomFilter.storeFilter(outfile);
}
