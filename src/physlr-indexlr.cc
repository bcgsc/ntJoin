// Convert linked-reads to minimizers using ntHash-2.0.0.
// Usage:  physlr-indexlr -k K -w W [-v] [-o file] file...
// Output: Each line of output is a barcode followed by a list of minimzers.
// Originally written for Physlr: (https://github.com/bcgsc/physlri
// Written by Vladimir Nikolic (schutzekatze) and Shaun Jackman (@sjackman)

#include "indexlr-workers.h"

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <getopt.h>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

// Read a FASTQ file and reduce each read to a set of minimizers
static void
minimizeReads(
    const std::string& ipath,
    const std::string& opath,
    const size_t k,
    const size_t w,
    const size_t t,
    const bool withPositions,
    const bool withStrands,
    const bool verbose)
{
	InputWorker inputWorker(ipath);
	OutputWorker outputWorker(opath, inputWorker);

	inputWorker.start();
	outputWorker.start();

	auto minimizeWorkers = std::vector<MinimizeWorker>(
	    t, MinimizeWorker(k, w, withPositions, withStrands, verbose, inputWorker, outputWorker));
	for (auto& worker : minimizeWorkers) {
		worker.start();
	}
	for (auto& worker : minimizeWorkers) {
		worker.join();
	}

	inputWorker.join();
	outputWorker.join();
}

static void
printErrorMsg(const std::string& progname, const std::string& msg)
{
	std::cerr << progname << ": " << msg << "\nTry 'physlr-indexlr --help' for more information.\n";
}

static void
printUsage(const std::string& progname)
{
	std::cout << "Usage:  " << progname
	          << "  -k K -w W [-v] [-o FILE] FILE...\n\n"
	             "  -k K       use K as k-mer size\n"
	             "  -w W       use W as sliding-window size\n"
	             "  --pos      include minimizer positions in the output\n"
	             "  --strand   include minimizer strand in the output\n"
	             "  -v         enable verbose output\n"
	             "  -o FILE    write output to FILE, default is stdout\n"
	             "  -t N       use N number of threads (default 1, max 5)\n"
	             "  --help     display this help and exit\n"
	             "  FILE       space separated list of FASTQ files\n";
}

int
main(int argc, char* argv[])
{
	auto progname = "physlr-indexlr";
	int c;
	int optindex = 0;
	static int help = 0;
	unsigned k = 0;
	unsigned w = 0;
	bool verbose = false;
	unsigned t = 1;
	bool failed = false;
	bool w_set = false;
	bool k_set = false;
	static int withPositions = 0;
	static int withStrands = 0;
	char* end = nullptr;
	std::string outfile("/dev/stdout");
	static const struct option longopts[] = { { "pos", no_argument, &withPositions, 1 },
		                                      { "strand", no_argument, &withStrands, 1 },
		                                      { "help", no_argument, &help, 1 },
		                                      { nullptr, 0, nullptr, 0 } };
	while ((c = getopt_long(argc, argv, "k:w:o:vt:", longopts, &optindex)) != -1) {
		switch (c) {
		case 0:
			break;
		case 'k':
			k_set = true;
			k = strtoul(optarg, &end, 10);
			break;
		case 'w':
			w_set = true;
			w = strtoul(optarg, &end, 10);
			break;
		case 'o':
			outfile.assign(optarg);
			break;
		case 'v':
			verbose = true;
			break;
		case 't':
			t = strtoul(optarg, &end, 10);
			if (t > 5) {
				t = 5;
				std::cerr << progname
				          << ": Using more than 5 threads does not scale, reverting to 5.\n";
			}
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
	} else if (!w_set) {
		printErrorMsg(progname, "missing option -- 'w'");
		failed = true;
	} else if (k == 0) {
		printErrorMsg(progname, "option has incorrect argument -- 'k'");
		failed = true;
	} else if (w == 0) {
		printErrorMsg(progname, "option has incorrect argument -- 'w'");
		failed = true;
	} else if (infiles.empty()) {
		printErrorMsg(progname, "missing file operand");
		failed = true;
	}
	if (failed) {
		exit(EXIT_FAILURE);
	}

	for (auto& infile : infiles) {
		minimizeReads(
		    infile == "-" ? "/dev/stdin" : infile,
		    outfile,
		    k,
		    w,
		    t,
		    withPositions,
		    withStrands,
		    verbose);
	}

	return 0;
}