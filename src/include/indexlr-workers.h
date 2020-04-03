#ifndef INDEXLR_WORKERS_H
#define INDEXLR_WORKERS_H

#include "btl_bloomfilter/BloomFilter.hpp"
#include "indexlr-buffer.h"
#include "indexlr-minimize.h"

#include "kseq.h" // NOLINT
#include <zlib.h>
KSEQ_INIT(gzFile, gzread) // NOLINT

#include <atomic>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

// A read from input file with a number associated with it
// where first read is 0, second 1, etc...
struct Read
{
	size_t num = 0;
	std::string id;
	std::string barcode;
	std::string sequence;
};

// Result stores barcodes and minimizers of a block of reads
// so that it can be written to output file in one go
struct Result
{
	size_t num = 0;
	std::string barcodesAndMinimizers;
	size_t lastNum = 0;

	Result() { barcodesAndMinimizers.reserve(BLOCK_SIZE * 1024); }
};

class InputWorker;
class MinimizeWorker;
class OutputWorker;

class InputWorker
{

  public:
	explicit InputWorker(std::string ipath)
	  : ipath(std::move(ipath))
	{}

	void start() { t = std::thread(doWork, this); }

	void join() { t.join(); }

	std::atomic<bool> allRead{ false };
	size_t inputNum = 0;
	std::atomic<bool> fasta{ false };
	InputBuffer<Block<Read>> buffer;

  private:
	const std::string ipath;

	inline void work();

	static void doWork(InputWorker* worker) { worker->work(); }

	std::thread t;
};

class MinimizeWorker
{

  public:
	MinimizeWorker(
	    size_t k,
	    size_t w,
	    bool withRepeat,
	    bool withSolid,
	    bool withPositions,
	    bool withStrands,
	    bool verbose,
	    const BloomFilter& repeatBF,
	    const BloomFilter& solidBF,
	    InputWorker& inputWorker,
	    OutputWorker& outputWorker)
	  : k(k)
	  , w(w)
	  , withRepeat(withRepeat)
	  , withSolid(withSolid)
	  , withPositions(withPositions)
	  , withStrands(withStrands)
	  , verbose(verbose)
	  , repeatBF(repeatBF)
	  , solidBF(solidBF)
	  , inputWorker(inputWorker)
	  , outputWorker(outputWorker)
	{}

	MinimizeWorker(const MinimizeWorker& worker)
	  : k(worker.k)
	  , w(worker.w)
	  , withRepeat(worker.withRepeat)
	  , withSolid(worker.withSolid)
	  , withPositions(worker.withPositions)
	  , withStrands(worker.withStrands)
	  , verbose(worker.verbose)
	  , repeatBF(worker.repeatBF)
	  , solidBF(worker.solidBF)
	  , inputWorker(worker.inputWorker)
	  , outputWorker(worker.outputWorker)
	{}

	MinimizeWorker(MinimizeWorker&& worker) noexcept
	  : k(worker.k)
	  , w(worker.w)
	  , withRepeat(worker.withRepeat)
	  , withSolid(worker.withRepeat)
	  , withPositions(worker.withPositions)
	  , withStrands(worker.withStrands)
	  , verbose(worker.verbose)
	  , repeatBF(worker.repeatBF)
	  , solidBF(worker.solidBF)
	  , inputWorker(worker.inputWorker)
	  , outputWorker(worker.outputWorker)
	{}

	MinimizeWorker& operator=(const MinimizeWorker& worker) = delete;
	MinimizeWorker& operator=(MinimizeWorker&& worker) = delete;

	~MinimizeWorker() = default;

	void start() { t = std::thread(doWork, this); }

	void join() { t.join(); }

  private:
	size_t k = 0;
	size_t w = 0;
	bool withRepeat = false;
	bool withSolid = false;
	bool withPositions = false;
	bool withStrands = false;
	bool verbose = false;
	const BloomFilter& repeatBF;
	const BloomFilter& solidBF;
	InputWorker& inputWorker;
	OutputWorker& outputWorker;

	inline void work();

	static void doWork(MinimizeWorker* worker) { worker->work(); }

	std::thread t;
};

class OutputWorker
{

  public:
	OutputWorker(std::string opath, const InputWorker& inputWorker)
	  : opath(std::move(opath))
	  , inputWorker(inputWorker)
	{
		ofs.open(this->opath);
	}

	void start() { t = std::thread(doWork, this); }

	void join() { t.join(); }

	OutputBuffer<Result> buffer;

  private:
	const std::string opath;
	std::ofstream ofs;
	const InputWorker& inputWorker;

	inline void work();

	static void doWork(OutputWorker* worker) { worker->work(); }

	std::thread t;
};

inline void
InputWorker::work()
{
	gzFile fp = gzopen(ipath.c_str(), "r");

	char c;
	c = gzgetc(fp);
	gzungetc(c, fp);
	if (c == '>') {
		fasta = true;
	} else {
		fasta = false;
	}

	if (gzeof(fp)) {
		std::cerr << "physlr-indexlr: error: Empty input file: " << ipath << '\n';
		exit(EXIT_FAILURE);
	} else {
		kseq_t* seq = kseq_init(fp);

		bool done = false;
		while (true) {
			size_t currentNum = inputNum;
			Block<Read>& reads = buffer.getWriteAccess(currentNum);

			reads.dataCounter = 0;
			for (auto& read : reads.data) {
				if (kseq_read(seq) < 0) {
					done = true;
					break;
				}

				read.id = seq->name.l > 0 ? seq->name.s : "";
				read.barcode = seq->comment.l > 0 ? seq->comment.s : "";
				read.sequence = seq->seq.l > 0 ? seq->seq.s : "";

				read.num = inputNum;
				++inputNum;
				++reads.dataCounter;
			}
			if (reads.dataCounter > 0) {
				reads.num = reads.data[0].num;
			} else {
				reads.num = inputNum - 1;
			}

			if (done) {
				allRead = true;
				buffer.releaseWriteAccess(currentNum);
				if (buffer.elements() == 0) {
					buffer.close();
				}
				break;
			}
			buffer.releaseWriteAccess(currentNum);
		}

		kseq_destroy(seq);
	}

	gzclose(fp);
}

inline void
MinimizeWorker::work()
{
	Block<Read> reads;
	std::stringstream ss;
	Result result;
	while (!inputWorker.allRead || inputWorker.buffer.elements() != 0) {
		inputWorker.buffer.read(reads);
		if (inputWorker.buffer.isClosed()) {
			break;
		}
		if (inputWorker.allRead && inputWorker.buffer.elements() == 0) {
			inputWorker.buffer.close();
		}

		ss.str("");
		for (size_t i = 0; i < reads.dataCounter; i++) {
			assert(i < sizeof(reads.data) / sizeof(reads.data[0]));
			Read& read = reads.data[i];

			if (startsWith(read.barcode, "BX:Z:")) {
				auto pos = read.barcode.find(' ');
				if (pos != std::string::npos) {
					read.barcode.erase(pos);
				}
				read.barcode.erase(0, 5);
			} else if (inputWorker.fasta) {
				// For FASTA, use the sequence ID.
				read.barcode = read.id;
			} else {
				// No barcode tag is present. Check for stLFR barcode within read.id.
				read.barcode = "NA";
				size_t sharpPos = read.id.find('#');
				if (sharpPos != std::string::npos) {
					size_t slashPos = read.id.rfind('/');
					if (slashPos > sharpPos) {
						read.barcode = read.id.substr(sharpPos + 1, slashPos - 1 - sharpPos);
					}
				}
			}

			if (read.sequence.size() < k) {
				if (verbose) {
					std::stringstream ss;
					ss << "physlr-indexlr: warning: Skip read " << (read.num + 1) << " on line "
					   << (read.num + 1) * 4 - 2 << "; k > read length "
					   << "(k = " << k << ", read length = " << read.sequence.size() << ")\n";
					std::cerr << ss.str();
				}
			}

			auto hashes = hashKmers(read.sequence, k);
			if (w > hashes.size()) {
				if (verbose) {
					std::stringstream ss;
					ss << "physlr-indexlr: warning: Skip read " << (read.num + 1) << " on line "
					   << (read.num + 1) * 4 - 2 << "; window size > #hashes (w = " << w
					   << ", #hashes = " << hashes.size() << ")\n";
					std::cerr << ss.str();
				}
			}

			if (withRepeat && withSolid) {
				for (auto hashesIt = hashes.begin(); hashesIt < hashes.end(); ++hashesIt) {
					vector<uint64_t> vect{ (*hashesIt).hash1 };
					if (repeatBF.contains(vect) || !solidBF.contains(vect)) {
						(*hashesIt).hash1 = UINT64_MAX;
					}
				}
			} else {
				if (withRepeat) {
					for (auto hashesIt = hashes.begin(); hashesIt < hashes.end(); ++hashesIt) {
						vector<uint64_t> vect{ (*hashesIt).hash1 };
						if (repeatBF.contains(vect)) {
							(*hashesIt).hash1 = UINT64_MAX;
						}
					}
				}

				if (withSolid) {
					for (auto hashesIt = hashes.begin(); hashesIt < hashes.end(); ++hashesIt) {
						vector<uint64_t> vect{ (*hashesIt).hash1 };
						if (!solidBF.contains(vect)) {
							(*hashesIt).hash1 = UINT64_MAX;
						}
					}
				}
			}

			auto minimizers = getMinimizers(hashes, w);

			ss << read.barcode;

			char sep = '\t';
			if (minimizers.empty()) {
				ss << sep;
			}
			for (auto& m : minimizers) {
				if (m.hash1 != UINT64_MAX) {
					ss << sep << m.hash2;
					if (withPositions) {
						ss << ':' << m.pos;
					}
					if (withStrands) {
						ss << ':' << m.strand;
					}
					sep = ' ';
				}
			}
			ss << '\n';
		}
		if (reads.dataCounter > 0) {
			result.num = reads.num;
			assert(reads.dataCounter - 1 < sizeof(reads.data) / sizeof(reads.data[0]));
			result.lastNum = reads.data[reads.dataCounter - 1].num;
			result.barcodesAndMinimizers = ss.str();
		} else {
			result.num = reads.num;
			result.lastNum = reads.num;
			result.barcodesAndMinimizers = "";
		}

		outputWorker.buffer.write(result);
	}
}

inline void
OutputWorker::work()
{
	Result result;
	size_t lastWritten;

	do {
		buffer.read(result);
		ofs << result.barcodesAndMinimizers;
		lastWritten = result.lastNum;
		assert_good(ofs, opath);
	} while (!inputWorker.allRead || lastWritten != inputWorker.inputNum - 1);
}

#endif
