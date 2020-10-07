#include "btllib/indexlr.hpp"
#include "btllib/bloom_filter.hpp"

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
#include <sstream>
#include <memory>
#include <mutex>
#include <queue>
#include <condition_variable>
#include <cstdio>

const static char* PROGNAME = "indexlr";
const static size_t OUTPUT_PERIOD = 512;

static void
print_error_msg(const std::string& msg)
{
  std::cerr << PROGNAME << ": " << msg << std::endl;
}

static void
print_usage()
{
  std::cerr << "Usage: " << PROGNAME
            << "  -k K -w W [-r repeat_bf_path] [-s solid_bf_path] [--id] [--bx] [--pos] [--seq] [-o FILE] FILE...\n\n"
               "  -k K        use K as k-mer size\n"
               "  -w W        use W as sliding-window size\n"
               "  --id        include read ids in the output\n"
               "  --bx        include read barcodes in the output\n"
               "  --pos       include minimizer positions in the output\n"
               "  --seq       include minimizer sequences in the output\n"
               "  -r repeat_bf_path  use a Bloom filter to filter out repetitive minimizers\n"
               "  -s solid_bf_path  use a Bloom filter to only select solid minimizers\n"
               "  -o FILE     write output to FILE, default is stdout\n"
               "  -t T        use T number of threads (default 5, max 5) per input file\n"
               "  --help      display this help and exit\n"
               "  FILE        space separated list of FASTA/Q files" << std::endl;
}

int
main(int argc, char* argv[])
{
  int c;
  int optindex = 0;
  int help = 0;
  unsigned k = 0;
  unsigned w = 0;
  bool w_set = false;
  bool k_set = false;	
  unsigned t = 5;
  int with_id = 0;
  int with_bx = 0;
  int with_pos = 0;
  int with_seq = 0;
  std::unique_ptr<btllib::BloomFilter> repeat_bf;
  std::unique_ptr<btllib::BloomFilter> solid_bf;
  bool with_repeat = false;
  bool with_solid = false;
  std::string outfile("-");
  bool failed = false;
  static const struct option longopts[] = {
                                          { "id", no_argument, &with_id, 1 },
                                          { "bx", no_argument, &with_bx, 1 },
                                          { "pos", no_argument, &with_pos, 1 },
                                          { "seq", no_argument, &with_seq, 1 },
                                          { "help", no_argument, &help, 1 },
                                          { nullptr, 0, nullptr, 0 } };
  while ((c = getopt_long(argc, argv, "k:w:o:t:r:s:", longopts, &optindex)) != -1) {
    switch (c) {
    case 0:
      break;
    case 'k':
      k_set = true;
      k = std::stoul(optarg);
      break;
    case 'w':
      w_set = true;
      w = std::stoul(optarg);
      break;
    case 'o':
      outfile = optarg;
      break;
    case 't':
      t = std::stoul(optarg);
      break;
    case 'r': {
      with_repeat = true;
      std::cerr << "Loading repeat Bloom filter from " << optarg << std::endl;
      try {
        repeat_bf = std::unique_ptr<btllib::BloomFilter>(new btllib::BloomFilter(optarg));
      } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
      }
      std::cerr << "Finished loading repeat Bloom filter" << std::endl;
      break;
    }
    case 's': {
      with_solid = true;
      std::cerr << "Loading solid Bloom filter from " << optarg << std::endl;
      try {
        solid_bf = std::unique_ptr<btllib::BloomFilter>(new btllib::BloomFilter(optarg));
      } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
      }
      std::cerr << "Finished loading solid Bloom filter" << std::endl;
      break;
    }
    default:
      std::exit(EXIT_FAILURE);
    }
  }
  if (t > 5) {
    t = 5;
    std::cerr
        << PROGNAME
        << ": Using more than 5 threads does not scale, reverting to 5."
        << std::endl;
  }
  std::vector<std::string> infiles(&argv[optind], &argv[argc]);
  if (argc < 2) {
    print_usage();
    std::exit(EXIT_FAILURE);
  }
  if (help != 0) {
    print_usage();
    std::exit(EXIT_SUCCESS);
  }
  if (!k_set) {
    print_error_msg("missing option -- 'k'");
    failed = true;
  } else if (k == 0) {
    print_error_msg("option has incorrect value -- 'k'");
    failed = true;
  }
  if (!w_set) {
    print_error_msg("missing option -- 'w'");
    failed = true;
  } else if (w == 0) {
    print_error_msg("option has incorrect value -- 'w'");
    failed = true;
  }
  if (infiles.empty()) {
    print_error_msg("missing file operand");
    failed = true;
  }
  if (failed) {
    std::cerr << "Try '" << PROGNAME << " --help' for more information.\n";
    std::exit(EXIT_FAILURE);
  }

  int flags = 0;
  if (with_id)  { flags |= btllib::Indexlr::Flag::ID;  }
  if (with_bx)  { flags |= btllib::Indexlr::Flag::BX;  }
  if (with_seq) { flags |= btllib::Indexlr::Flag::SEQ; }

  btllib::Indexlr::Record record;
  FILE* out;
  if (outfile == "-") {
    out = stdout;
  } else {
    out = fopen(outfile.c_str(), "w");
  }
  for (auto& infile : infiles) {
    std::unique_ptr<btllib::Indexlr> indexlr;
    if (with_repeat && with_solid) {
      flags |= btllib::Indexlr::Flag::FILTER_IN;
      flags |= btllib::Indexlr::Flag::FILTER_OUT;
      indexlr = std::unique_ptr<btllib::Indexlr>(new btllib::Indexlr(infile, k, w, flags, t, *solid_bf, *repeat_bf));
    } else if (with_repeat) {
      flags |= btllib::Indexlr::Flag::FILTER_OUT;
      indexlr = std::unique_ptr<btllib::Indexlr>(new btllib::Indexlr(infile, k, w, flags, t, *repeat_bf));
    } else if (with_solid) {
      flags |= btllib::Indexlr::Flag::FILTER_IN;
      indexlr = std::unique_ptr<btllib::Indexlr>(new btllib::Indexlr(infile, k, w, flags, t, *solid_bf));
    } else {
      indexlr = std::unique_ptr<btllib::Indexlr>(new btllib::Indexlr(infile, k, w, flags, t));
    }
    std::queue<std::string> output_queue;
    std::mutex output_queue_mutex;
    std::condition_variable queue_empty, queue_full;
    size_t max_seen_output_size = 100;
    size_t queue_max_size = 128;
    std::unique_ptr<std::thread> info_compiler(new std::thread([&]() {
      std::stringstream ss;
      while ((record = indexlr->get_minimizers())) {
        if (record.num % OUTPUT_PERIOD == OUTPUT_PERIOD - 1) {
          max_seen_output_size = std::max(max_seen_output_size, ss.str().size());
          std::unique_lock<std::mutex> lock(output_queue_mutex);
          while (output_queue.size() == queue_max_size) { queue_full.wait(lock); }
          output_queue.push(std::move(ss.str()));
          queue_empty.notify_one();
          lock.unlock();
          std::string newstring;
          newstring.reserve(max_seen_output_size);
          ss.str(std::move(newstring));
        } else {
          if (with_id || (!with_id && !with_bx)) { ss << record.id << '\t'; }
          if (with_bx) { ss << record.barcode << '\t'; }
          int j = 0;
          for (const auto& min : record.minimizers) {
            if (j > 0) { ss << ' '; } 
            ss << min.hash2;
            if (with_pos) { ss << ':' << min.pos; }
            if (with_seq) { ss << ':' << min.seq; }
            j++;
          }
          ss << '\n';
        }
      }
      {
        std::unique_lock<std::mutex> lock(output_queue_mutex);
        output_queue.push(std::move(ss.str()));
        output_queue.push(std::string());
        queue_empty.notify_one();
      }
    }));
    std::unique_ptr<std::thread> output_worker(new std::thread([&]() {
      std::string to_write;
      for (;;) {
        {
          std::unique_lock<std::mutex> lock(output_queue_mutex);
          while (output_queue.empty()) { queue_empty.wait(lock); } 
          to_write = std::move(output_queue.front());
          output_queue.pop();
          queue_full.notify_one();
        }
        if (to_write.empty()) { break; }
        fwrite(to_write.c_str(), 1, to_write.size(), out);
      }
    }));
    info_compiler->join();
    output_worker->join();
  }
  if (out != stdout) {
    fclose(out);
  }

  return 0;
}