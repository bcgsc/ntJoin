#include "../include/btllib/seq_reader.hpp"

#include <iostream>

int
main()
{
  btllib::SeqReader reader("my_reads.fq.gz");
  for (btllib::SeqReader::Record record; (record = reader.read());) {
    std::cout << record.seq << '\n' << record.qual << '\n';
  }

  return 0;
}