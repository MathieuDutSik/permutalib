// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Group.h"
#include "Permutation.h"
#include "gmpxx.h"
#include <fstream>

int main(int argc, char *argv[]) {
  try {
    std::vector<int64_t> ListIdx;
    for (int i=0; i<100; i++) {
      int val = (random() % 5) - 2;
      if (val != 0)
        ListIdx.push_back(val);
    }
    permutalib::SequenceType<true> seq1(ListIdx);
    std::cerr << "seq1 = " << seq1 << "\n";
    permutalib::SimplifySequence(ListIdx);
    permutalib::SequenceType<true> seq2(ListIdx);
    std::cerr << "seq2 = " << seq2 << "\n";
    std::cerr << "CPP Normal completion of the program\n";
  } catch (permutalib::PermutalibException const &e) {
    std::cerr << "Erroneous completion of the program\n";
    exit(e.eVal);
  }
}
