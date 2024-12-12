// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Permutation.h"
#include "gmpxx.h"
#include <fstream>

#include "Group.h"

int main(int argc, char *argv[]) {
  try {
    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    if (argc != 3 && argc != 4) {
      std::cerr << "GapAscendingChainPair [H] [G]\n";
      std::cerr << "or\n";
      std::cerr << "GapAscendingChainPair [H] [G] [OUT]\n";
      throw permutalib::PermutalibException{1};
    }
    std::string File_H = argv[1];
    std::string File_H = argv[2];
    //
    std::ifstream is1(File_H);
    Tgroup eH = permutalib::ReadGroupFromStream<Tgroup>(is1);
    //
    std::ifstream is2(File_G);
    Tgroup eG = permutalib::ReadGroupFromStream<Tgroup>(is2);
    //
    std::vector<Tgroup> ListGroup = eG.GetAscendingChainSubgroup(eH);
    auto prt = [&](std::ostream &os) -> void {
      os << "return [";
      for (size_t i = 0; i < ListGroup.size(); i++) {
        if (i > 0)
          os << ",\n";
        os << ListGroup[i].GapString();
      }
      os << "];\n";
    };
    if (argc == 2) {
      prt(std::cerr);
    } else {
      std::ofstream os(argv[2]);
      prt(os);
    }
    //
    std::cerr << "CPP Normal completion of the program\n";
  } catch (permutalib::PermutalibException const &e) {
    std::cerr << "Erroneous completion of the program\n";
    exit(e.eVal);
  }
  return 0;
}
