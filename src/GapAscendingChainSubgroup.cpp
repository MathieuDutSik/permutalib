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
    if (argc != 2 && argc != 3) {
      std::cerr << "GapAscendingChainPair [H_G]\n";
      std::cerr << "or\n";
      std::cerr << "GapAscendingChainPair [H_G] [OUT]\n";
      throw permutalib::PermutalibException{1};
    }
    std::string File_HG = argv[1];
    //
    std::ifstream is(File_HG);
    Tgroup eH = permutalib::ReadGroupFromStream<Tgroup>(is);
    Tgroup eG = permutalib::ReadGroupFromStream<Tgroup>(is);
    Tint size_G = eG.size();
    Tint size_H = eH.size();
    std::cerr << "size_G=" << size_G << " size_H=" << size_H << "\n";
    //
    std::vector<Tgroup> ListGroup = eG.GetAscendingChainSubgroup(eH);
    for (size_t i = 0; i < ListGroup.size(); i++) {
      std::cerr << "i=" << i << " |GRP|=" << ListGroup[i].size() << "\n";
    }
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
