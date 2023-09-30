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
      std::cerr << "GapComputeSequenceBlockSystems [EXMP]\n";
      std::cerr << "or\n";
      std::cerr << "GapComputeSequenceBlockSystems [EXMP] [OUT]\n";
      throw permutalib::PermutalibException{1};
    }
    std::string InputFile = argv[1];
    std::ifstream is(InputFile);
    Tgroup eG = permutalib::ReadGroupFromStream<Tgroup>(is);
    //
    std::vector<permutalib::BlockDecomposition<Tidx>> ListBlkDec =
        eG.GetSequenceBlockDecomposition();
    auto prt = [&](std::ostream &os) -> void {
      os << "return [";
      for (size_t i = 0; i < ListBlkDec.size(); i++) {
        if (i > 0)
          os << ",\n";
        os << "rec(map_vert_block:=[";
        bool IsFirst = true;
        for (auto &iBlock : ListBlkDec[i].map_vert_block) {
          if (!IsFirst)
            os << ",";
          IsFirst = false;
          os << iBlock;
        }
        os << "])";
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
