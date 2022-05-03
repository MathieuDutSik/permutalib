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
      std::cerr << "GapAscendingChain [EXMP]\n";
      std::cerr << "or\n";
      std::cerr << "GapAscendingChain [EXMP] [OUT]\n";
      throw permutalib::PermutalibException{1};
    }
    std::string InputFile = argv[1];
    //
    std::ifstream is(InputFile);
    size_t nbGen;
    int n_i;
    is >> nbGen;
    is >> n_i;
    Tidx n = Tidx(n_i);
    std::vector<Telt> LGen(nbGen);
    for (size_t iGen = 0; iGen < nbGen; iGen++) {
      std::vector<Tidx> ePermV(n);
      for (Tidx i = 0; i < n; i++) {
        int eVal_i;
        is >> eVal_i;
        ePermV[i] = Tidx(eVal_i);
      }
      Telt ePerm(ePermV);
      LGen[iGen] = ePerm;
    }
    //
    Telt id(n);
    Tgroup eG(LGen, id);
    //
    std::vector<Tgroup> ListGroup = eG.GetAscendingChain();
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
