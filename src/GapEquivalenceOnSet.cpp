// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Permutation.h"
#include "gmpxx.h"
#include <fstream>

#include "Group.h"

int main(int argc, char *argv[]) {
  try {
    using Tidx = uint16_t;
    using Telt = permutalib::DoubleSidedPerm<Tidx>;
    using Tint = mpz_class;
    if (argc != 3) {
      std::cerr << "GapEquivalenceOnSet [EXMP] |OutFile]\n";
      std::cerr << "with EXMP generated by GenerateExample.g\n";
      throw permutalib::PermutalibException{1};
    }
    std::ifstream is(argv[1]);
    size_t nbGen;
    int n_i;
    is >> n_i;
    is >> nbGen;
    Tidx n = Tidx(n_i);
    std::vector<Telt> LGen(nbGen);
    for (size_t iGen = 0; iGen < nbGen; iGen++) {
      std::vector<Tidx> ePermV(n);
      for (Tidx i = 0; i < n; i++) {
        int eVal_i;
        is >> eVal_i;
        Tidx eVal = Tidx(eVal_i);
        ePermV[i] = eVal;
      }
      Telt ePerm(ePermV);
      LGen[iGen] = ePerm;
    }
    std::cerr.setf(std::ios::boolalpha);
    //
    std::cerr << "CPP Before call to MinimalStabChain\n";
    //    permutalib::StabChain<Telt> eG =
    //    permutalib::MinimalStabChain<Telt,Tint>(LGen, n);
    Telt id(n);
    permutalib::Group<Telt, Tint> eG(LGen, id);
    std::cerr << "CPP After call to MinimalStabChain\n";
    //
    std::cerr << "CPP |eG|=" << eG.size() << "\n";
    //
    permutalib::Face f1(n);
    for (Tidx i = 0; i < n; i++) {
      int eVal;
      is >> eVal;
      f1[i] = eVal;
    }
    //
    permutalib::Face f2(n);
    for (Tidx i = 0; i < n; i++) {
      int eVal;
      is >> eVal;
      f2[i] = eVal;
    }
    std::optional<Telt> test = eG.RepresentativeAction_OnSets(f1, f2);
    //
    std::ofstream os(argv[2]);
    if (test) {
      os << "return " << *test << ";\n";
    } else {
      os << "return fail;\n";
    }
    std::cerr << "CPP Normal completion of the program\n";
  } catch (permutalib::PermutalibException const &e) {
    std::cerr << "Erroneous completion of the program\n";
    exit(e.eVal);
  }
  return 0;
}
