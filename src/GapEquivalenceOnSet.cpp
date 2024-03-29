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
    using Tgroup = permutalib::Group<Telt, Tint>;
    if (argc != 3) {
      std::cerr << "GapEquivalenceOnSet [EXMP] |OutFile]\n";
      std::cerr << "with EXMP generated by GenerateExample.g\n";
      throw permutalib::PermutalibException{1};
    }
    std::ifstream is(argv[1]);
    Tgroup eG = permutalib::ReadGroupFromStream<Tgroup>(is);
    Tidx n = eG.n_act();
    std::cerr.setf(std::ios::boolalpha);
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
