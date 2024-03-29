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
      std::cerr << "The program is used as\n";
      std::cerr << "GapSmallGeneratingSet [EXMP]\n";
      std::cerr << "or\n";
      std::cerr << "GapSmallGeneratingSet [EXMP] [OutFile]\n";
      std::cerr << "with EXMP generated by GenerateExample.g\n";
      throw permutalib::PermutalibException{1};
    }
    std::string InputFile = argv[1];
    //
    std::ifstream is(InputFile);
    Tgroup eG = permutalib::ReadGroupFromStream<Tgroup>(is);
    std::cerr.setf(std::ios::boolalpha);
    std::vector<Telt> Sma = eG.SmallGeneratingSet();
    //
    if (argc == 3) {
      std::string OutputFile = argv[2];
      std::ofstream os(OutputFile);
      os << "return " + GapStringTVector(Sma) + ";\n";
    } else {
      std::cerr << "CPP Sma=" << GapStringTVector(Sma) << "\n";
    }
    std::cerr << "CPP Normal completion of the program\n";
  } catch (permutalib::PermutalibException const &e) {
    std::cerr << "Erroneous completion of the program\n";
    exit(e.eVal);
  }
  return 0;
}
