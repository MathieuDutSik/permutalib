// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Permutation.h"
#include "gmpxx.h"
#include <fstream>

#include "Group.h"
#include "Resolution.h"

int main(int argc, char *argv[]) {
  try {
    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tres = permutalib::FiniteGroupResolution<Telt, Tint>;
    if (argc != 3 && argc != 4) {
      std::cerr << "GapResolution [EXMP] [K]\n";
      std::cerr << "or\n";
      std::cerr << "GapResolution [EXMP] [K] [OutFile]\n";
      std::cerr << "with EXMP a group file and K the length of the resolution.\n";
      std::cerr << "The output is a GAP record with the elements, the ranks and\n";
      std::cerr << "the boundaries in the word format of HAP, so that it can be\n";
      std::cerr << "compared with ResolutionFiniteGroup(GeneratorsOfGroup(G), K).\n";
      throw permutalib::PermutalibException{1};
    }
    std::string InputFile = argv[1];
    std::ifstream is(InputFile);
    std::pair<std::vector<Telt>, Telt> pair = permutalib::ReadListGenFromStream<Telt>(is);
    int K_i;
    (void)sscanf(argv[2], "%d", &K_i);
    size_t K = size_t(K_i);
    //
    Tres R(pair.first, pair.second, K);
    //
    if (argc == 4) {
      std::string OutputFile = argv[3];
      std::ofstream os(OutputFile);
      os << "return " << R.GapString() << ";\n";
    } else {
      std::cerr << "CPP |G|=" << R.group_order() << " ranks=[";
      for (size_t i = 0; i <= K; i++) {
        if (i > 0)
          std::cerr << ",";
        std::cerr << R.dimension(int(i));
      }
      std::cerr << "]\n";
      for (size_t n = 1; n < K; n++) {
        std::cerr << "CPP H_" << n << " = [";
        std::vector<Tint> hom = R.integral_homology(n);
        for (size_t u = 0; u < hom.size(); u++) {
          if (u > 0)
            std::cerr << ",";
          std::cerr << hom[u];
        }
        std::cerr << "]\n";
      }
    }
    std::cerr << "CPP Normal completion of the program\n";
  } catch (permutalib::PermutalibException const &e) {
    std::cerr << "Erroneous completion of the program\n";
    exit(e.eVal);
  }
  return 0;
}
