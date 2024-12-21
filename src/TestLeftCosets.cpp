// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "Group.h"
#include "Permutation.h"
#include "gmpxx.h"
#include <fstream>
// clang-format on

int main(int argc, char *argv[]) {
  try {
    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    using Tidx_label = typename Tgroup::Tidx_label;
    if (argc != 2 && argc != 3) {
      std::cerr << "TestLeftCosets [H_G]\n";
      std::cerr << "or\n";
      std::cerr << "TestLeftCosets [H_G] [OUT]\n";
      throw permutalib::PermutalibException{1};
    }
    std::string File_HG = argv[1];
    //
    std::ifstream is(File_HG);
    if (!is.good()) {
      std::cerr << "is stream is invalid, not possible to read H_G\n";
      throw permutalib::PermutalibException{1};
    }
    Tgroup eH = permutalib::ReadGroupFromStream<Tgroup>(is);
    Tgroup eG = permutalib::ReadGroupFromStream<Tgroup>(is);
    std::vector<Telt> l_cos = eG.get_all_left_cosets(eH);
    permutalib::KernelCheckLeftCosets<Telt,Tidx_label,Tint>(eG.stab_chain(), eH.stab_chain(), l_cos);
    //
    std::cerr << "CPP Normal completion of the program\n";
  } catch (permutalib::PermutalibException const &e) {
    std::cerr << "Erroneous completion of the program\n";
    exit(e.eVal);
  }
  return 0;
}
