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
    using RightCosets = typename Tgroup::RightCosets;
    using Tidx_label = typename Tgroup::Tidx_label;
    if (argc != 2 && argc != 3) {
      std::cerr << "TestRightCosets [H_G]\n";
      std::cerr << "or\n";
      std::cerr << "TestRightCosets [H_G] [OUT]\n";
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
    RightCosets rc = eG.right_cosets(eH);
    std::vector<Telt> l_cos;
    for (auto & eCos: rc) {
      l_cos.push_back(eCos);
    }
    permutalib::KernelCheckRightCosets<Telt,Tidx_label,Tint>(eG.stab_chain(), eH.stab_chain(), l_cos);
    //
    std::cerr << "CPP Normal completion of the program\n";
  } catch (permutalib::PermutalibException const &e) {
    std::cerr << "Erroneous completion of the program\n";
    exit(e.eVal);
  }
  return 0;
}
