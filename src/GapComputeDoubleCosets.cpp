// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "Permutation.h"
#include "gmpxx.h"
#include <fstream>
#include "Group.h"
// clang-format on

int main(int argc, char *argv[]) {
  try {
    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    using DoubleCosetComputer = typename Tgroup::DoubleCosetComputer;
    if (argc != 2 && argc != 3) {
      std::cerr << "GapAscendingChainPair [H_G]\n";
      std::cerr << "or\n";
      std::cerr << "GapAscendingChainPair [H_G] [OUT]\n";
      throw permutalib::PermutalibException{1};
    }
    std::string File_HG = argv[1];
    //
    std::ifstream is(File_HG);
    if (!is.good()) {
      std::cerr << "is stream is invalid, not possible to read H_G\n";
      throw permutalib::PermutalibException{1};
    }
    Tgroup eG = permutalib::ReadGroupFromStream<Tgroup>(is);
    Tgroup eU = permutalib::ReadGroupFromStream<Tgroup>(is);
    Tgroup eV = permutalib::ReadGroupFromStream<Tgroup>(is);
    Tint size_G = eG.size();
    Tint size_U = eU.size();
    Tint size_V = eV.size();
    std::cerr << "size_G=" << size_G << " size_U=" << size_U << " size_V=" << size_V << "\n";
    //
    permutalib::MicrosecondTime_perm time1;
    DoubleCosetComputer dcc_v = eG.double_coset_computer_v(eU);
    std::vector<Telt> list_dcc1 = dcc_v.double_cosets(eV);
    std::cerr << "We have list_dcc1, |list_dcc1|=" << list_dcc1.size() << " time=" << time1 << "\n";
    KernelCheckDoubleCosets(eG.stab_chain(), eU.stab_chain(), eV.stab_chain(), list_dcc1);
    //
    permutalib::MicrosecondTime_perm time2;
    DoubleCosetComputer dcc_u = eG.double_coset_computer_u(eV);
    std::vector<Telt> list_dcc2 = dcc_u.double_cosets(eU);
    std::cerr << "We have list_dcc2, |list_dcc2|=" << list_dcc2.size() << " time=" << time2 << "\n";
    KernelCheckDoubleCosets(eG.stab_chain(), eU.stab_chain(), eV.stab_chain(), list_dcc2);
    //
    permutalib::MicrosecondTime_perm time3;
    std::vector<Telt> list_dcc3 = eG.double_cosets(eU, eV);
    std::cerr << "We have list_dcc3, |list_dcc3|=" << list_dcc3.size() << " time=" << time3 << "\n";
    KernelCheckDoubleCosets(eG.stab_chain(), eU.stab_chain(), eV.stab_chain(), list_dcc3);
    //
    std::cerr << "CPP Normal completion of the program\n";
  } catch (permutalib::PermutalibException const &e) {
    std::cerr << "Erroneous completion of the program\n";
    exit(e.eVal);
  }
  return 0;
}
