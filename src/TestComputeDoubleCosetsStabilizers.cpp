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
    using DccEntry = typename Tgroup::DccEntry;
    if (argc != 2 && argc != 3) {
      std::cerr << "TestComputeDoubleCosetsStabilizer [G_UV]\n";
      std::cerr << "or\n";
      std::cerr << "TestComputeDoubleCosetsStabilizer [G_UV] [OUT]\n";
      throw permutalib::PermutalibException{1};
    }
    std::string File_G_UV = argv[1];
    //
    std::ifstream is(File_G_UV);
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
    permutalib::MicrosecondTime_perm time;
    std::vector<DccEntry> list_dccent = eG.double_cosets_and_stabilizers(eU, eV);
    std::cerr << "We have list_dccent, |list_dccent|=" << list_dccent.size() << " time=" << time << "\n";
    std::vector<Telt> list_dcc;
    for (auto& dccent : list_dccent) {
      list_dcc.push_back(dccent.cos);
      std::unordered_set<Telt> set_cos;
      for (auto & elt : eU) {
        Telt p = elt * dccent.cos;
        set_cos.insert(p);
      }
      auto is_stab=[&](Telt const& y) -> bool {
        for (auto & x : set_cos) {
          Telt x_y = x * y;
          if (set_cos.count(x_y) == 0) {
            return false;
          }
        }
        return true;
      };
      for (auto & gen : dccent.stab_gens) {
        if (!is_stab(gen)) {
          std::cerr << "The element does not stabilize set_cos\n";
          throw permutalib::PermutalibException{1};
        }
      }
    }
    std::cerr << "Pass the test of coset stabilization time=" << time << "\n";
    ExhaustiveCheck_DoubleCosets(eG.stab_chain(), eU.stab_chain(), eV.stab_chain(), list_dcc);
    std::cerr << "Pass KernelCheckDoubleCosets time=" << time << "\n";
    //
    std::cerr << "CPP Normal completion of the program\n";
  } catch (permutalib::PermutalibException const &e) {
    std::cerr << "Erroneous completion of the program\n";
    exit(e.eVal);
  }
  return 0;
}
