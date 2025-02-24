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
    if (argc != 2 && argc != 4) {
      std::cerr << "GapComputeDoubleCosets [H_UV]\n";
      std::cerr << "or\n";
      std::cerr << "GapComputeDoubleCosets [G_UV] [OutFormat] [FileO]\n";
      std::cerr << "\n";
      std::cerr << "OutFormat: Count (default) or GAP\n";
      std::cerr << "FileO: stderr (default), stdout or something else\n";
      throw permutalib::PermutalibException{1};
    }
    std::string File_G_UV = argv[1];
    std::string OutFormat = "Count";
    std::string FileO = "stderr";
    if (argc == 4) {
      OutFormat = argv[2];
      FileO = argv[3];
    }
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
    permutalib::MicrosecondTime_perm time1;
    std::vector<Telt> list_dcc = eG.double_cosets(eU, eV);
    std::cerr << "We have list_dcc, |list_dcc|=" << list_dcc.size() << " time=" << time1 << "\n";
    //
    auto f=[&](std::ostream& osf) -> void {
      if (OutFormat == "Count") {
        osf << "|list_dcc|=" << list_dcc.size() << "\n";
        return;
      }
      if (OutFormat == "GAP") {
        osf << "return [";
        bool IsFirst = true;
        for (auto &dcc : list_dcc) {
          if (!IsFirst) {
            osf << ",";
          }
          if (IsFirst) {
            IsFirst = false;
          }
          osf << dcc;
        }
        osf << "];\n";
        return;
      }
      std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat << "\n";
      throw permutalib::PermutalibException{1};
    };
    if (FileO == "stderr") {
      f(std::cerr);
    } else {
      if (FileO == "stdout") {
        f(std::cout);
      } else {
        std::ofstream osf(FileO);
        f(osf);
      }
    }
    std::cerr << "CPP Normal completion of the program\n";
  } catch (permutalib::PermutalibException const &e) {
    std::cerr << "Erroneous completion of the program\n";
    exit(e.eVal);
  }
  return 0;
}
