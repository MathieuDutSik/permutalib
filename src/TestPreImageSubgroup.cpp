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
    if (argc != 3) {
      std::cerr << "TestPreImageSubgroup [GRP_big] [GRP_sma]\n";
      throw permutalib::PermutalibException{1};
    }
    std::string File_GRP_big = argv[1];
    std::string File_GRP_sma = argv[2];
    //
    Tgroup GRP_big = permutalib::ReadGroupFromFile<Tgroup>(File_GRP_big);
    Tgroup GRP_sma = permutalib::ReadGroupFromFile<Tgroup>(File_GRP_sma);
    //
    Tint size_GRP_big = GRP_big.size();
    Tint size_GRP_sma = GRP_sma.size();
    bool test = GRP_big.IsSubgroup(GRP_sma);
    std::cerr << "|GRP_big|=" << size_GRP_big << " |GRP_sma|=" << size_GRP_sma << " test=" << test << "\n";

    std::vector<Telt> ListPermGens = GRP_big.SmallGeneratingSet();
    std::cerr << "|ListPermGens|=" << ListPermGens.size() << "\n";

    size_t n_act_GRP_big = GRP_big.n_act();
    size_t n_act_GRP_sma = GRP_sma.n_act();
    std::cerr << "n_act_GRP_big=" << n_act_GRP_big << " n_act_GRP_sma=" << n_act_GRP_sma << "\n";

    Telt id = Telt(n_act_GRP_big);
    std::vector<Telt> LGen_pre = permutalib::PreImageSubgroup<Tgroup,Telt>(ListPermGens, ListPermGens, id, GRP_sma);
    Tgroup PreGroup = Tgroup(LGen_pre, id);
    std::cerr << "|PreGroup|=" << PreGroup.size() << "\n";
    if (PreGroup == GRP_sma) {
      std::cerr << "The groups are equal\n";
    } else {
      std::cerr << "The groups are NOT equal\n";
    }
    //
    std::cerr << "CPP Normal completion of the program\n";
  } catch (permutalib::PermutalibException const &e) {
    std::cerr << "Erroneous completion of the program\n";
    exit(e.eVal);
  }
  return 0;
}
