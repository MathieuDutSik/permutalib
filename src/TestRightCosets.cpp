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
    using Tcosets = permutalib::RightCosets<Telt, Tint>;
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
    Tcosets rc = eG.right_cosets(eH);
    std::vector<Telt> H_elts = eH.get_all_element();
    std::unordered_set<Telt> set;
    std::vector<std::unordered_set<Telt>> l_cos;
    for (auto &eCos : rc) {
      std::cerr << "eCos=" << eCos << " eCos.siz=" << static_cast<size_t>(eCos.siz) << "\n";
      std::unordered_set<Telt> set;
      for (auto x_h: H_elts) {
        Telt p = x_h * eCos;
        set.insert(p);
      }
      l_cos.push_back(set);
    }
    auto size_int=[&](size_t const& idx1, size_t const& idx2) -> size_t {
      size_t the_int = 0;
      for (auto & elt: l_cos[idx1]) {
        if (l_cos[idx2].count(elt) == 1) {
          the_int += 1;
        }
      }
      return the_int;
    };
    std::cerr << "|l_cos|=" << l_cos.size() << "\n";
    std::cerr << "|G|=" << eG.size() << " |H|=" << eH.size() << "\n";
    for (size_t i_cos=0; i_cos<l_cos.size(); i_cos++) {
      for (size_t j_cos=i_cos + 1; j_cos<l_cos.size(); j_cos++) {
        size_t the_int = size_int(i_cos, j_cos);
        if (the_int > 0) {
          std::cerr << "Cosets i_cos=" << i_cos << " and j_cos=" << j_cos << " are intersecting, the_int=" << the_int << "\n";
          throw permutalib::PermutalibException{1};
        }
      }
    }
    //
    std::cerr << "CPP Normal completion of the program\n";
  } catch (permutalib::PermutalibException const &e) {
    std::cerr << "Erroneous completion of the program\n";
    exit(e.eVal);
  }
  return 0;
}
