// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "Permutation.h"
#include "gmpxx.h"
#include <algorithm>
#include <fstream>
#include <string>
#include <vector>
#include "Group.h"
// clang-format on

/*
  The double cosets can be obtained in two ways:
  * By computing the full list of them level by level (the V1 functions).
  * By iterating over them one by one (the iterators).
  Both have to return exactly the same list and this program checks that on
  a (G, U, V) triple read from a file.
 */

template <typename Telt>
void CheckSameList(std::vector<Telt> const &list_full,
                   std::vector<Telt> const &list_iter,
                   std::string const &context) {
  if (list_full.size() != list_iter.size()) {
    std::cerr << "For " << context << " |list_full|=" << list_full.size()
              << " |list_iter|=" << list_iter.size() << "\n";
    std::cerr << "The full list and the iterator disagree on the number of "
                 "double cosets\n";
    throw permutalib::PermutalibException{1};
  }
  std::vector<Telt> sort_full = list_full;
  std::vector<Telt> sort_iter = list_iter;
  std::sort(sort_full.begin(), sort_full.end());
  std::sort(sort_iter.begin(), sort_iter.end());
  if (sort_full != sort_iter) {
    std::cerr << "For " << context
              << " the full list and the iterator do not return the same "
                 "double coset representatives\n";
    throw permutalib::PermutalibException{1};
  }
  if (list_full != list_iter) {
    std::cerr << "For " << context
              << " the full list and the iterator return the same double coset "
                 "representatives but in a different order\n";
    throw permutalib::PermutalibException{1};
  }
}

int main(int argc, char *argv[]) {
  try {
    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tidx_label = uint16_t;
    using Tgroup = permutalib::Group<Telt, Tint>;
    using DccEntry = typename Tgroup::DccEntry;
    /*
      The iterators are not exposed by the Group API, so the test is done on
      the InnerDoubleCosetComputer that is behind it. Only the V side is
      considered, on the U side the representatives have to be inverted.
     */
    using InnerDcc = permutalib::InnerDoubleCosetComputer<Telt, Tidx_label, Tint>;
    if (argc != 3) {
      std::cerr << "TestDoubleCosetIterator [G_UV] [meth_check]\n";
      throw permutalib::PermutalibException{1};
    }
    std::string File_G_UV = argv[1];
    std::string meth_check = argv[2];
    //
    std::ifstream is(File_G_UV);
    if (!is.good()) {
      std::cerr << "is stream is invalid, not possible to read " << File_G_UV
                << "\n";
      throw permutalib::PermutalibException{1};
    }
    Tgroup eG = permutalib::ReadGroupFromStream<Tgroup>(is);
    Tgroup eU = permutalib::ReadGroupFromStream<Tgroup>(is);
    Tgroup eV = permutalib::ReadGroupFromStream<Tgroup>(is);
    auto check_dcc=[&](std::vector<Telt> const& list_dcc) -> void {
      if (meth_check == "exhaustive") {
        return permutalib::ExhaustiveCheck_DoubleCosets(eG.stab_chain(), eU.stab_chain(), eV.stab_chain(), list_dcc);
      }
      if (meth_check == "fast_check_sizes") {
        return permutalib::FastCheckSizes_DoubleCosets<Telt,Tidx_label,Tint>(eG.stab_chain(), eU.stab_chain(), eV.stab_chain(), list_dcc);
      }
      if (meth_check == "fast_check_intersection") {
        return permutalib::FastCheckIntersection_DoubleCosets<Telt,Tidx_label,Tint>(eG.stab_chain(), eU.stab_chain(), eV.stab_chain(), list_dcc);
      }
      if (meth_check == "nothing") {
        return;
      }
      std::cerr << "Failed to find a matching method for the check\n";
      throw permutalib::PermutalibException{1};
    };
    Tint size_G = eG.size();
    Tint size_U = eU.size();
    Tint size_V = eV.size();
    std::cerr << "size_G=" << size_G << " size_U=" << size_U << " size_V=" << size_V << "\n";
    //
    InnerDcc dcc_v(eG.stab_chain(), eU.stab_chain());
    //
    permutalib::MicrosecondTime_perm time1;
    std::vector<Telt> list_full = dcc_v.double_cosets_V1(eV.stab_chain());
    std::cerr << "We have list_full, |list_full|=" << list_full.size() << " time=" << time1 << "\n";
    //
    permutalib::MicrosecondTime_perm time2;
    std::vector<Telt> list_iter;
    for (auto iter = dcc_v.begin_elt(eV.stab_chain()); iter != dcc_v.end_elt(); ++iter) {
      list_iter.push_back(*iter);
    }
    std::cerr << "We have list_iter, |list_iter|=" << list_iter.size() << " time=" << time2 << "\n";
    CheckSameList(list_full, list_iter, "double_cosets");
    std::cerr << "Pass the check of the double coset iterator\n";
    //
    permutalib::MicrosecondTime_perm time3;
    std::vector<DccEntry> listent_full = dcc_v.double_cosets_and_stabilizers_V1(eV.stab_chain());
    std::cerr << "We have listent_full, |listent_full|=" << listent_full.size() << " time=" << time3 << "\n";
    //
    permutalib::MicrosecondTime_perm time4;
    std::vector<DccEntry> listent_iter;
    for (auto iter = dcc_v.begin_elt_stab(eV.stab_chain()); iter != dcc_v.end_elt_stab(); ++iter) {
      listent_iter.push_back(*iter);
    }
    std::cerr << "We have listent_iter, |listent_iter|=" << listent_iter.size() << " time=" << time4 << "\n";
    auto get_cosets=[](std::vector<DccEntry> const& l_ent) -> std::vector<Telt> {
      std::vector<Telt> l_cos;
      for (auto & ent : l_ent) {
        l_cos.push_back(ent.cos);
      }
      return l_cos;
    };
    CheckSameList(get_cosets(listent_full), get_cosets(listent_iter), "double_cosets_and_stabilizers");
    /*
      The generating sets themselves cannot be compared since
      Kernel_SmallGeneratingSet is randomized. What has to match is the group
      that they generate.
     */
    Telt id = eG.get_identity();
    for (size_t i_ent=0; i_ent<listent_full.size(); i_ent++) {
      Tgroup stab_full(listent_full[i_ent].stab_gens, id);
      Tgroup stab_iter(listent_iter[i_ent].stab_gens, id);
      if (stab_full != stab_iter) {
        std::cerr << "For i_ent=" << i_ent << " the stabilizer of the full list and of the iterator differ\n";
        std::cerr << "|stab_full|=" << stab_full.size() << " |stab_iter|=" << stab_iter.size() << "\n";
        throw permutalib::PermutalibException{1};
      }
    }
    std::cerr << "Pass the check of the double coset and stabilizer iterator\n";
    //
    check_dcc(list_iter);
    std::cerr << "Pass the check meth_check=" << meth_check << "\n";
    //
    std::cerr << "CPP Normal completion of the program\n";
  } catch (permutalib::PermutalibException const &e) {
    std::cerr << "Erroneous completion of the program\n";
    exit(e.eVal);
  }
  return 0;
}
