// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Group.h"
#include "Permutation.h"
#include "gmpxx.h"
#include <fstream>
#include "TestingFct.h"


template<typename Tgroup>
void full_check(Tgroup const& eG, std::string const& opt, int64_t const& n_iter, size_t & siz_control) {
  using Tint = typename Tgroup::Tint;
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  Tidx n = eG.n_act();
  auto random_face = [](const Tidx &len) -> permutalib::Face {
    permutalib::Face eFace(len);
    for (Tidx i = 0; i < len; i++) {
      int eVal = Tidx(random()) % 2;
      eFace[i] = eVal;
    }
    return eFace;
  };
  auto bench_canonical = [&]() -> void {
    for (int64_t iter = 0; iter < n_iter; iter++) {
      permutalib::Face eFace = random_face(n);
      permutalib::Face set_can = eG.CanonicalImage(eFace);
      siz_control += set_can.count();
    }
  };
  auto bench_stabilizer = [&]() -> void {
    for (int64_t iter = 0; iter < n_iter; iter++) {
      permutalib::Face eFace = random_face(n);
      Tgroup eG2 = eG.Stabilizer_OnSets(eFace);
      siz_control += eFace.count();
    }
  };
  auto bench_pointstabilizer = [&]() -> void {
    for (int64_t iter = 0; iter < n_iter; iter++) {
      Tidx pos = Tidx(random()) % n;
      Tgroup eG2 = eG.Stabilizer_OnPoints(pos);
      siz_control += size_t(pos);
    }
  };
  auto bench_pointrepresentative = [&]() -> void {
    for (int64_t iter = 0; iter < n_iter; iter++) {
      Tidx pos1 = Tidx(random()) % n;
      Tidx pos2 = Tidx(random()) % n;
      std::optional<Telt> eP = eG.RepresentativeAction_OnPoints(pos1, pos2);
      siz_control += size_t(pos1) + size_t(pos2);
      if (eP)
        siz_control++;
    }
  };
  auto check_canonical = [&]() -> void {
    for (int64_t iter = 0; iter < n_iter; iter++) {
      permutalib::Face eFace1 = random_face(n);
      permutalib::Face set_can1 = eG.CanonicalImage(eFace1);
      for (int i = 0; i < 4; i++) {
        Telt u = eG.random();
        permutalib::Face eFace2 = OnSets(eFace1, u);
        permutalib::Face set_can2 = eG.CanonicalImage(eFace2);
        if (set_can1 != set_can2) {
          std::cerr << "Canonicalization failed\n";
          throw permutalib::PermutalibException{1};
        }
        std::pair<permutalib::Face,Tgroup> pairCan = eG.PairCanonicalImageSubgroupStabilizer(eFace2);
        if (pairCan.first != set_can2) {
          std::cerr << "We fail at that obvious step\n";
          throw permutalib::PermutalibException{1};
        }
        Tgroup eStab = eG.Stabilizer_OnSets(pairCan.first);
        if (!eStab.IsSubgroup(pairCan.second)) {
          std::cerr << "We fail that pairCan.second is not a subgroup of the stabilizer\n";
          std::cerr << "|eStab|=" << eStab.size() << " |pairCan.second|=" << pairCan.second.size() << "\n";
          throw permutalib::PermutalibException{1};
        }
      }
      siz_control += eFace1.count();
    }
  };
  auto check_exhaustive_canonical = [&]() -> void {
    for (int64_t iter = 0; iter < n_iter; iter++) {
      permutalib::Face eFace1 = random_face(n);
      permutalib::Face set_can1 = eG.ExhaustiveCanonicalImage(eFace1);
      for (int i = 0; i < 4; i++) {
        Telt u = eG.random();
        permutalib::Face eFace2 = OnSets(eFace1, u);
        permutalib::Face set_can2 = eG.ExhaustiveCanonicalImage(eFace2);
        if (set_can1 != set_can2) {
          std::cerr << "ExhaustiveCanonicalization failed\n";
          std::cerr << "set_can1=" << set_can1 << "\n";
          std::cerr << "set_can2=" << set_can2 << "\n";
          throw permutalib::PermutalibException{1};
        }
      }
      siz_control += eFace1.count();
    }
  };
  auto check_store_canonical = [&]() -> void {
    for (int64_t iter = 0; iter < n_iter; iter++) {
      permutalib::Face eFace1 = random_face(n);
      permutalib::Face set_can1 = eG.StoreCanonicalImage(eFace1);
      for (int i = 0; i < 4; i++) {
        Telt u = eG.random();
        permutalib::Face eFace2 = OnSets(eFace1, u);
        permutalib::Face set_can2 = eG.StoreCanonicalImage(eFace2);
        if (set_can1 != set_can2) {
          std::cerr << "ExhaustiveCanonicalization failed\n";
          std::cerr << "set_can1=" << set_can1 << "\n";
          std::cerr << "set_can2=" << set_can2 << "\n";
          throw permutalib::PermutalibException{1};
        }
      }
      siz_control += eFace1.count();
    }
  };
  auto check_canonical_orbitsize = [&]() -> void {
    for (int64_t iter = 0; iter < n_iter; iter++) {
      permutalib::Face eFace1 = random_face(n);
      permutalib::Face set_can1 = eG.OptCanonicalImage(eFace1);
      for (int i = 0; i < 4; i++) {
        Telt u = eG.random();
        permutalib::Face eFace2 = OnSets(eFace1, u);
        std::pair<permutalib::Face,Tint> pair = eG.OptCanonicalImageOrbitSize(eFace2);
        if (set_can1 != pair.first) {
          std::cerr << "ExhaustiveCanonicalization failed\n";
          std::cerr << "set_can1   = " << set_can1 << "\n";
          std::cerr << "pair.first = " << pair.first << "\n";
          throw permutalib::PermutalibException{1};
        }
        Tgroup eStab = eG.Stabilizer_OnSets(eFace2);
        Tint eProd = eStab.size() * pair.second;
        if (eProd != eG.size()) {
          std::cerr << "Discrepancy in the order size\n";
          std::cerr << "eProd=" << eProd << "\n";
          std::cerr << " |eG|=" << eG.size() << "\n";
          throw permutalib::PermutalibException{1};
        }
      }
      siz_control += eFace1.count();
    }
  };
  auto timing_canonical_algorithms = [&]() -> void {
    std::vector<permutalib::Face> ListF;
    for (int64_t iter = 0; iter < n_iter; iter++) {
      permutalib::Face eFace1 = random_face(n);
      ListF.push_back(eFace1);
    }
    permutalib::NanosecondTime time;
    for (auto & f : ListF) {
      (void)eG.CanonicalImage(f);
    }
    double time_canonic = double(time.eval());
    for (auto & f : ListF) {
      (void)eG.ExhaustiveCanonicalImage(f);
    }
    double time_exhaustive_canonic = double(time.eval());
    for (auto & f : ListF) {
      (void)eG.StoreCanonicalImage(f);
    }
    double time_store_canonic = double(time.eval());
    double frac1 = time_canonic / time_exhaustive_canonic;
    double frac2 = time_exhaustive_canonic / time_store_canonic;
    double frac3 = time_canonic / time_store_canonic;
    std::cerr << "|eG|=" << eG.size() << " frac3=" << frac3 << " frac1=" << frac1 << " frac2=" << frac2 << " |canonic|=" << time_canonic << " |exhaust|=" << time_exhaustive_canonic << " |store|=" << time_store_canonic << "\n";
  };
  auto check_representative = [&]() -> void {
    for (int64_t iter = 0; iter < n_iter; iter++) {
      permutalib::Face eFace1 = random_face(n);
      Telt u = eG.random();
      permutalib::Face eFace2 = OnSets(eFace1, u);
      std::optional<Telt> test =
        eG.RepresentativeAction_OnSets(eFace1, eFace2);
      if (!test) {
        std::cerr << "RepresentativeAction_OnSets error\n";
        throw permutalib::PermutalibException{1};
      }
      siz_control += eFace1.count();
    }
  };
  auto check_stabilizer = [&]() -> void {
    for (int64_t iter = 0; iter < n_iter; iter++) {
      permutalib::Face eFace1 = random_face(n);
      Tgroup eG1 = eG.Stabilizer_OnSets(eFace1);
      for (int i = 0; i < 4; i++) {
        Telt u = eG.random();
        permutalib::Face eFace2 = OnSets(eFace1, u);
        Tgroup eG2 = eG.Stabilizer_OnSets(eFace2);
        if (eG1.size() != eG2.size()) {
          std::cerr << "Stabilizer_OnSets error\n";
          throw permutalib::PermutalibException{1};
        }
      }
      siz_control += eFace1.count();
    }
  };
  //
  if (opt == "canonical")
    bench_canonical();
  if (opt == "stabilizer")
    bench_stabilizer();
  if (opt == "pointstabilizer")
    bench_pointstabilizer();
  if (opt == "pointrepresentative")
    bench_pointrepresentative();
  if (opt == "check_canonical")
    check_canonical();
  if (opt == "check_exhaustive_canonical")
    check_exhaustive_canonical();
  if (opt == "check_canonical_orbitsize")
    check_canonical_orbitsize();
  if (opt == "check_store_canonical")
    check_store_canonical();
  if (opt == "timing_canonical_algorithms")
    timing_canonical_algorithms();
  if (opt == "check_representative")
    check_representative();
  if (opt == "check_stabilizer")
    check_stabilizer();
  //
  if (opt == "all") {
    bench_canonical();
    bench_stabilizer();
    bench_pointstabilizer();
    bench_pointrepresentative();
    check_canonical();
    check_exhaustive_canonical();
    check_canonical_orbitsize();
    check_store_canonical();
    check_representative();
    check_stabilizer();
  }
}


int main(int argc, char *argv[]) {
  try {
    //    using Tidx = int16_t;
    using Tidx = uint8_t;
    // using Telt = permutalib::DoubleSidedPerm<Tidx>;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    int64_t n_iter = 50;
    std::vector<std::string> ListOpts = {
        "canonical",        "stabilizer",
        "pointstabilizer",  "pointrepresentative",
        "check_canonical",   "check_exhaustive_canonical",
        "check_store_canonical", "check_canonical_orbitsize",
        "timing_canonical_algorithms",
        "check_representative", "check_stabilizer",
        "all"};
    if (argc != 4 && argc != 5) {
      std::cerr << "BenchmarkingPermutalib [InputFile] [limit] [opt]\n";
      std::cerr << "or\n";
      std::cerr << "BenchmarkingPermutalib [InputFile] [limit] [opt] [n_iter]\n";
      std::cerr << "\n";
      std::cerr << "    ----- possible options ----\n";
      std::cerr << "\n";
      std::cerr << "InputFile: The list of groups as an input file\n";
      std::cerr << "limit: the limit in group size. E.g. 1000 or -1 for no limit\n";
      std::cerr << "opt: The list of allowed options:\n";
      for (auto &opt : ListOpts)
        std::cerr << " " << opt;
      std::cerr << "\n";
      std::cerr << "n_iter: the number of random iteration done. if absent then 50\n";
      throw permutalib::PermutalibException{1};
    }
    std::string InputFile = argv[1];
    std::cerr << "InputFile=" << InputFile << "\n";
    std::string limit_s = argv[2];
    int limit = permutalib::ReadScalar<int64_t>(limit_s);
    std::cerr << "limit=" << limit << "\n";
    std::string opt = argv[3];
    std::cerr << "opt=" << opt << "\n";
    if (argc == 5) {
      n_iter = permutalib::ReadScalar<int64_t>(argv[4]);
      std::cerr << "Using input value n_iter=" << n_iter << "\n";
    } else {
      std::cerr << "Using default value of 50 on n_iter\n";
    }
    std::cerr << "n_iter=" << n_iter << "\n";
    bool IsMatch = false;
    for (auto &e_opt : ListOpts) {
      if (e_opt == opt)
        IsMatch = true;
    }
    if (!IsMatch) {
      std::cerr << "opt=" << opt << " does not match\n";
      std::cerr << "Please select an option that is allowed\n";
      throw permutalib::PermutalibException{1};
    }
    //
    std::ifstream is(InputFile);
    size_t siz_control = 0;
    int nGroup;
    is >> nGroup;
    int n_group_treated = 0;
    for (int iGroup = 0; iGroup < nGroup; iGroup++) {
      //      std::cerr << "iGroup=" << iGroup << " / " << nGroup << "\n";
      size_t nbGen;
      int n_i;
      is >> nbGen;
      is >> n_i;
      Tidx n = Tidx(n_i);
      std::cerr << "iGroup=" << iGroup << "/" << nGroup << " n=" << int(n)
                << " nbGen=" << nbGen << "\n";
      std::vector<Telt> LGen(nbGen);
      for (size_t iGen = 0; iGen < nbGen; iGen++) {
        //        std::cerr << "iGen=" << iGen << "/" << nbGen << " n=" << int(n) << "\n";
        std::vector<Tidx> ePermV(n);
        for (Tidx i = 0; i < n; i++) {
          int eVal_i;
          is >> eVal_i;
          Tidx eVal = Tidx(eVal_i);
          if (eVal >= n) {
            std::cerr << "Values is above range\n";
            std::cerr << "i=" << int(i) << " n=" << int(n)
                      << " eVal=" << int(eVal) << "\n";
            throw permutalib::PermutalibException{1};
          }
          ePermV[i] = eVal;
        }
        Telt ePerm(ePermV);
        LGen[iGen] = ePerm;
      }
      Telt id(n);
      Tgroup eG(LGen, id);
      std::cerr << "  |eG|=" << eG.size() << "\n";
      //
      if (eG.size() < limit || limit < 0) {
        full_check(eG, opt, n_iter, siz_control);
        n_group_treated++;
      }
    }
    //
    std::cerr << "Control values (should stay the same because we do "
                 "deterministic random sets)="
              << siz_control << "\n";
    std::cerr << "Number of groups treated = " << n_group_treated << "\n";
    std::cerr << "CPP Normal completion of the program\n";
  } catch (permutalib::PermutalibException const &e) {
    std::cerr << "Erroneous completion of the program\n";
    exit(e.eVal);
  }
}
