// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Permutation.h"
#include "gmpxx.h"
#include <fstream>

#include "Group.h"

//#define USE_LIBGAP
#undef USE_LIBGAP

#ifdef USE_LIBGAP

#include <libgap-api.h>

std::string GAP_string_face(const permutalib::Face &f) {
  std::string comm2 = "[";
  boost::dynamic_bitset<>::size_type ePt = f.find_first();
  for (size_t u = 0; u < f.count(); u++) {
    if (u > 0)
      comm2 += ",";
    comm2 += std::to_string(ePt + 1);
    ePt = f.find_next(ePt);
  }
  comm2 += "]";
  return comm2;
}

template <typename Tgroup>
bool GetRepresentativeAction_OnSets_libgap(const Tgroup &G,
                                           const permutalib::Face &eFace1,
                                           const permutalib::Face &eFace2) {
  std::string comm1 = "G:=" + G.GapString() + ";";
  std::string comm2 = "f1:=" + GAP_string_face(eFace1) + ";";
  std::string comm3 = "f2:=" + GAP_string_face(eFace2) + ";";
  std::string comm4 = "RepresentativeAction(G, f1, f2, OnSets)<>fail;";
  std::string commtot = comm1 + ";  " + comm2 + ";  " + comm3 + ";  " + comm4;
  //  std::cerr << "commtot=" << commtot << "\n";
  //
  Obj res = GAP_EvalString(commtot.c_str());
  Int rc = GAP_LenList(res);
  //  std::cerr << "rc=" << rc << "\n";
  for (Int i = 1; i <= rc; i++) {
    Obj ires = GAP_ElmList(res, 1);
    if (GAP_ElmList(ires, 1) == GAP_True) {
      Char *buffer = GAP_CSTR_STRING(GAP_ElmList(ires, 5));
      if (buffer) {
        //        printf("%s\n", buffer);
      }
      if (strlen(buffer) > 0) {
        std::string str = buffer;
        if (str == "true")
          return true;
        if (str == "false")
          return false;
        std::cerr << "Failed to find a matching entry for str=" << str << "\n";
        throw permutalib::PermutalibException{1};
      }
    }
  }
  return true;
  std::cerr << "We should not reach that stage\n";
  throw permutalib::PermutalibException{1};
}

template <typename Tgroup>
Tgroup GetStabilizer_libgap(const Tgroup &G, const permutalib::Face &f) {
  std::string comm1 = "G:=" + G.GapString() + ";";
  std::string comm2 = "f:=" + GAP_string_face(f) + ";";
  std::string comm3 = "GeneratorsOfGroup(Stabilizer(G, f, OnSets));";
  std::string commtot = comm1 + ";  " + comm2 + ";  " + comm3;
  std::cerr << "commtot=" << commtot << "\n";
  //
  Obj res = GAP_EvalString(commtot.c_str());
  Int rc = GAP_LenList(res);
  for (Int i = 1; i <= rc; i++) {
    Obj ires = GAP_ElmList(res, 1);
    if (GAP_ElmList(ires, 1) == GAP_True) {
      Char *buffer = GAP_CSTR_STRING(GAP_ElmList(ires, 5));
      if (buffer) {
        std::cerr << "Parsing code is needed\n";
        printf("%s\n", buffer);
      }
    }
  }
  std::cerr << "Code is incomplete\n";
  return G.Stabilizer_OnSets(f);
}

void set_argc_argv_gap(int *argc, char ***argv) {
  std::vector<std::string> LStr = {
      "./BenchmarkingPermutalib", "-A", "-l", ".", "-q", "-T", "--nointeract"};
  *argc = LStr.size();
  *argv = (char **)malloc(LStr.size() * sizeof(char **));
  for (size_t i_str = 0; i_str < LStr.size(); i_str++) {
    size_t len_str = LStr[i_str].size();
    (*argv)[i_str] = (char *)malloc((len_str + 1) * sizeof(char *));
    for (size_t i_c = 0; i_c < len_str; i_c++)
      (*argv)[i_str][i_c] = LStr[i_str][i_c];
    (*argv)[i_str][len_str] = '\0';
  }
}

void free_argc_argv_gap(int *argc, char ***argv) {
  size_t len = *argc;
  for (size_t i_str = 0; i_str < len; i_str++) {
    free((*argv)[i_str]);
  }
  free(*argv);
}

#endif

int main(int argc, char *argv[]) {
#ifdef USE_LIBGAP
  int argc_gap;
  char **argv_gap;
  std::cerr << "Before set_argc_argv_gap\n";
  set_argc_argv_gap(&argc_gap, &argv_gap);
  for (int i = 0; i < argc; i++) {
    std::cerr << "main : i=" << i << " argv=" << argv[i] << "\n";
  }
  for (int i = 0; i < argc_gap; i++) {
    std::cerr << "gap : i=" << i << " argv=" << argv_gap[i] << "\n";
  }

  std::cerr << "Before GAP_Initialize\n";
  GAP_Initialize(argc, argv, 0, 0, 1);
  std::cerr << "Before GAP_Enter\n";
  Int ok = GAP_Enter();
  std::cerr << "ok=" << ok << "\n";

#endif
  try {
    //    using Tidx = int16_t;
    using Tidx = uint8_t;
    // using Telt = permutalib::DoubleSidedPerm<Tidx>;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    long n_iter = 50;
#ifdef USE_LIBGAP
    std::string InputFile = "AllFileTest";
    std::string opt = "check_representative";
#else
    if (argc != 3 && argc != 4) {
      std::cerr << "BenchmarkingPermutalib [EXMP] [opt]\n";
      std::cerr << "or\n";
      std::cerr << "BenchmarkingPermutalib [EXMP] [opt] [n_iter]\n";
      throw permutalib::PermutalibException{1};
    }
    std::string InputFile = argv[1];
    std::string opt = argv[2];
    if (argc == 4) {
      sscanf(argv[3], "%ld", &n_iter);
      std::cerr << "Using input value n_iter=" << n_iter << "\n";
    } else {
      std::cerr << "Using default value of 50 on n_iter\n";
    }
#endif
    std::vector<std::string> ListOpts = {
        "canonical",        "stabilizer",
        "pointstabilizer",  "pointrepresentative",
        "check_canonical",  "check_representative",
        "check_stabilizer", "all"};
    bool IsMatch = false;
    for (auto &e_opt : ListOpts) {
      if (e_opt == opt)
        IsMatch = true;
    }
    if (!IsMatch) {
      std::cerr << "opt=" << opt << "\n";
      std::cerr << "ListOpts =";
      for (auto &e_opt : ListOpts)
        std::cerr << " " << e_opt;
      std::cerr << "\n";
      std::cerr << "Please select an option that is allowed\n";
      throw permutalib::PermutalibException{1};
    }
    //
    std::ifstream is(InputFile);
    int siz_control = 0;
    int nGroup;
    is >> nGroup;
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
        std::cerr << "iGen=" << iGen << "/" << nbGen << " n=" << int(n) << "\n";
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
      auto random_face = [](const Tidx &len) -> permutalib::Face {
        permutalib::Face eFace(len);
        for (Tidx i = 0; i < len; i++) {
          int eVal = Tidx(rand()) % 2;
          eFace[i] = eVal;
        }
        return eFace;
      };
      auto bench_canonical = [&]() -> void {
        for (long iter = 0; iter < n_iter; iter++) {
          permutalib::Face eFace = random_face(n);
          permutalib::Face set_can = eG.CanonicalImage(eFace);
          siz_control += set_can.count();
        }
      };
      auto bench_stabilizer = [&]() -> void {
        for (long iter = 0; iter < n_iter; iter++) {
          permutalib::Face eFace = random_face(n);
#ifdef USE_LIBGAP
          Tgroup eG2 = GetStabilizer_libgap(eG, eFace);
#else
          Tgroup eG2 = eG.Stabilizer_OnSets(eFace);
#endif
          siz_control += eG2.n_act();
        }
      };
      auto bench_pointstabilizer = [&]() -> void {
        for (long iter = 0; iter < n_iter; iter++) {
          Tidx pos = Tidx(rand()) % n;
          Tgroup eG2 = eG.Stabilizer_OnPoints(pos);
          siz_control += eG2.n_act();
        }
      };
      auto bench_pointrepresentative = [&]() -> void {
        for (long iter = 0; iter < n_iter; iter++) {
          Tidx pos1 = Tidx(rand()) % n;
          Tidx pos2 = Tidx(rand()) % n;
          std::pair<bool, Telt> eP =
              eG.RepresentativeAction_OnPoints(pos1, pos2);
          siz_control += int(eP.first);
        }
      };
      auto check_canonical = [&]() -> void {
        for (long iter = 0; iter < n_iter; iter++) {
          permutalib::Face eFace1 = random_face(n);
          permutalib::Face set_can1 = eG.CanonicalImage(eFace1);
          for (int i = 0; i < 4; i++) {
            Telt u = eG.rand();
            permutalib::Face eFace2 = OnSets(eFace1, u);
            permutalib::Face set_can2 = eG.CanonicalImage(eFace2);
            if (set_can1 != set_can2) {
              std::cerr << "Canonicalization failed\n";
              throw permutalib::PermutalibException{1};
            }
          }
        }
      };
      auto check_representative = [&]() -> void {
        for (long iter = 0; iter < n_iter; iter++) {
          permutalib::Face eFace1 = random_face(n);
          Telt u = eG.rand();
          permutalib::Face eFace2 = OnSets(eFace1, u);
#ifdef USE_LIBGAP
          bool test = GetRepresentativeAction_OnSets_libgap(eG, eFace1, eFace2);
#else
          bool test = eG.RepresentativeAction_OnSets(eFace1, eFace2).first;
#endif
          if (!test) {
            std::cerr << "RepresentativeAction_OnSets error\n";
            throw permutalib::PermutalibException{1};
          }
        }
      };
      auto check_stabilizer = [&]() -> void {
        for (long iter = 0; iter < n_iter; iter++) {
          permutalib::Face eFace1 = random_face(n);
          Tgroup eG1 = eG.Stabilizer_OnSets(eFace1);
          for (int i = 0; i < 4; i++) {
            Telt u = eG.rand();
            permutalib::Face eFace2 = OnSets(eFace1, u);
            Tgroup eG2 = eG.Stabilizer_OnSets(eFace2);
            if (eG1.size() != eG2.size()) {
              std::cerr << "Stabilizer_OnSets error\n";
              throw permutalib::PermutalibException{1};
            }
          }
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
        check_representative();
        check_stabilizer();
      }
    }
    //
    std::cerr << "CPP Normal completion of the program\n";
  } catch (permutalib::PermutalibException const &e) {
    std::cerr << "Erroneous completion of the program\n";
    exit(e.eVal);
  }
#ifdef USE_LIBGAP
  GAP_Leave();
  free_argc_argv_gap(&argc_gap, &argv_gap);
#endif
}
