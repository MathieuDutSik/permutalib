#include "Permutation.h"
#include "gmpxx.h"
#include <fstream>
#include "Group.h"





int main(int argc, char *argv[])
{
  try {
    //    using Tidx = int16_t;
    using Tidx = uint8_t;
    //using Telt = permutalib::DoubleSidedPerm<Tidx>;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt,Tint>;
    long n_iter = 50;
    if (argc != 3 && argc != 4) {
      std::cerr << "BenchmarkingPermutalib [EXMP] [opt]\n";
      std::cerr << "or\n";
      std::cerr << "BenchmarkingPermutalib [EXMP] [opt] [n_iter]\n";
      throw PermutalibException{1};
    }
    std::string InputFile = argv[1];
    std::string opt = argv[2];
    if (argc == 4) {
      sscanf(argv[3], "%ld", &n_iter);
      std::cerr << "Using input value n_iter=" << n_iter << "\n";
    } else {
      std::cerr << "Using default value of 50 on n_iter\n";
    }
    std::vector<std::string> ListOpts={"canonical", "stabilizer", "pointstabilizer", "pointrepresentative",
      "check_canonical", "check_representative", "check_stabilizer", "all"};
    bool IsMatch=false;
    for (auto & e_opt : ListOpts) {
      if (e_opt == opt)
        IsMatch=true;
    }
    if (!IsMatch) {
      std::cerr << "opt=" << opt << "\n";
      std::cerr << "ListOpts =";
      for (auto& e_opt : ListOpts)
        std::cerr << " " << e_opt;
      std::cerr << "\n";
      std::cerr << "Please select an option that is allowed\n";
      throw PermutalibException{1};
    }
    //
    std::ifstream is(InputFile);
    size_t siz_control = 0;
    int nGroup;
    is >> nGroup;
    for (int iGroup=0; iGroup<nGroup; iGroup++) {
      //      std::cerr << "iGroup=" << iGroup << " / " << nGroup << "\n";
      size_t nbGen;
      int n_i;
      is >> nbGen;
      is >> n_i;
      Tidx n=Tidx(n_i);
      std::cerr << "iGroup=" << iGroup << "/" << nGroup << " n=" << int(n) << " nbGen=" << nbGen << "\n";
      std::vector<Telt> LGen(nbGen);
      for (size_t iGen=0; iGen<nbGen; iGen++) {
        std::cerr << "iGen=" << iGen << "/" << nbGen << " n=" << int(n) << "\n";
        std::vector<Tidx> ePermV(n);
        for (Tidx i=0; i<n; i++) {
          int eVal_i;
          is >> eVal_i;
          Tidx eVal = Tidx(eVal_i);
          if (eVal >= n) {
            std::cerr << "Values is above range\n";
            std::cerr << "i=" << int(i) << " n=" << int(n) << " eVal=" << int(eVal) << "\n";
            throw PermutalibException{1};
          }
          ePermV[i]=eVal;
        }
        Telt ePerm(ePermV);
        LGen[iGen] = ePerm;
      }
      Telt id(n);
      Tgroup eG(LGen, id);
      std::cerr << "  |eG|=" << eG.size() << "\n";
      //
      auto random_face=[](const Tidx& len) -> permutalib::Face {
        permutalib::Face eFace(len);
        for (Tidx i=0; i<len; i++) {
          int eVal = Tidx(rand()) % 2;
          eFace[i] = eVal;
        }
        return eFace;
      };
      auto bench_canonical=[&]() -> void {
        for (long iter=0; iter<n_iter; iter++) {
          permutalib::Face eFace = random_face(n);
          permutalib::Face set_can = eG.CanonicalImage(eFace);
          siz_control += set_can.count();
        }
      };
      auto bench_stabilizer=[&]() -> void {
        for (long iter=0; iter<n_iter; iter++) {
          permutalib::Face eFace = random_face(n);
          Tgroup eG2 = eG.Stabilizer_OnSets(eFace);
          siz_control += eFace.count();
        }
      };
      auto bench_pointstabilizer=[&]() -> void {
        for (long iter=0; iter<n_iter; iter++) {
          Tidx pos = Tidx(rand()) % n;
          Tgroup eG2 = eG.Stabilizer_OnPoints(pos);
          siz_control += size_t(pos);
        }
      };
      auto bench_pointrepresentative=[&]() -> void {
        for (long iter=0; iter<n_iter; iter++) {
          Tidx pos1 = Tidx(rand()) % n;
          Tidx pos2 = Tidx(rand()) % n;
          std::optional<Telt> eP = eG.RepresentativeAction_OnPoints(pos1, pos2);
          siz_control += size_t(pos1) + size_t(pos2);
          if (eP)
            siz_control++;
        }
      };
      auto check_canonical=[&]() -> void {
        for (long iter=0; iter<n_iter; iter++) {
          permutalib::Face eFace1 = random_face(n);
          permutalib::Face set_can1 = eG.CanonicalImage(eFace1);
          for (int i=0; i<4; i++) {
            Telt u = eG.rand();
            permutalib::Face eFace2 = OnSets(eFace1, u);
            permutalib::Face set_can2 = eG.CanonicalImage(eFace2);
            if (set_can1 != set_can2) {
              std::cerr << "Canonicalization failed\n";
              throw PermutalibException{1};
            }
          }
          siz_control += eFace1.count();
        }
      };
      auto check_representative=[&]() -> void {
        for (long iter=0; iter<n_iter; iter++) {
          permutalib::Face eFace1 = random_face(n);
          Telt u = eG.rand();
          permutalib::Face eFace2 = OnSets(eFace1, u);
          std::optional<Telt> test = eG.RepresentativeAction_OnSets(eFace1, eFace2);
          if (!test) {
            std::cerr << "RepresentativeAction_OnSets error\n";
            throw PermutalibException{1};
          }
          siz_control += eFace1.count();
        }
      };
      auto check_stabilizer=[&]() -> void {
        for (long iter=0; iter<n_iter; iter++) {
          permutalib::Face eFace1 = random_face(n);
          Tgroup eG1 = eG.Stabilizer_OnSets(eFace1);
          for (int i=0; i<4; i++) {
            Telt u = eG.rand();
            permutalib::Face eFace2 = OnSets(eFace1, u);
            Tgroup eG2 = eG.Stabilizer_OnSets(eFace2);
            if (eG1.size() != eG2.size()) {
              std::cerr << "Stabilizer_OnSets error\n";
              throw PermutalibException{1};
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
    std::cerr << "Control values (should stay the same because we do deterministic random sets)=" << siz_control << "\n";
    std::cerr << "CPP Normal completion of the program\n";
  }
  catch (PermutalibException const& e) {
    std::cerr << "Erroneous completion of the program\n";
    exit(e.eVal);
  }
}
