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
    if (argc != 3 && argc != 4) {
      std::cerr << "We should have argc = 2\n";
      std::cerr << "BenchmarkingPermutalib [EXMP] [opt]\n";
      std::cerr << "BenchmarkingPermutalib [EXMP] [opt] [n_iter]\n";
      throw PermutalibException{1};
    }
    std::string InputFile = argv[1];
    std::string opt = argv[2];
    std::vector<std::string> ListOpts={"canonical", "stabilizer", "pointstabilizer", "pointrepresentative", "all"};
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
    long n_iter = 50;
    if (argc == 4) {
      sscanf(argv[3], "%ld", &n_iter);
      std::cerr << "Using input value n_iter=" << n_iter << "\n";
    } else {
      std::cerr << "Using default value of 50 on n_iter\n";
    }
    //
    std::ifstream is(InputFile);
    int siz_control = 0;
    int nGroup;
    is >> nGroup;
    for (int iGroup=0; iGroup<nGroup; iGroup++) {
      //      std::cerr << "iGroup=" << iGroup << " / " << nGroup << "\n";
      size_t nbGen;
      Tidx n;
      is >> nbGen;
      is >> n;
      std::cerr << "iGroup=" << iGroup << " n=" << n << " nbGen=" << nbGen << "\n";
      std::vector<Telt> LGen(nbGen);
      for (size_t iGen=0; iGen<nbGen; iGen++) {
        std::vector<Tidx> ePermV(n);
        for (Tidx i=0; i<n; i++) {
          Tidx eVal;
          is >> eVal;
          ePermV[i]=eVal;
        }
        Telt ePerm(ePermV);
        LGen[iGen] = ePerm;
      }
      Tgroup eG(LGen, n);
      //
      auto bench_canonical=[&]() -> void {
        for (long iter=0; iter<n_iter; iter++) {
          permutalib::Face eFace(n);
          for (Tidx i=0; i<n; i++) {
            int eVal = Tidx(rand()) % 2;
            eFace[i] = eVal;
          }
          permutalib::Face set_can = eG.CanonicalImage(eFace);
          siz_control += set_can.count();
        }
      };
      auto bench_stabilizer=[&]() -> void {
        for (long iter=0; iter<n_iter; iter++) {
          permutalib::Face eFace(n);
          for (Tidx i=0; i<n; i++) {
            int eVal = Tidx(rand()) % 2;
            eFace[i] = eVal;
          }
          Tgroup eG2 = eG.Stabilizer_OnSets(eFace);
          siz_control += eG2.n_act();
        }
      };
      auto bench_pointstabilizer=[&]() -> void {
        for (long iter=0; iter<n_iter; iter++) {
          Tidx pos = Tidx(rand()) % n;
          Tgroup eG2 = eG.Stabilizer_OnPoints(pos);
          siz_control += eG2.n_act();
        }
      };
      auto bench_pointrepresentative=[&]() -> void {
        for (long iter=0; iter<n_iter; iter++) {
          Tidx pos1 = Tidx(rand()) % n;
          Tidx pos2 = Tidx(rand()) % n;
          std::pair<bool,Telt> eP = eG.RepresentativeAction_OnPoints(pos1, pos2);
          siz_control += int(eP.first);
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
      //
      if (opt == "all") {
        bench_canonical();
        bench_stabilizer();
        bench_pointstabilizer();
        bench_pointrepresentative();
      }
    }
    //
    std::cerr << "CPP Normal completion of the program\n";
  }
  catch (PermutalibException const& e) {
    std::cerr << "Erroneous completion of the program\n";
    exit(e.eVal);
  }
  return 0;
}
