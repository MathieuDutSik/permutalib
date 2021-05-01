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
    if (opt != "canonical" && opt != "stabilizer") {
      std::cerr << "We should have opt=canonical or stabilizer\n";
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
      int nbGen, n;
      is >> nbGen;
      is >> n;
      std::cerr << "iGroup=" << iGroup << " n=" << n << " nbGen=" << nbGen << "\n";
      std::vector<Telt> LGen(nbGen);
      for (int iGen=0; iGen<nbGen; iGen++) {
        std::vector<Tidx> ePermV(n);
        for (int i=0; i<n; i++) {
          int eVal;
          is >> eVal;
          ePermV[i]=eVal;
        }
        Telt ePerm(ePermV);
        LGen[iGen] = ePerm;
      }
      Tgroup eG(LGen, n);
      //
      for (long iter=0; iter<n_iter; iter++) {
        permutalib::Face eFace(n);
        for (int i=0; i<n; i++) {
          int eVal = rand() % 2;
          eFace[i] = eVal;
        }
        //        std::cerr << "  iter=" << iter << " |eFace|=" << eFace.count() << " / " << eFace.size() << "\n";
        //
        if (opt == "canonical") {
          permutalib::Face set_can = eG.CanonicalImage(eFace);
          siz_control += set_can.count();
        }
        if (opt == "stabilizer") {
          Tgroup eG2 = eG.Stabilizer_OnSets(eFace);
          siz_control += eG2.n_act();
        }
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
