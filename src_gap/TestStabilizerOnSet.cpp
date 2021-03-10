#include "list.h"
#include "Permutation.h"
#include "StabChainMain.h"
#include "gmpxx.h"
#include "stbcbckt.h"
#include <fstream>

int main(int argc, char *argv[])
{
  try {
    using Telt = permutalib::DoubleSidedPerm;
    using Tint = mpz_class;
    if (argc != 2) {
      std::cerr << "We should have argc = 2\n";
      std::cerr << "TestStabilizerOnSet [EXMP]\n";
      std::cerr << "with EXMP generated by GenerateExample.g\n";
      throw PermutalibException{1};
    }
    std::ifstream is(argv[1]);
    int nbGen, n;
    is >> nbGen;
    is >> n;
    std::vector<Telt> LGen(nbGen);
    for (int iGen=0; iGen<nbGen; iGen++) {
      std::vector<int> ePermV(n);
      for (int i=0; i<n; i++) {
	int eVal;
	is >> eVal;
	ePermV[i]=eVal;
      }
      Telt ePerm(ePermV);
      LGen[iGen] = ePerm;
    }
    std::cerr.setf(std::ios::boolalpha);
    //
    std::cerr << "CPP Before call to MinimalStabChain\n";
    //    permutalib::StabChain<Telt> eG = permutalib::MinimalStabChain<Telt,Tint>(LGen, n);
    permutalib::StabChain<Telt> eG = permutalib::Group<Telt,Tint>(LGen, n);
    //    std::cerr << "CPP After call to MinimalStabChain\n";
    //    PrintStabChain(eG);
    //    std::cerr << "CPP eG=" << eG << "\n";
    //
    std::cerr << "CPP |eG|=" << permutalib::Order<Telt,Tint>(eG) << "\n";
    //
    Face eFace(n);
    for (int i=0; i<n; i++) {
      int eVal;
      is >> eVal;
      eFace[i] = eVal;
    }
    permutalib::StabChain<Telt> eG2 = permutalib::Stabilizer_OnSets<Telt,Tint>(eG, eFace);
    std::cerr << "CPP eG2=" << eG2 << "\n";
    std::cerr << "CPP |eG2|=" << permutalib::Order<Telt,Tint>(eG2) << "\n";
    std::cerr << "CPP Normal completion of the program\n";
  }
  catch (PermutalibException const& e) {
    std::cerr << "Erroneous completion of the program\n";
    exit(e.eVal);
  }
  return 0;
}
