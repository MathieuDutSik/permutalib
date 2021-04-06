#include "Permutation.h"
#include "gmpxx.h"
#include <fstream>

#include "Group.h"

int main(int argc, char *argv[])
{
  try {
    using Tidx = int16_t;
    using Telt = permutalib::DoubleSidedPerm<Tidx>;
    using Tint = mpz_class;
    if (argc != 3) {
      std::cerr << "We should have argc = 2\n";
      std::cerr << "GapEquivalenceOnSet [EXMP] |OutFile]\n";
      std::cerr << "with EXMP generated by GenerateExample.g\n";
      throw PermutalibException{1};
    }
    std::ifstream is(argv[1]);
    int nbGen, n;
    is >> nbGen;
    is >> n;
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
    std::cerr.setf(std::ios::boolalpha);
    //
    std::cerr << "CPP Before call to MinimalStabChain\n";
    //    permutalib::StabChain<Telt> eG = permutalib::MinimalStabChain<Telt,Tint>(LGen, n);
    permutalib::Group<Telt,Tint> eG = permutalib::Group<Telt,Tint>(LGen, n);
    std::cerr << "CPP After call to MinimalStabChain\n";
    //
    std::cerr << "CPP |eG|=" << eG.size() << "\n";
    //
    permutalib::Face f1(n);
    for (int i=0; i<n; i++) {
      int eVal;
      is >> eVal;
      f1[i] = eVal;
    }
    //
    permutalib::Face f2(n);
    for (int i=0; i<n; i++) {
      int eVal;
      is >> eVal;
      f2[i] = eVal;
    }
    std::pair<bool,Telt> epair = eG.RepresentativeAction_OnSets(f1, f2);
    //
    std::ofstream os(argv[2]);
    if (epair.first) {
      os << "return " << epair.second << ";\n";
    } else {
      os << "return fail;\n";
    }
    std::cerr << "CPP Normal completion of the program\n";
  }
  catch (PermutalibException const& e) {
    std::cerr << "Erroneous completion of the program\n";
    exit(e.eVal);
  }
  return 0;
}
