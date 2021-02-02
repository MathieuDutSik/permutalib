#include "list.h"
#include "Permutation.h"
#include "StabChainMain.h"
#include "NumberTheory.h"
#include "stbcbckt.h"
//#include "partition.h"

int main(int argc, char *argv[])
{
  try {
    using Telt = permutalib::DoubleSidedPerm;
    using Tint = mpz_class;
    if (argc != 3) {
      std::cerr << "We should have argc = 2\n";
      std::cerr << "GapEquivalenceOnSet [EXMP] |OutFile]\n";
      std::cerr << "with EXMP generated by GenerateExample.g\n";
      throw TerminalException{1};
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
    permutalib::StabChain<Telt> eG = permutalib::StabChainOp_group_options<Telt,Tint>(LGen, n);
    std::cerr << "CPP After call to MinimalStabChain\n";
    std::cerr << "CPP eG=" << eG << "\n";
    //
    std::cerr << "CPP |eG|=" << permutalib::Order<Telt,Tint>(eG) << "\n";
    //
    Face f1(n);
    for (int i=0; i<n; i++) {
      int eVal;
      is >> eVal;
      f1[i] = eVal;
    }
    //
    Face f2(n);
    for (int i=0; i<n; i++) {
      int eVal;
      is >> eVal;
      f2[i] = eVal;
    }
    std::pair<bool,Telt> epair = permutalib::RepresentativeAction_OnSets<Telt,Tint>(eG, f1, f2);
    //
    std::ofstream os(argv[2]);
    if (test.first) {
      os << "return " << test.second << ";\n";
    } else {
      os << "return fail;\n";
    }
    std::cerr << "CPP Normal completion of the program\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Erroneous completion of the program\n";
    exit(e.eVal);
  }
  return 0;
}
