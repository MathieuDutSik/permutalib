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
    if (argc != 2) {
      std::cerr << "We should have argc = 2\n";
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
    //
    permutalib::StabChain<Telt> eG = permutalib::MinimalStabChain<Telt,Tint>(LGen, n);
    std::cerr << "eG=" << eG << "\n";
    //
    std::cerr << "|eG|=" << permutalib::Order<Telt,Tint>(eG) << "\n";
    //
    Face eFace(n);
    for (int i=0; i<n; i++) {
      int eVal;
      is >> eVal;
      eFace[i] = eVal;
    }
    permutalib::StabChain<Telt> eG2 = permutalib::Stabilizer_OnSets<Telt,Tint>(eG, eFace);
    std::cerr << "eG2=" << eG2 << "\n";
    std::cerr << "|eG2|=" << permutalib::Order<Telt,Tint>(eG2) << "\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
  std::cerr << "Completion of the program\n";
}
