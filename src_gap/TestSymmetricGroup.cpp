#include "Permutation.h"
#include "StabChainMain.h"
#include "gmpxx.h"
#include <fstream>

int main(int argc, char *argv[])
{
  try {
    using Tidx = uint16_t;
    using Telt = permutalib::DoubleSidedPerm<Tidx>;
    using Tint = mpz_class;
    if (argc != 2) {
      std::cerr << "The program is used as\n";
      std::cerr << "TestSymmetricGroup [n]\n";
      throw PermutalibException{1};
    }
    int n;
    (void)sscanf(argv[1], "%d", &n);
    //
    std::vector<Tidx> ePermV1(n);
    for (int i=0; i<n; i++) {
      int iNext=i+1;
      if (iNext == n)
	iNext=0;
      ePermV1[i] = iNext;
    }
    Telt ePerm1(ePermV1);
    //
    std::vector<Tidx> ePermV2(n);
    for (int i=0; i<n; i++)
      ePermV2[i] = i;
    ePermV2[1]=0;
    ePermV2[0]=1;
    Telt ePerm2(ePermV2);
    //
    std::vector<Telt> LGen{ePerm1, ePerm2};
    //
    permutalib::StabChain<Telt> S = permutalib::MinimalStabChain<Telt,Tint>(LGen, n);
    std::cerr << "S=" << S << "\n";
  }
  catch (PermutalibException const& e) {
    exit(e.eVal);
  }
  std::cerr << "Completion of the program\n";
}
