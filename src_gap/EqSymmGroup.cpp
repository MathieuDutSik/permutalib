#include "Permutation.h"
#include "Group.h"
#include "gmpxx.h"
#include <fstream>

int main(int argc, char *argv[])
{
  try {
    using Tidx = int16_t;
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
    // First set of generators
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
    std::vector<Telt> LGen_A{ePerm1, ePerm2};
    permutalib::Group<Telt,Tint> eG_A = permutalib::Group<Telt,Tint>(LGen_A, n);
    //
    // Second set of generators
    //
    std::vector<Telt> LGen_B;
    for (size_t i=0; i<n-1; i++) {
      std::vector<Tidx> ePerm(n);
      for (int j=0; j<n; j++)
        ePerm[j] = j;
      ePerm[i] = i+1;
      ePerm[i+1] = i;
      Telt eElt(ePerm);
      LGen_B.push_back(eElt);
    }
    permutalib::Group<Telt,Tint> eG_B = permutalib::Group<Telt,Tint>(LGen_B, n);
    //
    // Now testing for equality
    //
    bool test = eG_A == eG_B;
    std::cerr << "test=" << test << "\n";
  }
  catch (PermutalibException const& e) {
    exit(e.eVal);
  }
  std::cerr << "Completion of the program\n";
}
