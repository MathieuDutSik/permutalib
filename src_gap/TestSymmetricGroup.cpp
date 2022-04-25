#include "Permutation.h"
#include "StabChainMain.h"
#include "gmpxx.h"
#include <fstream>

int main(int argc, char *argv[]) {
  try {
    using Tidx = uint16_t;
    using Tidx_label = uint16_t;
    using Telt = permutalib::DoubleSidedPerm<Tidx>;
    using Tint = mpz_class;
    if (argc != 2) {
      std::cerr << "The program is used as\n";
      std::cerr << "TestSymmetricGroup [n]\n";
      throw permutalib::PermutalibException{1};
    }
    int n;
    (void)sscanf(argv[1], "%d", &n);
    //
    std::vector<Tidx> ePermV1(n);
    for (int i = 0; i < n; i++) {
      int iNext = i + 1;
      if (iNext == n)
        iNext = 0;
      ePermV1[i] = iNext;
    }
    Telt ePerm1(ePermV1);
    //
    std::vector<Tidx> ePermV2(n);
    for (int i = 0; i < n; i++)
      ePermV2[i] = i;
    ePermV2[1] = 0;
    ePermV2[0] = 1;
    Telt ePerm2(ePermV2);
    //
    std::vector<Telt> LGen{ePerm1, ePerm2};
    //
    Telt id(n);
    permutalib::StabChain<Telt, Tidx_label> S =
        permutalib::MinimalStabChain<Telt, Tidx_label, Tint>(LGen, id);
    std::cerr << "S=" << S << "\n";
  } catch (permutalib::PermutalibException const &e) {
    exit(e.eVal);
  }
  std::cerr << "Completion of the program\n";
}
