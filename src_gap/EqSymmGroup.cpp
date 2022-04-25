#include "Group.h"
#include "Permutation.h"
#include "gmpxx.h"
#include <fstream>

int main(int argc, char *argv[]) {
  try {
    using Tidx = uint16_t;
    using Telt = permutalib::DoubleSidedPerm<Tidx>;
    using Tint = mpz_class;
    if (argc != 2) {
      std::cerr << "The program is used as\n";
      std::cerr << "TestSymmetricGroup [n]\n";
      throw permutalib::PermutalibException{1};
    }
    int n_pre;
    (void)sscanf(argv[1], "%d", &n_pre);
    size_t n = size_t(n_pre);
    //
    // First set of generators
    //
    std::vector<Tidx> ePermV1(n);
    for (size_t i = 0; i < n; i++) {
      size_t iNext = i + 1;
      if (iNext == n)
        iNext = 0;
      ePermV1[i] = Tidx(iNext);
    }
    Telt ePerm1(ePermV1);
    //
    std::vector<Tidx> ePermV2(n);
    for (size_t i = 0; i < n; i++)
      ePermV2[i] = Tidx(i);
    ePermV2[1] = 0;
    ePermV2[0] = 1;
    Telt ePerm2(ePermV2);
    //
    std::vector<Telt> LGen_A{ePerm1, ePerm2};
    Telt id(n);
    permutalib::Group<Telt, Tint> eG_A(LGen_A, id);
    //
    // Second set of generators
    //
    std::vector<Telt> LGen_B;
    for (size_t i = 0; i < n - 1; i++) {
      std::vector<Tidx> ePerm(n);
      for (size_t j = 0; j < n; j++)
        ePerm[j] = Tidx(j);
      ePerm[i] = Tidx(i + 1);
      ePerm[i + 1] = Tidx(i);
      Telt eElt(ePerm);
      LGen_B.push_back(eElt);
    }
    permutalib::Group<Telt, Tint> eG_B(LGen_B, id);
    //
    // Now testing for equality
    //
    bool test = eG_A == eG_B;
    std::cerr << "test=" << test << "\n";
  } catch (permutalib::PermutalibException const &e) {
    exit(e.eVal);
  }
  std::cerr << "Completion of the program\n";
}
