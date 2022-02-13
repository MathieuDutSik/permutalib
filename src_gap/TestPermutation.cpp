#include "Permutation.h"
#include "Group.h"
#include "gmpxx.h"
#include <fstream>

int main(int argc, char *argv[])
{
  try {
    using Tidx = uint16_t;
    //    using Telt = permutalib::DoubleSidedPerm<Tidx>;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    if (argc != 2) {
      std::cerr << "The program is used as\n";
      std::cerr << "TestPermutation [n]\n";
      throw PermutalibException{1};
    }
    int n_i;
    (void)sscanf(argv[1], "%d", &n_i);
    Tidx n = Tidx(n_i);
    // Testing initializer_list
    //
    Telt p = {1,0,2,3};
    //
    std::vector<Tidx> ePermV1(n);
    for (Tidx i=0; i<n; i++) {
      Tidx iNext=i+1;
      if (iNext == n)
	iNext=0;
      ePermV1[i] = iNext;
    }
    Telt ePerm1(ePermV1);
    //
    std::vector<Tidx> ePermV2(n);
    for (Tidx i=0; i<n; i++)
      ePermV2[i] = i;
    ePermV2[1]=0;
    ePermV2[0]=1;
    Telt ePerm2(ePermV2);
    //
    std::vector<Telt> LGen{ePerm1, ePerm2};
    Telt id(n);
    permutalib::Group<Telt,Tint> SymmGrp(LGen, id);
    //
    for (int i=0; i<400; i++) {
      std::cerr << "i=" << i << "\n";
      Telt elt1 = SymmGrp.rand();
      Telt elt2 = SymmGrp.rand();
      Telt prod = elt1 * elt2;
      elt1 *= elt2;
      if (elt1 != prod) {
        std::cerr << "Error in the product operator\n";
        throw PermutalibException{1};
      }
      for (Tidx i=0; i<n; i++) {
        Tidx img1 = SlashAct(i, elt1);
        Tidx img2 = PowAct(img1, elt1);
        if (img2 != i) {
          std::cerr << "We have i=" << i << " but img2=" << img2 << "\n";
          std::cerr << "Please correct the bug\n";
          throw PermutalibException{1};
        }
      }
    }

    std::cerr << "Normal completion of the program\n";
  }
  catch (PermutalibException const& e) {
    std::cerr << "Something wrong happened\n";
    exit(e.eVal);
  }
}
