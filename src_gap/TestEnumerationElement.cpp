#include "Permutation.h"
#include "gmpxx.h"
#include <fstream>

#include "Group.h"

int main(int argc, char *argv[])
{
  try {
    using Tidx = uint8_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt,Tint>;
    if (argc != 3) {
      std::cerr << "TestEnumerationElement [EXMP] [limit]\n";
      throw PermutalibException{1};
    }
    std::string InputFile = argv[1];
    int limit_i;
    (void)sscanf(argv[2], "%d", &limit_i);
    
    Tint limit = limit_i;
    //
    std::ifstream is(InputFile);
    int nGroup;
    is >> nGroup;
    size_t n_treat = 0;
    for (int iGroup=0; iGroup<nGroup; iGroup++) {
      size_t nbGen;
      int n_i;
      is >> nbGen;
      is >> n_i;
      Tidx n=Tidx(n_i);
      std::cerr << "iGroup=" << iGroup << "/" << nGroup << " n=" << int(n) << " nbGen=" << nbGen << "\n";
      std::vector<Telt> LGen(nbGen);
      for (size_t iGen=0; iGen<nbGen; iGen++) {
        std::cerr << "iGen=" << iGen << "/" << nbGen << " n=" << int(n) << "\n";
        std::vector<Tidx> ePermV(n);
        for (Tidx i=0; i<n; i++) {
          int eVal_i;
          is >> eVal_i;
          Tidx eVal = Tidx(eVal_i);
          if (eVal >= n) {
            std::cerr << "Values is above range\n";
            std::cerr << "i=" << int(i) << " n=" << int(n) << " eVal=" << int(eVal) << "\n";
            throw PermutalibException{1};
          }
          ePermV[i]=eVal;
        }
        Telt ePerm(ePermV);
        LGen[iGen] = ePerm;
      }
      Tgroup eG(LGen, n);
      std::cerr << "|eG|=" << eG.size() << " limit=" << limit << "\n";
      if (eG.size() <= limit) { // Some groups can be too large to iterate over their elements.
        std::unordered_set<Telt> ListElt;
        size_t n_iter=0;
        for (auto & elt : eG) {
          ListElt.insert(elt);
          if (Tint(ListElt.size()) > eG.size()) {
            std::cerr << "We found more elements than the size of the group. Clear bug\n";
            throw PermutalibException{1};
          }
          n_iter++;
        }
        if (n_iter != ListElt.size()) {
          std::cerr << "Some elements were found several times. Clear bug\n";
          throw PermutalibException{1};
        }
        if (Tint(n_iter) != eG.size()) {
          std::cerr << "Some elements were missed. Clear bug\n";
          throw PermutalibException{1};
        }
        n_treat++;
      }
    }
    std::cerr << "n_treat = " << n_treat << "\n";
    //
    std::cerr << "Normal completion of the program\n";
  }
  catch (PermutalibException const& e) {
    std::cerr << "Erroneous completion of the program\n";
    exit(e.eVal);
  }
  return 0;
}
