#include "Permutation.h"
#include "gmpxx.h"
#include <fstream>

#include "Group.h"

int main(int argc, char *argv[])
{
  try {
    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    if (argc != 2 && argc != 3) {
      std::cerr << "GapStabilizerOnSet [EXMP]\n";
      std::cerr << "or\n";
      std::cerr << "GapStabilizerOnSet [EXMP] [OutFile]\n";
      std::cerr << "with EXMP generated by GenerateExample.g\n";
      throw PermutalibException{1};
    }
    std::string InputFile = argv[1];
    //
    std::ifstream is(InputFile);
    size_t nbGen;
    int n_i;
    is >> nbGen;
    is >> n_i;
    Tidx n = Tidx(n_i);
    std::vector<Telt> LGen(nbGen);
    for (size_t iGen=0; iGen<nbGen; iGen++) {
      std::vector<Tidx> ePermV(n);
      for (Tidx i=0; i<n; i++) {
	int eVal_i;
	is >> eVal_i;
	Tidx eVal = Tidx(eVal_i);
	ePermV[i]=eVal;
      }
      Telt ePerm(ePermV);
      LGen[iGen] = ePerm;
    }
    std::cerr.setf(std::ios::boolalpha);
    //
    std::cerr << "CPP Before call to MinimalStabChain\n";
    //    permutalib::StabChain<Telt> eG = permutalib::MinimalStabChain<Telt,Tint>(LGen, n);
    Telt id(n);
    permutalib::Group<Telt,Tint> eG(LGen, id);
    std::cerr << "CPP After call to MinimalStabChain\n";
    //
    std::cerr << "CPP |eG|=" << eG.size() << "\n";
    //
    permutalib::Face eFace(n);
    for (int i=0; i<n; i++) {
      int eVal;
      is >> eVal;
      eFace[i] = eVal;
    }
    permutalib::Group<Telt,Tint> eG2 = eG.Stabilizer_OnSets(eFace);
    std::cerr << "CPP |eG2|=" << eG2.size() << "\n";
    std::cerr << "CPP Normal completion of the program\n";
    //
    if (argc == 3) {
      std::string OutputFile = argv[2];
      std::ofstream os(OutputFile);
      os << "return " + eG2.GapString() + ";\n";
    } else {
      std::cerr << "CPP |eG2|=" << eG2.size() << "\n";
    }
  }
  catch (PermutalibException const& e) {
    std::cerr << "Erroneous completion of the program\n";
    exit(e.eVal);
  }
  return 0;
}
