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
      std::cerr << "We should have argc = 2\n";
      std::cerr << "GroupProperties [EXMP]\n";
      std::cerr << "or\n";
      std::cerr << "GroupProperties [EXMP] [OutFile]\n";
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
    //
    permutalib::Group<Telt,Tint> eG = permutalib::Group<Telt,Tint>(LGen, n);
    bool IsPrimitive = eG.IsPrimitive();
    bool IsTransitive = eG.IsTransitive();
    bool IsCommutative = eG.IsCommutative();
    bool IsCyclic = eG.IsCyclic();
    //
    auto prt=[&](std::ostream & os) -> void {
      auto fct=[&](const bool& val) -> std::string {
        if (val)
          return "true";
        return "false";
      };
      os << "return rec(IsPrimitive:=" << fct(IsPrimitive) <<
        ", IsTransitive:=" << fct(IsTransitive) <<
        ", IsCommutative:=" << fct(IsCommutative) << 
        ", IsCyclic:=" << fct(IsCyclic) << ");\n\n";
    };

    if (argc == 3) {
      std::string OutputFile = argv[2];
      std::ofstream os(OutputFile);
      prt(os);
    } else {
      prt(std::cerr);
    }
  }
  catch (PermutalibException const& e) {
    std::cerr << "Erroneous completion of the program\n";
    exit(e.eVal);
  }
  return 0;
}
