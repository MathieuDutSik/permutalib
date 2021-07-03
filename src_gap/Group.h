#ifndef DEFINE_PERMUTALIB_GROUP_H
#define DEFINE_PERMUTALIB_GROUP_H


#include "StabChainMain.h"
#include "stbcbckt.h"
#include "nsi.h"
#include <map>


namespace permutalib {



template<typename Telt>
Telt RandomElement(const std::vector<Telt>& LGen, const int& n)
{
  size_t len = rand() % 100;
  size_t n_gen = LGen.size();
  Telt eElt(n);
  for (size_t iIter=0; iIter<len; iIter++) {
    size_t pos = size_t(rand()) % n_gen;
    eElt = eElt * LGen[pos];
  }
  return eElt;
}





template<typename Telt_inp, typename Tint_inp>
struct Group {
public:
  // dependent types
  using Telt = Telt_inp;
  using Tidx = typename Telt::Tidx;
  using Tint = Tint_inp;
  using Tidx_label = uint16_t;
  // constructors
  Group(const StabChain<Telt,Tidx_label>& _S) : S(_S), size_tint(Order<Telt,Tidx_label,Tint>(_S))
  {
  }
  Group(const std::vector<Telt>& LGen, const Tidx& n)
  {
#ifdef DEBUG_STABCHAINMAIN
    std::cerr << "CPP Beginning of MinimalStabChain\n";
#endif
    StabChainOptions<Tint,Tidx> options = GetStandardOptions<Tint,Tidx>(n);
    options.base = {};
#ifdef DEBUG_STABCHAINMAIN
    std::cerr << "CPP Before StabChainOp_listgen\n";
#endif
    S = StabChainOp_listgen<Telt,Tidx_label,Tint>(LGen, options);
    UnbindCycles(S);
    size_tint = Order<Telt,Tidx_label,Tint>(S);
  }
  Group(const Tidx& n) : Group({}, n)
  {
  }
  Group() : Group(0)
  {
  }
  // Basic getters
  std::vector<Telt> GeneratorsOfGroup() const
  {
    return Kernel_GeneratorsOfGroup(S);
  }
  Tint size() const
  {
    return size_tint;
  }
  std::map<Tidx, int> factor_size() const
  {
    return FactorsSizeStabChain(S);
  }
  Tidx n_act() const
  {
    return S->comm->n;
  }
  // operation
  Group<Telt,Tint> GroupConjugate(const Telt& x) const
  {
    std::vector<Telt> LGen;
    Telt xInv =~x;
    for (auto & eGen : Kernel_GeneratorsOfGroup(S)) {
      Telt eGenCj = xInv * eGen * x;
      LGen.emplace_back(eGenCj);
    }
    return Group<Telt,Tint>(LGen, S->comm->n);
  }
  // Action on points or sets
  Group<Telt,Tint> Stabilizer_OnPoints(const Tidx& x) const
  {
    return Group(Kernel_Stabilizer_OnPoints<Telt,Tidx_label,Tint>(S, x));
  }
  std::pair<bool,Telt> RepresentativeAction_OnPoints(const Tidx& x1, const Tidx& x2) const
  {
    return Kernel_RepresentativeAction_OnPoints<Telt,Tidx_label,Tint>(S, x1, x2);
  }
  Group<Telt,Tint> Stabilizer_OnSets(const Face& f) const
  {
    return Group(Kernel_Stabilizer_OnSets<Telt,Tidx_label,Tint>(S, f));
  }
  std::pair<bool,Telt> RepresentativeAction_OnSets(const Face& f1, const Face& f2) const
  {
    return Kernel_RepresentativeAction_OnSets<Telt,Tidx_label,Tint>(S, f1, f2);
  }
  bool operator==(const Group& g) const
  {
    return EqualityTest(S, g.S);
  }
  Face CanonicalImage(const Face& f) const
  {
    return Kernel_CanonicalImage<Telt,Tidx_label,Tint>(S, f);
  }
  Telt rand() const
  {
    return RandomElement(Kernel_GeneratorsOfGroup(S), S->comm->n);
  }
private:
  StabChain<Telt,Tidx_label> S;
  Tint size_tint;
};



}


#endif
