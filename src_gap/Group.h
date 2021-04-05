#ifndef DEFINE_PERMUTALIB_GROUP_H
#define DEFINE_PERMUTALIB_GROUP_H


#include "StabChainMain.h"
#include "stbcbckt.h"
#include "nsi.h"


namespace permutalib {



template<typename Telt>
Telt RandomElement(std::vector<Telt> const& LGen, int const& n)
{
  int len = rand() % 100;
  size_t n_gen = LGen.size();
  Telt eElt(n);
  for (int iIter=0; iIter<len; iIter++) {
    size_t pos = rand() % n_gen;
    eElt = eElt * LGen[pos];
  }
  return eElt;
}





template<typename Telt_inp, typename Tint_inp>
struct Group {
public:
  using Telt = Telt_inp;
  using Tidx = typename Telt::Tidx;
  using Tint = Tint_inp;
  Group(StabChain<Telt> const& _S) : S(_S), size_tint(Order<Telt,Tint>(_S))
  {
  }
  Group(std::vector<Telt> const& LGen, int const& n)
  {
#ifdef DEBUG_STABCHAINMAIN
    std::cerr << "CPP Beginning of MinimalStabChain\n";
#endif
    StabChainOptions<Tint,Tidx> options = GetStandardOptions<Tint,Tidx>(n);
    options.base = {};
#ifdef DEBUG_STABCHAINMAIN
    std::cerr << "CPP Before StabChainOp_listgen\n";
#endif
    S = StabChainOp_listgen(LGen, options);
    UnbindCycles(S);
    size_tint = Order<Telt,Tint>(S);
  }
  Group(Tidx const& n) : Group({}, n)
  {
  }
  Group() : Group(0)
  {
  }
  Group<Telt,Tint> Stabilizer_OnPoints(int const& x) const
  {
    return Group(Kernel_Stabilizer_OnPoints<Telt,Tint>(S, x));
  }
  std::pair<bool,Telt> RepresentativeAction_OnPoints(int const& x1, int const& x2) const
  {
    return Kernel_RepresentativeAction_OnPoints<Telt,Tint>(S, x1, x2);
  }
  Group<Telt,Tint> Stabilizer_OnSets(Face const& f) const
  {
    return Group(Kernel_Stabilizer_OnSets<Telt,Tint>(S, f));
  }
  std::pair<bool,Telt> RepresentativeAction_OnSets(Face const& f1, Face const& f2) const
  {
    return Kernel_RepresentativeAction_OnSets<Telt,Tint>(S, f1, f2);
  }
  std::vector<Telt> GeneratorsOfGroup() const
  {
    return Kernel_GeneratorsOfGroup(S);
  }
  Face CanonicalImage(Face const& f) const
  {
    return Kernel_CanonicalImage<Telt,Tint>(S, f);
  }
  Tint size() const
  {
    return size_tint;
  }
  Tidx n_act() const
  {
    return S->comm->n;
  }
  Telt rand() const
  {
    return RandomElement(Kernel_GeneratorsOfGroup(S), S->comm->n);
  }
private:
  StabChain<Telt> S;
  Tint size_tint;
};



}


#endif
