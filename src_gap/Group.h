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
  using Tint = Tint_inp;
  Group(StabChain<Telt_inp> const& _S) : S(_S), size_tint(Order<Telt_inp,Tint_inp>(_S))
  {
  }
  Group(std::vector<Telt_inp> const& LGen, int const& n)
  {
#ifdef DEBUG_STABCHAINMAIN
    std::cerr << "CPP Beginning of MinimalStabChain\n";
#endif
    StabChainOptions<Tint_inp> options = GetStandardOptions<Tint_inp>(n);
    options.base = {};
#ifdef DEBUG_STABCHAINMAIN
    std::cerr << "CPP Before StabChainOp_listgen\n";
#endif
    S = StabChainOp_listgen(LGen, options);
    UnbindCycles(S);
    size_tint = Order<Telt_inp,Tint_inp>(S);
  }
  Group(int const& n) : Group({}, n)
  {
  }
  Group() : Group(0)
  {
  }
  Group<Telt_inp,Tint_inp> Stabilizer_OnPoints(int const& x) const
  {
    return Group(Kernel_Stabilizer_OnPoints<Telt_inp,Tint_inp>(S, x));
  }
  std::pair<bool,Telt_inp> RepresentativeAction_OnPoints(int const& x1, int const& x2) const
  {
    return Kernel_RepresentativeAction_OnPoints<Telt_inp,Tint_inp>(S, x1, x2);
  }
  Group<Telt_inp,Tint_inp> Stabilizer_OnSets(Face const& f) const
  {
    return Group(Kernel_Stabilizer_OnSets<Telt_inp,Tint_inp>(S, f));
  }
  std::pair<bool,Telt_inp> RepresentativeAction_OnSets(Face const& f1, Face const& f2) const
  {
    return Kernel_RepresentativeAction_OnSets<Telt_inp,Tint_inp>(S, f1, f2);
  }
  std::vector<Telt_inp> GeneratorsOfGroup() const
  {
    return Kernel_GeneratorsOfGroup(S);
  }
  Face CanonicalImage(Face const& f) const
  {
    return Kernel_CanonicalImage<Telt_inp,Tint_inp>(S, f);
  }
  Tint_inp size() const
  {
    return size_tint;
  }
  int n_act() const
  {
    return S->comm->n;
  }
  Telt_inp rand() const
  {
    return RandomElement(Kernel_GeneratorsOfGroup(S), S->comm->n);
  }
private:
  StabChain<Telt_inp> S;
  Tint_inp size_tint;
};



}


#endif
