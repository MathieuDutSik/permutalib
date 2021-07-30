#ifndef PERMUTALIB_INCLUDE_NORMAL_STRUCTURE_H
#define PERMUTALIB_INCLUDE_NORMAL_STRUCTURE_H

#include "StabChain.h"


namespace permutalib {


template<typename Telt, typename Tidx_label>
bool Kernel_IsNormalSubgroup(const StabChain<Telt,Tidx_label>& G, const StabChain<Telt,Tidx_label>& U)
{
  std::vector<Telt> gens_G = Kernel_GeneratorsOfGroup(G);
  std::vector<Telt> gens_U = Kernel_GeneratorsOfGroup(U);

  for (auto & eGen_G : gens_G) {
    for (auto & eGen_U : gens_U) {
      Telt eConj = Conjugation(eGen_U, eGen_G);
      if (!IsElementInStabChain(U, eConj))
        return false;
    }
  }
  return true;
}



template<typename Telt, typename Tidx_label>
StabChain<Telt,Tidx_label> Kernel_NormalClosure(const StabChain<Telt,Tidx_label>& G, const StabChain<Telt,Tidx_label>& H)
{
  std::vector<Telt> gens_G = Kernel_GeneratorsOfGroup(G);


  StabChain<Telt,Tidx_label> Hret = H;

  auto test=[&]() -> std::optional<Telt> {
    std::vector<Telt> gens_Hret = Kernel_GeneratorsOfGroup(Hret);
    for (auto & eGen_G : gens_G) {
      for (auto & eGen_Hret : gens_Hret) {
        Telt eConj = Conjugation(eGen_Hret, eGen_G);
        if (!IsElementInStabChain(Hret, eConj))
          return eConj;
      }
    }
    return {};
  };
  while(true) {
    std::optional<Telt> ret = test();
    if (!ret)
      break;
    ClosureGroup(Hret, *ret);
  }
  return Hret;
}


template<typename Telt, typename Tidx_label>
StabChain<Telt,Tidx_label> Kernel_DerivedSubgroup(const StabChain<Telt,Tidx_label>& G)
{
  using Tidx = typename Telt::Tidx;
  Tidx n = G->comm->n;
  StabChain<Telt,Tidx_label> S = EmptyStabChain<Telt,Tidx_label>(n);
  std::vector<Telt> LGen = Kernel_GeneratorsOfGroup(S);
  size_t n_gen = LGen.size();
  for (size_t i=0; i<n_gen; i++) {
    const Telt& g1 = LGen[i];
    for (size_t j=i+1; j<n_gen; j++) {
      const Telt& g2 = LGen[j];
      Telt comm = Conjugation(g1, g2) * Inverse(g1);
      ClosureGroup(S, comm);
    }
  }
  return Kernel_NormalClosure(G, S);
}





}

#endif
