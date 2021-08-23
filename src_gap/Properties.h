#ifndef DEFINE_PERMUTALIB_PROPERTIES_H
#define DEFINE_PERMUTALIB_PROPERTIES_H



#include "StabChain.h"
#include "BlockSystem.h"

namespace permutalib {



template<typename Telt>
bool Kernel_IsCommutativeGenerators(const std::vector<Telt>& LGen)
{
  size_t len = LGen.size();
  for (size_t i=0; i<len; i++) {
    for (size_t j=i+1; j<len; j++) {
      Telt prod1 = LGen[i] * LGen[j];
      Telt prod2 = LGen[j] * LGen[i];
      if (prod1 != prod2)
        return false;
    }
  }
  return true;
}

template<typename Telt, typename Tidx_label>
bool Kernel_IsCommutative(const StabChain<Telt,Tidx_label>& S)
{
  std::vector<Telt> LGen = Kernel_GeneratorsOfGroup(S);
  return Kernel_IsCommutativeGenerators(LGen);
}

template<typename Telt, typename Tidx_label>
bool Kernel_IsTransitive(const StabChain<Telt,Tidx_label>& S)
{
  using Tidx = typename Telt::Tidx;
  Tidx len=Tidx(S->orbit.size());
  Tidx n = S->comm->n;
  return len == n;
}


template<typename Telt, typename Tidx_label>
bool Kernel_IsPrimitive(const StabChain<Telt,Tidx_label>& S)
{
  if (!Kernel_IsTransitive(S)) {
    return false;
  }
  using Tidx = typename Telt::Tidx;
  std::vector<Telt> acts = Kernel_GeneratorsOfGroup(S);
  Tidx n = S->comm->n;
  std::vector<std::vector<Tidx>> blocks = Blocks(acts, n);
  return blocks.size() == 1;
}



template<typename K, typename V>
bool EqualityMap(const std::map<K,V>& map1, const std::map<K,V>& map2) {
  if (map1.size() != map2.size())
    return false;
  for (auto & kv : map1) {
    if (map2.count(kv.first) != 1)
      return false;
    if (kv.second != map2.at(kv.first))
      return false;
  }
  return true;
}




// Adapted from GAP. Quite inefficient
template<typename Telt, typename Tidx_label, typename Tint>
bool Kernel_IsCyclic(const StabChain<Telt,Tidx_label>& S)
{
  using Tidx=typename Telt::Tidx;
  std::vector<Telt> acts = Kernel_GeneratorsOfGroup(S);
  if (!Kernel_IsCommutativeGenerators(acts))
    return false;
  Tidx n = S->comm->n;
  std::map<Tidx,int> LFact = FactorsSizeStabChain(S);
  for (auto & kv : LFact) {
    std::vector<Telt> GenPow;
    for (auto & eGen : acts)
      GenPow.push_back(ElementPower(eGen, kv.first));
    StabChainOptions<Tint,Tidx> options = GetStandardOptions<Tint,Tidx>(n);
#ifdef DEBUG_STABCHAINMAIN
    std::cerr << "CPP Before StabChainOp_listgen\n";
#endif
    StabChain<Telt,Tidx_label> Spow = StabChainOp_listgen<Telt,Tidx_label,Tint>(GenPow, options);
    std::map<Tidx,int> LFactPow = FactorsSizeStabChain(Spow);
    std::map<Tidx,int> Quot1 = QuotientMapMultiplicity(LFact, LFactPow);
    std::map<Tidx,int> Quot2;
    Quot2[kv.first] = 1;
    if (!EqualityMap(Quot1, Quot2))
      return false;
  }
  return true;
}


/*
  We follow here a different approach to GAP.
  Testing cyclicity should not be a computationally intensive
  operation.
*/
/*
template<typename Telt, typename Tidx_label>
std::optional<Telt> Kernel_IsCyclic(const StabChain<Telt,Tidx_label>& S_in)
{
  using Tidx=typename Telt::Tidx;
  Tidx n = S_in->comm->n;
  Tidx order = 1;
  auto transposition=[&](const Tidx& a, const Tidx& b) -> Telt {
    std::vector<Tidx> V(n);
    for (Tidx i=0; i<n; i++)
      V[i] = i;
    V[a] = b;
    V[b] = a;
    return Telt(V);
  };
  auto get_spanning_generator=[&](const StabChain<Telt,Tidx_label>& S) -> std::optional<Telt> {
    Tidx len = Tidx(S->orbit.size());
    Tidx bpt = S->orbit[0];
    for (Tidx i=1; i<len; i++) {
      const Tidx img = S->orbit[i];
      Tidx img_work = img;
      Telt g = S->comm->identity;
      while(true) {
        if (img_work == bpt)
          break;
        Tidx_label idx = S->transversal[img];
        g *= S->comm->labels[idx];
        img_work = PowAct(img, g);
      }
      Tidx ord = OrderElement(g);
      Tidx res = ord % len;
      if (res == 0) { // Finding an element of order len
        Face f(n);
        Tidx work = bpt;
        Tidx siz_match = 0;
        while(true) {
          f[workpt] = 1;
          siz_match++;
          workpt = PowAct(workpt, g);
          if (f[workpt] == 1)
            break;
        }
        if (siz_match == len)
          return g;
      }
    }
    return {};
  };
  StabChain<Telt,Tidx_label> Swork = S_in;
  Face status(n);
  while(true) {
    std::optional<Telt> test = get_spanning_generator(Swork);
    if (!test) {
      return {};
    }
    
  }
}
*/


}



#endif
