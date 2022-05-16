// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GAP_PROPERTIES_H_
#define SRC_GAP_PROPERTIES_H_

#include "BlockSystem.h"
#include "StabChain.h"
#include <map>
#include <vector>

namespace permutalib {

template <typename Telt>
bool Kernel_IsCommutativeGenerators(
    const std::vector<Telt> &LGen,
    const std::vector<typename Telt::Tidx> &bas) {
  using Tidx = typename Telt::Tidx;
  size_t len = LGen.size();
  for (size_t i = 0; i < len; i++) {
    for (size_t j = i + 1; j < len; j++) {
      for (auto &x : bas) {
        Tidx img1 = PowAct(PowAct(x, LGen[i]), LGen[j]);
        Tidx img2 = PowAct(PowAct(x, LGen[j]), LGen[i]);
        if (img1 != img2)
          return false;
      }
    }
  }
  return true;
}

template <typename Telt, typename Tidx_label>
bool Kernel_IsCommutative(const StabChain<Telt, Tidx_label> &S) {
  using Tidx = typename Telt::Tidx;
  std::vector<Tidx> bas = BaseStabChain(S);
  std::vector<Telt> LGen = Kernel_GeneratorsOfGroup(S);
  return Kernel_IsCommutativeGenerators(LGen, bas);
}

template <typename Telt, typename Tidx_label>
bool Kernel_IsTransitive(const StabChain<Telt, Tidx_label> &S) {
  using Tidx = typename Telt::Tidx;
  Tidx len = Tidx(S->orbit.size());
  Tidx n = S->comm->n;
  return len == n;
}

template <typename Telt, typename Tidx_label>
bool Kernel_IsPrimitive(const StabChain<Telt, Tidx_label> &S) {
  if (!Kernel_IsTransitive(S)) {
    return false;
  }
  using Tidx = typename Telt::Tidx;
  std::vector<Telt> acts = Kernel_GeneratorsOfGroup(S);
  Tidx n = S->comm->n;
  std::vector<std::vector<Tidx>> blocks = Blocks(acts, n);
  return blocks.size() == 1;
}

template <typename K, typename V>
bool EqualityMap(const std::map<K, V> &map1, const std::map<K, V> &map2) {
  if (map1.size() != map2.size())
    return false;
  for (auto &kv : map1) {
    if (map2.count(kv.first) != 1)
      return false;
    if (kv.second != map2.at(kv.first))
      return false;
  }
  return true;
}

// Adapted from GAP. Quite inefficient
template <typename Telt, typename Tidx_label, typename Tint>
bool Kernel_IsCyclic(const StabChain<Telt, Tidx_label> &S) {
  using Tidx = typename Telt::Tidx;
  std::vector<Tidx> bas = BaseStabChain(S);
  std::vector<Telt> acts = Kernel_GeneratorsOfGroup(S);
  if (!Kernel_IsCommutativeGenerators(acts, bas))
    return false;
  Telt id = S->comm->identity;
  std::map<Tidx, int> LFact = FactorsSizeStabChain(S);
  for (auto &kv : LFact) {
    std::vector<Telt> GenPow;
    for (auto &eGen : acts)
      GenPow.push_back(ElementPower(eGen, kv.first));
    StabChainOptions<Tint, Telt> options = GetStandardOptions<Tint, Telt>(id);
    StabChain<Telt, Tidx_label> Spow =
        StabChainOp_listgen<Telt, Tidx_label, Tint>(GenPow, options);
    std::map<Tidx, int> LFactPow = FactorsSizeStabChain(Spow);
    std::map<Tidx, int> Quot1 = QuotientMapMultiplicity(LFact, LFactPow);
    std::map<Tidx, int> Quot2;
    Quot2[kv.first] = 1;
    if (!EqualityMap(Quot1, Quot2))
      return false;
  }
  return true;
}

// clang-format off
}  // namespace permutalib
#endif  // SRC_GAP_PROPERTIES_H_
// clang-format on
