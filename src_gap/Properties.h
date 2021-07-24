#ifndef DEFINE_PERMUTALIB_PROPERTIES_H
#define DEFINE_PERMUTALIB_PROPERTIES_H


#include "StabChain.h"

namespace permutalib {


template<typename Telt, typename Tidx_label>
bool Kernel_IsCommutative(const StabChain<Telt,Tidx_label>& S)
{
  std::vector<Telt> LGen = Kernel_GeneratorsOfGroup(S);
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
bool Kernel_IsTransitive(const StabChain<Telt,Tidx_label>& S)
{
  using Tidx = typename Telt::Tidx;
  Tidx len=Tidx(S->orbit.size());
  Tidx n = S->comm->n;
  return len == n;
}


}



#endif
