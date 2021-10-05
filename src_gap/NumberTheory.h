#ifndef INCLUDE_PERMUTALIB_NUMBER_THEORY
#define INCLUDE_PERMUTALIB_NUMBER_THEORY

#include "gmpxx.h"


namespace permutalib {

  template <typename T>
  struct is_mpz_class {
    static const bool value = false;
  };

  template <>
  struct is_mpz_class<mpz_class> {
    static const bool value = true;
  };

  template<typename Tint, typename Tidx>
  inline typename std::enable_if<(not is_mpz_class<Tint>::value),Tint>::type ComputePower(const Tint& x, const Tidx& expo)
  {
    Tint ret = 1;
    for (Tidx i=0; i<expo; i++)
      ret *= x;
    return ret;
  }

  template<typename Tint, typename Tidx>
  inline typename std::enable_if<is_mpz_class<Tint>::value,Tint>::type ComputePower(const Tint& x, const Tidx& expo)
  {
    Tint ret = 1;
    for (Tidx i=0; i<expo; i++)
      ret *= x;
    return ret;
    /*
    mpz_class ret;
    unsigned long int expo_uli = expo;
    mpz_pow_ui(ret.get_mpz_t(), x.get_mpz_t(), expo_uli);
    return ret;
    */
  }
}






#endif
