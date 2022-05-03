// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GAP_NUMBERTHEORY_H_
#define SRC_GAP_NUMBERTHEORY_H_

#include "gmpxx.h"

namespace permutalib {

template <typename T> struct is_mpz_class { static const bool value = false; };

template <> struct is_mpz_class<mpz_class> { static const bool value = true; };

template <typename Tint, typename Tidx>
inline typename std::enable_if<(not is_mpz_class<Tint>::value), Tint>::type
ComputePower(const Tint &x, const Tidx &expo) {
  Tint ret = 1;
  for (Tidx i = 0; i < expo; i++)
    ret *= x;
  return ret;
}

template <typename Tint, typename Tidx>
inline typename std::enable_if<is_mpz_class<Tint>::value, Tint>::type
ComputePower(const Tint &x, const Tidx &expo) {
  Tint ret = 1;
  for (Tidx i = 0; i < expo; i++)
    ret *= x;
  return ret;
  /*
  mpz_class ret;
  unsigned long int expo_uli = expo;
  mpz_pow_ui(ret.get_mpz_t(), x.get_mpz_t(), expo_uli);
  return ret;
  */
}

// clang-format off
}  // namespace permutalib
#endif  // SRC_GAP_NUMBERTHEORY_H_
// clang-format on
