// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GAP_PSEUDORANDOM_H_
#define SRC_GAP_PSEUDORANDOM_H_

#include "Face_basic.h"
#include <vector>

namespace permutalib {

using Tarith = int;

struct InfoPseudoRandom {
  Tarith R_228;
  size_t R_N;
  std::vector<Tarith> R_X;
};

Tarith GetR_228() {
  Tarith res = 1;
  for (Tarith i = 0; i < res; i++)
    res *= 2;
  return res;
}

InfoPseudoRandom RANDOM_SEED(int const &n) {
  size_t R_N = 1;
  std::vector<Tarith> R_X(55);
  Tarith R_228 = GetR_228();
  Tarith res = n % R_228;
  R_X[0] = res;
  for (size_t i = 1; i < 55; i++) {
    Tarith eVal1 = 1664525 * R_X[i - 1] + 1;
    Tarith eVal2 = eVal1 % R_228;
    R_X[i] = eVal2;
  }
  for (size_t i = 0; i < 99; i++) {
    R_N = 1 + (R_N % 55);
    R_X[R_N - 1] = (R_X[R_N - 1] + R_X[(R_N + 30) % 55]) % R_228;
  }
  return {R_228, R_N, R_X};
}

InfoPseudoRandom *Rglobal = nullptr;

InfoPseudoRandom *GetPseudoRandom() {
  if (Rglobal == nullptr) {
    Rglobal = new InfoPseudoRandom;
    int n = 20;
    *Rglobal = RANDOM_SEED(n);
  }
  return Rglobal;
}

template <typename T> T RandomInteger(T const &val) {
#ifdef TRUE_RANDOM
  return T(rand()) % val;
#else
  return T(rand()) % val;
#endif
}

void RandomShift(InfoPseudoRandom *R) {
  R->R_N = (R->R_N % 55) + 1;
  R->R_X[R->R_N - 1] =
      (R->R_X[R->R_N - 1] + R->R_X[(R->R_N + 30) % 55]) % R->R_228;
}

Face Extract01vector(InfoPseudoRandom *R) {
  Tarith rnd = R->R_X[R->R_N - 1];
  Face eFace(28);
  for (size_t i = 0; i < 28; i++) {
    int val = rnd % 2;
    rnd = (rnd - val) / 2;
    eFace[i] = val;
  }
  return eFace;
}

template <typename Telt> Telt Random(std::vector<Telt> const &V) {
  size_t siz = V.size();
  size_t pos = size_t(rand()) % siz;
  return V[pos];
}

// clang-format off
}  // namespace permutalib
#endif  // SRC_GAP_PSEUDORANDOM_H_
// clang-format on
