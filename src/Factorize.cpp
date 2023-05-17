// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "GapPrint.h"
#include "factorize.h"
#include <iostream>
#include <stdio.h>

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "Program is used as\n";
    std::cerr << "Factorize [nb]\n";
    return 1;
  }
  using T = int64_t;
  T val = permutalib::ReadScalar<T>(argv[1]);
  std::cerr << "val=" << val << "\n";
  std::vector<T> V = permutalib::factorize(val);
  std::cerr << "V =";
  for (auto &eVal : V)
    std::cerr << " " << eVal;
  std::cerr << "\n";
}
