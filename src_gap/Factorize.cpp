#include <iostream>
#include "factorize.h"
#include <stdio.h>

int main(int argc, char *argv[])
{
  if (argc != 2) {
    std::cerr << "Program is used as\n";
    std::cerr << "Factorize [nb]\n";
    return 1;
  }
  using T = long;
  T val;
  sscanf(argv[1], "%ld", &val);
  std::cerr << "val=" << val << "\n";
  std::vector<T> V = permutalib::factorize(val);
  std::cerr << "V =";
  for (auto & eVal : V)
    std::cerr << " " << eVal;
  std::cerr << "\n";
}
