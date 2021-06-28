#ifndef PERMUTALIB_INCLUDE_FACTORIZE_H
#define PERMUTALIB_INCLUDE_FACTORIZE_H


#include <vector>

namespace permutalib {


  template<typename T>
  T gcd(T a, T b)
  {
    T remainder;
    while (b != 0) {
      remainder = a % b;
      a = b;
      b = remainder;
    }
    return a;
  }

  template<typename T>
  std::pair<bool, T> rho_pollard_factorize(T const& number)
  {
    T count;
    T x_fixed = 2, x = 2, size = 2, factor, diff;
    do {
      count = size;
      do {
        x = (x * x + 1) % number;
        diff = x - x_fixed;
        if (diff < 0)
          diff = -diff;
        factor = gcd(diff, number);
      } while (--count && (factor == 1));
      size *= 2;
      x_fixed = x;
    } while (factor == 1);
    if (factor == number) {
      return {false, -1};
    } else {
      return {true, factor};
    }
  }


  template<typename T>
  std::vector<T> successive_division_factorize(T const& N)
  {
    T pos = 2;
    while(true) {
      T res = N % pos;
      if (res == 0) {
        T quot = N / pos;
        if (quot > 1) {
          std::vector<T> eVect = successive_division_factorize(quot);
          eVect.push_back(pos);
          return eVect;
        }
        return {pos};
      }
      pos++;
      if (pos * pos > N)
        break;
    }
    return {N};
  }


  template<typename T>
  std::vector<T> factorize(T const& N)
  {
    std::pair<bool, T> epair = rho_pollard_factorize(N);
    if (epair.first) {
      T fact1 = epair.second;
      T fact2 = N / fact1;
      std::vector<T> ListPrime;
      std::vector<T> V1 = factorize(fact1);
      std::vector<T> V2 = factorize(fact2);
      ListPrime.insert(ListPrime.end(), V1.begin(), V1.end());
      ListPrime.insert(ListPrime.end(), V2.begin(), V2.end());
      return ListPrime;
    } else {
      return successive_division_factorize(N);
    }
  }


}
#endif
