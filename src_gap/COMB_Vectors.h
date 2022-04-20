#ifndef SRC_GAP_COMB_VECTORS_H_
#define SRC_GAP_COMB_VECTORS_H_

#include <algorithm>
#include <limits>
#include <map>
#include <set>
#include <unordered_set>
#include <vector>

namespace permutalib {

/*
template<typename T>
int PositionVect(std::vector<T> const& V, T const& eVal)
{
  size_t len=V.size();
  for (size_t i=0; i<len; i++)
    if (V[i] == eVal)
      return i;
  return -1;
}
*/

template <typename T, typename Tui>
Tui PositionVect_ui(std::vector<T> const &V, T const &eVal) {
  Tui len = Tui(V.size());
  for (Tui i = 0; i < len; i++)
    if (V[i] == eVal)
      return i;
  return std::numeric_limits<Tui>::max();
}

template <typename T> T VectorMax(std::vector<T> const &eVect) {
  T eMax = eVect[0];
  for (T eVal : eVect)
    if (eVal > eMax)
      eMax = eVal;
  return eMax;
}

template <typename T>
std::vector<T> ConcatenationTwo(std::vector<T> const &L1,
                                std::vector<T> const &L2) {
  std::vector<T> ret = L1;
  ret.insert(ret.end(), L2.begin(), L2.end());
  return ret;
}

template <typename T> std::vector<T> Concatenation(std::vector<T> const &L) {
  return L;
}

template <typename T, typename... Args>
std::vector<T> Concatenation(std::vector<T> const &first, Args... args) {
  std::vector<T> ret = first;
  std::vector<T> part = Concatenation(args...);
  ret.insert(ret.end(), part.begin(), part.end());
  return ret;
}

template <typename T> struct PreAllocatedVector {
  PreAllocatedVector(size_t n) {
    V = std::vector<T>(n);
    pos = 0;
  }
  size_t size() const { return pos; }
  void push_back(T val) {
    V[pos] = val;
    pos++;
  }
  void clear() { pos = 0; }
  T &operator[](size_t idx) { return V[idx]; }
  T const &operator[](size_t idx) const { return V[idx]; }

private:
  std::vector<T> V;
  size_t pos;
};

template <typename T, typename Tidx>
std::map<T, Tidx> Collected(const std::vector<T> &eVect) {
  std::map<T, Tidx> map;
  for (auto &eVal : eVect)
    map[eVal] += 1;
  return map;
}

template <typename T, class UnaryPredicate>
bool PML_ForAll(std::vector<T> const &V, UnaryPredicate const &f) {
  for (auto &eVal : V)
    if (!f(eVal))
      return false;
  return true;
}

template <typename T, class UnaryPredicate>
std::vector<T> PML_Filtered(std::vector<T> const &V, UnaryPredicate const &f) {
  std::vector<T> LRet;
  for (auto &eVal : V)
    if (f(eVal))
      LRet.push_back(eVal);
  return LRet;
}

template <typename Tidx, typename T, class UnaryPredicate>
Tidx PositionProperty(std::vector<T> const &V, UnaryPredicate const &f) {
  Tidx len = Tidx(V.size());
  for (Tidx i = 0; i < len; i++)
    if (f(V[i]))
      return i;
  return std::numeric_limits<Tidx>::max();
}

template <typename T, class UnaryPredicate>
std::vector<T> PML_ListT(std::vector<T> const &V, UnaryPredicate const &f) {
  size_t len = V.size();
  std::vector<T> retV(len);
  for (size_t i = 0; i < len; i++)
    retV[i] = f(V[i]);
  return retV;
}

// The difference V1 - V2
template <typename T>
std::vector<T> DifferenceVect(std::vector<T> const &V1,
                              std::vector<T> const &V2) {
  std::unordered_set<T> eSet2;
  for (auto &eVal : V2)
    eSet2.insert(eVal);
  std::vector<T> eV;
  for (auto &eVal : V1) {
    typename std::unordered_set<T>::iterator iter = eSet2.find(eVal);
    if (iter == eSet2.end())
      eV.push_back(eVal);
  }
  return eV;
}

template <typename T>
bool IsSubset(std::vector<T> const &S1, std::vector<T> const &S2) {
  std::unordered_set<T> eSet;
  for (auto &eVal : S1)
    eSet.insert(eVal);
  for (auto &eVal : S2)
    if (eSet.count(eVal) == 0)
      return false;
  return true;
}

template <typename T> std::vector<T> VectorAsSet(std::vector<T> const &V) {
  std::set<T> eSet;
  for (auto &eVal : V)
    eSet.insert(eVal);
  std::vector<T> eV;
  for (auto &eVal : eSet)
    eV.push_back(eVal);
  return eV;
}

template <typename T> std::vector<T> SortVector(std::vector<T> const &f) {
  std::vector<T> RetF = f;
  sort(RetF.begin(), RetF.end(),
       [](T const &x, T const &y) -> bool { return x < y; });
  return RetF;
}

} // namespace permutalib
#endif  // SRC_GAP_COMB_VECTORS_H_
