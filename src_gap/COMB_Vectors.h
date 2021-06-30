#ifndef DEFINE_PERMUTALIB_COMB_VECTORS_H
#define DEFINE_PERMUTALIB_COMB_VECTORS_H


#include <vector>
#include <set>
#include <unordered_set>

namespace permutalib {

template<typename T>
int PositionVect(std::vector<T> const& V, T const& eVal)
{
  size_t len=V.size();
  for (size_t i=0; i<len; i++)
    if (V[i] == eVal)
      return i;
  return -1;
}

template<typename T>
size_t PositionVect_sz(std::vector<T> const& V, T const& eVal)
{
  size_t len=V.size();
  for (size_t i=0; i<len; i++)
    if (V[i] == eVal)
      return i;
  return std::numeric_limits<size_t>::max();
}

template<typename T>
T VectorMax(std::vector<T> const& eVect)
{
  T eMax=eVect[0];
  for (T eVal : eVect)
    if (eVal > eMax)
      eMax=eVal;
  return eMax;
}


template<typename T>
std::vector<T> ConcatenationTwo(std::vector<T> const& L1, std::vector<T> const& L2)
{
  std::vector<T> ret = L1;
  for (auto & eVal : L2)
    ret.push_back(eVal);
  return ret;
}


template<typename T>
std::vector<T> Concatenation(std::vector<T> const& L)
{
  return L;
}


template<typename T, typename... Args>
std::vector<T> Concatenation(std::vector<T> const& first, Args... args)
{
  std::vector<T> ret = first;
  for (auto & eVal : Concatenation(args...))
    ret.push_back(eVal);
  return ret;
}


template<typename T>
struct CollectedResult {
  std::vector<T> LVal;
  std::vector<int> LMult;
};




template<typename T>
struct PreAllocatedVector {
  PreAllocatedVector(size_t n)
  {
    V = std::vector<T>(n);
    pos = 0;
  }
  size_t size() const
  {
    return pos;
  }
  void push_back(T val)
  {
    V[pos] = val;
    pos++;
  }
  void clear()
  {
    pos = 0;
  }
  T& operator[](size_t idx)
  {
    return V[idx];
  }
  T const& operator[](size_t idx) const
  {
    return V[idx];
  }

private:
  std::vector<T> V;
  size_t pos;
};




template<typename T>
CollectedResult<T> Collected(std::vector<T> const& eVect)
{
  std::set<T> SetVal;
  for (auto & eVal : eVect)
    SetVal.insert(eVal);
  std::vector<T> LVal;
  for (auto & eVal : SetVal)
    LVal.push_back(eVal);
  size_t eSize=LVal.size();
  std::vector<int> LMult(eSize,0);
  auto UpPosition=[&](T const& eVal) -> void {
    for (size_t i=0; i<eSize; i++)
      if (LVal[i] == eVal) {
        LMult[i] += 1;
        return;
      }
  };
  for (auto & eVal : eVect)
    UpPosition(eVal);
  return {std::move(LVal), std::move(LMult)};
}


template<typename T, class UnaryPredicate>
bool ForAll(std::vector<T> const& V, UnaryPredicate const& f)
{
  for (auto & eVal : V)
    if (!f(eVal))
      return false;
  return true;
}

template<typename T, class UnaryPredicate>
std::vector<T> Filtered(std::vector<T> const& V, UnaryPredicate const& f)
{
  std::vector<T> LRet;
  for (auto & eVal : V)
    if (f(eVal))
      LRet.push_back(eVal);
  return LRet;
}


template<typename Tidx, typename T, class UnaryPredicate>
Tidx PositionProperty(std::vector<T> const& V, UnaryPredicate const& f)
{
  Tidx len=Tidx(V.size());
  for (Tidx i=0; i<len; i++)
    if (f(V[i]))
      return i;
  return std::numeric_limits<Tidx>::max();
}

template<typename T, class UnaryPredicate>
std::vector<T> ListT(std::vector<T> const& V, UnaryPredicate const& f)
{
  std::vector<T> retV;
  for (auto & eVal : V)
    retV.push_back(f(eVal));
  return retV;
}


// The difference V1 - V2
template<typename T>
std::vector<T> DifferenceVect(std::vector<T> const& V1, std::vector<T> const& V2)
{
  std::unordered_set<T> eSet2;
  for (auto & eVal : V2)
    eSet2.insert(eVal);
  std::vector<T> eV;
  for (auto & eVal : V1) {
    typename std::unordered_set<T>::iterator iter=eSet2.find(eVal);
    if (iter == eSet2.end())
      eV.push_back(eVal);
  }
  return eV;
}

template<typename T>
bool IsSubset(std::vector<T> const& S1, std::vector<T> const& S2)
{
  std::unordered_set<T> eSet;
  for (auto & eVal : S1)
    eSet.insert(eVal);
  for (auto & eVal : S2)
    if (eSet.count(eVal) == 0)
      return false;
  return true;
}

template<typename T>
std::vector<T> VectorAsSet(std::vector<T> const& V)
{
  std::set<T> eSet;
  for (auto & eVal : V)
    eSet.insert(eVal);
  std::vector<T> eV;
  for (auto & eVal : eSet)
    eV.push_back(eVal);
  return eV;
}

}
#endif
