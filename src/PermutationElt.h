// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GAP_PERMUTATIONELT_H_
#define SRC_GAP_PERMUTATIONELT_H_

#include "Face_basic.h"
#include "exception.h"
#include "hash_fct.h"
#include <iostream>
#include <limits>
#include <stdlib.h>
#include <string>
#include <utility>
#include <vector>

namespace permutalib {

  /*
    This is for containing words.
   */

template<typename T>
void SimplifySequence(std::vector<T> & V)
{
  size_t len = V.size();
  std::vector<uint8_t> status(len,1);
  size_t miss_val = std::numeric_limits<size_t>::max();
  std::vector<size_t> next(len,1);
  std::vector<size_t> prev(len,1);
  for (size_t i=0; i<len-1; i++)
    next[i] = i+1;
  next[len-1] = miss_val;
  for (size_t i=1; i<len; i++)
    prev[i] = i-1;
  prev[0] = miss_val;
  size_t n_ent = len;
  while(true) {
    size_t n_change = 0;
    for (size_t i=0; i<len-1; i++) {
      if (status[i] == 1) {
        size_t iNext = next[i];
        if (iNext != miss_val) {
          if (V[i] == -V[iNext]) {
            n_change++;
            n_ent -= 2;
            status[i] = 0;
            status[iNext] = 0;
            size_t iP = prev[i];
            size_t iN = next[iNext];
            next[i] = miss_val;
            prev[i] = miss_val;
            next[iNext] = miss_val;
            prev[iNext] = miss_val;
            //
            if (iP != miss_val) {
              // We can assign, but note that the iN can be missing
              next[iP] = iN;
            }
            if (iN != miss_val) {
              // We can assign, but note that the iP can be missing
              prev[iN] = iP;
            }
          }
        }
      }
    }
    if (n_change == 0)
      break;
  }
  if (n_ent == len)
    return;
  size_t pos = 0;
  for (size_t i=0; i<len; i++) {
    if (status[i] == 1) {
      V[pos] = V[i];
      pos++;
    }
  }
  V.resize(n_ent);
}




#ifdef DEBUG_PERMUTATION_ELT
void PrintListIdx(std::string const& mesg, std::vector<int64_t> const& ListIdx) {
  size_t len = ListIdx.size();
  std::cerr << mesg << " " << len << ":[";
  bool IsFirst = true;
  for (auto & eval : ListIdx) {
    if (!IsFirst)
      std::cerr << ",";
    IsFirst = false;
    std::cerr << eval;
  }
  std::cerr << "]\n";
}
#endif

template<typename T>
bool IsSimplifiable(std::vector<T> const& V)
{
  size_t len = V.size();
  for (size_t i=1; i<len; i++) {
    if (V[i-1] == -V[i])
      return true;
  }
  return false;
}

template<bool always_equal>
struct SequenceType {
  SequenceType() : ListIdx() {
#ifdef DEBUG_PERMUTATION_ELT
    PrintListIdx("default constructor", ListIdx);
#endif
  }
  SequenceType(std::vector<int64_t> &&v) {
    ListIdx = v;
#ifdef DEBUG_PERMUTATION_ELT
    PrintListIdx("constructor 1", ListIdx);
#endif
  }
  SequenceType(std::vector<int64_t> const &v) {
    ListIdx = v;
#ifdef DEBUG_PERMUTATION_ELT
    PrintListIdx("constructor 2", ListIdx);
#endif
  }
  SequenceType(SequenceType<always_equal> const &seq) {
    ListIdx = seq.ListIdx;
#ifdef DEBUG_PERMUTATION_ELT
    PrintListIdx("constructor 3", ListIdx);
#endif
  }
  SequenceType(SequenceType<always_equal> &&seq) {
    ListIdx = std::move(seq.ListIdx);
#ifdef DEBUG_PERMUTATION_ELT
    PrintListIdx("constructor 4", ListIdx);
#endif
  }
  //
  // Copy operator
  //
  SequenceType operator=(SequenceType<always_equal> const &seq) {
    ListIdx = seq.ListIdx;
    return *this;
  }
  SequenceType operator=(SequenceType<always_equal> &&seq) {
    ListIdx = std::move(seq.ListIdx);
    return *this;
  }
  //
  // The destructor
  //
  ~SequenceType() {}
  //
  // Other stuff
  //
  bool isIdentity() const {
    return ListIdx.size() == 0;
  }
  const std::vector<int64_t>& getVect() const {
    return ListIdx;
  }
  std::vector<int64_t>& getVect() {
    return ListIdx;
  }
  size_t complexity() const {
    return ListIdx.size();
  }
private:
  std::vector<int64_t> ListIdx;
};

template<bool always_equal>
SequenceType<always_equal> operator*(SequenceType<always_equal> const& v1, SequenceType<always_equal> const& v2) {
  std::vector<int64_t> ListIdx1 = v1.getVect();
  const std::vector<int64_t> &ListIdx2 = v2.getVect();
  ListIdx1.insert(ListIdx1.end(), ListIdx2.begin(), ListIdx2.end());
  SimplifySequence(ListIdx1);
#ifdef DEBUG_PERMUTATION_ELT
  PrintListIdx("operator*", ListIdx1);
#endif
  return SequenceType<always_equal>(std::move(ListIdx1));
}

template<bool always_equal>
bool operator<(SequenceType<always_equal> const& v1, SequenceType<always_equal> const& v2) {
  std::vector<int64_t> const& ListIdx1 = v1.getVect();
  std::vector<int64_t> const& ListIdx2 = v2.getVect();
  size_t len1 = ListIdx1.size();
  size_t len2 = ListIdx2.size();
  if (len1 != len2) {
    // Shorter entries are preferred.
    return len1 < len2;
  }
  for (size_t i=0; i<len1; i++) {
    int64_t val1 = ListIdx1[i];
    int64_t val2 = ListIdx2[i];
    if (val1 != val2) {
      return val1 < val2;
    }
  }
  return false;
}

template<bool always_equal>
void operator*=(SequenceType<always_equal> &v1,
                SequenceType<always_equal> const &v2) {
  std::vector<int64_t> &ListIdx1 = v1.getVect();
  const std::vector<int64_t> &ListIdx2 = v2.getVect();
  ListIdx1.insert(ListIdx1.end(), ListIdx2.begin(), ListIdx2.end());
  SimplifySequence(ListIdx1);
#ifdef DEBUG_PERMUTATION_ELT
  PrintListIdx("operator*=", ListIdx1);
#endif
}


template<bool always_equal>
SequenceType<always_equal> Conjugation(SequenceType<always_equal> const &v1,
                                       SequenceType<always_equal> const &v2) {
  const std::vector<int64_t> &ListIdx1 = v1.getVect();
  const std::vector<int64_t> &ListIdx2 = v2.getVect();
  size_t siz1 = ListIdx1.size();
  size_t siz2 = ListIdx2.size();
  std::vector<int64_t> ListIdx(siz2 + siz1 + siz2);
  for (size_t i=0; i<siz2; i++)
    ListIdx[i] = - ListIdx2[siz2 - 1 - i];
  for (size_t i=0; i<siz1; i++)
    ListIdx[siz2 + i] = ListIdx1[i];
  for (size_t i=0; i<siz2; i++)
    ListIdx[siz2 + siz1 + i] = ListIdx2[i];
  SimplifySequence(ListIdx);
#ifdef DEBUG_PERMUTATION_ELT
  PrintListIdx("Conjugation", ListIdx);
#endif
  return SequenceType<always_equal>(std::move(ListIdx));
}


// LeftQuotient(a,b) = a^{-1} * b in the list.gi file
template<bool always_equal>
SequenceType<always_equal> LeftQuotient(SequenceType<always_equal> const &a, SequenceType<always_equal> const &b) {
  const std::vector<int64_t> &Val_A = a.getVect();
  const std::vector<int64_t> &Val_B = b.getVect();
  size_t siz_a = Val_A.size();
  size_t siz_b = Val_B.size();
  std::vector<int64_t> ListIdx(siz_a + siz_b);
  for (size_t i=0; i<siz_a; i++)
    ListIdx[i] = - Val_A[siz_a - 1 - i];
  for (size_t i=0; i<siz_b; i++)
    ListIdx[siz_a + i] = Val_B[i];
  SimplifySequence(ListIdx);
#ifdef DEBUG_PERMUTATION_ELT
  PrintListIdx("LeftQuotient", ListIdx);
#endif
  return SequenceType<always_equal>(std::move(ListIdx));
}


template<bool always_equal>
SequenceType<always_equal> operator~(SequenceType<always_equal> const &seq) {
  const std::vector<int64_t> & ListIdx = seq.getVect();
  size_t len = ListIdx.size();
  std::vector<int64_t> vret(len);
  for (size_t i=0; i<len; i++)
    vret[len - 1 - i] = - ListIdx[i];
#ifdef DEBUG_PERMUTATION_ELT
  PrintListIdx("operator~", vret);
#endif
  return SequenceType<always_equal>(std::move(vret));
}


template<bool always_equal>
SequenceType<always_equal> Inverse(SequenceType<always_equal> const &seq) {
  return ~seq;
}


template<bool always_equal>
bool operator==(SequenceType<always_equal> const &v1,
                SequenceType<always_equal> const &v2) {
  if (always_equal) {
    return true;
  } else {
    const std::vector<int64_t> & LIdx1 = v1.getVect();
    const std::vector<int64_t> & LIdx2 = v2.getVect();
    size_t siz = LIdx1.size();
    if (siz != LIdx2.size())
      return false;
    for (size_t i = 0; i < siz; i++)
      if (LIdx1.at(i) != LIdx2.at(i))
        return false;
    return true;
  }
}


template<bool always_equal>
bool operator!=(SequenceType<always_equal> const &v1,
                SequenceType<always_equal> const &v2) {
  if (always_equal) {
    return false;
  } else {
    const std::vector<int64_t> & LIdx1 = v1.getVect();
    const std::vector<int64_t> & LIdx2 = v2.getVect();
    size_t siz = LIdx1.size();
    if (siz != LIdx2.size())
      return true;
    for (size_t i = 0; i < siz; i++)
      if (LIdx1.at(i) != LIdx2.at(i))
        return true;
    return false;
  }
}


template<bool always_equal>
std::ostream &operator<<(std::ostream &os, SequenceType<always_equal> const& seq)
{
  const std::vector<int64_t>& ListIdx = seq.getVect();
  size_t len = ListIdx.size();
  os << len << ":[";
  bool IsFirst = true;
  for (auto & eval : ListIdx) {
    if (!IsFirst)
      os << ",";
    IsFirst = false;
    os << eval;
  }
  os << "]";
  return os;
}


  /*
    This is for containing pairs of Element and Permutation.
   */


template <typename Tidx> void CheckSize(size_t siz) {
  if (siz >= std::numeric_limits<Tidx>::max() - 1) {
    std::cerr << "siz=" << siz << "\n";
    std::cerr << "but std::numeric_limits<Tidx>::max() = "
              << std::numeric_limits<Tidx>::max() << "\n";
    std::cerr << "Tidx is too small for representing the vector\n";
    throw PermutalibException{1};
  }
}

template <typename Tidx_inp, typename Telt_impl> struct PermutationElt;

template <typename Tidx, typename Telt>
void NicePrint(std::string const &name,
               PermutationElt<Tidx, Telt> const &ePermElt) {
  //  Tidx siz = ePermElt.size();
  const std::vector<Tidx> &LVal = ePermElt.getListVal();
  std::cerr << "name=" << name << " : V =";
  for (auto &val : LVal)
    std::cerr << " " << val;
  std::cerr << "\n";
  //
  std::cerr << "M =";
  WriteMatrix(std::cerr, ePermElt.getElt());
}

template <typename Tidx_inp, typename Telt_impl> struct PermutationElt {
public:
  using Tidx = Tidx_inp;
  using Telt = Telt_impl;
  //
  // The constructors
  //
  PermutationElt() : siz(0), ListVal(), elt() {
  }
  PermutationElt(std::vector<Tidx> &&v, Telt &&_elt) {
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
    CheckSize<Tidx>(v.size());
#endif
    ListVal = v;
    siz = Tidx(v.size());
    elt = std::move(_elt);
  }
  PermutationElt(std::vector<Tidx> const &v, Telt const &_elt) {
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
    CheckSize<Tidx>(v.size());
#endif
    ListVal = v;
    siz = Tidx(v.size());
    elt = _elt;
  }
  PermutationElt(PermutationElt const &ePermElt) {
    siz = ePermElt.siz;
    ListVal = ePermElt.ListVal;
    elt = ePermElt.elt;
  }
  PermutationElt(PermutationElt &&ePermElt) {
    siz = ePermElt.siz;
    ListVal = std::move(ePermElt.ListVal);
    elt = std::move(ePermElt.elt);
    ePermElt.siz = 0;
  }
  //
  // Copy operator
  //
  PermutationElt<Tidx, Telt>
  operator=(PermutationElt<Tidx, Telt> const &ePermElt) {
    siz = ePermElt.siz;
    ListVal = ePermElt.ListVal;
    elt = ePermElt.elt;
    return *this;
  }
  PermutationElt<Tidx, Telt> operator=(PermutationElt<Tidx, Telt> &&ePermElt) {
    siz = ePermElt.siz;
    ListVal = std::move(ePermElt.ListVal);
    elt = std::move(ePermElt.elt);
    ePermElt.siz = 0;
    return *this;
  }
  //
  // The destructor
  //
  ~PermutationElt() {}
  //
  // The other functionalities
  //
  bool isIdentity() const {
    for (Tidx i = 0; i < siz; i++)
      if (ListVal[i] != i)
        return false;
    return true;
  }
  Tidx at(Tidx const &i) const { return ListVal[i]; }
  Tidx atRev(Tidx const &i) const {
    Tidx i1 = i, i2 = ListVal[i];
    while (true) {
      if (i2 == i)
        return i1;
      i1 = i2;
      i2 = ListVal[i2];
    }
  }
  const Tidx *getPtr() const { return ListVal.data(); }
  const std::vector<Tidx> &getListVal() const { return ListVal; }
  const Telt &getElt() const { return elt; }
  Tidx operator[](Tidx const &i) const { return ListVal[i]; }
  Tidx size() const { return siz; }
  //
public:
  // Should be private in a more classic construction
  Tidx siz;
  std::vector<Tidx> ListVal;
  Telt elt;
};

template <typename Tidx, typename Telt>
bool operator==(PermutationElt<Tidx, Telt> const &v1,
                PermutationElt<Tidx, Telt> const &v2) {
  Tidx siz = v1.size();
  if (siz != v2.size())
    return false;
  for (Tidx i = 0; i < siz; i++)
    if (v1.at(i) != v2.at(i))
      return false;
  return v1.elt == v2.elt;
}

template <typename Tidx, typename Telt>
bool operator!=(PermutationElt<Tidx, Telt> const &v1,
                PermutationElt<Tidx, Telt> const &v2) {
  Tidx siz = v1.size();
  if (siz != v2.size())
    return true;
  for (Tidx i = 0; i < siz; i++)
    if (v1.at(i) != v2.at(i))
      return true;
  return v1.elt != v2.elt;
}

template <typename Tidx, typename Telt>
bool operator<(PermutationElt<Tidx, Telt> const &v1,
               PermutationElt<Tidx, Telt> const &v2) {
  Tidx siz1 = v1.size();
  Tidx siz2 = v2.size();
  if (siz1 != siz2)
    return siz1 < siz2;
  Tidx siz = siz1;
  for (Tidx i = 0; i < siz; i++) {
    if (v1.at(i) != v2.at(i))
      return v1.at(i) < v2.at(i);
  }
  return v1.elt < v2.elt;
}

//
// Operations
//

template <typename Tidx, typename Telt>
bool CoherencyCheck(PermutationElt<Tidx, Telt> const &ePermElt) {
  Tidx siz = ePermElt.size();
  const std::vector<Tidx> &LVal = ePermElt.getListVal();
  Telt elt = ePermElt.getElt();
  bool is_ok = true;
  if (static_cast<int>(siz) == elt.rows() &&
      static_cast<int>(siz) == elt.cols()) {
    for (Tidx iLine = 0; iLine < siz; iLine++) {
      Tidx pos = LVal[iLine];
      if (elt(iLine, pos) != 1)
        is_ok = false;
    }
  }
  return is_ok;
}

template <typename Tidx, typename Telt>
PermutationElt<Tidx, Telt>
operator~(PermutationElt<Tidx, Telt> const &ePermElt) {
  //  NicePrint("operatorInverse : ePermElt", ePermElt);
  Tidx siz = ePermElt.size();
  const std::vector<Tidx> &LVal = ePermElt.getListVal();
  std::vector<Tidx> v(siz);
  for (Tidx i = 0; i < siz; i++)
    v[LVal[i]] = i;
  Telt Minv = Inverse(ePermElt.elt);
  //  NicePrint("eInv", PermutationElt<Tidx,Telt>(v, Minv));
  //  std::cerr << "operatorInverse status=" << CoherencyCheck(ePermElt) <<
  //  "\n";
  return PermutationElt<Tidx, Telt>(std::move(v), std::move(Minv));
}

// Form the product v1 * v2
template <typename Tidx, typename Telt>
PermutationElt<Tidx, Telt> operator*(PermutationElt<Tidx, Telt> const &v1,
                                     PermutationElt<Tidx, Telt> const &v2) {
  //  NicePrint("operator* : v1", v1);
  //  NicePrint("operator* : v2", v2);
  Tidx siz = v1.size();
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
  if (siz != v2.size()) {
    std::cerr << "Error in the PermutationElt product\n";
    throw PermutalibException{1};
  }
#endif
  const std::vector<Tidx> &LVal1 = v1.getListVal();
  const std::vector<Tidx> &LVal2 = v2.getListVal();
  std::vector<Tidx> vVal(siz);
  for (Tidx i = 0; i < siz; i++) {
    Tidx j = LVal1[i];
    Tidx k = LVal2[j];
    vVal[i] = k;
  }
  //  return SingleSidedPerm<Tidx>(vVal);
  Telt Mprod = v1.elt * v2.elt;
  //  NicePrint("prod", PermutationElt<Tidx,Telt>(vVal, Mprod));
  //  std::cerr << "operator* status=" <<
  //  CoherencyCheck(PermutationElt<Tidx,Telt>(vVal, Mprod)) << "\n";
  return PermutationElt<Tidx, Telt>(std::move(vVal), std::move(Mprod));
}

template <typename Tidx, typename Telt>
void operator*=(PermutationElt<Tidx, Telt> &v1,
                PermutationElt<Tidx, Telt> const &v2) {
  Tidx siz_idx = Tidx(v1.size());
  for (Tidx i = 0; i < siz_idx; i++) {
    Tidx j = v1.at(i);
    Tidx k = v2.at(j);
    v1.ListVal[i] = k;
  }
  v1.elt *= v2.elt;
}

template <typename Tidx, typename Telt>
PermutationElt<Tidx, Telt> Conjugation(PermutationElt<Tidx, Telt> const &v1,
                                       PermutationElt<Tidx, Telt> const &v2) {
  //  NicePrint("Conjugation : v1", v1);
  //  NicePrint("Conjugation : v2", v2);
  Tidx siz = v1.size();
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
  if (siz != v2.size()) {
    std::cerr << "Error in the DoubleSidedPerm conjugation\n";
    throw PermutalibException{1};
  }
#endif
  std::vector<Tidx> v(siz);
  for (Tidx i = 0; i < siz; i++) {
    Tidx j = v1[i];
    Tidx i2 = v2[i];
    Tidx j2 = v2[j];
    v[i2] = j2;
  }
  // That formula is the correct one
  Telt Mret = Inverse(v2.elt) * v1.elt * v2.elt;
  //  NicePrint("res", PermutationElt<Tidx,Telt>(v, Mret));
  //  std::cerr << "Conjugation status=" <<
  //  CoherencyCheck(PermutationElt<Tidx,Telt>(v, Mret)) << "\n";
  return PermutationElt<Tidx, Telt>(std::move(v), std::move(Mret));
}

template <typename Tidx, typename Telt>
Tidx PowAct(Tidx const &i, PermutationElt<Tidx, Telt> const &g) {
  return g.at(i);
}

template <typename Tidx, typename Telt>
Tidx SlashAct(Tidx const &i, PermutationElt<Tidx, Telt> const &g) {
  return g.atRev(i);
}

// LeftQuotient(a,b) = a^{-1} * b in the list.gi file
template <typename Tidx, typename Telt>
PermutationElt<Tidx, Telt> LeftQuotient(PermutationElt<Tidx, Telt> const &a,
                                        PermutationElt<Tidx, Telt> const &b) {
  //  NicePrint("LeftQuotient : a", a);
  //  NicePrint("LeftQuotient : b", b);
  Tidx siz = Tidx(a.size());
  const std::vector<Tidx> &Val_A = a.getListVal();
  const std::vector<Tidx> &Val_B = b.getListVal();
  std::vector<Tidx> ListVal(siz);
  for (Tidx i = 0; i < siz; i++)
    ListVal[Val_A[i]] = Val_B[i];
  Telt Mret = Inverse(a.elt) * b.elt;
  //  NicePrint("res", PermutationElt<Tidx,Telt>(ListVal, Mret));
  //  std::cerr << "LeftQuotient status=" <<
  //  CoherencyCheck(PermutationElt<Tidx,Telt>(ListVal, Mret)) << "\n";
  return PermutationElt<Tidx, Telt>(std::move(ListVal), std::move(Mret));
}

template <typename Tidx, typename Telt>
PermutationElt<Tidx, Telt> Inverse(PermutationElt<Tidx, Telt> const &ePermElt) {
  return ~ePermElt;
}

//
// Prinouts and other operations
//

template <typename Tidx, typename Telt>
std::string GapStyleStringShift(PermutationElt<Tidx, Telt> const &ePermElt,
                                int const &eShift) {
  Tidx n = ePermElt.size();
  Face ListStat(n);
  std::string eRet;

  for (Tidx i = 0; i < n; i++) {
    if (ListStat[i] == 0) {
      Tidx eFirst = i;
      Tidx eCurr = i;
      std::string ePart = "(";
      bool IsFirst = true;
      Tidx len = 0;
      while (true) {
        if (!IsFirst)
          ePart += ",";
        IsFirst = false;
        ePart += std::to_string(eCurr + eShift);
        ListStat[eCurr] = 1;
        Tidx eNext = ePermElt.at(eCurr);
        len++;
        if (eNext == eFirst)
          break;
        eCurr = eNext;
      }
      ePart += ")";
      if (len > 1)
        eRet += ePart;
    }
  }
  if (eRet.size() > 0)
    return eRet;
  return "()";
}

template <typename Tidx, typename Telt>
std::string GapStyleString(PermutationElt<Tidx, Telt> const &ePermElt) {
  return GapStyleStringShift(ePermElt, 1);
}

template <typename Tidx, typename T>
std::ostream &operator<<(std::ostream &os,
                         PermutationElt<Tidx, T> const &ePermElt) {
  os << GapStyleStringShift(ePermElt, 1);
  //  os << "M =";
  //  os << ePermElt.getElt();
  return os;
}

// clang-format off
}  // namespace permutalib
// clang-format on

namespace std {
template <typename Tidx, typename Telt>
struct hash<permutalib::PermutationElt<Tidx, Telt>> {
  std::size_t
  operator()(const permutalib::PermutationElt<Tidx, Telt> &e_val) const {
    uint32_t seed = 0x1b873540;
    const Tidx *ptr_tidx = e_val.getPtr();
    const uint8_t *ptr_i = (const uint8_t *)ptr_tidx;
    size_t len = sizeof(Tidx) * e_val.size();
    size_t hash1 = permutalib::robin_hood_hash_bytes(ptr_i, len, seed);
    size_t hash2 = std::hash<Telt>()(e_val.elt);
    hash1 ^= hash2 + 0x9e3779b9 + (hash1 << 6) + (hash1 >> 2);
    return hash1;
  }
};

// clang-format off
}  // namespace std
#endif  // SRC_GAP_PERMUTATIONELT_H_
// clang-format on
