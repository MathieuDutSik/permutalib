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
    //    std::cerr << "PermutationElt 1\n";
    //    NicePrint("Constructor 1", *this);
  }
  PermutationElt(std::vector<Tidx> &&v, Telt &&_elt) {
    //    std::cerr << "PermutationElt 2\n";
    //    std::cerr << "_elt=";
    //    WriteMatrix(std::cerr, _elt);
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
    if (v.size() >= std::numeric_limits<Tidx>::max() - 1) {
      std::cerr << "Tidx is too small for representing the vector\n";
      throw PermutalibException{1};
    }
#endif
    ListVal = v;
    siz = Tidx(v.size());
    elt = std::move(_elt);
    //    NicePrint("Constructor 2", *this);
  }
  PermutationElt(std::vector<Tidx> const &v, Telt const &_elt) {
    //    std::cerr << "PermutationElt 3\n";
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
    if (v.size() >= std::numeric_limits<Tidx>::max() - 1) {
      std::cerr << "Tidx is too small for representing the vector\n";
      throw PermutalibException{1};
    }
#endif
    ListVal = v;
    siz = Tidx(v.size());
    elt = _elt;
    //    NicePrint("Constructor 3", *this);
  }
  PermutationElt(PermutationElt const &ePermElt) {
    //    std::cerr << "PermutationElt 4\n";
    siz = ePermElt.siz;
    ListVal = ePermElt.ListVal;
    elt = ePermElt.elt;
    //    NicePrint("Constructor 4", *this);
  }
  PermutationElt(PermutationElt &&ePermElt) {
    //    std::cerr << "PermutationElt 5\n";
    siz = ePermElt.siz;
    ListVal = std::move(ePermElt.ListVal);
    elt = std::move(ePermElt.elt);
    ePermElt.siz = 0;
    //    NicePrint("Constructor 5", *this);
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
  // The destructor
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
    /*
    for (Tidx j=0; j<siz; j++)
      if (ListVal[j] == i)
        return j;
    return std::numeric_limits<Tidx>::max();
    */
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
    std::cerr << "Error in the DoubleSidedPerm product\n";
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
