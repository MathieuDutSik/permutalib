// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GAP_PERMUTATION_H_
#define SRC_GAP_PERMUTATION_H_

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

template <typename Tidx> bool CheckList(const std::vector<Tidx> &l) {
  size_t len = l.size();
  if (len >= size_t(std::numeric_limits<Tidx>::max())) {
    std::cerr << "Number of rows is higher than expected\n";
    return false;
  }
  Face f(len);
  for (size_t i = 0; i < len; i++) {
    if (size_t(l[i]) >= len) {
      std::cerr << "The value of l is too high\n";
      return false;
    }
    f[l[i]] = 1;
  }
  size_t covered = 0;
  for (size_t i = 0; i < len; i++)
    covered += f[i];
  if (covered != len) {
    std::cerr << "covering errors\n";
    return false;
  }
  return true;
}

template <typename Tidx>
std::pair<std::vector<Tidx>, std::vector<Tidx>>
GetListValRev(std::string const &estr) {
  size_t maxlen = 0;
  std::vector<Tidx> ListVal;
  std::vector<Tidx> ListRev;
  auto insertLVal = [&](std::vector<Tidx> const &LVal) -> void {
    for (auto &eVal : LVal)
      if (eVal + 1 >= static_cast<int>(maxlen))
        maxlen = eVal + 1;
    for (size_t pos = ListVal.size(); pos < maxlen; pos++) {
      ListVal[pos] = pos;
      ListRev[pos] = pos;
    }
    size_t len = LVal.size();
    for (size_t i = 0; i < len; i++) {
      size_t j = i + 1;
      if (j == len)
        j = 0;
      Tidx val1 = LVal[i];
      Tidx val2 = LVal[j];
      ListVal[val1] = val2;
      ListRev[val2] = val1;
    }
  };
  auto ParseStringByComma = [&](std::string const &estr) -> std::vector<Tidx> {
    size_t n_char = estr.size();
    size_t pos_start = 0;
    std::vector<Tidx> LVal;
    auto insert = [&](size_t const &pos1, size_t const &pos2) -> void {
      size_t len = pos2 - pos1;
      std::string ustr = estr.substr(pos_start, len);
      Tidx eVal = std::stoi(std::string(ustr)) - 1;
      LVal.push_back(eVal);
      pos_start = pos2 + 1;
    };
    for (size_t i_char = 0; i_char < n_char; i_char++) {
      std::string echar = estr.substr(i_char, 1);
      if (echar == ",")
        insert(pos_start, i_char);
    }
    insert(pos_start, n_char);
    return LVal;
  };
  //
  size_t n_char = estr.size();
  size_t pos_start = 0;
  for (size_t i_char = 0; i_char < n_char; i_char++) {
    if (estr.substr(i_char, 1) == "(") {
      pos_start = i_char + 1;
    }
    if (estr.substr(i_char, 1) == ")") {
      size_t pos_end = i_char;
      size_t len = pos_end - pos_start;
      std::string sstr = estr.substr(pos_start, len);
      std::vector<Tidx> LVal = ParseStringByComma(sstr);
      insertLVal(LVal);
    }
  }
  return {std::move(ListVal), std::move(ListRev)};
}

template <typename Tidx_inp> struct DoubleSidedPerm {
public:
  using Tidx = Tidx_inp;
  //
  // The constructors
  //
  DoubleSidedPerm(std::string const &estr) {
    std::pair<std::vector<Tidx>, std::vector<Tidx>> epair =
        GetListValRev<Tidx>(estr);
    ListVal = epair.first;
    ListRev = epair.second;
    siz = ListVal.size();
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
    if (!CheckList(ListVal) || !CheckList(ListVal)) {
      std::cerr << "ListVal or ListRev are not rightly organized\n";
      throw PermutalibException{1};
    }
#endif
  }
  DoubleSidedPerm(DoubleSidedPerm const &ePerm, Tidx const &n) {
    if (ePerm.size() > n) {
      std::cerr << "ePerm.size()=" << ePerm.size() << " n=" << n << "\n";
      std::cerr << "ExtendPermutation to a size that is lower than the current "
                   "size\n";
      throw PermutalibException{1};
    }
    ListVal = ePerm.getListVal();
    ListRev = ePerm.getListRev();
    for (Tidx pos = ePerm.size(); pos < n; pos++) {
      ListVal.push_back(pos);
      ListRev.push_back(pos);
    }
    siz = n;
  }
  DoubleSidedPerm() {
    siz = 0;
    ListVal = {};
    ListRev = {};
  }
  DoubleSidedPerm(Tidx const &n) {
    siz = n;
    ListVal = std::vector<Tidx>(n);
    ListRev = std::vector<Tidx>(n);
    for (Tidx i = 0; i < n; i++) {
      ListVal[i] = i;
      ListRev[i] = i;
    }
  }
  DoubleSidedPerm(std::vector<Tidx> const &v) {
    ListVal = v;
    siz = Tidx(v.size());
    ListRev.resize(siz);
    for (Tidx i = 0; i < siz; i++)
      ListRev[v[i]] = i;
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
    if (!CheckList(ListVal) || !CheckList(ListRev)) {
      std::cerr << "ListVal or ListRev are not rightly organized\n";
      throw PermutalibException{1};
    }
#endif
  }
  DoubleSidedPerm(std::vector<Tidx> const &v1, std::vector<Tidx> const &v2) {
    siz = Tidx(v1.size());
    ListVal = v1;
    ListRev = v2;
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
    if (!CheckList(ListVal) || !CheckList(ListRev)) {
      std::cerr << "ListVal or ListRev are not rightly organized\n";
      throw PermutalibException{1};
    }
#endif
  }
  DoubleSidedPerm(DoubleSidedPerm const &ePerm) {
    siz = ePerm.siz;
    ListVal = ePerm.ListVal;
    ListRev = ePerm.ListRev;
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
    if (!CheckList(ListVal) || !CheckList(ListRev)) {
      std::cerr << "ListVal or ListRev are not rightly organized\n";
      throw PermutalibException{1};
    }
#endif
  }
  DoubleSidedPerm(DoubleSidedPerm &&ePerm) {
    siz = ePerm.siz;
    ListVal = std::move(ePerm.ListVal);
    ListRev = std::move(ePerm.ListRev);
    ePerm.siz = 0;
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
    if (!CheckList(ListVal) || !CheckList(ListRev)) {
      std::cerr << "ListVal or ListRev are not rightly organized\n";
      throw PermutalibException{1};
    }
#endif
  }
  //
  // Copy operator
  //
  DoubleSidedPerm<Tidx> operator=(DoubleSidedPerm const &ePerm) {
    siz = ePerm.siz;
    ListVal = ePerm.ListVal;
    ListRev = ePerm.ListRev;
    return *this;
  }
  DoubleSidedPerm<Tidx> operator=(DoubleSidedPerm &&ePerm) {
    siz = ePerm.siz;
    ListVal = std::move(ePerm.ListVal);
    ListRev = std::move(ePerm.ListRev);
    ePerm.siz = 0;
    return *this;
  }
  //
  // The destructor
  //
  ~DoubleSidedPerm() {}
  //
  // The functionality
  //
  bool isIdentity() const {
    for (Tidx i = 0; i < siz; i++)
      if (ListVal[i] != i)
        return false;
    return true;
  }
  Tidx at(Tidx const &i) const { return ListVal[i]; }
  Tidx atRev(Tidx const &i) const { return ListRev[i]; }
  const Tidx *getPtr() const { return ListVal.data(); }
  const std::vector<Tidx> &getListVal() const { return ListVal; }
  const std::vector<Tidx> &getListRev() const { return ListRev; }
  Tidx operator[](Tidx const &i) const { return ListVal[i]; }
  Tidx size() const { return siz; }
  //
public: // Should be private but simpler to do it like that.
  Tidx siz;
  std::vector<Tidx> ListVal;
  std::vector<Tidx> ListRev;
};

template <typename Tidx>
bool operator==(DoubleSidedPerm<Tidx> const &v1,
                DoubleSidedPerm<Tidx> const &v2) {
  Tidx siz = v1.size();
  if (siz != v2.size())
    return false;
  for (Tidx i = 0; i < siz; i++)
    if (v1.at(i) != v2.at(i))
      return false;
  return true;
}

template <typename Tidx>
bool operator!=(DoubleSidedPerm<Tidx> const &v1,
                DoubleSidedPerm<Tidx> const &v2) {
  Tidx siz = v1.size();
  if (siz != v2.size())
    return true;
  for (Tidx i = 0; i < siz; i++)
    if (v1.at(i) != v2.at(i))
      return true;
  return false;
}

template <typename Tidx>
bool operator<(DoubleSidedPerm<Tidx> const &v1,
               DoubleSidedPerm<Tidx> const &v2) {
  Tidx siz1 = v1.size();
  Tidx siz2 = v2.size();
  if (siz1 != siz2)
    return siz1 < siz2;
  Tidx siz = siz1;
  for (Tidx i = 0; i < siz; i++) {
    if (v1.at(i) != v2.at(i))
      return v1.at(i) < v2.at(i);
  }
  return false;
}

template <typename Tidx>
DoubleSidedPerm<Tidx> operator~(DoubleSidedPerm<Tidx> const &ePerm) {
  return DoubleSidedPerm<Tidx>(ePerm.getListRev(), ePerm.getListVal());
}

// Form the product v1 * v2
template <typename Tidx>
DoubleSidedPerm<Tidx> operator*(DoubleSidedPerm<Tidx> const &v1,
                                DoubleSidedPerm<Tidx> const &v2) {
  size_t siz = v1.size();
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
  if (siz != v2.size()) {
    std::cerr << "Error in the DoubleSidedPerm product\n";
    throw PermutalibException{1};
  }
#endif
  std::vector<Tidx> vVal(siz), vRev(siz);
  Tidx siz_idx = Tidx(siz);
  for (Tidx i = 0; i < siz_idx; i++) {
    Tidx j = v1.at(i);
    Tidx k = v2.at(j);
    vVal[i] = k;
    //
    Tidx j2 = v2.atRev(i);
    Tidx k2 = v1.atRev(j2);
    vRev[i] = k2;
  }
  return DoubleSidedPerm<Tidx>(vVal, vRev);
}

template <typename Tidx>
void operator*=(DoubleSidedPerm<Tidx> &v1, DoubleSidedPerm<Tidx> const &v2) {
  size_t siz = v1.size();
  std::vector<Tidx> vRev(siz);
  Tidx siz_idx = Tidx(siz);
  for (Tidx i = 0; i < siz_idx; i++) {
    Tidx j = v1.at(i);
    Tidx k = v2.at(j);
    v1.ListVal[i] = k;
    //
    Tidx j2 = v2.atRev(i);
    Tidx k2 = v1.atRev(j2);
    vRev[i] = k2;
  }
  v1.ListRev = vRev;
}

template <typename Tidx>
DoubleSidedPerm<Tidx> Conjugation(DoubleSidedPerm<Tidx> const &v1,
                                  DoubleSidedPerm<Tidx> const &v2) {
  size_t siz = v1.size();
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
  if (siz != v2.size()) {
    std::cerr << "Error in the DoubleSidedPerm conjugation\n";
    throw PermutalibException{1};
  }
#endif
  std::vector<Tidx> v(siz);
  for (size_t i = 0; i < siz; i++) {
    int j = v1[i];
    int i2 = v2[i];
    int j2 = v2[j];
    v[i2] = j2;
  }
  return DoubleSidedPerm<Tidx>(v);
}

template <typename Tidx>
Tidx PowAct(Tidx const &i, DoubleSidedPerm<Tidx> const &g) {
  return g.at(i);
}

template <typename Tidx>
Tidx SlashAct(Tidx const &i, DoubleSidedPerm<Tidx> const &g) {
  return g.atRev(i);
}

// LeftQuotient(x,y) = x^{-1}*y in the list.gi file
template <typename Tidx>
DoubleSidedPerm<Tidx> LeftQuotient(DoubleSidedPerm<Tidx> const &a,
                                   DoubleSidedPerm<Tidx> const &b) {
  size_t siz = a.size();
  std::vector<Tidx> ListVal(siz), ListRev(siz);
  Tidx siz_idx = Tidx(siz);
  for (Tidx i = 0; i < siz_idx; i++) {
    Tidx i1 = a.atRev(i);
    Tidx j1 = b.at(i1);
    ListVal[i] = j1;
    Tidx i2 = b.atRev(i);
    Tidx j2 = a.at(i2);
    ListRev[i] = j2;
  }
  return DoubleSidedPerm<Tidx>(ListVal, ListRev);
}

template <typename Tidx> DoubleSidedPerm<Tidx> SCRandomPerm(int const &d) {
  std::vector<Tidx> rnd(d);
  for (int i = 0; i < d; i++)
    rnd[i] = i;
  for (int i = 0; i < d; i++) {
    int idx = d - i;
    int res = d - i;
    int k = random() % res;
    if (k != idx) {
      int tmp = rnd[idx];
      rnd[idx] = rnd[k];
      rnd[k] = tmp;
    }
  }
  return DoubleSidedPerm<Tidx>(rnd);
}

template <typename Tidx>
DoubleSidedPerm<Tidx> Inverse(DoubleSidedPerm<Tidx> const &ePerm) {
  return ~ePerm;
}

// Input / Output

template <typename Tidx>
std::string GapStyleStringShift(DoubleSidedPerm<Tidx> const &ePerm,
                                int const &eShift) {
  Tidx n = ePerm.size();
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
        Tidx eNext = ePerm.at(eCurr);
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

template <typename Tidx>
std::string GapStyleString(DoubleSidedPerm<Tidx> const &ePerm) {
  return GapStyleStringShift(ePerm, 1);
}

template <typename Tidx>
std::ostream &operator<<(std::ostream &os, DoubleSidedPerm<Tidx> const &ePerm) {
  os << GapStyleStringShift(ePerm, 1);
  return os;
}

template <typename Tidx_inp> struct SingleSidedPerm {
public:
  using Tidx = Tidx_inp;
  //
  // The constructors
  //
  SingleSidedPerm(std::string const &estr) {
    std::pair<std::vector<Tidx>, std::vector<Tidx>> epair =
        GetListValRev<Tidx>(estr);
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
    if (epair.first.size() >= std::numeric_limits<Tidx>::max() - 1) {
      std::cerr << "Tidx is too small for representing the vector\n";
      throw PermutalibException{1};
    }
#endif
    ListVal = std::move(epair.first);
    siz = ListVal.size();
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
    if (!CheckList(ListVal)) {
      std::cerr << "ListVal does not define a permutation\n";
      throw PermutalibException{1};
    }
#endif
  }
  SingleSidedPerm(SingleSidedPerm const &ePerm, int const &n) {
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
    if (n >= std::numeric_limits<Tidx>::max() - 1) {
      std::cerr << "Tidx is too small for representing the vector\n";
      throw PermutalibException{1};
    }
#endif
    if (ePerm.size() > n) {
      std::cerr << "ePerm.size()=" << ePerm.size() << " n=" << n << "\n";
      std::cerr << "ExtendPermutation to a size that is lower than the current "
                   "size\n";
    }
    ListVal = ePerm.getListVal();
    for (int pos = ePerm.size(); pos < n; pos++)
      ListVal.push_back(pos);
    siz = n;
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
    if (!CheckList(ListVal)) {
      std::cerr << "ListVal does not define a permutation\n";
      throw PermutalibException{1};
    }
#endif
  }
  SingleSidedPerm() {
    siz = 0;
    ListVal = {};
  }
  SingleSidedPerm(Tidx const &n) {
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
    if (n >= std::numeric_limits<Tidx>::max() - 1) {
      std::cerr << "Tidx is too small for representing the vector\n";
      throw PermutalibException{1};
    }
#endif
    siz = n;
    ListVal = std::vector<Tidx>(n);
    for (Tidx i = 0; i < n; i++)
      ListVal[i] = i;
  }
  SingleSidedPerm(std::vector<Tidx> &&v) {
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
    if (v.size() >= std::numeric_limits<Tidx>::max() - 1) {
      std::cerr << "Tidx is too small for representing the vector\n";
      throw PermutalibException{1};
    }
#endif
    ListVal = v;
    siz = Tidx(v.size());
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
    if (!CheckList(ListVal)) {
      std::cerr << "ListVal does not define a permutation\n";
      throw PermutalibException{1};
    }
#endif
  }
  SingleSidedPerm(std::vector<Tidx> const &v) {
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
    if (v.size() >= std::numeric_limits<Tidx>::max() - 1) {
      std::cerr << "Tidx is too small for representing the vector\n";
      throw PermutalibException{1};
    }
#endif
    ListVal = v;
    siz = Tidx(v.size());
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
    if (!CheckList(ListVal)) {
      std::cerr << "ListVal does not define a permutation\n";
      throw PermutalibException{1};
    }
#endif
  }
  SingleSidedPerm(std::vector<Tidx> const &v1,
                  [[maybe_unused]] std::vector<Tidx> const &v2) {
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
    if (v1.size() >= std::numeric_limits<Tidx>::max() - 1) {
      std::cerr << "Tidx is too small for representing the vector\n";
      throw PermutalibException{1};
    }
#endif
    siz = v1.size();
    ListVal = v1;
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
    if (!CheckList(ListVal)) {
      std::cerr << "ListVal does not define a permutation\n";
      throw PermutalibException{1};
    }
#endif
  }
  SingleSidedPerm(SingleSidedPerm const &ePerm) {
    siz = ePerm.siz;
    ListVal = ePerm.ListVal;
  }
  SingleSidedPerm(SingleSidedPerm &&ePerm) {
    siz = ePerm.siz;
    ListVal = std::move(ePerm.ListVal);
    ePerm.siz = 0;
  }
  SingleSidedPerm(std::initializer_list<Tidx> l)
      : siz(Tidx(l.size())), ListVal(l) {
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
    if (ListVal.size() >= std::numeric_limits<Tidx>::max() - 1) {
      std::cerr << "Tidx is too small for representing the vector\n";
      throw PermutalibException{1};
    }
#endif
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
    if (!CheckList(ListVal)) {
      std::cerr << "ListVal does not define a permutation\n";
      throw PermutalibException{1};
    }
#endif
  }
  //
  // Copy operator
  //
  SingleSidedPerm<Tidx> operator=(SingleSidedPerm const &ePerm) {
    siz = ePerm.siz;
    ListVal = ePerm.ListVal;
    return *this;
  }
  SingleSidedPerm<Tidx> operator=(SingleSidedPerm &&ePerm) {
    siz = ePerm.siz;
    ListVal = std::move(ePerm.ListVal);
    ePerm.siz = 0;
    return *this;
  }
  //
  // The destructor
  //
  ~SingleSidedPerm() {}
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
  }
  const Tidx *getPtr() const { return ListVal.data(); }
  const std::vector<Tidx> &getListVal() const { return ListVal; }
  Tidx operator[](Tidx const &i) const { return ListVal[i]; }
  Tidx size() const { return siz; }
  //
public:  // Should be private in a more classic construction
  Tidx siz;
  std::vector<Tidx> ListVal;
};

template <typename Tidx>
bool operator==(SingleSidedPerm<Tidx> const &v1,
                SingleSidedPerm<Tidx> const &v2) {
  Tidx siz = v1.size();
  if (siz != v2.size())
    return false;
  for (Tidx i = 0; i < siz; i++)
    if (v1.at(i) != v2.at(i))
      return false;
  return true;
}

template <typename Tidx>
bool operator!=(SingleSidedPerm<Tidx> const &v1,
                SingleSidedPerm<Tidx> const &v2) {
  Tidx siz = v1.size();
  if (siz != v2.size())
    return true;
  for (Tidx i = 0; i < siz; i++)
    if (v1.at(i) != v2.at(i))
      return true;
  return false;
}

template <typename Tidx>
bool operator<(SingleSidedPerm<Tidx> const &v1,
               SingleSidedPerm<Tidx> const &v2) {
  Tidx siz1 = v1.size();
  Tidx siz2 = v2.size();
  if (siz1 != siz2)
    return siz1 < siz2;
  Tidx siz = siz1;
  for (Tidx i = 0; i < siz; i++) {
    if (v1.at(i) != v2.at(i))
      return v1.at(i) < v2.at(i);
  }
  return false;
}

template <typename Tidx>
SingleSidedPerm<Tidx> operator~(SingleSidedPerm<Tidx> const &ePerm) {
  Tidx siz = ePerm.size();
  const std::vector<Tidx> &LVal = ePerm.getListVal();
  std::vector<Tidx> v(siz);
  for (Tidx i = 0; i < siz; i++)
    v[LVal[i]] = i;
  return SingleSidedPerm<Tidx>(std::move(v));
}

// Form the product v1 * v2
template <typename Tidx>
SingleSidedPerm<Tidx> operator*(SingleSidedPerm<Tidx> const &v1,
                                SingleSidedPerm<Tidx> const &v2) {
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
  return SingleSidedPerm<Tidx>(
      std::move(vVal)); // Does not seem to get us a speedup
}

template <typename Tidx>
void operator*=(SingleSidedPerm<Tidx> &v1, SingleSidedPerm<Tidx> const &v2) {
  Tidx siz_idx = Tidx(v1.size());
  for (Tidx i = 0; i < siz_idx; i++) {
    Tidx j = v1.at(i);
    Tidx k = v2.at(j);
    v1.ListVal[i] = k;
  }
}

template <typename Tidx>
SingleSidedPerm<Tidx> Conjugation(SingleSidedPerm<Tidx> const &v1,
                                  SingleSidedPerm<Tidx> const &v2) {
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
  return SingleSidedPerm<Tidx>(std::move(v));
}

template <typename Tidx>
Tidx PowAct(Tidx const &i, SingleSidedPerm<Tidx> const &g) {
  return g.at(i);
}

template <typename Tidx>
Tidx SlashAct(Tidx const &i, SingleSidedPerm<Tidx> const &g) {
  return g.atRev(i);
}

// LeftQuotient(x,y) = x^{-1}*y in the list.gi file
template <typename Tidx>
SingleSidedPerm<Tidx> LeftQuotient(SingleSidedPerm<Tidx> const &a,
                                   SingleSidedPerm<Tidx> const &b) {
  Tidx siz = Tidx(a.size());
  const std::vector<Tidx> &Val_A = a.getListVal();
  const std::vector<Tidx> &Val_B = b.getListVal();
  std::vector<Tidx> ListVal(siz);
  for (Tidx i = 0; i < siz; i++)
    ListVal[Val_A[i]] = Val_B[i];
  return SingleSidedPerm<Tidx>(std::move(ListVal));
}

template <typename Tidx> SingleSidedPerm<Tidx> SCRandomPerm(Tidx const &d) {
  std::vector<Tidx> rnd(d);
  for (Tidx i = 0; i < d; i++)
    rnd[i] = i;
  for (Tidx i = 0; i < d; i++) {
    Tidx idx = d - i;
    Tidx res = d - i;
    Tidx k = random() % res;
    if (k != idx) {
      Tidx tmp = rnd[idx];
      rnd[idx] = rnd[k];
      rnd[k] = tmp;
    }
  }
  return DoubleSidedPerm<Tidx>(rnd);
}

template <typename Tidx>
SingleSidedPerm<Tidx> Inverse(SingleSidedPerm<Tidx> const &ePerm) {
  return ~ePerm;
}

template <typename Tidx>
std::string GapStyleStringShift(SingleSidedPerm<Tidx> const &ePerm,
                                int const &eShift) {
  Tidx n = ePerm.size();
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
        Tidx eNext = ePerm.at(eCurr);
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

template <typename Tidx>
std::string GapStyleString(SingleSidedPerm<Tidx> const &ePerm) {
  return GapStyleStringShift(ePerm, 1);
}

template <typename Tidx>
std::ostream &operator<<(std::ostream &os, SingleSidedPerm<Tidx> const &ePerm) {
  os << GapStyleStringShift(ePerm, 1);
  return os;
}

}  // namespace permutalib

namespace std {
template <typename Tidx>
std::string to_string(permutalib::DoubleSidedPerm<Tidx> const &ePerm) {
  return GapStyleStringShift(ePerm, 1);
}
template <typename Tidx>
std::string to_string(permutalib::SingleSidedPerm<Tidx> const &ePerm) {
  return GapStyleStringShift(ePerm, 1);
}
}  // namespace std

namespace std {
template <typename Tidx> struct hash<permutalib::SingleSidedPerm<Tidx>> {
  std::size_t operator()(const permutalib::SingleSidedPerm<Tidx> &e_val) const {
    uint32_t seed = 0x1b873540;
    const Tidx *ptr_tidx = e_val.getPtr();
    const uint8_t *ptr_i = (const uint8_t *)ptr_tidx;
    size_t len = sizeof(Tidx) * e_val.size();
    return permutalib::robin_hood_hash_bytes(ptr_i, len, seed);
  }
};
template <typename Tidx> struct hash<permutalib::DoubleSidedPerm<Tidx>> {
  std::size_t operator()(const permutalib::DoubleSidedPerm<Tidx> &e_val) const {
    uint32_t seed = 0x1b873540;
    const Tidx *ptr_tidx = e_val.getPtr();
    const uint8_t *ptr_i = (const uint8_t *)ptr_tidx;
    size_t len = sizeof(Tidx) * e_val.size();
    return permutalib::robin_hood_hash_bytes(ptr_i, len, seed);
  }
};

// clang-format off
}  // namespace std
#endif  // SRC_GAP_PERMUTATION_H_
// clang-format on
