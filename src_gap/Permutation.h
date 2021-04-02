#ifndef DEFINE_PERMUTALIB_PERMUTATION_H
#define DEFINE_PERMUTALIB_PERMUTATION_H

#include <vector>
#include <stdlib.h>
#include <string>
#include <iostream>

namespace permutalib {


template<typename Tidx>
std::pair<std::vector<Tidx>, std::vector<Tidx>> GetListValRev(std::string const& estr)
{
  size_t maxlen = 0;
  std::vector<Tidx> ListVal;
  std::vector<Tidx> ListRev;
  auto insertLVal=[&](std::vector<Tidx> const& LVal) -> void {
    for (auto & eVal : LVal)
      if (eVal+1 >= int(maxlen))
        maxlen = eVal + 1;
    for (size_t pos=ListVal.size(); pos<maxlen; pos++) {
      ListVal[pos] = pos;
      ListRev[pos] = pos;
    }
    size_t len = LVal.size();
    for (size_t i=0; i<len; i++) {
      size_t j = i+1;
      if (j == len)
        j = 0;
      Tidx val1 = LVal[i];
      Tidx val2 = LVal[j];
      ListVal[val1] = val2;
      ListRev[val2] = val1;
    }
  };
  auto ParseStringByComma=[&](std::string const& estr) -> std::vector<Tidx> {
    size_t n_char=estr.size();
    size_t pos_start = 0;
    std::vector<Tidx> LVal;
    auto insert=[&](size_t const& pos1, size_t const& pos2) -> void {
      size_t len = pos2 - pos1;
      std::string ustr = estr.substr(pos_start, len);
      Tidx eVal = std::stoi(std::string(ustr)) - 1;
      LVal.push_back(eVal);
      pos_start = pos2 + 1;
    };
    for (size_t i_char=0; i_char<n_char; i_char++) {
      std::string echar = estr.substr(i_char, 1);
      if (echar == ",")
        insert(pos_start, i_char);
    }
    insert(pos_start, n_char);
    return LVal;
  };
  //
  size_t n_char = estr.size();
  size_t pos_start=0;
  for (size_t i_char=0; i_char<n_char; i_char++) {
    if (estr.substr(i_char,1) == "(") {
      pos_start = i_char + 1;
    }
    if (estr.substr(i_char,1) == ")") {
      size_t pos_end = i_char;
      size_t len = pos_end - pos_start;
      std::string sstr = estr.substr(pos_start, len);
      std::vector<Tidx> LVal = ParseStringByComma(sstr);
      insertLVal(LVal);
    }
  }
  return {ListVal,ListRev};
}




template<typename Tidx_inp>
struct DoubleSidedPerm {
public:
  using Tidx = Tidx_inp;
  //
  // The constructors
  //
  DoubleSidedPerm(std::string const& estr)
  {
    std::pair<std::vector<Tidx>, std::vector<Tidx>> epair = GetListValRev<Tidx>(estr);
    ListVal = epair.first;
    ListRev = epair.second;
    siz = ListVal.size();
  }
  DoubleSidedPerm(DoubleSidedPerm const& ePerm, int const& n)
  {
    if (ePerm.size() > n) {
      std::cerr << "ePerm.size()=" << ePerm.size() << " n=" << n << "\n";
      std::cerr << "ExtendPermutation to a size that is lower than the current size\n";
    }
    ListVal = ePerm.getListVal();
    ListRev = ePerm.getListRev();
    for (int pos=ePerm.size(); pos<n; pos++) {
      ListVal.push_back(pos);
      ListRev.push_back(pos);
    }
    siz = n;
  }
  DoubleSidedPerm ()
  {
    siz=0;
    ListVal = {};
    ListRev = {};
  }
  DoubleSidedPerm (int const& n)
  {
    siz=n;
    ListVal = std::vector<Tidx>(n);
    ListRev = std::vector<Tidx>(n);
    for (int i=0; i<n; i++) {
      ListVal[i]=i;
      ListRev[i]=i;
    }
  }
  DoubleSidedPerm(std::vector<Tidx> const& v)
  {
    ListVal=v;
    siz=v.size();
    ListRev.resize(siz);
    for (int i=0; i<siz; i++)
      ListRev[v[i]]=i;
  }
  DoubleSidedPerm(std::vector<Tidx> const& v1, std::vector<Tidx> const& v2)
  {
    siz=v1.size();
    ListVal=v1;
    ListRev=v2;
  }
  DoubleSidedPerm(DoubleSidedPerm const& ePerm)
  {
    siz     = ePerm.siz;
    ListVal = ePerm.ListVal;
    ListRev = ePerm.ListRev;
  }
  DoubleSidedPerm(DoubleSidedPerm&& ePerm)
  {
    siz = ePerm.siz;
    ListVal = std::move(ePerm.ListVal);
    ListRev = std::move(ePerm.ListRev);
    ePerm.siz = 0;
  }
  //
  // Copy operator
  //
  DoubleSidedPerm<Tidx> operator=(DoubleSidedPerm const& ePerm)
  {
    siz     = ePerm.siz;
    ListVal = ePerm.ListVal;
    ListRev = ePerm.ListRev;
    return *this;
  }
  DoubleSidedPerm<Tidx> operator=(DoubleSidedPerm&& ePerm)
  {
    siz = ePerm.siz;
    ListVal = std::move(ePerm.ListVal);
    ListRev = std::move(ePerm.ListRev);
    ePerm.siz = 0;
    return *this;
  }
  //
  // The destructor
  //
  ~DoubleSidedPerm()
  {
  }
  //
  // The destructor
  //
  bool isIdentity() const
  {
    for (int i=0; i<siz; i++)
      if (ListVal[i] != i)
	return false;
    return true;
  }
  int at(int const& i) const
  {
    return ListVal[i];
  }
  int atRev(int const& i) const
  {
    return ListRev[i];
  }
  std::vector<Tidx> getListVal() const
  {
    return ListVal;
  }
  std::vector<Tidx> getListRev() const
  {
    return ListRev;
  }
  int operator[](int const& i) const
  {
    return ListVal[i];
  }
  int size() const
  {
    return siz;
  }
  //
private:
  int siz;
  std::vector<Tidx> ListVal;
  std::vector<Tidx> ListRev;
};



template<typename Tidx>
bool operator==(DoubleSidedPerm<Tidx> const& v1, DoubleSidedPerm<Tidx> const& v2)
{
  int siz=v1.size();
  if (siz != v2.size() )
    return false;
  for (int i=0; i<siz; i++)
    if (v1.at(i) != v2.at(i))
      return false;
  return true;
}


template<typename Tidx>
bool operator!=(DoubleSidedPerm<Tidx> const& v1, DoubleSidedPerm<Tidx> const& v2)
{
  int siz=v1.size();
  if (siz != v2.size() )
    return true;
  for (int i=0; i<siz; i++)
    if (v1.at(i) != v2.at(i))
      return true;
  return false;
}


template<typename Tidx>
bool operator<(DoubleSidedPerm<Tidx> const& v1, DoubleSidedPerm<Tidx> const& v2)
{
  int siz1=v1.size();
  int siz2=v2.size();
  if (siz1 != siz2)
    return siz1<siz2;
  int siz=siz1;
  for (int i=0; i<siz; i++) {
    if (v1.at(i) != v2.at(i))
      return v1.at(i) < v2.at(i);
  }
  return false;
}

template<typename Tidx>
DoubleSidedPerm<Tidx> operator~(DoubleSidedPerm<Tidx> const& ePerm)
{
  return DoubleSidedPerm<Tidx>(ePerm.getListRev(), ePerm.getListVal());
}







// Form the product v1 * v2
template<typename Tidx>
DoubleSidedPerm<Tidx> operator*(DoubleSidedPerm<Tidx> const& v1, DoubleSidedPerm<Tidx> const& v2)
{
  int siz=v1.size();
#ifdef DEBUG
  if (siz != v2.size() ) {
    std::cerr << "Error in the DoubleSidedPerm product\n";
    throw PermutalibException{1};
  }
#endif
  std::vector<Tidx> vVal(siz), vRev(siz);
  for (int i=0; i<siz; i++) {
    int j=v1.at(i);
    int k=v2.at(j);
    vVal[i]=k;
    //
    int j2=v2.atRev(i);
    int k2=v1.atRev(j2);
    vRev[i]=k2;
  }
  return DoubleSidedPerm<Tidx>(vVal, vRev);
}



template<typename Tidx>
DoubleSidedPerm<Tidx> Conjugation(DoubleSidedPerm<Tidx> const& v1, DoubleSidedPerm<Tidx> const& v2)
{
  int siz=v1.size();
#ifdef DEBUG
  if (siz != v2.size() ) {
    std::cerr << "Error in the DoubleSidedPerm conjugation\n";
    throw PermutalibException{1};
  }
#endif
  std::vector<Tidx> v(siz);
  for (int i=0; i<siz; i++) {
    int j=v1[i];
    int i2=v2[i];
    int j2=v2[j];
    v[i2]=j2;
  }
  return DoubleSidedPerm<Tidx>(v);
}



template<typename Tidx>
int PowAct(int const& i, DoubleSidedPerm<Tidx> const& g)
{
  return g.at(i);
}



template<typename Tidx>
int SlashAct(int const& i, DoubleSidedPerm<Tidx> const& g)
{
  return g.atRev(i);
}


// LeftQuotient(x,y) = x^{-1}*y in the list.gi file
template<typename Tidx>
DoubleSidedPerm<Tidx> LeftQuotient(DoubleSidedPerm<Tidx> const& a, DoubleSidedPerm<Tidx> const& b)
{
  int siz=a.size();
  std::vector<Tidx> ListVal(siz), ListRev(siz);
  for (int i=0; i<siz; i++) {
    int i1=a.atRev(i);
    int j1=b.at(i1);
    ListVal[i]=j1;
    int i2=b.atRev(i);
    int j2=a.at(i2);
    ListRev[i]=j2;
  }
  return DoubleSidedPerm<Tidx>(ListVal, ListRev);
}



template<typename Tidx>
DoubleSidedPerm<Tidx> SCRandomPerm(int const& d)
{
  std::vector<Tidx> rnd(d);
  for (int i=0; i<d; i++)
    rnd[i]=i;
  for (int i=0; i<d; i++) {
    int idx=d-i;
    int res=d-i;
    int k=rand() % res;
    if (k != idx) {
      int tmp=rnd[idx];
      rnd[idx]=rnd[k];
      rnd[k]=tmp;
    }
  }
  return DoubleSidedPerm<Tidx>(rnd);
}



template<typename Tidx>
DoubleSidedPerm<Tidx> Inverse(DoubleSidedPerm<Tidx> const& ePerm)
{
  return ~ePerm;
}



// Input / Output

template<typename Tidx>
std::string GapStyleStringShift(DoubleSidedPerm<Tidx> const& ePerm, int const& eShift)
{
  int n=ePerm.size();
  std::vector<int> ListStat(n,1);
  std::string eRet;

  for (int i=0; i<n; i++) {
    if (ListStat[i] == 1) {
      int eFirst=i;
      int eCurr=i;
      std::string ePart = "(";
      bool IsFirst=true;
      int len=0;
      while(true) {
	if (!IsFirst)
	  ePart += ",";
	IsFirst=false;
	ePart += std::to_string(eCurr + eShift);
	ListStat[eCurr]=0;
	int eNext = ePerm.at(eCurr);
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

template<typename Tidx>
std::string GapStyleString(DoubleSidedPerm<Tidx> const& ePerm)
{
  return GapStyleStringShift(ePerm, 1);
}


template<typename Tidx>
std::ostream& operator<<(std::ostream& os, DoubleSidedPerm<Tidx> const& ePerm)
{
  os << GapStyleStringShift(ePerm,1);
  return os;
}




template<typename Tidx_inp>
struct SingleSidedPerm {
public:
  using Tidx = Tidx_inp;
  //
  // The constructors
  //
  SingleSidedPerm(std::string const& estr)
  {
    std::pair<std::vector<Tidx>, std::vector<Tidx>> epair = GetListValRev<Tidx>(estr);
    ListVal = epair.first;
    siz = ListVal.size();
  }
  SingleSidedPerm(SingleSidedPerm const& ePerm, int const& n)
  {
    if (ePerm.size() > n) {
      std::cerr << "ePerm.size()=" << ePerm.size() << " n=" << n << "\n";
      std::cerr << "ExtendPermutation to a size that is lower than the current size\n";
    }
    ListVal = ePerm.getListVal();
    for (int pos=ePerm.size(); pos<n; pos++)
      ListVal.push_back(pos);
    siz = n;
  }
  SingleSidedPerm ()
  {
    siz=0;
    ListVal = {};
  }
  SingleSidedPerm (int const& n)
  {
    siz=n;
    ListVal = std::vector<Tidx>(n);
    for (int i=0; i<n; i++)
      ListVal[i]=i;
  }
  SingleSidedPerm(std::vector<Tidx> const& v)
  {
    ListVal=v;
    siz=v.size();
  }
  SingleSidedPerm(std::vector<Tidx> const& v1, std::vector<Tidx> const& v2)
  {
    siz=v1.size();
    ListVal=v1;
  }
  SingleSidedPerm(SingleSidedPerm const& ePerm)
  {
    siz     = ePerm.siz;
    ListVal = ePerm.ListVal;
  }
  SingleSidedPerm(SingleSidedPerm&& ePerm)
  {
    siz = ePerm.siz;
    ListVal = std::move(ePerm.ListVal);
    ePerm.siz = 0;
  }
  //
  // Copy operator
  //
  SingleSidedPerm<Tidx> operator=(SingleSidedPerm const& ePerm)
  {
    siz     = ePerm.siz;
    ListVal = ePerm.ListVal;
    return *this;
  }
  SingleSidedPerm<Tidx> operator=(SingleSidedPerm&& ePerm)
  {
    siz = ePerm.siz;
    ListVal = std::move(ePerm.ListVal);
    ePerm.siz = 0;
    return *this;
  }
  //
  // The destructor
  //
  ~SingleSidedPerm()
  {
  }
  //
  // The destructor
  //
  bool isIdentity() const
  {
    for (int i=0; i<siz; i++)
      if (ListVal[i] != i)
	return false;
    return true;
  }
  int at(int const& i) const
  {
    return ListVal[i];
  }
  int atRev(int const& i) const
  {
    for (int j=0; j<siz; j++)
      if (ListVal[j] == i)
        return j;
    return -1;
  }
  std::vector<Tidx> getListVal() const
  {
    return ListVal;
  }
  int operator[](int const& i) const
  {
    return ListVal[i];
  }
  int size() const
  {
    return siz;
  }
  //
private:
  int siz;
  std::vector<Tidx> ListVal;
};


template<typename Tidx>
bool operator==(SingleSidedPerm<Tidx> const& v1, SingleSidedPerm<Tidx> const& v2)
{
  int siz=v1.size();
  if (siz != v2.size() )
    return false;
  for (int i=0; i<siz; i++)
    if (v1.at(i) != v2.at(i))
      return false;
  return true;
}


template<typename Tidx>
bool operator!=(SingleSidedPerm<Tidx> const& v1, SingleSidedPerm<Tidx> const& v2)
{
  int siz=v1.size();
  if (siz != v2.size() )
    return true;
  for (int i=0; i<siz; i++)
    if (v1.at(i) != v2.at(i))
      return true;
  return false;
}


template<typename Tidx>
bool operator<(SingleSidedPerm<Tidx> const& v1, SingleSidedPerm<Tidx> const& v2)
{
  int siz1=v1.size();
  int siz2=v2.size();
  if (siz1 != siz2)
    return siz1<siz2;
  int siz=siz1;
  for (int i=0; i<siz; i++) {
    if (v1.at(i) != v2.at(i))
      return v1.at(i) < v2.at(i);
  }
  return false;
}

template<typename Tidx>
SingleSidedPerm<Tidx> operator~(SingleSidedPerm<Tidx> const& ePerm)
{
  int siz = ePerm.size();
  std::vector<Tidx> LVal = ePerm.getListVal();
  std::vector<Tidx> v(siz);
  for (int i=0; i<siz; i++)
    v[LVal[i]] = i;
  return SingleSidedPerm<Tidx>(v);
}







// Form the product v1 * v2
template<typename Tidx>
SingleSidedPerm<Tidx> operator*(SingleSidedPerm<Tidx> const& v1, SingleSidedPerm<Tidx> const& v2)
{
  int siz=v1.size();
#ifdef DEBUG
  if (siz != v2.size() ) {
    std::cerr << "Error in the DoubleSidedPerm product\n";
    throw PermutalibException{1};
  }
#endif
  std::vector<Tidx> vVal(siz);
  for (int i=0; i<siz; i++) {
    int j=v1.at(i);
    int k=v2.at(j);
    vVal[i]=k;
  }
  return SingleSidedPerm<Tidx>(vVal);
}



template<typename Tidx>
SingleSidedPerm<Tidx> Conjugation(SingleSidedPerm<Tidx> const& v1, SingleSidedPerm<Tidx> const& v2)
{
  int siz=v1.size();
#ifdef DEBUG
  if (siz != v2.size() ) {
    std::cerr << "Error in the DoubleSidedPerm conjugation\n";
    throw PermutalibException{1};
  }
#endif
  std::vector<Tidx> v(siz);
  for (int i=0; i<siz; i++) {
    int j=v1[i];
    int i2=v2[i];
    int j2=v2[j];
    v[i2]=j2;
  }
  return SingleSidedPerm<Tidx>(v);
}



template<typename Tidx>
int PowAct(int const& i, SingleSidedPerm<Tidx> const& g)
{
  return g.at(i);
}



template<typename Tidx>
int SlashAct(int const& i, SingleSidedPerm<Tidx> const& g)
{
  return g.atRev(i);
}


// LeftQuotient(x,y) = x^{-1}*y in the list.gi file
template<typename Tidx>
SingleSidedPerm<Tidx> LeftQuotient(SingleSidedPerm<Tidx> const& a, SingleSidedPerm<Tidx> const& b)
{
  int siz=a.size();
  std::vector<Tidx> ListVal(siz);
  for (int i=0; i<siz; i++) {
    int i1=a.atRev(i);
    int j1=b.at(i1);
    ListVal[i]=j1;
  }
  return SingleSidedPerm<Tidx>(ListVal);
}



template<typename Tidx>
SingleSidedPerm<Tidx> SCRandomPerm(int const& d)
{
  std::vector<Tidx> rnd(d);
  for (int i=0; i<d; i++)
    rnd[i]=i;
  for (int i=0; i<d; i++) {
    int idx=d-i;
    int res=d-i;
    int k=rand() % res;
    if (k != idx) {
      int tmp=rnd[idx];
      rnd[idx]=rnd[k];
      rnd[k]=tmp;
    }
  }
  return DoubleSidedPerm<Tidx>(rnd);
}



template<typename Tidx>
SingleSidedPerm<Tidx> Inverse(SingleSidedPerm<Tidx> const& ePerm)
{
  return ~ePerm;
}







}





#endif
