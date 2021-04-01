#ifndef DEFINE_PERMUTALIB_PERMUTATION_H
#define DEFINE_PERMUTALIB_PERMUTATION_H

#include <vector>
#include <stdlib.h>
#include <string>
#include <iostream>

namespace permutalib {



std::pair<std::vector<int>, std::vector<int>> GetListValRev(std::string_view const& estr)
{
  size_t maxlen = 0;
  std::vector<int> ListVal;
  std::vector<int> ListRev;
  auto insertLVal=[&](std::vector<int> const& LVal) -> void {
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
      int val1 = LVal[i];
      int val2 = LVal[j];
      ListVal[val1] = val2;
      ListRev[val2] = val1;
    }
  };
  auto ParseStringByComma=[&](std::string_view const& estr) -> std::vector<int> {
    size_t n_char=estr.size();
    size_t pos_start = 0;
    std::vector<int> LVal;
    auto insert=[&](size_t const& pos1, size_t const& pos2) -> void {
      size_t len = pos2 - pos1;
      std::string_view ustr = estr.substr(pos_start, len);
      int eVal = std::stoi(std::string(ustr)) - 1;
      LVal.push_back(eVal);
      pos_start = pos2 + 1;
    };
    for (size_t i_char=0; i_char<n_char; i_char++) {
      std::string_view echar = estr.substr(i_char, 1);
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
      std::string_view sstr = estr.substr(pos_start, len);
      std::vector<int> LVal = ParseStringByComma(sstr);
      insertLVal(LVal);
    }
  }
  return {ListVal,ListRev};
}




struct DoubleSidedPerm {
public:
  //
  // The constructors
  //
  DoubleSidedPerm(std::string_view const& estr)
  {
    std::pair<std::vector<int>, std::vector<int>> epair = GetListValRev(estr);
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
    ListVal = std::vector<int>(n);
    ListRev = std::vector<int>(n);
    for (int i=0; i<n; i++) {
      ListVal[i]=i;
      ListRev[i]=i;
    }
  }
  DoubleSidedPerm(std::vector<int> const& v)
  {
    ListVal=v;
    siz=v.size();
    ListRev.resize(siz);
    for (int i=0; i<siz; i++)
      ListRev[v[i]]=i;
  }
  DoubleSidedPerm(std::vector<int> const& v1, std::vector<int> const& v2)
  {
    siz=v1.size();
    ListVal=v1;
    ListRev=v2;
  }
  DoubleSidedPerm(DoubleSidedPerm const& ePerm)
  {
    siz=ePerm.siz;
    ListVal=ePerm.ListVal;
    ListRev=ePerm.ListRev;
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
  DoubleSidedPerm operator=(DoubleSidedPerm const& ePerm)
  {
    siz = ePerm.siz;
    ListVal = ePerm.ListVal;
    ListRev = ePerm.ListRev;
    return *this;
  }
  DoubleSidedPerm operator=(DoubleSidedPerm&& ePerm)
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
  std::vector<int> getListVal() const
  {
    return ListVal;
  }
  std::vector<int> getListRev() const
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
  std::vector<int> ListVal;
  std::vector<int> ListRev;
};


bool operator==(DoubleSidedPerm const& v1, DoubleSidedPerm const& v2)
{
  int siz=v1.size();
  if (siz != v2.size() )
    return false;
  for (int i=0; i<siz; i++)
    if (v1.at(i) != v2.at(i))
      return false;
  return true;
}


bool operator!=(DoubleSidedPerm const& v1, DoubleSidedPerm const& v2)
{
  int siz=v1.size();
  if (siz != v2.size() )
    return true;
  for (int i=0; i<siz; i++)
    if (v1.at(i) != v2.at(i))
      return true;
  return false;
}


bool operator<(DoubleSidedPerm const& v1, DoubleSidedPerm const& v2)
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

DoubleSidedPerm operator~(DoubleSidedPerm const& ePerm)
{
  return DoubleSidedPerm(ePerm.getListRev(), ePerm.getListVal());
}







// Form the product v1 * v2
DoubleSidedPerm operator*(DoubleSidedPerm const& v1, DoubleSidedPerm const& v2)
{
  int siz=v1.size();
#ifdef DEBUG
  if (siz != v2.size() ) {
    std::cerr << "Error in the DoubleSidedPerm product\n";
    throw PermutalibException{1};
  }
#endif
  std::vector<int> vVal(siz), vRev(siz);
  for (int i=0; i<siz; i++) {
    int j=v1.at(i);
    int k=v2.at(j);
    vVal[i]=k;
    //
    int j2=v2.atRev(i);
    int k2=v1.atRev(j2);
    vRev[i]=k2;
  }
  return DoubleSidedPerm(vVal, vRev);
}

DoubleSidedPerm Conjugation(DoubleSidedPerm const& v1, DoubleSidedPerm const& v2)
{
  int siz=v1.size();
#ifdef DEBUG
  if (siz != v2.size() ) {
    std::cerr << "Error in the DoubleSidedPerm conjugation\n";
    throw PermutalibException{1};
  }
#endif
  std::vector<int> v(siz);
  for (int i=0; i<siz; i++) {
    int j=v1[i];
    int i2=v2[i];
    int j2=v2[j];
    v[i2]=j2;
  }
  return DoubleSidedPerm(v);
}

int PowAct(int const& i, DoubleSidedPerm const& g)
{
  return g.at(i);
}

int SlashAct(int const& i, DoubleSidedPerm const& g)
{
  return g.atRev(i);
}


// LeftQuotient(x,y) = x^{-1}*y in the list.gi file
DoubleSidedPerm LeftQuotient(DoubleSidedPerm const& a, DoubleSidedPerm const& b)
{
  int siz=a.size();
  std::vector<int> ListVal(siz), ListRev(siz);
  for (int i=0; i<siz; i++) {
    int i1=a.atRev(i);
    int j1=b.at(i1);
    ListVal[i]=j1;
    int i2=b.atRev(i);
    int j2=a.at(i2);
    ListRev[i]=j2;
  }
  return DoubleSidedPerm(ListVal, ListRev);
}


DoubleSidedPerm SCRandomPerm(int const& d)
{
  std::vector<int> rnd(d);
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
  return DoubleSidedPerm(rnd);
}

DoubleSidedPerm Inverse(DoubleSidedPerm const& ePerm)
{
  return ~ePerm;
}



// Input / Output

std::string GapStyleStringShift(DoubleSidedPerm const& ePerm, int const& eShift)
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

std::string GapStyleString(DoubleSidedPerm const& ePerm)
{
  return GapStyleStringShift(ePerm, 1);
}


std::ostream& operator<<(std::ostream& os, DoubleSidedPerm const& ePerm)
{
  os << GapStyleStringShift(ePerm,1);
  return os;
}






}





#endif
