#ifndef DEFINE_PERMUTALIB_LIST_H
#define DEFINE_PERMUTALIB_LIST_H


#include <vector>
#include "Face_basic.h"
#include "COMB_Vectors.h"

namespace permutalib {

template<typename T>
std::vector<T> Reversed(std::vector<T> const& eList)
{
  int len=eList.size();
  std::vector<T> retList(len);
  for (int u=0; u<len; u++) {
    int pos=len - 1 - u;
    retList[pos] = eList[u];
  }
  return retList;
}



template<typename Tidx>
std::vector<Tidx> ClosedInterval(Tidx const& begin, Tidx const& end)
{
  std::vector<Tidx> clInt;
  for (Tidx u=begin; u<end; u++)
    clInt.push_back(u);
  return clInt;
}


template<typename Tidx>
Face BlistList(std::vector<Tidx> const& list, std::vector<Tidx> const& sub)
{
  int len=list.size();
#ifdef SYNCHRONIZED_DEBUG_GAP478
  std::cerr << "DEBUG list=" << GapStringTVectorB(list) << "\n";
  std::cerr << "DEBUG sub=" << GapStringTVectorB(sub) << "\n";
#endif
  Face ret(len);
  for (auto & eVal : sub) {
    //    std::cerr << "eVal=" << eVal << "\n";
    int pos=PositionVect(list, eVal);
    //    std::cerr << "pos=" << pos << "\n";
    ret[pos]=true;
  }
  return ret;
}



boost::dynamic_bitset<>::size_type PositionNthTrueBlist(Face const& blist, int const& hpos)
{
  boost::dynamic_bitset<>::size_type b=blist.find_first();
  int epos=0;
  while (true) {
    if (b == boost::dynamic_bitset<>::npos)
      return b;
    if (epos == hpos)
      return b;
    // We increase
    b=blist.find_next(b);
    epos++;
  }
}



template<typename Tidx>
std::vector<Tidx> ListBlist(std::vector<Tidx> const& list, std::vector<int8_t> const& blist)
{
  std::vector<Tidx> ret;
  Tidx len=list.size();
  for (Tidx i=0; i<len; i++) {
    if (blist[i] == 1)
      ret.push_back(list[i]);
  }
  return ret;
}

int SizeBlist(Face const& blist)
{
  return blist.count();
}

bool IsSubsetBlist(Face const& a, Face const& b)
{
  int siz=a.size();
  for (int i=0; i<siz; i++) {
    if (b[i] == 1 && a[i] == 0)
      return false;
  }
  return true;
}



void UniteBlist(Face & a, Face const& b)
{
  int siz=a.size();
  for (int i=0; i<siz; i++) {
    if (b[i] == 1)
      a[i]=1;
  }
}



void IntersectBlist(Face & a, Face const& b)
{
  int siz=a.size();
  for (int i=0; i<siz; i++) {
    int val=0;
    if (a[i] == 1 && b[i] == 1)
      val=1;
    a[i]=val;
  }
}


void SubtractBlist(Face & a, Face const& b)
{
  int siz=a.size();
  for (int i=0; i<siz; i++) {
    if (b[i] == 1)
      a[i]=0;
  }
}


}
#endif
