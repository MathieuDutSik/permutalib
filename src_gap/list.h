#ifndef DEFINE_PERMUTALIB_LIST_H
#define DEFINE_PERMUTALIB_LIST_H


#include <vector>
#include "Face_basic.h"
#include "COMB_Vectors.h"

namespace permutalib {

template<typename T>
std::vector<T> Reversed(std::vector<T> const& eList)
{
  size_t len=eList.size();
  std::vector<T> retList(len);
  for (size_t u=0; u<len; u++) {
    size_t pos=len - 1 - u;
    retList[pos] = eList[u];
  }
  return retList;
}



template<typename Tidx>
std::vector<Tidx> ClosedInterval(Tidx const& begin, Tidx const& end)
{
  size_t len = end - begin;
  std::vector<Tidx> clInt(len);
  for (size_t u=0; u<len; u++)
    clInt[u] = Tidx(u + begin);
  return clInt;
}


template<typename Tidx>
Face BlistList(std::vector<Tidx> const& list, std::vector<Tidx> const& sub)
{
  size_t len=list.size();
#ifdef SYNCHRONIZED_DEBUG_GAP478
  std::cerr << "DEBUG list=" << GapStringTVectorB(list) << "\n";
  std::cerr << "DEBUG sub=" << GapStringTVectorB(sub) << "\n";
#endif
  Face ret(len);
  for (auto & eVal : sub) {
    size_t pos=PositionVect_ui<Tidx,size_t>(list, eVal);
    ret[pos]=true;
  }
  return ret;
}




template<typename Tidx>
Face BlistList_direct(size_t const& len, std::vector<Tidx> const& sub)
{
  Face ret(len);
  for (auto & eVal : sub)
    ret[eVal]=true;
  return ret;
}


boost::dynamic_bitset<>::size_type PositionNthTrueBlist(Face const& blist, size_t const& hpos)
{
  boost::dynamic_bitset<>::size_type b=blist.find_first();
  size_t epos=0;
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
  size_t len=list.size();
  for (size_t i=0; i<len; i++)
    if (blist[i] == 1)
      ret.push_back(list[i]);
  return ret;
}

size_t SizeBlist(Face const& blist)
{
  return blist.count();
}

bool IsSubsetBlist(Face const& a, Face const& b)
{
  boost::dynamic_bitset<>::size_type pos=b.find_first();
  while (pos != boost::dynamic_bitset<>::npos) {
    if (a[pos] == 0)
      return false;
    pos = b.find_next(pos);
  }
  return true;
}



void UniteBlist(Face & a, Face const& b)
{
  boost::dynamic_bitset<>::size_type pos=b.find_first();
  while (pos != boost::dynamic_bitset<>::npos) {
    a[pos] = 1;
    pos = b.find_next(pos);
  }
}



void IntersectBlist(Face & a, Face const& b)
{
  a &= b;
}


void SubtractBlist(Face & a, Face const& b)
{
  boost::dynamic_bitset<>::size_type pos=b.find_first();
  while (pos != boost::dynamic_bitset<>::npos) {
    a[pos] = 0;
    pos = b.find_next(pos);
  }
}


}
#endif
