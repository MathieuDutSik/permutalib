#ifndef DEFINE_PERMUTALIB_PERM_GROUP_H
#define DEFINE_PERMUTALIB_PERM_GROUP_H

#include "Face_basic.h"
#include <set>

namespace permutalib {

template<typename Telt>
typename Telt::Tidx SmallestMovedPoint(std::vector<Telt> const& LGen)
{
  using Tidx = typename Telt::Tidx;
  if (LGen.size() == 0)
    return std::numeric_limits<Tidx>::max();
  Tidx n=LGen[0].size();
  for (Tidx u=0; u<n; u++) {
    bool IsOK=false;
    for (auto & eGen : LGen)
      if (eGen.at(u) != u)
	IsOK=true;
    if (IsOK)
      return u;
  }
  return std::numeric_limits<Tidx>::max();
}


template<typename Telt>
std::vector<typename Telt::Tidx> OrbitPerms(std::vector<Telt> const& gens, typename Telt::Tidx const& n, typename Telt::Tidx const& d)
{
  using Tidx=typename Telt::Tidx;
  std::vector<Tidx> orb;
  Face eFace(n);
  auto InsertValue=[&](Tidx const& val) -> void {
    orb.push_back(val);
    eFace[val]=1;
  };
  InsertValue(d);
  size_t posDone=0;
  while(true) {
    size_t posTot=orb.size();
    if (posTot == posDone)
      break;
    for (size_t u=posDone; u<posTot; u++) {
      Tidx pnt=orb[u];
      for (auto & eGen : gens) {
	Tidx img=PowAct(pnt, eGen);
	if (eFace[img] == 0)
	  InsertValue(img);
      }
    }
    posDone = posTot;
  }
  return orb;
}


template<typename Telt>
std::vector<std::vector<typename Telt::Tidx>> OrbitsPerms(std::vector<Telt> const& gens, typename Telt::Tidx const&n, std::vector<typename Telt::Tidx> const& D)
{
  using Tidx=typename Telt::Tidx;
  std::vector<std::vector<Tidx>> orbs;
  Face dom(n);
  for (auto & eV : D)
    dom[eV]=1;
  while(true) {
    if (dom.count() == 0)
      break;
    std::vector<Tidx> orb;
    auto insert=[&](Tidx const& eV) -> void {
      orb.push_back(eV);
      dom[eV]=0;
    };
    boost::dynamic_bitset<>::size_type fst=dom.find_first();
    insert(fst);
    size_t posDone=0;
    while(true) {
      size_t posTot=orb.size();
      if (posTot == posDone)
	break;
      for (size_t u=posDone; u<posTot; u++) {
	Tidx pnt=orb[u];
	for (auto & eGen : gens) {
	  Tidx img=PowAct(pnt, eGen);
	  if (dom[img] == 1)
	    insert(img);
	}
      }
      posDone = posTot;
    }
    orbs.push_back(orb);
  }
  return orbs;
}

template<typename Telt>
typename Telt::Tidx SmallestMovedPointsPerms(std::vector<Telt> const& gens)
{
  using Tidx=typename Telt::Tidx;
  Tidx siz=0;
  for (Tidx i=0; i<siz; i++) {
    for (auto & eGen : gens)
      if (PowAct(i, eGen) != i)
	return i;
  }
  return std::numeric_limits<Tidx>::max();
}



template<typename Telt>
std::vector<int> MovedPointsPerms(std::vector<Telt> const& gens)
{
  using Tidx=typename Telt::Tidx;
  if (gens.size() == 0)
    return {};
  std::vector<int> ListMoved;
  Tidx siz=Tidx(gens[0].size());
  for (Tidx i=0; i<siz; i++) {
    auto IsMoved=[&](int const& ePt) -> bool {
      for (auto & eGen : gens)
	if (PowAct(ePt, eGen) != ePt)
	  return true;
      return false;
    };
    if (IsMoved(i))
      ListMoved.push_back(i);
  }
  return ListMoved;
}




template<typename Telt>
Telt RestrictedPermNC(Telt const& x, std::vector<int> const& listRes)
{
  using Tidx=typename Telt::Tidx;
  int n=x.size();
  std::vector<Tidx> MapRev(n);
  size_t nbRes=listRes.size();
  for (size_t iRes=0; iRes<nbRes; iRes++) {
    Tidx ePt=listRes[iRes];
    MapRev[ePt] = iRes;
  }
  std::vector<Tidx> eList(nbRes);
  for (size_t iRes=0; iRes<nbRes; iRes++) {
    Tidx ePt=listRes[iRes];
    Tidx ePtImg=x.at(ePt);
    Tidx iResImg=MapRev[ePtImg];
    eList[iRes]=iResImg;
  }
  Telt eRet(eList);
  return eRet;
}

template<typename Telt, typename Tobj>
std::vector<Tobj> Orbit(std::vector<Telt> const& ListGen, Tobj const& x, std::function<Tobj(Tobj const&,Telt const&)> const& act)
{
  std::vector<Tobj> ListObj{x};
  std::vector<uint8_t> ListStat(0);
  std::set<Tobj> SetObj{x};
  while(true) {
    bool IsFinished=true;
    size_t len=ListObj.size();
    for (size_t u=0; u<len; u++) {
      if (ListStat[u] == 0) {
	IsFinished=false;
	ListStat[u]=1;
	for (auto & eElt : ListGen) {
	  Tobj eImg = act(ListObj[u], eElt);
	  if (SetObj.find(eImg) == SetObj.end()) {
	    ListObj.push_back(eImg);
	    ListStat.push_back(0);
	    SetObj.insert(eImg);
	  }
	}
      }
    }
    if (IsFinished)
      break;
  }
  return ListObj;
}

/*
template<typename Telt>
std::vector<int> OrbitPerms(std::vector<Telt> const& ListGen, int const & ePt)
{
  auto act=[](Telt const& u, int const& eVal) -> int {
    return u.at(eVal);
  };
  return Orbit(ListGen, ePt, act);
}
*/


template<typename Telt>
int CycleLength(Telt const& u, int const& x)
{
  int CycleLen=0;
  int xFirst=x;
  int xWork=x;
  while(true) {
    int xImg=u.at(xWork);
    CycleLen++;
    if (xImg == xFirst)
      break;
    xWork = xImg;
  }
  return CycleLen;
}


template<typename Telt>
std::vector<int> OnSets(std::vector<int> const& V, Telt const& u)
{
  std::vector<int> Vret;
  for (auto & ePt : V)
    Vret.push_back(u.at(ePt));
  sort(Vret.begin(), Vret.end());
  return Vret;
}

template<typename Telt>
Telt PowerGroupElement(Telt const& u, int const& n)
{
  if (n <= 0) {
    std::cerr << "We should have n >= 1\n";
    throw PermutalibException{1};
  }
  Telt pow = u;
  for (int i=1; i<n; i++)
    pow *= u;
  return pow;
}



}

#endif
