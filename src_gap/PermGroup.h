#ifndef DEFINE_PERMUTALIB_PERM_GROUP_H
#define DEFINE_PERMUTALIB_PERM_GROUP_H

#include "Face_basic.h"
#include <unordered_set>
#include <set>

namespace permutalib {

//
// Actions
//




template<typename Tidx, typename Telt>
std::vector<Tidx> OnSets(std::vector<Tidx> const& V, Telt const& u)
{
  std::vector<Tidx> Vret;
  Vret.reserve(V.size());
  for (auto & ePt : V)
    Vret.push_back(u.at(ePt));
  sort(Vret.begin(), Vret.end());
  return Vret;
}

//
// One single element being used
//



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
Telt PowerGroupElement(Telt const& u, int const& n)
{
  if (n < 0) {
    std::cerr << "We should have n >= 1\n";
    throw PermutalibException{1};
  }
  Telt pow = u;
  for (int i=1; i<n; i++)
    pow *= u;
  return pow;
}

template<typename Telt>
typename Telt::Tidx OrderElement(const Telt& x)
{
  using Tidx = typename Telt::Tidx;
  if (x.isIdentity())
    return 1;
  Telt xw = x;
  Tidx ord = 1;
  while(true) {
    xw *= x;
    ord++;
    if (xw.isIdentity())
      return ord;
  }
}


/* Return:
   --- <Tidx>::max() if no such power exist
   --- e if a^e = b.
*/
template<typename Telt>
typename Telt::Tidx LogPerm(const Telt& a, const Telt& b)
{
  using Tidx = typename Telt::Tidx;
  if (a.isIdentity()) {
    if (b.isIdentity())
      return 1;
    else
      return std::numeric_limits<Tidx>::max();
  }
  Telt apow = a;
  Tidx ord = 1;
  while(true) {
    apow *= a;
    ord++;
    if (apow == b)
      return ord;
    if (apow.isIdentity())
      return std::numeric_limits<Tidx>::max();
  }
}


template<typename Telt>
Telt ElementPower(const Telt& x, const typename Telt::Tidx& n) {
  using Tidx=typename Telt::Tidx;
  std::function<Telt(const Tidx& n)> pow=[&](const Tidx& u) -> Telt {
    if (u == 0)
      return Telt(x.size());
    Tidx res = u % 2;
    if (res == 0) {
      Telt xp = pow(u / 2);
      return xp * xp;
    }
    Telt xp = pow(u / 2);
    return x * (xp * xp);
  };
  return pow(n);
}


//
// vector<Telt> taking functions
//

template<typename Telt>
typename Telt::Tidx SmallestMovedPoint(std::vector<Telt> const& LGen)
{
  using Tidx = typename Telt::Tidx;
  if (LGen.size() == 0)
    return std::numeric_limits<Tidx>::max();
  Tidx n=LGen[0].size();
  for (Tidx u=0; u<n; u++) {
    for (auto & eGen : LGen)
      if (eGen.at(u) != u)
        return u;
  }
  return std::numeric_limits<Tidx>::max();
}


template<typename Telt>
std::pair<std::vector<typename Telt::Tidx>,Face> OrbitPerms(std::vector<Telt> const& gens, typename Telt::Tidx const& n, typename Telt::Tidx const& d)
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
  return {std::move(orb), std::move(eFace)};
}


template<typename Telt>
std::vector<std::vector<typename Telt::Tidx>> Kernel_OrbitsPerms(const std::vector<Telt>& gens, Face & dom)
{
  using Tidx=typename Telt::Tidx;
  std::vector<std::vector<Tidx>> orbs;
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
    orbs.emplace_back(std::move(orb));
  }
  return orbs;
}

template<typename Telt>
std::vector<std::vector<typename Telt::Tidx>> OrbitsPerms(const std::vector<Telt>& gens, const typename Telt::Tidx&n, const std::vector<typename Telt::Tidx>& D)
{
  Face dom(n);
  for (auto & eV : D)
    dom[eV]=1;
  return Kernel_OrbitsPerms(gens, dom);
}

template<typename Telt>
std::vector<std::vector<typename Telt::Tidx>> OrbitsPerms(const std::vector<Telt>& gens, const typename Telt::Tidx&n)
{
  using Tidx=typename Telt::Tidx;
  Face dom(n);
  for (Tidx i=0; i<n; i++)
    dom[i] = 1;
  return Kernel_OrbitsPerms(gens, dom);
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
std::vector<typename Telt::Tidx> MovedPointsPerms(std::vector<Telt> const& gens)
{
  using Tidx=typename Telt::Tidx;
  if (gens.size() == 0)
    return {};
  std::vector<Tidx> ListMoved;
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

template<typename Telt, typename Tobj, typename Tact>
std::vector<Tobj> Orbit(std::vector<Telt> const& ListGen, Tobj const& x, Tact act)
{
  std::vector<Tobj> ListObj{x};
  std::unordered_set<Tobj> SetObj{x};
  size_t curr_pos = 0;
  while(true) {
    size_t len=ListObj.size();
    if (curr_pos == len)
      break;
    for (size_t u=curr_pos; u<len; u++) {
      for (auto & eElt : ListGen) {
        Tobj eImg = act(ListObj[u], eElt);
        if (SetObj.count(eImg) == 0) {
          ListObj.push_back(eImg);
          SetObj.insert(eImg);
        }
      }
    }
    curr_pos = len;
  }
  return ListObj;
}



template<typename Telt, typename Tobj, typename Tact>
std::vector<std::pair<Tobj,Telt>> OrbitPairEltRepr(std::vector<Telt> const& ListGen, Telt const& id, Tobj const& x, Tact act)
{
  std::vector<std::pair<Tobj,Telt>> ListPair{ {x,id} };
  std::unordered_set<Tobj> SetObj{x};
  size_t curr_pos = 0;
  while(true) {
    size_t len=ListPair.size();
    if (curr_pos == len)
      break;
    for (size_t u=curr_pos; u<len; u++) {
      for (auto & eElt : ListGen) {
        Tobj eImg = act(ListPair[u].first, eElt);
        if (SetObj.count(eImg) == 0) {
          Telt eProd = ListPair[u].second * eElt;
          ListPair.push_back({eImg,eProd});
          SetObj.insert(eImg);
        }
      }
    }
    curr_pos = len;
  }
  return ListPair;
}





template<typename Telt>
std::vector<typename Telt::Tidx> MovedPoints(const std::vector<Telt>& LGen, const typename Telt::Tidx& n)
{
  using Tidx = typename Telt::Tidx;
  auto IsMoved=[&](Tidx const& ePt) -> bool {
    for (auto & eGen : LGen)
      if (eGen.at(ePt) != ePt)
        return true;
    return false;
  };
  std::vector<Tidx> LMoved;
  LMoved.reserve(n);
  for (Tidx i=0; i<n; i++)
    if (IsMoved(i))
      LMoved.push_back(i);
  return LMoved;
}




}

#endif
