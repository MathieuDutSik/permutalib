#ifndef DEFINE_PERMUTALIB_BLOCK_SYSTEMS_H
#define DEFINE_PERMUTALIB_BLOCK_SYSTEMS_H

#include "Face_basic.h"
#include "GraphicFunctionality.h"
#include "PermGroup.h"


namespace permutalib {

template<typename Telt>
std::pair<std::vector<std::vector<size_t>>,std::vector<Face>> Blocks_Kernel(std::vector<Telt> const& ListGen, std::vector<int> const& Omega, std::vector<int> const& seed)
{
  using Tidx=typename Telt::Tidx;
  int nbMax=VectorMax(Omega);
  size_t n=Omega.size();
  std::vector<Tidx> OmegaRev(nbMax);
  for (size_t u=0; u<n; u++) {
    size_t ePt=Omega[u];
    OmegaRev[ePt] = Tidx(u);
  }
  std::vector<Telt> ListGenRed;
  for (auto & eGen : ListGen) {
    std::vector<int> eList(n);
    for (size_t u=0; u<n; u++) {
      Tidx ePt=Omega[u];
      Tidx ePtImg=eGen.at(ePt);
      Tidx uImg=OmegaRev[ePtImg];
      eList[u] = uImg;
    }
    Telt eGenRed(eList);
    ListGenRed.push_back(eGenRed);
  }
  // Building the orbit
  Face eFace(n);
  for (auto & ePt : seed) {
    Tidx u=OmegaRev[ePt];
    eFace[u]=1;
  }
  std::function<Face(Face const&,Telt const&)> act=[&](Face const& x, Telt const& eElt) -> Face {
    Face eImg(n);
    int siz=x.count();
    boost::dynamic_bitset<>::size_type ePt=x.find_first();
    for (int u=0; u<siz; u++) {
      Tidx ePtImg=eElt.at(ePt);
      eImg[ePtImg]=1;
      ePt = x.find_next(ePt);
    }
    return eImg;
  };
  std::vector<Face> SeedOrbit = Orbit(ListGenRed, eFace, act);
  // Building the adjacency graph
  auto IsIntersecting=[&](Face const& face1, Face const& face2) -> bool {
    int siz=face1.count();
    boost::dynamic_bitset<>::size_type ePt=face1.find_first();
    for (int u=0; u<siz; u++) {
      if (face2[ePt] == 1)
	return true;
      ePt = face1.find_next(ePt);
    }
    return false;
  };
  size_t sizOrb=SeedOrbit.size();
  std::vector<size_t> ListEdge;
  for (size_t u=0; u<sizOrb; u++)
    for (size_t v=u+1; v<sizOrb; v++) {
      if (IsIntersecting(SeedOrbit[u], SeedOrbit[v])) {
	ListEdge.push_back(u);
	ListEdge.push_back(v);
      }
    }
  GraphSparseImmutable eGR(ListEdge, sizOrb);
  std::vector<std::vector<size_t>> ListConn = ConnectedComponents_set(eGR);
  std::vector<std::vector<size_t>> ListBlocks;
  for (auto & eConn : ListConn) {
    Face eFaceComb(n);
    for (auto & iElt : eConn) {
      Face eFace = SeedOrbit[iElt];
      size_t siz=eFace.count();
      boost::dynamic_bitset<>::size_type ePt=eFace.find_first();
      for (size_t u=0; u<siz; u++) {
	eFaceComb[ePt]=1;
	ePt = eFace.find_next(ePt);
      }
    }
    std::vector<size_t> eBlock;
    size_t blkSiz=eFaceComb.count();
    boost::dynamic_bitset<>::size_type ePt=eFaceComb.find_first();
    for (size_t u=0; u<blkSiz; u++) {
      eBlock.push_back(Omega[ePt]);
      ePt = eFaceComb.find_next(ePt);
    }
    ListBlocks.push_back(eBlock);
  }
  return {std::move(ListBlocks),std::move(SeedOrbit)};
}


template<typename Telt>
std::vector<std::vector<int>> Blocks_from_seed(std::vector<Telt> const& ListGen, std::vector<int> const& Omega, std::vector<int> const& seed)
{
  return Blocks_Kernel(ListGen, Omega, seed).first;
}



// Maybe we can do better by using std::vector for the edges, but
// that is a secondary optimization
template<typename Telt>
std::vector<std::vector<int>> Blocks_without_seed(std::vector<Telt> const& ListGen, std::vector<int> const& Omega)
{
  size_t n=Omega.size();
  size_t nbTotalEdge=n * (n-1) / 2;
  Face UsedEdge(nbTotalEdge);
  size_t ThePos=0;
  auto UpdatePosition=[&]() -> bool {
    while(true) {
      if (UsedEdge[ThePos] == 0)
	return true;
      ThePos++;
      if (ThePos == nbTotalEdge)
	break;
    }
    return false;
  };
  auto GetEdgeFromPosition=[&](int const& pos) -> std::vector<int> {
    int posTot=pos;
    int firstElement=0;
    while(true) {
      int allowedSiz = n - 1 - firstElement;
      if (posTot < allowedSiz) {
	int secondElement=firstElement + 1 + posTot;
	return {firstElement,secondElement};
      }
      posTot -= allowedSiz;
      firstElement++;
    }
  };
  auto GetPosFromEdge=[&](std::vector<int> const& eEdge) -> int {
    int eFirst=eEdge[0];
    int eSecond=eEdge[1];
    int pos=0;
    for (int u=0; u<eFirst; u++) {
      int allowedSiz = n - 1 - eFirst;
      pos += allowedSiz;
    }
    return pos + eSecond - 1 - eFirst;
  };
  while(true) {
    std::vector<int> eEdge=GetEdgeFromPosition(ThePos);
    std::pair<std::vector<std::vector<int>>,std::vector<Face>> ePair = Blocks_Kernel(ListGen, Omega, eEdge);
    if (ePair.first.size() > 1) {
      return ePair.first;
    }
    for (auto & eEdgeFace : ePair.second) {
      std::vector<int> eEdge = FaceToVector<int>(eEdgeFace);
      int pos=GetPosFromEdge(eEdge);
      UsedEdge[pos]=1;
    }
    bool test = UpdatePosition();
    if (!test)
      break;
  }
  return {Omega};
}

template<typename Telt>
std::function<Telt(Telt const&)> MapElementToSetPlusBlocks(std::vector<std::vector<int>> const& blks, int const& n)
{
  using Tidx=typename Telt::Tidx;
  int TheMax=0;
  for (auto & eBlock : blks)
    for (auto & ePt : eBlock)
      if (ePt > TheMax)
	TheMax = ePt;
  std::vector<int> VectStatus(TheMax);
  size_t nbBlock=blks.size();
  for (size_t iBlk=0; iBlk<nbBlock; iBlk++)
    for (auto & ePt : blks[iBlk])
      VectStatus[ePt] = iBlk;
  std::function<Telt(Telt const&)> f=[=](Telt const& u) -> Telt {
    std::vector<Tidx> eList(n + nbBlock);
    for (size_t i=0; i<n; i++)
      eList[i] = u.at(i);
    for (size_t iBlk=0; iBlk<nbBlock; iBlk++) {
      int ePt=blks[iBlk][0];
      Tidx ePtImg=u.at(ePt);
      size_t iBlkImg=VectStatus[ePtImg];
      eList[n+iBlk] = Tidx(n+iBlkImg);
    }
    return Telt(eList);
  };
  return f;
}

std::vector<int> GetBlock(std::vector<std::vector<int>> const& ListBlock, int const& ePt)
{
  for (auto & eBlock : ListBlock)
    if (PositionVect(eBlock, ePt) != -1)
      return eBlock;
  return {-1};
}


int GetNonTrivialPointInBlock(std::vector<int> const& eBlock, int const& ePt)
{
  for (auto & uPt : eBlock)
    if (uPt != ePt)
      return uPt;
  return -1;
}

}

#endif

