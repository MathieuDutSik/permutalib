#ifndef DEFINE_PERMUTALIB_BLOCK_SYSTEMS_H
#define DEFINE_PERMUTALIB_BLOCK_SYSTEMS_H

#include "Face_basic.h"
#include "GraphicFunctionality.h"
#include "PermGroup.h"


namespace permutalib {


template<typename Telt>
std::vector<std::vector<typename Telt::Tidx>> Blocks(const std::vector<Telt>& acts, const typename Telt::Tidx& n)
{
  using Tidx = typename Telt::Tidx;
  if (n == 1) {
    std::vector<Tidx> V = {0};
    return {V};
  }
  if (acts.size() == 0) {
    std::cerr << "We need at least one generator for Blocks computation\n";
    throw PermutalibException{1};
  }
  if (IsPrime(n)) {
    std::vector<Tidx> eVect(n);;
    for (Tidx i=0; i<n; i++)
      eVect[i] = i;
    return {eVect};
  }
  //
  // Now no easy case. Need to do a complex computation
  //
  Telt one(n);
  std::vector<Tidx> orbit{0};
  Face status(n);
  std::vector<Telt> trans(n, Telt());
  trans[0] = one;
  while(true) {
    size_t len = orbit.size();
    bool is_finished = true;
    for (size_t i=0; i<len; i++) {
      if (status[i] == 0) {
        status[i] = 1;
        Tidx pnt = orbit[i];
        for (auto & gen : acts) {
          Tidx pnt_s_gen = SlashAct(pnt, gen);
          if (trans[pnt_s_gen].size() == 0) {
            orbit.push_back(pnt_s_gen);
            trans[pnt_s_gen] = gen;
            is_finished = false;
          }
        }
      }
    }
    if (is_finished)
      break;
  }

  // check that the group is transitive
  if (Tidx(orbit.size()) != n) {
    std::cerr << "G must operate transitively\n";
    throw PermutalibException{1};
  }

  // since $i \in k^{G_1}$ implies $\beta(i)=\beta(k)$,  we initialize <eql>
  // so that the connected components are orbits of some subgroup  $H < G_1$
  std::vector<Tidx> eql(n);
  std::vector<Tidx> leq(n);
  std::vector<Tidx> next(n);
  std::vector<Tidx> last(n);
  for (Tidx i=0; i<n; i++) {
    eql[i] = i;
    leq[i] = i;
    next[i] = 0; // not completely sure
    last[i] = i;
  }

  // repeat until we have a block system
  int changed = 0;
  Tidx cur = orbit[1];
  Telt rnd = one;
  std::vector<std::vector<Tidx>> blocks;
  while(true) {
    std::cerr << "Passing by the while loop changed=" << changed << "\n";
    // compute such an $H$ by taking random  Schreier generators  of $G_1$
    // and stop if 2 successive generators dont change the orbits any more
    while (changed < 2) {

      // compute a random Schreier generator of $G_1$
      size_t i_siz = orbit.size();
      while (1 <= i_siz) {
        rnd *= Random(acts);
        i_siz = i_siz / 2;
      }
      Telt gen = rnd;
      std::cerr << "gen=" << gen << "\n";
      Tidx d1g = PowAct(Tidx(0), gen);
      while (d1g != 0) {
        Telt tr = trans[d1g];
        gen *= tr;
        d1g = PowAct(d1g, tr);
      }
      changed++;

      // compute the image of every point under <gen>
      for (auto & pre_pnt : orbit) { // not sure about this
        Tidx pnt = pre_pnt;
        Tidx img = PowAct(pnt, gen);

        // find the representative of the orbit of <pnt>
        while (eql[pnt] != pnt) {
          pnt = eql[pnt];
        }
        // find the representative of the orbit of <img>
        while (eql[img] != img) {
          img = eql[img];
        }

        // if the do not agree merge their orbits
        if (pnt < img) {
          eql[img] = pnt;
          next[ last[pnt] ] = img;
          last[pnt] = last[img];
          changed = 0;
        } else {
          if (img < pnt) {
            eql[pnt] = img;
            next[ last[img] ] = pnt;
            last[img] = last[pnt];
            changed = 0;
          }
        }
      }
      std::cerr << "changed=" << changed << "\n";
    }

    // take arbitrary point <cur>,  and an element <gen> taking 1 to <cur>
    while (eql[cur] != cur) {
      cur = eql[cur];
    }
    std::vector<Telt> gen_list;
    Tidx img = cur;
    while (img != 0) {
      gen_list.push_back(trans[img]);
      img = PowAct(img, trans[img]);
    }
    gen_list = Reversed(gen_list);
    std::cerr << "|gen_list|=" << gen_list.size() << "\n";

    // compute an alleged block as orbit of 1 under $< H, gen >$
    Tidx pnt = cur;
    std::cerr << "cur=" << cur << "\n";
    while (pnt != 0) {

      // compute the representative of the block containing the image
      img = pnt;
      for (auto & e_gen : gen_list) {
        img = SlashAct(img, e_gen);
      }
      while (eql[img] != img) {
        img = eql[img];
      }

      // if it is not our current block but a minimal block
      if (img != 0 && img != cur && leq[img] == img) {

        // then try <img> as a new start
        leq[cur] = img;
        cur = img;
        gen_list.clear();
        img = cur;
        while (img != 0) {
          gen_list.push_back(trans[img]);
          img = PowAct(img, trans[img]);
        }
        gen_list = Reversed( gen_list );
        pnt = cur;


        // otherwise if it is not our current block but contains it
        // by construction a nonminimal block contains the current block
      } else {
        if (img != 0 && img != cur && leq[img] != img) {

          // then merge all blocks it contains with <cur>
          while (img != cur) {
            eql[img] = cur;
            next[ last[cur] ] = img;
            last[ cur ] = last[ img ];
            img = leq[img];
            while (img != eql[img]) {
              img = eql[img];
            }
          }
          pnt = next[pnt];

          // go on to the next point in the orbit
        } else {
          pnt = next[pnt];
        }
      }
    }

    // make the alleged block
    std::vector<Tidx> block{0};
    pnt = cur;
    while (pnt != 0) {
      block.push_back(pnt);
      pnt = next[pnt];
    }
    block = SortVector(block);
    blocks = {block};

    // quick test to see if the group is primitive
    if (block.size() == orbit.size()) {
      return blocks;
    }

    std::cerr << "|block|=" << block.size() << "\n";
    // quick test to see if the orbit can be a block
    if (orbit.size() % block.size() != 0) {
      changed = -1000;
    }

    // '<rep>[<i>]' is the representative of the block containing <i>
    Tidx miss_val = std::numeric_limits<Tidx>::max();
    std::vector<Tidx> rep(n,miss_val);
    for (auto & pnt : block)
      rep[pnt] = 0;

    // compute the block system with an orbit algorithm
    int i = 0;
    while (0 <= changed &&  i < int(blocks.size())) {

      // loop over the generators
      for (auto & gen : acts) {

        // compute the image of the block under the generator
        std::vector<Tidx> img = OnSets(blocks[i], gen);

        // if this block is new
        if (rep[ img[0] ] == miss_val) {

          // add the new block to the list of blocks
          blocks.push_back(img);

          // check that all points in the image are new
          for (auto & pnt : img) {
            if (rep[pnt] != miss_val) {
              changed = -1000;
            }
            rep[pnt] = img[0];
          }

          // if this block is old
        } else {

          // check that all points in the image lie in the block
          for (auto & pnt : img) {
            if (rep[pnt] != rep[img[0]]) {
              changed = -1000;
            }
          }
        }
      }

      // on to the next block in the orbit
      i++;
    }
    std::cerr << "Before until changed=" << changed << "\n";
    if (changed >= 0) {
      break;
    }
  }
  std::cerr << "|blocks|=" << blocks.size() << "\n";
  // return the block system
  return blocks;
}


template<typename Telt>
std::pair<std::vector<std::vector<typename Telt::Tidx>>,std::vector<Face>> Blocks_Kernel(std::vector<Telt> const& ListGen, std::vector<typename Telt::Tidx> const& Omega, std::vector<typename Telt::Tidx> const& seed)
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
    ListGenRed.emplace_back(std::move(eGenRed));
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
    std::vector<Tidx> eBlock;
    size_t blkSiz=eFaceComb.count();
    boost::dynamic_bitset<>::size_type ePt=eFaceComb.find_first();
    for (size_t u=0; u<blkSiz; u++) {
      eBlock.push_back(Omega[ePt]);
      ePt = eFaceComb.find_next(ePt);
    }
    ListBlocks.emplace_back(std::move(eBlock));
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
  auto GetEdgeFromPosition=[&](size_t const& pos) -> std::pair<size_t,size_t> {
    size_t posTot = pos;
    size_t firstElement=0;
    while(true) {
      size_t allowedSiz = n - 1 - firstElement;
      if (posTot < allowedSiz) {
	size_t secondElement=firstElement + 1 + posTot;
	return {firstElement,secondElement};
      }
      posTot -= allowedSiz;
      firstElement++;
    }
  };
  auto GetPosFromEdge=[&](std::pair<size_t,size_t> const& eEdge) -> size_t {
    size_t eFirst=eEdge.first;
    size_t eSecond=eEdge.second;
    size_t pos=0;
    for (size_t u=0; u<eFirst; u++) {
      size_t allowedSiz = n - 1 - eFirst;
      pos += allowedSiz;
    }
    return pos + eSecond - 1 - eFirst;
  };
  while(true) {
    std::vector<size_t> eEdge=GetEdgeFromPosition(ThePos);
    std::pair<std::vector<std::vector<int>>,std::vector<Face>> ePair = Blocks_Kernel(ListGen, Omega, eEdge);
    if (ePair.first.size() > 1) {
      return ePair.first;
    }
    for (auto & eEdgeFace : ePair.second) {
      std::vector<size_t> eEdge = FaceToVector<size_t>(eEdgeFace);
      size_t pos=GetPosFromEdge(eEdge);
      UsedEdge[pos]=1;
    }
    bool test = UpdatePosition();
    if (!test)
      break;
  }
  return {Omega};
}

template<typename Telt>
std::function<Telt(Telt const&)> MapElementToSetPlusBlocks(std::vector<std::vector<typename Telt::Tidx>> const& blks, typename Telt::Tidx const& n)
{
  using Tidx=typename Telt::Tidx;
  Tidx TheMax=0;
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
      Tidx ePt=blks[iBlk][0];
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
    if (PositionVect_ui<int,size_t>(eBlock, ePt) != std::numeric_limits<size_t>::max())
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

