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
    std::vector<Tidx> V(size_t(n));;
    for (Tidx i=0; i<n; i++)
      V[i] = i;
    return {V};
  }
  //
  // Now no easy case. Need to do a complex computation
  //

    one:= One( G );
    orbit := [ D[1] ];
    trans := [];
    trans[ D[1] ] := one;
    for pnt  in orbit  do
        for gen  in acts  do
            if not IsBound( trans[ pnt / gen ] )  then
                Add( orbit, pnt / gen );
                trans[ pnt / gen ] := gen;
            fi;
        od;
    od;

    # check that the group is transitive
    if Length( orbit ) <> Length( D )  then
        Error("<G> must operate transitively on <D>");
    fi;
    nrorbs := Length( orbit );

    # since $i \in k^{G_1}$ implies $\beta(i)=\beta(k)$,  we initialize <eql>
    # so that the connected components are orbits of some subgroup  $H < G_1$
    eql := [];
    leq := [];
    next := [];
    last := [];
    for pnt  in orbit  do
        eql[pnt]  := pnt;
        leq[pnt]  := pnt;
        next[pnt] := 0;
        last[pnt] := pnt;
    od;

    # repeat until we have a block system
    changed := 0;
    cur := orbit[2];
    rnd := one;
    repeat

        # compute such an $H$ by taking random  Schreier generators  of $G_1$
        # and stop if 2 successive generators dont change the orbits any more
	while changed < 2  do

            # compute a random Schreier generator of $G_1$
            i := Length( orbit );
            while 1 <= i  do
                rnd := rnd * Random( acts );
                i   := QuoInt( i, 2 );
            od;
            gen := rnd;
            d1g:=D[1]^gen;
            while d1g <> D[1]  do
                tr:=trans[ d1g ];
		gen := gen * tr;
                d1g:=d1g^tr;
            od;
            changed := changed + 1;

            # compute the image of every point under <gen>
            for pnt  in orbit  do
                img := pnt ^ gen;

                # find the representative of the orbit of <pnt>
                while eql[pnt] <> pnt  do
                    pnt := eql[pnt];
		od;

                # find the representative of the orbit of <img>
		while eql[img] <> img  do
                    img := eql[img];
                od;

                # if the don't agree merge their orbits
                if   pnt < img  then
                    eql[img] := pnt;
                    next[ last[pnt] ] := img;
                    last[pnt] := last[img];
                    nrorbs := nrorbs - 1;
                    changed := 0;
                elif img < pnt  then
                    eql[pnt] := img;
                    next[ last[img] ] := pnt;
                    last[img] := last[pnt];
                    nrorbs := nrorbs - 1;
                    changed := 0;
                fi;
            od;

        od;

        # take arbitrary point <cur>,  and an element <gen> taking 1 to <cur>
        while eql[cur] <> cur  do
            cur := eql[cur];
        od;
        gen := [];
        img := cur;
        while img <> D[1]  do
            Add( gen, trans[img] );
            img := img ^ trans[img];
        od;
        gen := Reversed( gen );

        # compute an alleged block as orbit of 1 under $< H, gen >$
        pnt := cur;
        while pnt <> 0  do

            # compute the representative of the block containing the image
            img := pnt;
            for i  in gen  do
                img := img / i;
            od;
            while eql[img] <> img  do
                img := eql[img];
            od;

            # if it is not our current block but a minimal block
            if   img <> D[1]  and img <> cur  and leq[img] = img  then

                # then try <img> as a new start
                leq[cur] := img;
                cur := img;
                gen := [];
                img := cur;
                while img <> D[1]  do
                    Add( gen, trans[img] );
                    img := img ^ trans[img];
                od;
                gen := Reversed( gen );
                pnt := cur;


            # otherwise if it is not our current block but contains it
            # by construction a nonminimal block contains the current block
            elif img <> D[1]  and img <> cur  and leq[img] <> img  then

                # then merge all blocks it contains with <cur>
                while img <> cur  do
                    eql[img] := cur;
                    next[ last[cur] ] := img;
                    last[ cur ] := last[ img ];
                    img := leq[img];
                    while img <> eql[img]  do
                        img := eql[img];
                    od;
                od;
                pnt := next[pnt];

            # go on to the next point in the orbit
            else
                pnt := next[pnt];
            fi;
        od;

        # make the alleged block
        block := [ D[1] ];
        pnt := cur;
        while pnt <> 0  do
            Add( block, pnt );
            pnt := next[pnt];
        od;
        block := Set( block );
        blocks := [ block ];

        # quick test to see if the group is primitive
        if Length( block ) = Length( orbit )  then
            return Immutable( [ D ] );
        fi;

        # quick test to see if the orbit can be a block
        if Length( orbit ) mod Length( block ) <> 0  then
            changed := -1000;
        fi;

        # '<rep>[<i>]' is the representative of the block containing <i>
        rep := [];
        for pnt  in orbit  do
            rep[pnt] := 0;
        od;
        for pnt  in block  do
            rep[pnt] := 1;
        od;
        
        # compute the block system with an orbit algorithm
        i := 1;
        while 0 <= changed  and i <= Length( blocks )  do

            # loop over the generators
            for gen  in acts  do

                # compute the image of the block under the generator
                img := OnSets( blocks[i], gen );

                # if this block is new
                if rep[ img[1] ] = 0  then

                    # add the new block to the list of blocks
                    Add( blocks, img );

                    # check that all points in the image are new
                    for pnt  in img  do
                        if rep[pnt] <> 0  then
                            changed := -1000;
                        fi;
                        rep[pnt] := img[1];
                    od;

                # if this block is old
                else

                    # check that all points in the image lie in the block
                    for pnt  in img  do
                        if rep[pnt] <> rep[img[1]]  then
                            changed := -1000;
                        fi;
                    od;

                fi;

            od;

            # on to the next block in the orbit
            i := i + 1;
        od;

    until 0 <= changed;
    # force sortedness
    if Length(blocks[1])>0 and CanEasilySortElements(blocks[1][1]) then
      blocks:=AsSSortedList(List(blocks,i->Immutable(Set(i))));
      IsSSortedList(blocks);
    fi;
    # return the block system
    return Immutable( blocks );
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

