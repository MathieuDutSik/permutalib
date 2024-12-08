// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GAP_ASCENDINGCHAINS_AND_COSETS_H_
#define SRC_GAP_ASCENDINGCHAINS_AND_COSETS_H_

#include <limits>
#include <list>
#include <unordered_set>
#include <utility>
#include <vector>

/*
  Computation of ascending chain:
  ---A lot of tricks are applied in the generic GAP code like Normalizer, centralizer, etc.
  ---But this may due to the fact that the groups are general ones. For permutation
     groups, we can look for different stuff.
  ---By orbit structure, we look at the structure of the orbits of the group on the point set.
  ---If the orbit structure of G and H are the same, then we can take a stabilizer of one
     point (belonging to the smallest non-trivial orbit) and compute the stabilizer in G and H.
     ---Denote the stabilizers G_s, H_s and compute ascending chain for them by recursion.
     ---Find the transversal u_1, ..., u_k such that H = <H_s, u_1, ..., u_k>
     ---Then the the ascending chain for <H_s, G_s> can be lifted to an ascending chain
     by adding those elements.
  ---If the orbit structures are not the same, then their restriction on G leaves an interesting
     subgroup of H < K < G.
  ---And so we can go from that. By combining those tricks, we can get ascending chains
  from H to G.

  Uses of ascending chains for computing right cosets:
  ---We want to compute the right cosets of H in G.
  ---Suppose we have a chain
  H = G0 < G1 < G2 < .... < Gm = G
  ---The right cosets of G_i in G_{i+1} can be computed.
  This allows to build an iterator for the right cosets from H to G.

  Uses of ascending chains for computing with double cosets:
  ---We have two subgroups H and K of G.
  ---We want to decompose G as
  H g_0 K \cup H g_1 K \cup ... \cup H g_m K
  ---Those are conjugacy invariants. If we replace H and/or K by a conjugate
  then the decomposition can be computed.
  ---We can select a direction for the decomposition:
  H = G0 < G1 < G2 < ... < Gm = G
  ---Potentially, we could have decompositions for (Gi, K) and then decompose
  them into decomposition for (G_{i-1}, K).
  ---That should allow to go from a few double cosets to more double cosets.
  ---The first step of the algorithm is also clear:
  Write Gi = G_{i-1} u0 \cup ... \cup G_{i-1} uK
  Then we have for Gi h K the equality
  G_{i-1} u0 h K \cup G_{i-1} u1 K \cup ... \cup G_{i-1} uK h K
  ---Then the question is to identify when two cosets H u K and H v K are
  equal.
  ---OBJECTION: This is all very clever, but it collapses for the subgroup M24
  of S24:
     ---M24 is a maximal subgroup of S24 (besides A24). Therefore no ascending
     chain to speak of.
     ---The index is 2534272925184000 therefore single coset enumeration is not
     feasible.
     ---But actually GAP fails for the triple (G, U, V) = (S24, M24, M24) OOM
     ---But actually GAP fails for the triple (G, U, V) = (S22, M22, M22) Too long
     ---For (S12, M12, M12) the number of double cosets is 8. So, maybe some
     algorithm using cosets is used.
  ---If H and K generates a strict subgroup of G, can we have some room for
  improvements? NO: The above examples kill it: a common example is when H = K
  and those examples are clearly very large.
  ---We can count on the intersection working? Maybe? We have the functionality
  but we need to be sure that it works correctly.
  ---If so, that give us a useful tool: The size of a double coset
  H u K = u H^u K with H^u = u^{-1} H u being the conjugate subgroup.
  The size of the double coset should thus be |H| x |K| / |K \cap H^u|
  ---If we had the cosets, then we could have another enumeration algorithm.
  ---S: If H is a normal subgroup then HK is a subgroup and it becomes a coset
  algorithm.
  ---Q: If we compute the normalizer of H, how can that help? Not clear.
  ---S: The homomorphism method that we use for the orbit splitting
  in the polyhedral code is a workable method that is actually used in the
  "CalcDoubleCosets" function. So, it makes sense to look for such actions.
  ---S: The TryMaximalSubgroupClassReps cannot be realistically used for us.
  ---S: If H1 is a normal subgroup in H2 then we have with the cosets
  H2 = H1 g1 \cup ... \cup H1 gK.
  So, if we have H2 g K = H1 g1 g K \cup ... \cup H1 gK g K.
  Since the group is normal, we have H1 g1 g K = g1 H1 g K.
  Therefore all the double cosets H1 gI g K have the same size.
  However, we can have two phenomenons (separate or at once):
     ---A single coset H1 g1 g K could be equal to a H2 g K. Just that the
     intersection decreases in size.
     ---The coset H1 g1 g K is smaller than H2 g K.
     Both situation are visible with a commutative group G.
  ---Q: What is the algorithm used by the GAP in case of no homomorphism?
  Not very clear.


  Obtaining homomorphism:
  ---If H is a subgroup of G, we want to find an action of G on a set X
  such that H is the stabilizer of a point.
  ---One way to do that is to use the cosets, that is quite workable.
  The threshold in the GAP code is 1000000. So fairly high.
  ---The groups we work with are permutation groups acting on a space X.
  ---We can take decompose X under the action of H:
     ---If there are several orbits, we can compute the stabilizers of
     each of them.
     ---If any one of the stabilizrs is intermediate between H and G.
     then we should have detected that at the construction of ascending
     chains.
     ---So, if we have an orbit not completely stabilized, we can use
     that for having the homomorphism. Things will proceed just as in the
     polyhedral code.
     ---Therefore the failing scenario is when the stabilizer of the
     orbits under H are also orbits under G.
  ---So, together, that gets us two reasonable strategies to get such
  group actions on sets X.

 */



namespace permutalib {

/*
  From the stabilizer chain, we get an ascending chain almost automatically
  which is already very good.
  ---
  So, we face the situation of a group G acting on an orbit O.
  The stabilizer of O[0] is the subgroup H.
  We want to find a pyramid of block decompositions so as not just to test
  primitivity but get a sequence of groups.

 */
template <typename Telt, typename Tidx_label, typename Tint>
std::vector<StabChain<Telt, Tidx_label>>
Kernel_AscendingChain(StabChain<Telt, Tidx_label> const &G) {
  using Tidx = typename Telt::Tidx;
  using Tstab = StabChain<Telt, Tidx_label>;
  Tidx miss_val = std::numeric_limits<Tidx>::max();
  Tidx n_vert = G->comm->identity.size();
  std::list<Tstab> ListStab = StdListStabChain(G);
  ListStab.reverse();
  /*
  size_t pos = 0;
  for (auto & eStab : ListStab) {
    Tint ord = Order<Telt,Tidx_label,Tint>(eStab);
    std::cerr << "pos=" << pos << " ord=" << ord << "\n";
    pos++;
  }
  */
  auto iter = ListStab.begin();
  size_t len = ListStab.size();
  std::vector<Tstab> ListGroup;
  ListGroup.push_back(*iter);
  for (size_t i = 1; i < len; i++) {
    iter++;
    Tstab const &TheStab = *iter;
    std::vector<Telt> LGenBig = Kernel_GeneratorsOfGroup(TheStab);
    std::vector<Tidx> const &orbit = TheStab->orbit;
    Tidx pt_stab = orbit[0];
    Tidx len = Tidx(orbit.size());
    std::vector<Tidx> orbit_rev(n_vert, miss_val);
    for (Tidx i = 0; i < len; i++) {
      orbit_rev[orbit[i]] = i;
    }
    std::vector<Telt> LGenSma;
    for (auto &eGen : LGenBig) {
      std::vector<Tidx> eList(len);
      for (Tidx i = 0; i < len; i++) {
        Tidx pt1 = orbit[i];
        Tidx pt2 = OnPoints(pt1, eGen);
        Tidx pt3 = orbit_rev[pt2];
        eList[i] = pt3;
      }
      Telt eGenSma(std::move(eList));
      LGenSma.emplace_back(std::move(eGenSma));
    }
    std::vector<BlockDecomposition<Tidx>> l_blkdec =
        ComputeSequenceBlockDecomposition(LGenSma, len);
    for (size_t j = 1; j < l_blkdec.size() - 1; j++) {
      BlockDecomposition<Tidx> const &BlkDec = l_blkdec[j];
      Tidx pt_stab2 = orbit_rev[pt_stab];
      Tidx iBlock = BlkDec.map_vert_block[pt_stab2];
      Face Phi(n_vert);
      for (auto &ePt : BlkDec.ListBlocks[iBlock]) {
        Phi[orbit[ePt]] = 1;
      }
      Tstab eStab =
          Kernel_Stabilizer_OnSets<Telt, Tidx_label, Tint>(TheStab, Phi);
      ListGroup.push_back(eStab);
    }
    ListGroup.push_back(TheStab);
  }
  return ListGroup;
}


/*
  U is a subgroup of G.
  We compute the right transversals H g
  If the function f is triggered, then the enumeration ends prematurely.
*/
template <typename Telt, typename Tidx_label, typename Tint, typename Fterminate>
std::optional<std::vector<Telt>>
Kernel_RightTransversal_Direct_f(StabChain<Telt, Tidx_label> const &G,
                                 StabChain<Telt, Tidx_label> const &H,
                                 Fterminate f_terminate) {
  std::vector<Telt> ListTransversal;
  std::unordered_map<Telt, uint8_t> map;
  auto fInsert = [&](Telt const &x) -> void {
    Telt x_can = MinimalElementCosetStabChain(H, x);
    uint8_t &pos = map[x_can];
    if (pos == 0) {
      pos = 1;
      ListTransversal.push_back(x_can);
    }
  };
  Telt id = G->comm->identity;
  std::vector<Telt> LGen = Kernel_SmallGeneratingSet<Telt,Tidx_label,Tint>(G);
  Tint size_G = Order<Telt,Tidx_label,Tint>(G);
  Tint size_H = Order<Telt,Tidx_label,Tint>(H);
  Tint index = size_G / size_H;
  if (f_terminate(id)) {
    return {};
  }
  fInsert(id);
  size_t pos = 0;
  while (true) {
    size_t len = ListTransversal.size();
#ifndef CHECK_COSET_ENUMERATION
    Tint n_coset = len;
    if (n_coset == index) {
      // Early termination is possible
      break;
    }
#endif
    if (pos == len)
      break;
    for (size_t idx = pos; idx < len; idx++) {
      // Copy is needed because the insertion into ListTransversal
      // invalidates the reference.
      Telt eTrans = ListTransversal[idx];
      for (auto &eGen : LGen) {
        Telt eProd = eTrans * eGen;
        if (f_terminate(eProd)) {
          return {};
        }
        fInsert(eProd);
      }
    }
    pos = len;
  }
#ifdef CHECK_COSET_ENUMERATION
  Tint n_coset = ListTransversal.size();
  if (n_coset != index) {
    std::cerr << "n_coset=" << n_coset << " index=" << index << "\n";
    std::cerr << "G=[";
    for (auto & eElt : get_all_elements(G))
      std::cerr << " " << eElt;
    std::cerr << " ]\n";
    std::cerr << "H=[";
    for (auto & eElt : get_all_elements(H))
      std::cerr << " " << eElt;
    std::cerr << " ]\n";
    std::cerr << "ListTransversal=[";
    for (auto & eElt : ListTransversal)
      std::cerr << " " << eElt;
    std::cerr << " ]\n";
    std::cerr << "The enumeration found a wrong number of cosets\n";
    throw PermutalibException{1};
  }
#endif
  return ListTransversal;
}


/*
  U is a subgroup of G.
  We compute the right transversals H g
*/
template <typename Telt, typename Tidx_label, typename Tint>
std::vector<Telt>
Kernel_RightTransversal_Direct(StabChain<Telt, Tidx_label> const &G,
                               StabChain<Telt, Tidx_label> const &H) {
  auto f_terminate=[]([[maybe_unused]] Telt const& x) -> bool {
    return false;
  };
  std::optional<std::vector<Telt>> opt =
    Kernel_RightTransversal_Direct_f<Telt,Tidx_label,Tint,decltype(f_terminate)>(G, H, f_terminate);
  return *opt;
}

/*
  U is a subgroup of G.
  We compute the left transversals g H
*/
template <typename Telt, typename Tidx_label, typename Tint, typename Fterminate>
void Kernel_LeftTransversal_Direct_f(StabChain<Telt, Tidx_label> const &G,
                                     StabChain<Telt, Tidx_label> const &H,
                                     Fterminate f_terminate) {
  auto f_terminate_right=[&](Telt const& x) -> bool {
    Telt x_inv = Inverse(x);
    return f_terminate(x_inv);
  };
  (void)Kernel_RightTransversal_Direct_f<Telt,Tidx_label,Tint,decltype(f_terminate_right)>(G, H, f_terminate_right);
}


/*
  U is a subgroup of G.
  We compute the left transversals g H
*/
template <typename Telt, typename Tidx_label, typename Tint>
std::vector<Telt>
Kernel_LeftTransversal_Direct(StabChain<Telt, Tidx_label> const &G,
                              StabChain<Telt, Tidx_label> const &H) {
  std::vector<Telt> ListTransversal = Kernel_RightTransversal_Direct<Telt,Tidx_label,Tint>(G, H);
  size_t len = ListTransversal.size();
  std::vector<Telt> ListRet(len);
  for (size_t i = 0; i < len; i++)
    ListRet[i] = Inverse(ListTransversal[i]);
  return ListRet;
}

template <typename Telt, typename Tidx_label, typename Tint>
void CheckRightCosets(StabChain<Telt, Tidx_label> const &G,
                     StabChain<Telt, Tidx_label> const &H,
                     std::vector<Telt> const& ListRightTransversal) {
  std::vector<Telt> l_elt_g = get_all_elements(G);
  std::vector<Telt> l_elt_h = get_all_elements(H);
  std::unordered_set<Telt> set_elt;
  size_t ProdSize = l_elt_h.size() * ListRightTransversal.size();
  if (ProdSize != l_elt_g.size()) {
    std::cerr << "|l_elt_g|=" << l_elt_g.size() << "\n";
    std::cerr << "|l_elt_h|=" << l_elt_h.size() << "\n";
    std::cerr << "|ListRightTransversal|=" << ListRightTransversal.size() << "\n";
    std::cerr << "ProdSize=" << ProdSize << "\n";
    std::cerr << "Discrepancy at the order level\n";
    throw PermutalibException{1};
  }
  for (auto & eElt : ListRightTransversal) {
    for (auto & e_h : l_elt_h) {
      Telt eProd = e_h * eElt;
      if (set_elt.count(eProd) == 1) {
        std::cerr << "The element eProd is already present\n";
        throw PermutalibException{1};
      }
      set_elt.insert(eProd);
    }
  }
}

/*
  Functions depending on special data types.

  Those are
  --- lenflock / flock / lenblock / index / respStab
  --- group / subgroup / stabChainGroup / stabChainSubgroup

 */

/*
#############################################################################
##
#F  MinimizeExplicitTransversal( <U>, <maxmoved> )  . . . . . . . . . . local
##
*/
/*
InstallGlobalFunction( MinimizeExplicitTransversal, function( U, maxmoved )
  local   explicit,  lenflock,  flock,  lenblock,  index,  s;

  if     IsBound( U.explicit )
     and IsBound( U.stabilizer )  then
      explicit := U.explicit;
      lenflock := U.stabilizer.index * U.lenblock / Length( U.orbit );
      flock    := U.flock;
      lenblock := U.lenblock;
      index    := U.index;
      ChangeStabChain( U, [ 1 .. maxmoved ] );
      for s  in [ 1 .. Length( explicit ) ]  do
          explicit[ s ] := MinimalElementCosetStabChain( U, explicit[ s ] );
      od;
      Sort( explicit );
      U.explicit := explicit;
      U.lenflock := lenflock;
      U.flock    := flock;
      U.lenblock := lenblock;
      U.index    := index;
  fi;
end );
*/

/*
InstallGlobalFunction( AddCosetInfoStabChain, function( G, U, maxmoved )
  local   orb,  pimg,  img,  vert,  s,  t,  index,
          block,  B,  blist,  pos,  sliced,  lenflock,  i,  j,
          ss,  tt,t1,t1lim;

  if IsEmpty( G.genlabels )  then
      U.index    := 1;
      U.explicit := [ U.identity ];
      U.lenflock := 1;
      U.flock    := U.explicit;
  else
      AddCosetInfoStabChain( G.stabilizer, U.stabilizer, maxmoved );

      // U.index := [G_1:U_1];
      U.index := U.stabilizer.index * Length( G.orbit ) / Length( U.orbit );

      // block := 1 ^ <U,G_1>; is a block for G.
      block := OrbitPerms( Concatenation( U.generators, G.stabilizer.generators
), G.orbit[ 1 ] ); U.lenblock := Length( block ); lenflock := Length( G.orbit )
/ U.lenblock;

      # For small indices,  permutations   are multiplied,  so  we  need  a
      # multiplied transversal.
      if     IsBound( U.stabilizer.explicit )
         and U.lenblock * maxmoved <= MAX_SIZE_TRANSVERSAL
         and U.index    * maxmoved <= MAX_SIZE_TRANSVERSAL * lenflock  then
          U.explicit := [  ];
          U.flock    := [ G.identity ];
          tt := [  ];  tt[ G.orbit[ 1 ] ] := G.identity;
          for t  in G.orbit  do
              tt[ t ] := tt[ t ^ G.transversal[ t ] ] /
                         G.transversal[ t ];
          od;
      fi;

      // flock := { G.transversal[ B[1] ] | B in block system };
      blist := BlistList( G.orbit, block );
      pos := Position( blist, false );
      while pos <> fail  do
          img := G.orbit[ pos ];
          B := block{ [ 1 .. U.lenblock ] };
          sliced := [  ];
          while img <> G.orbit[ 1 ]  do
              Add( sliced, G.transversal[ img ] );
              img := img ^ G.transversal[ img ];
          od;
          for i  in Reversed( [ 1 .. Length( sliced ) ] )  do
              for j  in [ 1 .. Length( B ) ]  do
                  B[ j ] := B[ j ] / sliced[ i ];
              od;
          od;
          Append( block, B );
          if IsBound( U.explicit )  then
              Add( U.flock, tt[ B[ 1 ] ] );
          fi;
          UniteBlistList(G.orbit, blist, B );
          pos := Position( blist, false, pos );
      od;
      G.orbit := block;

      # Let <s> loop over the transversal elements in the stabilizer.
      U.repsStab := List( [ 1 .. U.lenblock ], x ->
                         BlistList( [ 1 .. U.stabilizer.index ], [  ] ) );
      U.repsStab[ 1 ] := BlistList( [ 1 .. U.stabilizer.index ],
                                    [ 1 .. U.stabilizer.index ] );
      index := U.stabilizer.index * lenflock;
      s := 1;

      # For  large  indices, store only   the  numbers of  the  transversal
      # elements needed.
      if not IsBound( U.explicit )  then

          # If  the   stabilizer   is the   topmost  level   with  explicit
          # transversal, this must contain minimal coset representatives.
          MinimizeExplicitTransversal( U.stabilizer, maxmoved );

          # if there are over 200 points, do a cheap test first.
          t1lim:=Length(G.orbit);
          if t1lim>200 then
            t1lim:=50;
          fi;

          orb := G.orbit{ [ 1 .. U.lenblock ] };
          pimg := [  ];
          while index < U.index  do
              pimg{ orb } := CosetNumber( G.stabilizer, U.stabilizer, s,
                                     orb );
              t := 2;
              while t <= U.lenblock  and  index < U.index  do

                  # do not test all points first if not necessary
                  # (test only at most t1lim points, if the test succeeds,
                  # test the rest)
                  # this gives a major speedup.
                  t1:=Minimum(t-1,t1lim);
                  # For this point  in the  block,  find the images  of the
                  # earlier points under the representative.
                  vert := G.orbit{ [ 1 .. t1 ] };
                  img := G.orbit[ t ];
                  while img <> G.orbit[ 1 ]  do
                      vert := OnTuples( vert, G.transversal[ img ] );
                      img  := img           ^ G.transversal[ img ];
                  od;

                  # If $Ust = Us't'$ then $1t'/t/s in 1U$. Also if $1t'/t/s
                  # in 1U$ then $st/t' =  u.g_1$ with $u  in U, g_1 in G_1$
                  # and $g_1  =  u_1.s'$ with $u_1  in U_1,  s' in S_1$, so
                  # $Ust = Us't'$.
                  if ForAll( [ 1 .. t1 ], i -> not IsBound
                     ( U.translabels[ pimg[ vert[ i ] ] ] ) )  then

                    # do all points
                    if t1<t-1 then
                      vert := G.orbit{ [ 1 .. t - 1 ] };
                      img := G.orbit[ t ];
                      while img <> G.orbit[ 1 ]  do
                          vert := OnTuples( vert, G.transversal[ img ] );
                          img  := img           ^ G.transversal[ img ];
                      od;
                      if ForAll( [ t1+1 .. t - 1 ], i -> not IsBound
                        ( U.translabels[ pimg[ vert[ i ] ] ] ) )  then
                          U.repsStab[ t ][ s ] := true;
                          index := index + lenflock;
                      fi;
                    else
                      U.repsStab[ t ][ s ] := true;
                      index := index + lenflock;
                    fi;
                  fi;

                  t := t + 1;
              od;
              s := s + 1;
          od;

      // For small indices, store a transversal explicitly.
      else
          for ss  in U.stabilizer.flock  do
              Append( U.explicit, U.stabilizer.explicit * ss );
          od;
          while index < U.index  do
              t := 2;
              while t <= U.lenblock  and  index < U.index  do
                  ss := U.explicit[ s ] * tt[ G.orbit[ t ] ];
                  if ForAll( [ 1 .. t - 1 ], i -> not IsBound
                         ( U.translabels[ G.orbit[ i ] / ss ] ) )  then
                      U.repsStab[ t ][ s ] := true;
                      Add( U.explicit, ss );
                      index := index + lenflock;
                  fi;
                  t := t + 1;
              od;
              s := s + 1;
          od;
          Unbind( U.stabilizer.explicit );
          Unbind( U.stabilizer.flock    );
      fi;

  fi;
end );
*/

/*
#############################################################################
##
#F  RightTransversalPermGroupConstructor( <filter>, <G>, <U> )  . constructor
##
*/
// MAX_SIZE_TRANSVERSAL := 100000;

/*
// If the option "noascendingchain" is selected then
ValueOption("noascendingchain") = true
// and thus noyet = false and so most of the text below does not apply
BindGlobal( "RightTransversalPermGroupConstructor", function( filter, G, U )
local GC, UC, noyet, orbs, domain, GCC, UCC, ac, nc, bpt, enum, i;

  GC := CopyStabChain( StabChainImmutable( G ) );
  UC := CopyStabChain( StabChainImmutable( U ) );
  noyet:=ValueOption("noascendingchain")<>true;
  if not IsTrivial( G )  then
      orbs := ShallowCopy( OrbitsDomain( U, MovedPoints( G ) ) );
      Sort( orbs, function( o1, o2 )
          return Length( o1 ) < Length( o2 ); end );
      domain := Concatenation( orbs );
      GCC:=GC;
      UCC:=UC;
      while    Length( GCC.genlabels ) <> 0
            or Length( UCC.genlabels ) <> 0  do
        if noyet and (
        (SizeStabChain(GCC)/SizeStabChain(UCC)*10 > MAX_SIZE_TRANSVERSAL) ||
        (Length(UCC.genlabels)=0 && SizeStabChain(GCC) > MAX_SIZE_TRANSVERSAL)
          ) then
          // we potentially go through many steps, making it expensive
          ac:=AscendingChain(G,U:cheap);
          // go in biggish steps through the chain
          nc:=[ac[1]];
          for i in [3..Length(ac)] do
            if Size(ac[i])/Size(nc[Length(nc)])>MAX_SIZE_TRANSVERSAL then
              Add(nc,ac[i-1]);
            fi;
          od;
          Add(nc,ac[Length(ac)]);
          if Length(nc)>2 then
            ac:=[];
            for i in [Length(nc),Length(nc)-1..2] do
              // do not try to factor again
              Add(ac,RightTransversal(nc[i],nc[i-1]:noascendingchain));
            od;
            return FactoredTransversal(G,U,ac);
          fi;
          noyet:=false;

        fi;
        bpt := First( domain, p -> not IsFixedStabilizer( GCC, p ) );
        ChangeStabChain( GCC, [ bpt ], true  );  GCC := GCC.stabilizer;
        ChangeStabChain( UCC, [ bpt ], false );  UCC := UCC.stabilizer;
      od;
  fi;

  AddCosetInfoStabChain(GC,UC,LargestMovedPoint(G));
  MinimizeExplicitTransversal(UC,LargestMovedPoint(G));

  enum := Objectify( NewType( FamilyObj( G ),
                         filter and IsList and IsDuplicateFreeList
                         and IsAttributeStoringRep ),
        rec( group := G,
          subgroup := U,
    stabChainGroup := GC,
 stabChainSubgroup := UC ) );

  return enum;
end );
*/

/*
#############################################################################
##
##  IntermediateGroup(<G>,<U>)  . . . . . . . . . subgroup of G containing U
##
##  This routine tries to find a subgroup E of G, such that G>E>U. If U is
##  maximal, it returns fail. This is done by using the maximal subgroups
machinery or
##  finding minimal blocks for
##  the operation of G on the Right Cosets of U.
##
*/
/*
template<typename Telt, typename Tidx_label>
std::optional<StabChain<Telt,Tidx_label>>
IntermediateGroup(StabChain<Telt,Tidx_label> const& G_in,
StabChain<Telt,Tidx_label> const& U)
{
  if (U.size() == G.size())
    return {};

  intersize:=Size(G);

  // use maximals, use `Try` as we call with limiting options
  // We disable the code using the TryMaximalSubgroupClassReps
  // as it requires a lot of prerequisite.

  Tidx_label hardlimit=1000000;

  if (Index(G,U) > hardlimit) {
    return {};
  }

  if Length(GeneratorsOfGroup(G))>3 then
    G1:=Group(SmallGeneratingSet(G));
    if HasSize(G) then
      SetSize(G1,Size(G));
    fi;
    G:=G1;
  fi;
  o:=ActionHomomorphism(G,RightTransversal(G,U:noascendingchain),
OnRight,"surjective"); img:=Range(o); b:=Blocks(img,MovedPoints(img)); if
Length(b)=1 then return fail; else b:=StabilizerOfBlockNC(img,First(b,i->1 in
i)); b:=PreImage(o,b); return b; fi;
}
*/

/*
#############################################################################
##
#F  RefinedChain(<G>,<c>) . . . . . . . . . . . . . . . .  refine chain links
##
*/
/*
template<typename Telt, typename Tidx_label>
std::vector<StabChain<Telt,Tidx_label>> RefinedChain(StabChain<Telt,Tidx_label>
const& G, std::vector<StabChain<Telt,Tidx_label>> const& cc)
{
  bound:=(10*LogInt(Size(G),10)+1)*Maximum(Factors(Size(G)));
  bound:=Minimum(bound,20000);
  cheap:=ValueOption("cheap")=true;
  c:=ValueOption("refineIndex");
  if IsInt(c) then
    bound:=c;
  fi;

  c:=[];
  for i in [2..Length(cc)] do
    Add(c,cc[i-1]);
    if Index(cc[i],cc[i-1]) > bound then
      a:=AsSubgroup(Parent(cc[i]),cc[i-1]);
      olda:=TrivialSubgroup(a);
      while Index(cc[i],a)>bound and Size(a)>Size(olda) do
        olda:=a;
        // try extension via normalizer
        b:=Normalizer(cc[i],a);
        if Size(b)>Size(a) then
           // extension by normalizer surely is a normal step
          normalStep:=true;
          bb:=b;
        else
          bb:=cc[i];
          normalStep:=false;
          b:=Centralizer(cc[i],Centre(a));
        fi;
        if Size(b)=Size(a) or Index(b,a)>bound then
          cnt:=8+2^(LogInt(Index(bb,a),9));
                      // if cheap then cnt:=Minimum(cnt,50);fi;
          cnt:=Minimum(cnt,40); # as we have better intermediate
          repeat
            if cnt<20 and not cheap then
                      // if random failed: do hard work
              b:=IntermediateGroup(bb,a);
              if b=fail then
                b:=bb;
              fi;
              cnt:=0;
            else
              // larger indices may take more tests...
              repeat
                r:=Random(bb);
              until not(r in a);
              if normalStep then
                # NC is safe
                b:=ClosureSubgroupNC(a,r);
              else
                // self normalizing subgroup: thus every element not in <a>
                // will surely map one generator out
                j:=0;
                gens:=GeneratorsOfGroup(a);
                repeat
                  j:=j+1;
                until not(gens[j]^r in a);
                r:=gens[j]^r;

                # NC is safe
                b:=ClosureSubgroupNC(a,r);
              fi;
              if Size(b)<Size(bb) then
                bb:=b;
              fi;
              cnt:=cnt-1;
            fi;
          until Index(bb,a)<=bound or cnt<1;
        fi;
        if Index(b,a)>bound and Length(c)>1 then
          bb:=IntermediateGroup(b,c[Length(c)-1]);
          if bb<>fail and Size(bb)>Size(c[Length(c)]) then
            c:=Concatenation(c{[1..Length(c)-1]},[bb],Filtered(cc,x->Size(x)>=Size(b)));
            return RefinedChain(G,c);
          fi;
        fi;

        a:=b;
        if a<>cc[i] then #not upper level
          Add(c,a);
        fi;

      od;
    fi;
  od;
  Add(c,cc[Length(cc)]);
  a:=c[Length(c)];
  for i in [Length(c)-1,Length(c)-2..1] do
          //  enforce parent relations
    if not HasParent(c[i]) then
      SetParent(c[i],a);
      a:=c[i];
    else
      a:=AsSubgroup(a,c[i]);
      c[i]:=a;
    fi;
  od;
  return c;
}
*/

/*
#############################################################################
##
#M  AscendingChainOp(<G>,<pnt>) . . . approximation of
##
*/
/*
template<typename Telt, typename Tidx_label>
std::vector<StabChain<Telt,typename Tidx_label>>
AscendingChain(StabChain<Telt,typename Tidx_label> const& G,
StabChain<Telt,typename Tidx_label> const& U)
{
  s:=G;
  c:=[G];
  repeat
    mp:=MovedPoints(s);
    o:=ShallowCopy(OrbitsDomain(s,mp));
    Sort(o,function(a,b) return Length(a)<Length(b);end);
    i:=1;
    step:=false;
    while i<=Length(o) and step=false do
      if not IsTransitive(U,o[i]) then
        o:=ShallowCopy(OrbitsDomain(U,o[i]));
        Sort(o,function(a,b) return Length(a)<Length(b);end);
        # union of same length -- smaller index
        a:=Union(Filtered(o,x->Length(x)=Length(o[1])));
        if Length(a)=Sum(o,Length) then
          a:=Set(o[1]);
        fi;
        s:=Stabilizer(s,a,OnSets);
        step:=true;
      elif Index(G,U)>NrMovedPoints(U)
          and IsPrimitive(s,o[i]) and not IsPrimitive(U,o[i]) then
        s:=Stabilizer(s,Set(List(MaximalBlocks(U,o[i]),Set)),
                      OnSetsDisjointSets);
        step:=true;
      else
        i:=i+1;
      fi;
    od;
    if step then
      Add(c,s);
    fi;
  until step=false or Index(s,U)=1; # we could not refine better
  if Index(s,U)>1 then
    Add(c,U);
  fi;
  return RefinedChain(G,Reversed(c));
}
*/

/*

#############################################################################
##
#F  CalcDoubleCosets( <G>, <A>, <B> ) . . . . . . . . .  double cosets: A\G/B
##
##  DoubleCosets routine using an
##  ascending chain of subgroups from A to G, using the fact, that a
##  double coset is an union of right cosets
##
BindGlobal("CalcDoubleCosets",function(G,a,b)
local c, flip, maxidx, refineChainActionLimit, cano, tryfct, p, r, t,
      stabs, dcs, homs, tra, a1, a2, indx, normal, hom, omi, omiz,c1,
      unten, compst, s, nr, nstab, lst, sifa, pinv, blist, bsz, cnt,
      ps, e, mop, mo, lstgens, lstgensop, rep, st, o, oi, i, img, ep,
      siz, rt, j, canrep,stab,step,nu,doneidx,orbcnt,posi,
      sizes,cluster,sel,lr,lstabs,ssizes,num,
      actlimit, uplimit, badlimit,avoidlimit;

  Print("Beginning of CalcDoubleCosets\n");
  actlimit:=300000; # maximal degree on which we try blocks
  uplimit:=10000; # maximal index for up step
  avoidlimit:=200000; # beyond this index we want to get smaller
  badlimit:=1000000; # beyond this index things might break down

  # if a is small and b large, compute cosets b\G/a and take inverses of the
  # representatives: Since we compute stabilizers in b and a chain down to
  # a, this is notably faster
  if ValueOption("noflip")<>true and 3*Size(a)<2*Size(b) then
    c:=b;
    b:=a;
    a:=c;
    flip:=true;
  else
    flip:=false;
  fi;

  if Index(G,a)=1 then
    return [[One(G),Size(G)]];
  fi;

  # maximal index of a series
  maxidx:=function(ser)
    return Maximum(List([1..Length(ser)-1],x->Size(ser[x+1])/Size(ser[x])));
  end;

  # compute ascending chain and refine if necessarily (we anyhow need action
  # on cosets).

  c:=AscendingChain(G,a:refineChainActionLimit:=actlimit,indoublecoset);

  # cano indicates whether there is a final up step (and thus we need to
  # form canonical representatives). ```Canonical'' means that on each
  # transversal level the orbit representative is chosen to be minimal (in
  # the transversal position).
  cano:=false;

  doneidx:=[]; # indices done already -- avoid duplicate
  if maxidx(c)>avoidlimit then
    # try to do better

    # what about flipping (back)?
    c1:=AscendingChain(G,b:refineChainActionLimit:=actlimit,indoublecoset);
    if maxidx(c1)<=avoidlimit then
      c:=b;
      b:=a;
      a:=c;
      flip:=not flip;
      c:=c1;

    elif IsPermGroup(G) then

      actlimit:=Maximum(actlimit,NrMovedPoints(G));
      avoidlimit:=Maximum(avoidlimit,NrMovedPoints(G));

      tryfct:=function(obj,act)
        local G1,a1,c1;
        if IsList(act) and Length(act)=2 then
          G1:=act[1];
          a1:=act[2];
        else
          G1:=Stabilizer(G,obj,act);
          if Index(G,G1)<maxidx(c) then
            a1:=Stabilizer(a,obj,act);
          else
            a1:=G;
          fi;
        fi;
        if Index(G,G1)<maxidx(c) and Index(a,a1)<=uplimit and (
          maxidx(c)>avoidlimit or Size(a1)>Size(c[1])) then
          c1:=AscendingChain(G1,a1:refineIndex:=avoidlimit,
                                   refineChainActionLimit:=actlimit,
                                   indoublecoset);
          if maxidx(c1)<maxidx(c) then
            c:=Concatenation(c1,[G]);
            cano:=true;
          fi;
        fi;
      end;

      r:=Filtered(TryMaximalSubgroupClassReps(G:cheap),
        x->Index(G,x)<=5*avoidlimit);
      SortBy(r,a->-Size(a));
      for i in r do
        if Index(G,i)<maxidx(c) then
          p:=Intersection(a,i);
          AddSet(doneidx,Index(a,p));
          if Index(a,p)<=uplimit then
            tryfct("max",[i,p]);
          fi;
        fi;
      od;

      p:=LargestMovedPoint(a);
      tryfct(p,OnPoints);

      for i in Orbits(Stabilizer(a,p),Difference(MovedPoints(a),[p])) do
        tryfct(Set([i[1],p]),OnSets);
      od;

    fi;

    if maxidx(c)>badlimit then

      r:=ShallowCopy(TryMaximalSubgroupClassReps(a:cheap));
      r:=Filtered(r,x->Index(a,x)<uplimit and not Index(a,x) in doneidx);

      SortBy(r,a->-Size(a));
      for j in r do
        t:=AscendingChain(G,j:refineIndex:=avoidlimit,
                              refineChainActionLimit:=actlimit,indoublecoset);
        if maxidx(t)<maxidx(c) and (maxidx(c)>badlimit or
          # only increase up-step if index gets better by extra index
          (maxidx(c)>maxidx(t)*Size(c[1])/Size(t[1])) ) then
          c:=t;
          cano:=true;
        fi;

      od;

    fi;

  elif ValueOption("sisyphus")=true then
    # purely to allow for tests of up-step mechanism in smaller examples.
    # This is creating unneccessary extra work and thus should never be used
    # in practice, but will force some code to be run through.
    c:=Concatenation([TrivialSubgroup(G)],c);
    cano:=true;
  fi;

  r:=[One(G)];
  stabs:=[b];
  dcs:=[];

  # Do we want to keep result for a smaller group (as cheaper fuse is possible
  # outside function at a later stage)?
  if ValueOption("noupfuse")=true then cano:=false; fi;

  # calculate setup for once
  homs:=[];
  tra:=[];
  for step in [1..Length(c)-1] do
    a1:=c[Length(c)-step+1];
    a2:=c[Length(c)-step];
    indx:=Index(a1,a2);
    normal:=IsNormal(a1,a2);
    # don't try to refine again for transversal, we've done so already.
    t:=RightTransversal(a1,a2:noascendingchain);
    tra[step]:=t;

    # is it worth using a permutation representation?
    if step>1 and Length(t)<badlimit and IsPermGroup(G) and
      not normal then
      # in this case, we can beneficially compute the action once and then use
      # homomorphism methods to obtain the permutation image
      hom:=Subgroup(G,SmallGeneratingSet(a1));
      hom:=ActionHomomorphism(hom,t,OnRight,"surjective");
    else
      hom:=fail;
    fi;
    homs[step]:=hom;
  od;

  for step in [1..Length(c)-1] do
    a1:=c[Length(c)-step+1];
    a2:=c[Length(c)-step];
    normal:=IsNormal(a1,a2);
    indx:=Index(a1,a2);

    # is this the last step?
    unten:=step=Length(c)-1 and cano=false;

    # shall we compute stabilizers?
    compst:=(not unten) or normal;

    t:=tra[step];
    hom:=homs[step];

    s:=[];
    nr:=[];
    nstab:=[];
    for nu in [1..Length(r)] do
      lst:=stabs[nu];
      sifa:=Size(a2)*Size(b)/Size(lst);
      p:=r[nu];
      pinv:=p^-1;
      blist:=BlistList([1..indx],[]);
      bsz:=indx;
      orbcnt:=0;

      # if a2 is normal in a1, the stabilizer is the same for all Orbits of
      # right cosets. Thus we need to compute only one, and will receive all
      # others by simple calculations afterwards

      if normal then
        cnt:=1;
      else
        cnt:=indx;
      fi;

      lstgens:=GeneratorsOfGroup(lst);
      if Length(lstgens)>2 then
        lstgens:=SmallGeneratingSet(lst);
      fi;
      lstgensop:=List(lstgens,i->i^pinv); # conjugate generators: operation
      # is on cosets a.p; we keep original cosets: Ua.p.g/p, this
      # corresponds to conjugate operation

      if hom<>fail then
        lstgensop:=List(lstgensop,i->Image(hom,i));
      fi;

      posi:=0;
      while bsz>0 and cnt>0 do
        cnt:=cnt - 1;

        # compute orbit and stabilizers for the next step
        # own Orbitalgorithm and stabilizer computation

        posi:=Position(blist,false,posi);
        ps:=posi;
        blist[ps]:=true;
        bsz:=bsz - 1;
        e:=t[ps];
        mop:=1;
        mo:=ps;

        rep := [ One(b) ];
        st := TrivialSubgroup(lst);

        o:=[ps];
        if compst then
          oi:=[];
          oi[ps]:=1; # reverse index
        fi;
        orbcnt:=orbcnt + 1;

        i:=1;
        while i<=Length(o)
          # will not grab if nonreg,. orbit and stabilizer not computed,
          # but comparatively low cost and huge help if hom=fail
          and Size(st)*Length(o)<Size(lst) do
          for j in [1..Length(lstgens)] do
            if hom=fail then
              img:=t[o[i]]*lstgensop[j];
              ps:=PositionCanonical(t,img);
            else
              ps:=o[i]^lstgensop[j];
            fi;
            if blist[ps] then
              if compst then
                # known image
                #NC is safe (initializing as TrivialSubgroup(G)
                st := ClosureSubgroupNC(st,rep[i]*lstgens[j]/rep[oi[ps]]);
              fi;
            else
              # new image
              blist[ps]:=true;
              bsz:=bsz - 1;
              Add(o,ps);
              if compst then
                Add(rep, rep[i] * lstgens[j]);
                oi[ps]:=Length(o);
              fi;
            fi;
          od;
          i:=i+1;
        od;

        ep:=e * rep[mop] * p;
        Add(nr, ep);

        if compst then
          st:=st^rep[mop];
          Add(nstab, st);
        fi;

        siz:=sifa*Length(o);

        if unten then
          if flip then
            Add(dcs,[ep^(-1),siz]);
          else
            Add(dcs,[ep,siz]);
          fi;
        fi;

      od;

      if normal then
        # in the normal case, we can obtain the other orbits easily via
        # the orbit theorem (same stabilizer)
        rt:=RightTransversal(lst,st);
        Assert(1, Length(rt)=Length(o));

        while bsz > 0 do
          ps:=Position(blist, false);
          e:=t[ps];
          blist[ps]:=true;

          ep:=e*p;
          mo:=ep;
          mop:=ps;
          # tick off the orbit
          for i in rt do
            j:=ep*i/p;
            ps:=PositionCanonical(t,ep*i/p);
            blist[ps]:=true;
          od;
          bsz:=bsz-Length(rt);

          Add(nr,mo);
          Add(nstab,st);

          if unten then
            if flip then
              Add(dcs,[ep^(-1),siz]);
            else
              Add(dcs,[ep,siz]);
            fi;
          fi;

        od;

      fi;

    od;
    stabs:=nstab;
    r:=nr;
  od;

  # NOTE FROM MDS: The cano code is deleted as a priori not needed for us

  return dcs;
end);

*/

// clang-format off
}  // namespace permutalib
#endif  // SRC_GAP_ASCENDINGCHAINS_AND_COSETS_H_
// clang-format on
