#ifndef PERMUTALIB_INCLUDE_NORMAL_STRUCTURE_H
#define PERMUTALIB_INCLUDE_NORMAL_STRUCTURE_H

#include "StabChain.h"


namespace permutalib {


template<typename Telt, typename Tidx_label>
bool Kernel_IsNormalSubgroup(const StabChain<Telt,Tidx_label>& G, const StabChain<Telt,Tidx_label>& U)
{
  std::vector<Telt> gens_G = Kernel_GeneratorsOfGroup(G);
  std::vector<Telt> gens_U = Kernel_GeneratorsOfGroup(U);

  for (auto & eGen_G : gens_G) {
    for (auto & eGen_U : gens_U) {
      Telt eConj = Conjugation(eGen_U, eGen_G);
      if (!IsElementInStabChain(U, eConj))
        return false;
    }
  }
  return true;
}



template<typename Telt, typename Tidx_label, typename Tint>
StabChain<Telt,Tidx_label> Kernel_NormalClosure(const StabChain<Telt,Tidx_label>& G, const StabChain<Telt,Tidx_label>& H)
{
  std::vector<Telt> gens_G = Kernel_GeneratorsOfGroup(G);
  StabChain<Telt,Tidx_label> Hret = H;
  auto test=[&]() -> std::optional<Telt> {
    std::vector<Telt> gens_Hret = Kernel_GeneratorsOfGroup(Hret);
    for (auto & eGen_G : gens_G) {
      for (auto & eGen_Hret : gens_Hret) {
        Telt eConj = Conjugation(eGen_Hret, eGen_G);
        if (!IsElementInStabChain(Hret, eConj))
          return eConj;
      }
    }
    return {};
  };
  while(true) {
    std::optional<Telt> ret = test();
    if (!ret)
      break;
    ClosureGroup<Telt,Tidx_label,Tint>(Hret, *ret);
  }
  return Hret;
}


template<typename Telt, typename Tidx_label, typename Tint>
StabChain<Telt,Tidx_label> Kernel_DerivedSubgroup(const StabChain<Telt,Tidx_label>& G)
{
  using Tidx = typename Telt::Tidx;
  Tidx n = G->comm->n;
  StabChain<Telt,Tidx_label> S = EmptyStabChain<Telt,Tidx_label>(n);
  std::vector<Telt> LGen = Kernel_GeneratorsOfGroup(S);
  size_t n_gen = LGen.size();
  for (size_t i=0; i<n_gen; i++) {
    const Telt& g1 = LGen[i];
    for (size_t j=i+1; j<n_gen; j++) {
      const Telt& g2 = LGen[j];
      Telt comm = Conjugation(g1, g2) * Inverse(g1);
      ClosureGroup<Telt,Tidx_label,Tint>(S, comm);
    }
  }
  return Kernel_NormalClosure<Telt,Tidx_label,Tint>(G, S);
}

  /*

template<typename Telt, typename Tidx_label, typename Tint>
std::vector<Telt> Kernel_SmallGeneratingSet(const StabChain<Telt,Tidx_label>& G)
{
  using Tidx = typename Telt::Tidx;
  std::unordered_set<Telt> gens_set;
  for (auto & eGen : Kernel_GeneratorsOfGroup(G)) {
    if (!eGen.isIdentity())
      gens.insert(eGen);
  }
  std::vector<Telt> gens;
  for (auto & eGen : gens_set)
    gens.push_back(eGen);

  size_t len=gens.size();
  Face status_remove(len);
  for (size_t i=0; i<len; i++) {
    if (status_remove[i] == 0) {
      if (size_t j=0; j<len; j++) {
        if (i != j && status_remove[j] == 0) {
          Tidx val = LogPerm(gens[i], gens[j]); // test if gens[i]^e = gens[j]
          if (val != std::numeric_limits<Tidx>::max()) {
            status_remove[j]=1;
          }
        }
      }
    }
  }
  std::vector<Telt> gens2;
  for (size_t i=0; i<len; i++)
    if (status_remove[i] == 0)
      gens2.push_back(gens[i]);

  orb:=Set(List(Orbits(G,MovedPoints(G)),Set));
  orp:=Filtered([1..Length(orb)],x->IsPrimitive(Action(G,orb[x])));

  min:=2;
  if Length(gens)>2 then
  # minimal: AbelianInvariants
    min:=Maximum(List(Collected(Factors(Size(G)/Size(DerivedSubgroup(G)))),x->x[2]));
    min:=Maximum(min,2);
    if min=Length(GeneratorsOfGroup(G)) then return GeneratorsOfGroup(G);fi;
    i:=Maximum(2,min);
    while i<=min+1 and i<Length(gens) do
      # try to find a small generating system by random search
      j:=1;
      while j<=5 and i<Length(gens) do
        U:=Subgroup(G,List([1..i],j->Random(G)));
        ok:=true;
        # first test orbits
        if ok then
          ok:=Length(orb)=Length(Orbits(U,MovedPoints(U))) and
              ForAll(orp,x->IsPrimitive(U,orb[x]));
        fi;
	StabChainOptions(U).random:=100; # randomized size
        if ok and Size(U)=Size(G) then
          gens:=Set(GeneratorsOfGroup(U));
        fi;
        j:=j+1;
      od;
      i:=i+1;
    od;
  fi;

  i := 1;
  if not IsAbelian(G) then
    i:=i+1;
  fi;
  while i <= Length(gens) and Length(gens)>min do
    # random did not improve much, try subsets
    U:=Subgroup(G,gens{Difference([1..Length(gens)],[i])});

    ok:=true;
    # first test orbits
    if ok then
      ok:=Length(orb)=Length(Orbits(U,MovedPoints(U))) and
          ForAll(orp,x->IsPrimitive(U,orb[x]));
    fi;

    StabChainOptions(U).random:=100; # randomized size
    if Size(U)<Size(G) then
      i:=i+1;
    else
      gens:=Set(GeneratorsOfGroup(U));
    fi;
  od;
  return gens;
end);




#############################################################################
##
#F  CentralizerTransSymmCSPG()  . . . . .  centralizer of transitive G in S_n
##
##  computes the centralizer of a transitive group G in S_n
##
InstallGlobalFunction( CentralizerTransSymmCSPG, function( G, chainG )
    local   n,          # the degree of G
            x,          # the first base point
            L,          # the set of fixed points of stabgroup
            orbitx,     # the orbit of x in the centralizer;
                        # eventually, orbitx=L
            y,          # a point in L
            z,          # loop variable running through permutation domain
            h,          # a coset representative of G, written as word in the
                        # generators
            gens,       # list of generators for the centralizer
            gen,        # an element of gens
            Ggens,      # generators of G
            Ginverses,  # list of inverses for the generators of G
            H;          # output group
    if IsTrivial(G)  then
        return TrivialSubgroup( Parent(G) );
    fi;

    if IsBound( chainG.stabFxdPnts ) then
       x := chainG.stabFxdPnts[1];
       L := chainG.stabFxdPnts[2];
       n := LargestMovedPoint(G);
       if not IsBound( chainG.orbit ) or chainG.orbit[1] <> x then
          chainG := EmptyStabChain( [  ], (), x );
          AddGeneratorsExtendSchreierTree( chainG, GeneratorsOfGroup(G) );
       fi;
    else
       n := LargestMovedPoint(G);
       x := chainG.orbit[1];
       L := Difference( [ 1 .. n ],
                        MovedPoints( chainG.stabilizer.generators ) );
    fi;
    Ginverses := GInverses( chainG );
    Ggens := chainG.generators;

    # the centralizer of G is semiregular, acting transitively on L
    orbitx := [x];
    gens := [];
    while Length(orbitx) < Length(L) do

        # construct element of centralizer which carries x to new point in L
        gen := [];
        y := Difference(L,orbitx)[1];
        for z in [1..n] do
            h := CosetRepAsWord( x, z, chainG.transversal );
            h := InverseAsWord(h,Ggens,Ginverses);
            gen[z] := ImageInWord(y,h);
        od;
        Add(gens,PermList(gen));
        orbitx := OrbitPerms(gens,x);
    od;

    H := SubgroupNC( G, gens );
    SetSize( H, Length( L ) );
    return H;
end );

#############################################################################
##
#F  IntersectionNormalClosurePermGroup(<G>,<H>[,order]) . . . intersection of
#F                                   normal closure of <H> under <G> with <G>
##
##  computes $H^G \cap G$ as subgroup of Parent(G)
##
InstallGlobalFunction( IntersectionNormalClosurePermGroup,
    function(arg)
    local   G,H,        # the groups to be handled
            n,          # maximum of degrees of G,H
            i,j,        # loop variables
            conperm,    # perm exchanging first and second n points
            newgens,    # set of extended generators
            options,    # options record for stabilizer computation
            group;      # the group generated by newgens
                        # stabilizing the second n points, we get H^G \cap G

    G := arg[1];
    H := arg[2];

    if IsTrivial(G) or IsTrivial(H)  then
        return TrivialSubgroup( Parent(G) );
    fi;

    n := Maximum(LargestMovedPoint(G),
                 LargestMovedPoint(H));
    conperm := PermList( Concatenation( [n+1 .. 2*n] , [1 .. n] ) );
    # extend the generators of G acting on [n+1..2n] exactly as on [1..n]
    newgens := List( StabChainMutable( G ).generators,
                     g -> g * ( g^conperm ) );

    # from the generators of H, create permutations which act on [n+1..2n]
    # as the original generator on [1..n] and which act trivially on [1..n]
    for i in StabChainMutable( H ).generators do
      Add( newgens, i^conperm );
    od;

    group := GroupByGenerators(newgens,());


    # create options record for stabilizer chain computation
    options := rec( base := [n+1..2*n] );
    #if size of group is part of input, use it
    if Length(arg) = 3 then
       options.size := arg[3];
       # if H is normalized by G and G,H already have stabilizer chains
       # then compute base for group
       #if ( IsBound(G.size) or IsBound(G.stabChain) ) and
       #   ( IsBound(H.size) or IsBound(H.stabChain) )  then
       #   if Size(G) * Size(H) = arg[3] then
       #      options.knownBase :=
       #      Concatenation( List( Base(H), x -> n + x ), Base(G) ) ;
       #   fi;
       #fi;
    fi;
    StabChain(group,options);
#T is this meaningful ??
    group := Stabilizer(group,[n+1 .. 2*n],OnTuples);
    return AsSubgroup( Parent(G),group);
end );




#############################################################################
##
#M  Centre( <G> ) . . . . . . . . . . . . . . . center of a permutation group
##
##  constructs the center of G.
##  Reference: Beals-Seress, 24th Symp. on Theory of Computing 1992, sect. 9
##
InstallMethod( Centre,
    "for a permutation group",
    [ IsPermGroup ],
    function(G)
    local   n,          # degree of G
            orbits,     # list of orbits of G
            base,       # lexicographically smallest (in list) base of G
            i,j,        # loop variables
            reps,       # array recording which orbit of G the points in
                        # perm. domain belong to
            domain,     # union of G orbits which contain base points
            significant,# indices of orbits in "orbits" that belong to domain
            max,        # loop variable, used at definition of significant
            len,        # length of domain
            tchom,      # trans. const. homom, restricting G to domain
            GG,         # the image of tchom
            chainGG,    # stabilizer chain of `GG'
            chainGGG,   # stabilizer chain of `GGG'
            orbit,      # an orbit of GG
            tchom2,     # trans. const. homom, restricting GG to orbit
            GGG,        # the image of GG at tchom2
            hgens,      # list of generators for the direct product of
                        # centralizers of GG in Sym(orbit), for orbits of GG
            order,      # order of `GroupByGenerators( hgens, () )'
            centr,      # the centralizer of GG in Sym(orbit)
            inverse2,   # inverse of the conjugating permutation of tchom2
            g,          # generator of centr
            cent;       # center of GG

    if IsTrivial(G)  then
       return TrivialSubgroup(G);
    fi;

    base := BaseStabChain(StabChainMutable(G));
    n := Maximum( Maximum( base ), LargestMovedPoint(G) );
    orbits := OrbitsDomain(G,[1..n]);
    # orbits := List( orbits, Set );

    # handle case of transitive G directly
    if Length(orbits) = 1  then
        centr := CentralizerTransSymmCSPG( G, StabChainMutable( G ) );
        if IsEmpty( GeneratorsOfGroup( centr ) ) then
           return TrivialSubgroup( G );
        else
           order := Size(centr);
           cent := IntersectionNormalClosurePermGroup(G,centr,order*Size(G));
           Assert( 1, IsAbelian( cent ) );
           SetIsAbelian( cent, true );
           return cent;
        fi;
    fi;
    # for intransitive G, find which orbit contains which
    # points of permutation domain
    reps := [];
    for i in [1..Length(orbits)] do
        for j in [1..Length(orbits[i])] do
            reps[orbits[i][j]] := i;
        od;
    od;

    # take union of significant orbits which contain base points
    max := reps[base[1]];
    significant := [max];
    domain := ShallowCopy(orbits[max]);
    for i in [2..Length(base)] do
        if not (reps[base[i]] in significant)  then
            max := reps[base[i]];
            Append(domain,orbits[max]);
            Add(significant,max);
        fi;
    od;
    len := Length(domain);

    # restrict G to significant orbits
    if n = len then
       GG := G;
    else
       tchom := ActionHomomorphism(G,domain,"surjective");
       GG := Image(tchom,G);
    fi;

    # handle case of transitive GG directly
    if Length(significant) = 1  then
        centr := CentralizerTransSymmCSPG( GG, StabChainMutable( GG ) );
        if IsEmpty( GeneratorsOfGroup( centr ) ) then
           return TrivialSubgroup( G );
        else
           order := Size( centr );
           cent := IntersectionNormalClosurePermGroup(GG,centr,order*Size(GG));
           cent:= PreImages(tchom,cent);
           Assert( 1, IsAbelian( cent ) );
           SetIsAbelian( cent, true );
           return cent;
        fi;
    fi;

    # case of intransitive GG
    # for each orbit of GG, construct generators of centralizer of GG in
    # Sym(orbit).  hgens is a list of generators for the direct product of
    # these centralizers.
    # the group generated by hgens contains the center of GG
    hgens := [];
    order := 1;
    for i in significant do
        if n = len then
           orbit := orbits[i];
        else
           orbit := OnTuples(orbits[i],tchom!.conperm);
        fi;
        tchom2 := ActionHomomorphism(GG,orbit,"surjective");
        GGG := Image(tchom2,GG);
        chainGG:= StabChainOp( GG, [ orbit[1] ] );
        chainGGG:= StabChainMutable( GGG );
        chainGGG.stabFxdPnts:=[ orbit[1]^tchom2!.conperm,
            OnTuples( Difference(orbit,
                      MovedPoints( chainGG.stabilizer.generators ) ),
                      tchom2!.conperm ) ];
        centr := CentralizerTransSymmCSPG( GGG, chainGGG );
        if not IsEmpty( GeneratorsOfGroup( centr ) ) then
           order := order * Size( centr );
           inverse2 := tchom2!.conperm^(-1);
           for g in StabChainMutable( centr ).generators do
               Add(hgens,g^inverse2);
           od;
        fi;
    od;
    if order = 1 then
        return TrivialSubgroup( G );
    else
        cent := IntersectionNormalClosurePermGroup
                 ( GG, GroupByGenerators(hgens,()), order*Size(GG) );
        if n <> len then
          cent:= PreImages( tchom, cent );
        fi;
        Assert( 1, IsAbelian( cent ) );
        SetIsAbelian( cent, true );
        return cent;
    fi;
end );




  */


}

#endif
