#ifndef PERMUTALIB_INCLUDE_NORMAL_STRUCTURE_H
#define PERMUTALIB_INCLUDE_NORMAL_STRUCTURE_H

#include "StabChainMain.h"
#include "BlockSystem.h"


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



template<typename Telt>
bool IsPrimitive_Subset(const std::vector<Telt>& LGen, const std::vector<typename Telt::Tidx>& subset, const typename Telt::Tidx& n) {
  using Tidx=typename Telt::Tidx;
  std::vector<Tidx> subset_rev(n);
  Tidx len = Tidx(subset.size());
  for (Tidx i=0; i<len; i++)
    subset_rev[subset[i]] = i;
  std::vector<Telt> gensB;
  for (auto & eGen : LGen) {
    std::vector<Tidx> eList(len);
    for (Tidx i=0; i<len; i++) {
      Tidx val1 = subset[i];
      Tidx val2 = PowAct(val1, eGen);
      Tidx val3 = subset_rev[val2];
      eList[i] = val3;
    }
    Telt eElt(std::move(eList));
    gensB.emplace_back(eElt);
  }
  std::vector<std::vector<Tidx>> blocks = Blocks(gensB, len);
  return blocks.size() == 1;
}


template<typename Telt, typename Tidx_label, typename Tint>
std::vector<Telt> Kernel_SmallGeneratingSet(const StabChain<Telt,Tidx_label>& G)
{
  using Tidx = typename Telt::Tidx;
  Tidx n = G->comm->n;
  std::unordered_set<Telt> gens_set;
  for (auto & eGen : Kernel_GeneratorsOfGroup(G)) {
    if (!eGen.isIdentity())
      gens_set.insert(eGen);
  }
  std::vector<Telt> gens;
  for (auto & eGen : gens_set)
    gens.push_back(eGen);

  size_t len=gens.size();
  Face status_remove(len);
  for (size_t i=0; i<len; i++) {
    if (status_remove[i] == 0) {
      for (size_t j=0; j<len; j++) {
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

  std::vector<Tidx> LMoved = MovedPoints(gens2, n);
  std::vector<std::vector<Tidx>> orb = OrbitsPerms(gens2, n, LMoved);
  size_t n_orb = orb.size();
  std::vector<size_t> orp;
  for (size_t i_orb=0; i_orb<n_orb; i_orb++)
    if (IsPrimitive_Subset(gens2, orb[i_orb], n))
      orp.push_back(i_orb);

  Tint order_G = Order<Telt,Tidx_label,Tint>(G);

  std::map<Tidx,int> LFact = FactorsSizeStabChain(G);
  auto check_correctness_gens=[&](const std::vector<Telt>& LGen) -> bool {
    if (LMoved.size() != MovedPoints(LGen, n).size())
      return false;
    if (orb.size() != OrbitsPerms(LGen, n, LMoved).size())
      return false;
    for (auto & i_orb : orp)
      if (!IsPrimitive_Subset(LGen, orb[i_orb], n))
        return false;
    StabChainOptions<Tint,Tidx> options = GetStandardOptions<Tint,Tidx>(n);
    StabChain<Telt,Tidx_label> U = StabChainOp_listgen<Telt,Tidx_label,Tint>(LGen, options);
    Tint order_U = Order<Telt,Tidx_label,Tint>(U);
    return order_G == order_U;
  };

  // Computing the lower bound on the number of generators
  size_t min = 1;
  if (!Kernel_IsCommutativeGenerators(gens2))
    min = 2;

  // Generating elements at
  auto get_and_test=[&](const size_t& i) -> bool {
    std::vector<Telt> gensB;
    for (size_t u=0; u<i; u++) {
      Telt g = RandomElement(gens2, n);
      gensB.push_back(std::move(g));
    }
    if (check_correctness_gens(gensB)) {
      gens2 = gensB;
      return true;
    }
    return false;
  };
  auto update_iife=[&]() -> void {
    size_t len = gens.size();
    for (size_t i=min; i<len; i++) {
      for (size_t j=0; j<5; j++) {
        bool test = get_and_test(i);
        if (test)
          return;
      }
    }
  };
  update_iife();

  size_t i = 1;

  while (i <= gens2.size() && gens2.size() > min) {
    // random did not improve much, try subsets

    std::vector<Telt> gensB;
    for (size_t i_orb=0; i_orb<gens2.size(); i_orb++) {
      if (i_orb != i-1)
        gensB.push_back(gens2[i_orb]);
    }
    if (check_correctness_gens(gensB)) {
      gens2 = gensB;
    } else {
      i++;
    }
  }
  return gens2;
}









template<typename Telt, typename Tidx_label, typename Tint>
StabChain<Telt,Tidx_label> SubsetStabChain(const StabChain<Telt,Tidx_label>& S, const std::vector<typename Telt::Tidx>& subset)
{
  using Tidx=typename Telt::Tidx;
  Tidx n = S->comm->n;
  Tidx miss_val = std::numeric_limits<Tidx>::max();
  std::vector<Tidx> subset_rev(n, miss_val);
  Tidx len = Tidx(subset.size());
  for (Tidx i=0; i<len; i++)
    subset_rev[subset[i]] = i;
  auto map=[&](const Telt& eGen) -> Telt {
    std::vector<Tidx> eList(len);
    for (Tidx i=0; i<len; i++) {
      Tidx i_big = subset[i];
      Tidx j_big = PowAct(i_big, eGen);
      Tidx j = subset_rev[j_big];
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
      if (j == miss_val) {
        std::cerr << "The subset is not stabilized. Clear bug\n";
        throw TerminalException{1};
      }
#endif
      eList[i] = j;
    }
    return Telt(std::move(eList));
  };

  // We cannot do an operation of subsetting directly on the stab chain
  // because the base might be outside of the subset
  std::vector<Telt> LGen = Kernel_GeneratorsOfGroup(S);
  size_t n_gen = LGen.size();
  std::vector<Telt> LGenRed(n_gen);
  for (size_t i=0; i<n_gen; i++)
    LGenRed[i] = map(LGen[i]);
  StabChainOptions<Tint,Tidx> options = GetStandardOptions<Tint,Tidx>(len);
  return StabChainOp_listgen<Telt,Tidx_label,Tint>(LGenRed, options);
}



template<typename Telt, typename Tidx_label>
StabChain<Telt,Tidx_label> Stabilizer_OnTuples_CorrectStabChain(const StabChain<Telt,Tidx_label>& S, const std::vector<typename Telt::Tidx>& subset)
{
  using Tidx=typename Telt::Tidx;
  Tidx n = S->comm->n;
  Tidx miss_val = std::numeric_limits<Tidx>::max();
  std::vector<Tidx> subset_rev(n, miss_val);
  Tidx len = Tidx(subset.size());
  for (Tidx i=0; i<len; i++)
    subset_rev[subset[i]] = i;
  StabChain<Telt,Tidx_label> Sptr = S;
  while (true) {
    if (Sptr == nullptr) {
      break;
    }
    if (Sptr->orbit.size() == 0) {
      break;
    }
    if (subset_rev[Sptr->orbit[0]] == miss_val) {
      break;
    }
    Sptr = Sptr->stabilizer;
  }
  return Sptr;
}





  /*
#############################################################################
##
#F  IntersectionNormalClosurePermGroup(<G>,<H>[,order]) . . . intersection of
#F                                   normal closure of <H> under <G> with <G>
##
##  computes $H^G \cap G$ as subgroup of Parent(G)
##
  */
template<typename Telt, typename Tidx_label, typename Tint>
StabChain<Telt,Tidx_label> Kernel_IntersectionNormalClosurePermGroup(const StabChain<Telt,Tidx_label>& G, const StabChain<Telt,Tidx_label>& H, const Tint& size)
{
  using Tidx = typename Telt::Tidx;
  Tidx n = G->comm->n;
  if (Order<Telt,Tidx_label,Tint>(G) == 1 || Order<Telt,Tidx_label,Tint>(H) == 1) {
    StabChainOptions<Tint,Tidx> options1 = GetStandardOptions<Tint,Tidx>(n);
    return StabChainOp_listgen<Telt,Tidx_label,Tint>({}, options1);
  }
  std::vector<Telt> newgens;
  for (auto & eGen : Kernel_GeneratorsOfGroup(G)) {
    std::vector<Tidx> eList(2*n);
    for (Tidx i=0; i<n; i++) {
      Tidx img = PowAct(i, eGen);
      eList[i] = img;
      eList[i + n] = img + n;
    }
    newgens.emplace_back(std::move(Telt(std::move(eList))));
  }
  for (auto & eGen : Kernel_GeneratorsOfGroup(H)) {
    std::vector<Tidx> eList(2*n);
    for (Tidx i=0; i<n; i++) {
      Tidx img = PowAct(i, eGen);
      eList[i] = i;
      eList[i + n] = img + n;
    }
    newgens.emplace_back(std::move(Telt(std::move(eList))));
  }
  StabChainOptions<Tint,Tidx> options2 = GetStandardOptions<Tint,Tidx>(2 * n);
  options2.size = size;
  options2.base.reserve(n);
  for (Tidx i=0; i<n; i++)
    options2.base.push_back(i + n);
  const std::vector<Telt> & tuple = options2.base;
  StabChain<Telt,Tidx_label> S = StabChainOp_listgen<Telt,Tidx_label,Tint>(newgens, options2);
  StabChain<Telt,Tidx_label> S_red = Stabilizer_OnTuples_CorrectStabChain(S, tuple);
  return SubsedtStabChain(S_red, tuple);
}










  /*

If g commutes with h then we have
g ( h(x) ) = h ( g(x) )




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
