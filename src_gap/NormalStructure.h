#ifndef SRC_GAP_NORMALSTRUCTURE_H_
#define SRC_GAP_NORMALSTRUCTURE_H_

#include "BlockSystem.h"
#include "StabChainMain.h"
#include <limits>
#include <map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace permutalib {

template <typename Telt, typename Tidx_label>
bool Kernel_IsNormalSubgroup(const StabChain<Telt, Tidx_label> &G,
                             const StabChain<Telt, Tidx_label> &U) {
  std::vector<Telt> gens_G = Kernel_GeneratorsOfGroup(G);
  std::vector<Telt> gens_U = Kernel_GeneratorsOfGroup(U);

  for (auto &eGen_G : gens_G) {
    for (auto &eGen_U : gens_U) {
      Telt eConj = Conjugation(eGen_U, eGen_G);
      if (!IsElementInStabChain(U, eConj))
        return false;
    }
  }
  return true;
}

template <typename Telt, typename Tidx_label, typename Tint>
StabChain<Telt, Tidx_label>
Kernel_NormalClosure(const StabChain<Telt, Tidx_label> &G,
                     const StabChain<Telt, Tidx_label> &H) {
  std::vector<Telt> gens_G = Kernel_GeneratorsOfGroup(G);
  StabChain<Telt, Tidx_label> Hret = H;
  auto test = [&]() -> std::optional<Telt> {
    std::vector<Telt> gens_Hret = Kernel_GeneratorsOfGroup(Hret);
    for (auto &eGen_G : gens_G) {
      for (auto &eGen_Hret : gens_Hret) {
        Telt eConj = Conjugation(eGen_Hret, eGen_G);
        if (!IsElementInStabChain(Hret, eConj))
          return eConj;
      }
    }
    return {};
  };
  while (true) {
    std::optional<Telt> ret = test();
    if (!ret)
      break;
    ClosureGroup<Telt, Tidx_label, Tint>(Hret, *ret);
  }
  return Hret;
}

template <typename Telt, typename Tidx_label, typename Tint>
StabChain<Telt, Tidx_label>
Kernel_DerivedSubgroup(const StabChain<Telt, Tidx_label> &G) {
  using Tidx = typename Telt::Tidx;
  Tidx n = G->comm->n;
  StabChain<Telt, Tidx_label> S = EmptyStabChain<Telt, Tidx_label>(n);
  std::vector<Telt> LGen = Kernel_GeneratorsOfGroup(G);
  size_t n_gen = LGen.size();
  for (size_t i = 0; i < n_gen; i++) {
    const Telt &g1 = LGen[i];
    for (size_t j = i + 1; j < n_gen; j++) {
      const Telt &g2 = LGen[j];
      Telt comm = Conjugation(g1, g2) * Inverse(g1);
      ClosureGroup<Telt, Tidx_label, Tint>(S, comm);
    }
  }
  return Kernel_NormalClosure<Telt, Tidx_label, Tint>(G, S);
}

template <typename Telt>
bool IsPrimitive_Subset(const std::vector<Telt> &LGen,
                        const std::vector<typename Telt::Tidx> &subset,
                        const typename Telt::Tidx &n) {
  using Tidx = typename Telt::Tidx;
  std::vector<Tidx> subset_rev(n);
  Tidx len = Tidx(subset.size());
  for (Tidx i = 0; i < len; i++)
    subset_rev[subset[i]] = i;
  std::vector<Telt> gensB;
  for (auto &eGen : LGen) {
    std::vector<Tidx> eList(len);
    for (Tidx i = 0; i < len; i++) {
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

template <typename Telt, typename Tidx_label, typename Tint>
std::vector<Telt>
Kernel_SmallGeneratingSet(const StabChain<Telt, Tidx_label> &G) {
  using Tidx = typename Telt::Tidx;
  Telt id = G->comm->identity;
  Tidx n = id.size();
  std::unordered_set<Telt> gens_set;
  for (auto &eGen : Kernel_GeneratorsOfGroup(G))
    if (!eGen.isIdentity())
      gens_set.insert(eGen);
  std::vector<Tidx> bas = BaseStabChain(G);
  std::vector<Telt> gens;
  for (auto &eGen : gens_set)
    gens.push_back(eGen);

  size_t len = gens.size();
  Face status_remove(len);
  for (size_t i = 0; i < len; i++) {
    if (status_remove[i] == 0) {
      for (size_t j = 0; j < len; j++) {
        if (i != j && status_remove[j] == 0) {
          Tidx val = LogPerm(gens[i], gens[j]); // test if gens[i]^e = gens[j]
          if (val != std::numeric_limits<Tidx>::max()) {
            status_remove[j] = 1;
          }
        }
      }
    }
  }
  std::vector<Telt> gens2;
  for (size_t i = 0; i < len; i++)
    if (status_remove[i] == 0)
      gens2.push_back(gens[i]);

  std::vector<Tidx> LMoved = MovedPoints(gens2, n);
  std::vector<std::vector<Tidx>> orb = OrbitsPerms(gens2, n, LMoved);
  size_t n_orb = orb.size();
  std::vector<size_t> orp;
  for (size_t i_orb = 0; i_orb < n_orb; i_orb++)
    if (IsPrimitive_Subset(gens2, orb[i_orb], n))
      orp.push_back(i_orb);

  Tint order_G = Order<Telt, Tidx_label, Tint>(G);

  std::map<Tidx, int> LFact = FactorsSizeStabChain(G);
  auto check_correctness_gens = [&](const std::vector<Telt> &LGen) -> bool {
    if (LMoved.size() != MovedPoints(LGen, n).size())
      return false;
    if (orb.size() != OrbitsPerms(LGen, n, LMoved).size())
      return false;
    for (auto &i_orb : orp)
      if (!IsPrimitive_Subset(LGen, orb[i_orb], n))
        return false;
    StabChainOptions<Tint, Telt> options = GetStandardOptions<Tint, Telt>(id);
    StabChain<Telt, Tidx_label> U =
        StabChainOp_listgen<Telt, Tidx_label, Tint>(LGen, options);
    Tint order_U = Order<Telt, Tidx_label, Tint>(U);
    return order_G == order_U;
  };

  // Computing the lower bound on the number of generators
  size_t min = 1;
  if (!Kernel_IsCommutativeGenerators(gens2, bas))
    min = 2;

  // Generating elements at
  auto get_and_test = [&](const size_t &i) -> bool {
    std::vector<Telt> gensB;
    for (size_t u = 0; u < i; u++) {
      Telt g = RandomElement(gens2, id);
      gensB.push_back(std::move(g));
    }
    if (check_correctness_gens(gensB)) {
      gens2 = gensB;
      return true;
    }
    return false;
  };
  auto update_iife = [&]() -> void {
    size_t len = gens.size();
    for (size_t i = min; i < len; i++) {
      for (size_t j = 0; j < 5; j++) {
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
    for (size_t i_orb = 0; i_orb < gens2.size(); i_orb++) {
      if (i_orb != i - 1)
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

template <typename Telt, typename Tidx_label, typename Tint>
StabChain<Telt, Tidx_label>
SubsetStabChain(const StabChain<Telt, Tidx_label> &S,
                const std::vector<typename Telt::Tidx> &subset) {
  using Tidx = typename Telt::Tidx;
  Tidx n = S->comm->n;
  Tidx miss_val = std::numeric_limits<Tidx>::max();
  std::vector<Tidx> subset_rev(n, miss_val);
  Tidx len = Tidx(subset.size());
  for (Tidx i = 0; i < len; i++)
    subset_rev[subset[i]] = i;
  auto map = [&](const Telt &eGen) -> Telt {
    std::vector<Tidx> eList(len);
    for (Tidx i = 0; i < len; i++) {
      Tidx i_big = subset[i];
      Tidx j_big = PowAct(i_big, eGen);
      Tidx j = subset_rev[j_big];
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
      if (j == miss_val) {
        std::cerr << "The subset is not stabilized. Clear bug\n";
        throw PermutalibException{1};
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
  for (size_t i = 0; i < n_gen; i++)
    LGenRed[i] = map(LGen[i]);
  StabChainOptions<Tint, Telt> options =
      GetStandardOptions<Tint, Telt>(Telt(len));
  return StabChainOp_listgen<Telt, Tidx_label, Tint>(LGenRed, options);
}

template <typename Telt, typename Tidx_label>
StabChain<Telt, Tidx_label> Stabilizer_OnTuples_CorrectStabChain(
    const StabChain<Telt, Tidx_label> &S,
    const std::vector<typename Telt::Tidx> &subset) {
  using Tidx = typename Telt::Tidx;
  Tidx n = S->comm->n;
  Tidx miss_val = std::numeric_limits<Tidx>::max();
  std::vector<Tidx> subset_rev(n, miss_val);
  Tidx len = Tidx(subset.size());
  for (Tidx i = 0; i < len; i++)
    subset_rev[subset[i]] = i;
  StabChain<Telt, Tidx_label> Sptr = S;
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

template <typename Telt> bool IsTrivialListGen(const std::vector<Telt> &LGen) {
  for (auto &eGen : LGen)
    if (!eGen.isIdentity())
      return false;
  return true;
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
template <typename Telt, typename Tidx_label, typename Tint>
StabChain<Telt, Tidx_label> Kernel_IntersectionNormalClosurePermGroup_LGen(
    const typename Telt::Tidx &n, const std::vector<Telt> &LGen_G,
    const std::vector<Telt> &LGen_H, const Tint &size) {
  using Tidx = typename Telt::Tidx;
  if (IsTrivialListGen(LGen_G) || IsTrivialListGen(LGen_H)) {
    StabChainOptions<Tint, Telt> options1 =
        GetStandardOptions<Tint, Telt>(Telt(n));
    return StabChainOp_listgen<Telt, Tidx_label, Tint>({}, options1);
  }
  std::vector<Telt> newgens;
  for (auto &eGen : LGen_G) {
    std::vector<Tidx> eList(2 * n);
    for (Tidx i = 0; i < n; i++) {
      Tidx img = PowAct(i, eGen);
      eList[i] = img;
      eList[i + n] = img + n;
    }
    newgens.emplace_back(std::move(Telt(std::move(eList))));
  }
  for (auto &eGen : LGen_H) {
    std::vector<Tidx> eList(2 * n);
    for (Tidx i = 0; i < n; i++) {
      Tidx img = PowAct(i, eGen);
      eList[i] = i;
      eList[i + n] = img + n;
    }
    newgens.emplace_back(std::move(Telt(std::move(eList))));
  }
  StabChainOptions<Tint, Telt> options2 =
      GetStandardOptions<Tint, Telt>(Telt(2 * n));
  options2.size = size;
  options2.base.reserve(n);
  for (Tidx i = 0; i < n; i++)
    options2.base.push_back(i + n);
  const std::vector<Tidx> &tuple = options2.base;
  StabChain<Telt, Tidx_label> S =
      StabChainOp_listgen<Telt, Tidx_label, Tint>(newgens, options2);
  StabChain<Telt, Tidx_label> S_red =
      Stabilizer_OnTuples_CorrectStabChain(S, tuple);
  std::vector<Tidx> subset(n);
  for (Tidx i = 0; i < n; i++)
    subset[i] = i;
  return SubsetStabChain<Telt, Tidx_label, Tint>(S_red, subset);
}

template <typename Telt, typename Tidx_label, typename Tint>
StabChain<Telt, Tidx_label>
Kernel_IntersectionNormalClosurePermGroup(const StabChain<Telt, Tidx_label> &G,
                                          const StabChain<Telt, Tidx_label> &H,
                                          const Tint &size) {
  using Tidx = typename Telt::Tidx;
  Tidx n = G->comm->n;
  std::vector<Telt> LGen_G = Kernel_GeneratorsOfGroup(G);
  std::vector<Telt> LGen_H = Kernel_GeneratorsOfGroup(H);
  return Kernel_IntersectionNormalClosurePermGroup_LGen<Telt, Tidx_label, Tint>(
      n, LGen_G, LGen_H, size);
}

template <typename Telt, typename Tidx_label>
struct InjectiveRestrictionHomomorphism_base {
  using Tidx = typename Telt::Tidx;
  const StabChain<Telt, Tidx_label> &G;
  const std::vector<Tidx> &subset;
  std::vector<Tidx> subset_rev;
  Tidx miss_val = std::numeric_limits<Tidx>::max();
  StabChain<Telt, Tidx_label> S_restr;
  InjectiveRestrictionHomomorphism_base(
      const StabChain<Telt, Tidx_label> &G,
      const std::vector<typename Telt::Tidx> &subset)
      : G(G), subset(subset) {
    Tidx n = G->comm->n;
    Tidx len = Tidx(subset.size());
    subset_rev = std::vector<Tidx>(n, miss_val);
    for (Tidx i = 0; i < len; i++)
      subset_rev[subset[i]] = i;
    auto f_pt = [&](const Tidx &x) -> Tidx { return subset_rev[x]; };
    auto f_elt = [&](const Telt &u) -> Telt {
      std::vector<Tidx> eList(len);
      for (Tidx i = 0; i < len; i++) {
        Tidx i_big = subset[i];
        Tidx j_big = PowAct(i_big, u);
        Tidx j = subset_rev[j_big];
        eList[i] = j;
      }
      return Telt(std::move(eList));
    };
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
    StabChain<Telt, Tidx_label> G_ptr = G;
    while (G_ptr != nullptr) {
      if (subset_rev[G_ptr->orbit[0]] == miss_val) {
        std::cerr
            << "The value is a missing one. Wrong base for the restriction\n";
        throw PermutalibException{1};
      }
    }
#endif
    S_restr = HomomorphismMapping<Telt, Telt, Tidx_label, decltype(f_pt),
                                  decltype(f_elt)>(G, f_pt, f_elt);
  }
  Telt PreImage_elt(const Telt &eElt) const {
    Tidx n = G->comm->n;
    std::vector<Tidx> eList(n);
    for (Tidx i = 0; i < n; i++)
      eList[i] = i;
    Tidx n_restr = S_restr->comm->n;
    for (Tidx i = 0; i < n_restr; i++) {
      eList[subset[i]] = subset[PowAct(i, eElt)];
    }
    Telt g(eList); // Build g, but keep eList separately for further use.
    Telt res = SiftedPermutation(G, g);
    Telt res_inv = Inverse(res);
    for (Tidx i = 0; i < n; i++)
      if (subset_rev[i] == miss_val)
        eList[i] = PowAct(i, res_inv);
    return Telt(std::move(eList));
  }
  StabChain<Telt, Tidx_label>
  PreImage_grp(const StabChain<Telt, Tidx_label> &U) const {
    auto f_pt = [&](const Tidx &x) -> Tidx { return subset[x]; };
    auto f_elt = [&](const Telt &u) -> Telt { return PreImage_elt(u); };
    return HomomorphismMapping<Telt, Telt, Tidx_label, decltype(f_pt),
                               decltype(f_elt)>(U, f_pt, f_elt);
  }
};

/*
#############################################################################
##
#F  CentralizerTransSymmCSPG()  . . . . .  centralizer of transitive G in S_n
##
##  computes the centralizer of a transitive group G in S_n
##
*/
template <typename Telt, typename Tidx_label>
std::pair<std::vector<Telt>, typename Telt::Tidx>
CentralizerTransSymmCSPG(const StabChain<Telt, Tidx_label> &S,
                         const typename Telt::Tidx &x, Face &L) {
  using Tidx = typename Telt::Tidx;
  size_t L_size = L.count();
  Tidx n = S->comm->n;
  std::vector<Telt> LGen = Kernel_GeneratorsOfGroup(S);
  if (IsTrivialListGen(LGen)) {
    return {};
  }

  // the centralizer of G is semiregular, acting transitively on L
  std::vector<Telt> gens;
  std::pair<std::vector<Tidx>, Face> epair = OrbitPerms(gens, n, x);
  while (epair.first.size() < L_size) {
    // construct element of centralizer which carries x to new point in L
    std::vector<Tidx> gen(n);
    Tidx y = Tidx(L.find_first());
    for (Tidx z = 0; z < n; z++) {
      Telt h = InverseRepresentative(S, z);
      gen[z] = SlashAct(y, h);
    }
    gens.emplace_back(std::move(Telt(std::move(gen))));
    epair = OrbitPerms(gens, n, x);
    for (auto &eVal : epair.first)
      L[eVal] = 0;
  }
  return {std::move(gens), Tidx(L_size)};
}

template <typename Telt, typename Tidx_label>
std::pair<std::vector<Telt>, typename Telt::Tidx>
CentralizerTransSymmCSPG_direct(const StabChain<Telt, Tidx_label> &S) {
  using Tidx = typename Telt::Tidx;
  Tidx n = S->comm->n;
  Tidx x = S->orbit[0];
  std::vector<Tidx> LMoved =
      MovedPoints(Kernel_GeneratorsOfGroup(S->stabilizer), n);
  Face L(n);
  for (Tidx i = 0; i < n; i++)
    L[i] = 1;
  for (auto &eVal : LMoved)
    L[eVal] = 0;
  return CentralizerTransSymmCSPG(S, x, L);
}

/*
#############################################################################
##
#M  Centre( <G> ) . . . . . . . . . . . . . . . center of a permutation group
##
##  constructs the center of G.
##  Reference: Beals-Seress, 24th Symp. on Theory of Computing 1992, sect. 9
##
*/
template <typename Telt, typename Tidx_label, typename Tint>
StabChain<Telt, Tidx_label>
Kernel_CentreSubgroup(const StabChain<Telt, Tidx_label> &G) {
  using Tidx = typename Telt::Tidx;
  Telt id = G->comm->identity;
  Tidx n = G->comm->n;
  Tidx miss_val = std::numeric_limits<Tidx>::max();
  std::vector<Telt> LGen = Kernel_GeneratorsOfGroup(G);
  if (IsTrivialListGen(LGen)) {
    StabChainOptions<Tint, Telt> options1 = GetStandardOptions<Tint, Telt>(n);
    return StabChainOp_listgen<Telt, Tidx_label, Tint>({}, options1);
  }
  std::vector<Tidx> base = BaseStabChain(G);
  std::vector<std::vector<Tidx>> orbits = OrbitsPerms(LGen, n);

  auto get_centralizer_transitive_case =
      [](const StabChain<Telt, Tidx_label> &U) -> StabChain<Telt, Tidx_label> {
    Tidx len = U->comm->n;
    std::pair<std::vector<Telt>, Tidx> pair_centr =
        CentralizerTransSymmCSPG_direct(U);
    if (IsTrivialListGen(pair_centr.first)) {
      StabChainOptions<Tint, Telt> options1 =
          GetStandardOptions<Tint, Telt>(Telt(len));
      return StabChainOp_listgen<Telt, Tidx_label, Tint>({}, options1);
    } else {
      Tint order_p_size = pair_centr.second * Order<Telt, Tidx_label, Tint>(U);
      std::vector<Telt> LGen_U = Kernel_GeneratorsOfGroup(U);
      return Kernel_IntersectionNormalClosurePermGroup_LGen<Telt, Tidx_label,
                                                            Tint>(
          len, LGen_U, pair_centr.first, order_p_size);
    }
  };
  // handle case of transitive G directly
  if (orbits.size() == 1)
    return get_centralizer_transitive_case(G);

  // for intransitive G, find which orbit contains which
  // points of permutation domain
  std::vector<Tidx> reps(n);
  Tidx n_orbit = Tidx(orbits.size());
  for (Tidx i_orbit = 0; i_orbit < n_orbit; i_orbit++)
    for (auto &eVal : orbits[i_orbit])
      reps[eVal] = i_orbit;

  // take union of significant orbits which contain base points
  Face significant_v(n_orbit);
  std::vector<Tidx> significant;
  Tidx len_base = Tidx(base.size());
  std::vector<Tidx> domain;
  for (Tidx i_base = 0; i_base < len_base; i_base++) {
    Tidx i_orbit = reps[base[i_base]];
    if (significant_v[i_orbit] == 0) {
      significant_v[i_orbit] = 1;
      significant.push_back(i_orbit);
      domain.insert(domain.end(), orbits[i_orbit].begin(),
                    orbits[i_orbit].end());
    }
  }
  Tidx len = Tidx(domain.size());

  auto compute_centr = [&](const StabChain<Telt, Tidx_label> &GG,
                           const std::vector<std::vector<Tidx>> &orbits)
      -> StabChain<Telt, Tidx_label> {
    // handle case of transitive GG directly
    if (orbits.size() == 1) {
      return get_centralizer_transitive_case(GG);
    }
    Tidx len = GG->comm->n;
    std::vector<Telt> hgens;
    Tint order = 1;

    /* case of intransitive GG
       for each orbit of GG, construct generators of centralizer of GG in
       Sym(orbit).  hgens is a list of generators for the direct product of
       these centralizers.
       the group generated by hgens contains the center of GG */
    std::vector<Telt> LGen_GG = Kernel_GeneratorsOfGroup(GG);
    for (auto &orbit : orbits) {
      Tidx len_o = Tidx(orbit.size());
      std::vector<Tidx> orbit_rev(len, miss_val);
      for (Tidx i = 0; i < len_o; i++)
        orbit_rev[orbit[i]] = i;
      StabChainOptions<Tint, Telt> options1 =
          GetStandardOptions<Tint, Telt>(Telt(len_o));
      options1.base.push_back(0);
      std::vector<Telt> LGen_o;
      for (auto &eGen : LGen_GG) {
        std::vector<Tidx> eList(len_o);
        for (Tidx i = 0; i < len_o; i++) {
          Tidx i_big = orbit[i];
          Tidx j_big = PowAct(i_big, eGen);
          Tidx j = orbit_rev[j_big];
          eList[i] = j;
        }
        Telt x = Telt(std::move(eList));
        LGen_o.emplace_back(std::move(x));
      }
      StabChain<Telt, Tidx_label> GGG =
          StabChainOp_listgen<Telt, Tidx_label, Tint>(LGen_o, options1);
      std::pair<std::vector<Telt>, Tidx> pair_centr =
          CentralizerTransSymmCSPG_direct(GGG);
      if (!IsTrivialListGen(pair_centr.first)) {
        order *= Tint(pair_centr.second);
        for (auto &eGen : pair_centr.first) {
          std::vector<Tidx> eList(len);
          for (Tidx i = 0; i < len; i++)
            eList[i] = i;
          for (Tidx i = 0; i < len_o; i++)
            eList[orbit[i]] = orbit[PowAct(i, eGen)];
          Telt x = Telt(std::move(eList));
          hgens.emplace_back(std::move(x));
        }
      }
    }
    if (order == 1) {
      StabChainOptions<Tint, Telt> options1 =
          GetStandardOptions<Tint, Telt>(Telt(len));
      return StabChainOp_listgen<Telt, Tidx_label, Tint>({}, options1);
    } else {
      Tint order_p_size = order * Order<Telt, Tidx_label, Tint>(GG);
      return Kernel_IntersectionNormalClosurePermGroup_LGen<Telt, Tidx_label,
                                                            Tint>(
          len, LGen_GG, hgens, order_p_size);
    }
  };

  // We restrict G to significant orbits
  std::vector<std::vector<Tidx>> orbits_red;
  if (n == len) {
    for (auto &eVal : significant)
      orbits_red.emplace_back(std::move(orbits[eVal]));
    return compute_centr(G, orbits_red);
  } else {
    InjectiveRestrictionHomomorphism_base<Telt, Tidx_label> InjResHom(G,
                                                                      domain);
    for (auto &eVal : significant) {
      std::vector<Tidx> &V = orbits[eVal];
      std::vector<Tidx> orbit_red;
      orbit_red.reserve(V.size());
      for (auto &eIdx : V)
        orbit_red.push_back(InjResHom.subset_rev[eIdx]);
      orbits_red.emplace_back(std::move(orbit_red));
    }
    StabChain<Telt, Tidx_label> cent =
        compute_centr(InjResHom.S_restr, orbits_red);
    return InjResHom.PreImage_grp(cent);
  }
}

// clang-format off
}  // namespace permutalib
#endif  // SRC_GAP_NORMALSTRUCTURE_H_
// clang-format on
