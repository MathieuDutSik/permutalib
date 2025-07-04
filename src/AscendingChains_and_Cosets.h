// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GAP_ASCENDINGCHAINS_AND_COSETS_H_
#define SRC_GAP_ASCENDINGCHAINS_AND_COSETS_H_

// clang-format off
#include "TestingFct.h"
#include <limits>
#include <list>
#include <unordered_set>
#include <utility>
#include <vector>
// clang-format on

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
     We compute the transversals, so we use homomorphisms in any case.
  ---S: Splitting G Id K into smaller double cosets is the same as decomposing G.
  ---S: The ascending chain being created is indeed ascending and the algorithm
  starts with the larger group and then downward.
  ---So, that is all about right cosets.

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

#ifdef DEBUG
#define DEBUG_ASCENDING_CHAINS_COSETS
#define DEBUG_SPAN_DOUBLE_COSETS
#endif

#ifdef TIMINGS
#define TIMINGS_ASCENDING_CHAINS_COSETS
#endif

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
Kernel_AscendingChainSingle(StabChain<Telt, Tidx_label> const &G) {
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
#ifdef TIMINGS_ASCENDING_CHAINS_COSETS
  MicrosecondTime_perm time;
#endif
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
#ifdef TIMINGS_ASCENDING_CHAINS_COSETS
  std::cerr << "|ACC: Kernel_RightTransversal_Direct_f, id|=" << time << "\n";
#endif
  std::vector<Telt> LGen = Kernel_SmallGeneratingSet<Telt,Tidx_label,Tint>(G);
#ifdef TIMINGS_ASCENDING_CHAINS_COSETS
  std::cerr << "|ACC: Kernel_RightTransversal_Direct_f, LGen|=" << time << "\n";
#endif
  Tint size_G = Order<Telt,Tidx_label,Tint>(G);
  Tint size_H = Order<Telt,Tidx_label,Tint>(H);
  Tint index = size_G / size_H;
#ifdef TIMINGS_ASCENDING_CHAINS_COSETS
  std::cerr << "|ACC: Kernel_RightTransversal_Direct_f, size_G / size_H / index|=" << time << "\n";
#endif
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
    std::cerr << "ACC: n_coset=" << n_coset << " index=" << index << "\n";
    std::cerr << "ACC: G=[";
    for (auto & eElt : get_all_elements(G))
      std::cerr << " " << eElt;
    std::cerr << " ]\n";
    std::cerr << "ACC: H=[";
    for (auto & eElt : get_all_elements(H))
      std::cerr << " " << eElt;
    std::cerr << " ]\n";
    std::cerr << "ACC: ListTransversal=[";
    for (auto & eElt : ListTransversal)
      std::cerr << " " << eElt;
    std::cerr << " ]\n";
    std::cerr << "ACC: The enumeration found a wrong number of cosets\n";
    throw PermutalibException{1};
  }
#endif
#ifdef TIMINGS_ASCENDING_CHAINS_COSETS
  std::cerr << "|ACC: Kernel_RightTransversal_Direct_f, ListTransversal|=" << time << "\n";
#endif
  return ListTransversal;
}



template<typename Tidx>
std::vector<std::pair<Tidx,Tidx>> get_belonging_vector(std::vector<std::vector<Tidx>> const& orbs, Tidx const& n) {
  std::vector<std::pair<Tidx,Tidx>> Vbelong(n);
  Tidx i_block = 0;
  for (auto & orb : orbs) {
    Tidx len_orb = orb.size();
    for (Tidx i_elt=0; i_elt<len_orb; i_elt++) {
      Tidx val = orb[i_elt];
      Vbelong[val] = {i_block, i_elt};
    }
    i_block += 1;
  }
  return Vbelong;
}

template<typename Telt, typename Tidx_label, typename Tint>
struct AscendingEntry {
  StabChain<Telt, Tidx_label> g;
  std::vector<std::vector<typename Telt::Tidx>> orbs;
  std::vector<std::pair<typename Telt::Tidx,typename Telt::Tidx>> Vbelong;
  std::vector<Telt> l_gens_small;
  Tint ord;
};

template<typename Telt, typename Tidx_label, typename Tint>
AscendingEntry<Telt,Tidx_label,Tint> get_ascending_entry(StabChain<Telt, Tidx_label> const& g) {
  using Tidx = typename Telt::Tidx;
  Tidx n_act = g->comm->n;
#ifdef DEBUG_ASCENDING_CHAINS_COSETS_DISABLE
  std::cerr << "ACC: get_ascending_entry, step 1\n";
#endif
  std::vector<Telt> l_gens_small = Kernel_SmallGeneratingSet<Telt, Tidx_label, Tint>(g);
#ifdef DEBUG_ASCENDING_CHAINS_COSETS_DISABLE
  std::cerr << "ACC: get_ascending_entry, step 2\n";
#endif
  std::vector<std::vector<typename Telt::Tidx>> orbs = OrbitsPerms(l_gens_small, n_act);
#ifdef DEBUG_ASCENDING_CHAINS_COSETS_DISABLE
  std::cerr << "ACC: get_ascending_entry, step 3\n";
#endif
  std::vector<std::pair<Tidx,Tidx>> Vbelong = get_belonging_vector(orbs, n_act);
#ifdef DEBUG_ASCENDING_CHAINS_COSETS_DISABLE
  std::cerr << "ACC: get_ascending_entry, step 4\n";
#endif
  Tint ord = Order<Telt,Tidx_label,Tint>(g);
#ifdef DEBUG_ASCENDING_CHAINS_COSETS_DISABLE
  std::cerr << "ACC: get_ascending_entry, step 5\n";
#endif
  return {std::move(g), std::move(orbs), std::move(Vbelong), std::move(l_gens_small), std::move(ord)};
}

template<typename Telt, typename Tidx_label, typename Tint>
void print_orb_sizes(AscendingEntry<Telt,Tidx_label,Tint> const& u, std::string const& name) {
  std::map<size_t, size_t> map_orbsize_mult;
  for (auto & orb: u.orbs) {
    map_orbsize_mult[orb.size()] += 1;
  }
  std::cerr << name << " |orbs| =";
  for (auto & kv: map_orbsize_mult) {
    std::cerr << " (" << kv.first << "," << kv.second << ")";
  }
  std::cerr << "\n";
}

template<typename Telt, typename Tidx_label, typename Tint>
void print_orbits(AscendingEntry<Telt,Tidx_label,Tint> const& u, std::string const& name) {
  std::map<size_t, std::vector<std::vector<size_t>>> map_mult_orbit;
  for (auto & orb: u.orbs) {
    std::vector<size_t> orb_s;
    for (auto & val: orb) {
      orb_s.push_back(val);
    }
    map_mult_orbit[orb.size()].push_back(orb_s);
  }
  std::cerr << name << " orbs =";
  for (auto & kv : map_mult_orbit) {
    for (auto & orb_s : kv.second) {
      std::cerr << " [";
      bool is_first = true;
      for (auto & val : orb_s) {
        if (!is_first)
          std::cerr << ",";
        is_first = false;
        std::cerr << val;
      }
      std::cerr << "]";
    }
  }
  std::cerr << "\n";
}

template<typename Telt, typename Tidx_label, typename Tint>
BlockDecomposition<typename Telt::Tidx> get_block_decomposition(AscendingEntry<Telt,Tidx_label,Tint> const& ent_H,
                                                                AscendingEntry<Telt,Tidx_label,Tint> const& ent_G, std::vector<typename Telt::Tidx> const& orb_G) {
  using Tidx = typename Telt::Tidx;
  Tidx miss_val = std::numeric_limits<Tidx>::max();
  size_t len = orb_G.size();
  Tidx n_act = ent_H.g->comm->n;
  std::vector<std::vector<Tidx>> orbs_H;
  Face f(len);
  std::vector<Tidx> map_vert_block_H(n_act, miss_val);
  Tidx i_block_dec = 0;
  for (size_t i=0; i<len; i++) {
    Tidx val = orb_G[i];
    if (f[i] == 0) {
      Tidx i_block_H = ent_H.Vbelong[val].first;
      std::vector<Tidx> const& orb_H = ent_H.orbs[i_block_H];
      for (auto & val : orb_H) {
        Tidx pos = ent_G.Vbelong[val].second;
        map_vert_block_H[val] = i_block_dec;
        f[pos] = 1;
      }
      i_block_dec += 1;
      orbs_H.push_back(orb_H);
    }
  }
  return {orbs_H, map_vert_block_H};
}

template<typename Telt, typename Tidx_label, typename Tint>
BlockDecomposition<typename Telt::Tidx> get_supercoarse_block_decomposition(AscendingEntry<Telt,Tidx_label,Tint> const& ent_H, std::vector<typename Telt::Tidx> const& orb_G) {
  using Tidx = typename Telt::Tidx;
  Tidx miss_val = std::numeric_limits<Tidx>::max();
  Tidx n_act = ent_H.g->comm->n;
  std::vector<std::vector<Tidx>> orbs_G{orb_G};
  std::vector<Tidx> map_vert_block_G(n_act, miss_val);
  for (auto & val : orb_G) {
    map_vert_block_G[val] = 0;
  }
  return {orbs_G, map_vert_block_G};
}

template<typename Telt>
bool is_alternating(std::vector<typename Telt::Tidx> const& v, Telt const& elt, typename Telt::Tidx const& n_act) {
  using Tidx = typename Telt::Tidx;
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
  //  std::cerr << "n_act=" << static_cast<size_t>(n_act) << " elt=" << elt << "\n";
#endif
  Tidx miss_val = std::numeric_limits<Tidx>::max();
  std::vector<Tidx> Vmap(n_act, miss_val);
  Tidx pos = 0;
  for (auto & val : v) {
    Vmap[val] = pos;
    pos += 1;
  }
  size_t len = v.size();
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
  //  std::cerr << "len=" << len << "\n";
#endif
  Face f(len);
  int sign = 1;
  for (size_t i = 0; i<len; i++) {
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
    //    std::cerr << "i=" << i << "\n";
#endif
    if (f[i] == 0) {
      Tidx val_first = v[i];
      size_t len_cycle = 0;
      Tidx val_curr = val_first;
      while(true) {
        Tidx pos = Vmap[val_curr];
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
        //        std::cerr << "pos=" << static_cast<size_t>(pos) << "\n";
#endif
        f[pos] = 1;
        val_curr = elt.at(val_curr);
        len_cycle += 1;
        if (val_curr == val_first) {
          break;
        }
      }
      size_t res = len_cycle % 2;
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
      //      std::cerr << "len_cycle=" << len_cycle << " res=" << res << " sign=" << sign << "\n";
#endif
      if (res == 0) {
        sign *= -1;
      }
    }
  }
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
  //  std::cerr << "Final sign=" << sign << "\n";
#endif
  if (sign == 1) {
    return true;
  }
  return false;
}

/*
  Try to find an intermediate group by finding some alternating group.
  We use a variant of Schreier's lemma in order to get the generators of the
  subgroup.
 */
template<typename Telt, typename Tidx_label, typename Tint>
std::optional<StabChain<Telt, Tidx_label>> Kernel_AscendingChain_Alt(AscendingEntry<Telt,Tidx_label,Tint> const& ent_H,
                                                                     AscendingEntry<Telt,Tidx_label,Tint> const& ent_G) {
  using Tidx = typename Telt::Tidx;
  Tidx n_act = ent_H.g->comm->n;
  auto is_grp_symmetric=[&](std::vector<Telt> const& LGen, std::vector<Tidx> const& orb) -> std::optional<Telt> {
    for (auto & eGen : LGen) {
      if (!is_alternating(orb, eGen, n_act)) {
        return eGen;
      }
    }
    return {};
  };
  Telt id = ent_H.g->comm->identity;
  for (auto & orb : ent_G.orbs) {
    std::optional<Telt> opt_G = is_grp_symmetric(ent_G.l_gens_small, orb);
    if (opt_G) { // We have a element that is not alternating.
      Telt const& eGenRef = *opt_G;
      Telt eGenRefInv = Inverse(eGenRef);
      std::optional<Telt> opt_H = is_grp_symmetric(ent_H.l_gens_small, orb);
      if (!opt_H) { // The subgroup is an alternating subgroup
        Tint size = ent_H.ord * 2;
        if (size != ent_G.ord) { // If equal then no chance to find an intermediate subgroup
          size_t n_gen = ent_G.l_gens_small.size();
          std::vector<int> l_sign;
          for (auto & u: ent_G.l_gens_small) {
            bool test = is_alternating(orb, u, n_act);
            if (test) {
              l_sign.push_back(1);
            } else {
              l_sign.push_back(-1);
            }
          }
          std::vector<Telt> LCos{id, eGenRef};
          std::vector<int> l_cos_sign{1, -1};
          std::vector<Telt> LGenAlt;
          for (size_t i_cos=0; i_cos<2; i_cos++) {
            Telt const& eCos = LCos[i_cos];
            int cos_sign = l_cos_sign[i_cos];
            for (size_t i_gen=0; i_gen<n_gen; i_gen++) {
              Telt u = ent_G.l_gens_small[i_gen];
              int e_sign = l_sign[i_gen];
              int tot_sign = cos_sign * e_sign;
              if (tot_sign == 1) {
                Telt eProd = eCos * u;
                LGenAlt.push_back(eProd);
              } else {
                Telt eProd = eCos * u * eGenRefInv;
                LGenAlt.push_back(eProd);
              }
            }
          }
          StabChainOptions<Tint, Telt> options = GetStandardOptions<Tint, Telt>(id);
          StabChain<Telt,Tidx_label> gAlt = StabChainOp_listgen<Telt, Tidx_label, Tint>(LGenAlt, options);
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
          Tint size_target = 2 * Order<Telt,Tidx_label,Tint>(gAlt);
          if (size_target != ent_G.ord) {
            std::cerr << "ACC: size_target is not of the right size\n";
            throw PermutalibException{1};
          }
#endif
          return gAlt;
        }
      }
    }
  }
  return {};
}

/*
 Try to find an intermediate group by using block decomposition.
 That is if some orbit can be joined then we will find an intermediate subgroup.
 */
template <typename Telt, typename Tidx_label, typename Tint>
std::optional<StabChain<Telt, Tidx_label>> Kernel_AscendingChain_Block(AscendingEntry<Telt,Tidx_label,Tint> const& ent_H,
                                                                       AscendingEntry<Telt,Tidx_label,Tint> const& ent_G) {
  using Tidx = typename Telt::Tidx;
  Tidx n_act = ent_H.g->comm->n;
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
  //  print_orbits(ent_H, "H");
  //  print_orbits(ent_G, "G");
  print_orb_sizes(ent_H, "H");
  print_orb_sizes(ent_G, "G");
#endif
  for (auto & orb_G : ent_G.orbs) {
    BlockDecomposition<Tidx> blk1 = get_block_decomposition(ent_H, ent_G, orb_G);
    //
    BlockDecomposition<Tidx> blk2 = get_supercoarse_block_decomposition(ent_H, orb_G);
    std::optional<BlockDecomposition<Tidx>> opt =
      FindIntermediateBlockDecomposition(ent_G.l_gens_small, blk1, blk2);
    if (opt) {
      std::vector<Tidx> orb_found = opt->ListBlocks[0];
      Face f(n_act);
      for (auto & val : orb_found) {
        f[val] = 1;
      }
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
      std::cerr << "ACC: orb_found = [";
      bool is_first = true;
      for (auto & val : orb_found) {
        if (!is_first)
          std::cerr << ",";
        is_first = false;
        std::cerr << static_cast<size_t>(val);
      }
      std::cerr << "]\n";
#endif
      StabChain<Telt,Tidx_label> gBlk = Kernel_Stabilizer_OnSets<Telt,Tidx_label,Tint>(ent_G.g, f);
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
      Tint size_blk = Order<Telt,Tidx_label,Tint>(gBlk);
      Tint size_G = Order<Telt,Tidx_label,Tint>(ent_G.g);
      Tint size_H = Order<Telt,Tidx_label,Tint>(ent_H.g);
      if (size_H == size_blk || size_blk == size_G) {
        std::cerr << "ACC: The sizes are not as they should be\n";
        throw PermutalibException{1};
      }
      if (!Kernel_IsSubgroup(gBlk, ent_H.g)) {
        std::cerr << "ACC: size_blk=" << size_blk << " size_G=" << size_G << " size_H=" << size_H << "\n";
        std::cerr << "ACC: |orb_G|=" << orb_G.size() << " |orb_found|=" << orb_found.size() << "\n";
        std::cerr << "ACC: H should be a subgroup of gBlk\n";
        throw PermutalibException{1};
      }
#endif
      return gBlk;
    }
  }
  return {};
}



/*
  Try to find an intermediate group by having some generators inserted
 */
template <typename Telt, typename Tidx_label, typename Tint>
std::optional<StabChain<Telt, Tidx_label>> Kernel_AscendingChain_Gens(AscendingEntry<Telt,Tidx_label,Tint> const& ent_H,
                                                                      AscendingEntry<Telt,Tidx_label,Tint> const& ent_G) {
  Tint size_G = ent_G.ord;
  Tint size_H = ent_H.ord;
  Telt id = ent_H.g->comm->identity;
  StabChainOptions<Tint, Telt> options = GetStandardOptions<Tint, Telt>(id);
  for (auto & eGen: Kernel_GeneratorsOfGroup(ent_G.g)) {
    if (!IsElementInStabChain(ent_H.g, eGen)) {
      std::vector<Telt> LGen = ent_H.l_gens_small;
      LGen.push_back(eGen);
      StabChain<Telt,Tidx_label> gExt = StabChainOp_listgen<Telt, Tidx_label, Tint>(LGen, options);
      Tint size_gExt = Order<Telt,Tidx_label,Tint>(gExt);
      if (size_H < size_gExt && size_gExt < size_G) {
        return gExt;
      }
    }
  }
  return {};
}




  /*
    Computes the orbits on the points.
   */
template <typename Telt, typename Tidx_label, typename Tint>
std::optional<StabChain<Telt, Tidx_label>> Kernel_AscendingChain_Subset(AscendingEntry<Telt,Tidx_label,Tint> const& ent_H,
                                                                        AscendingEntry<Telt,Tidx_label,Tint> const& ent_G) {
  using Tidx = typename Telt::Tidx;
  Tidx n_act = ent_G.g->comm->n;
  Tint const& size_G = ent_G.ord;
  Tint const& size_H = ent_H.ord;
  std::vector<std::vector<Tidx>> const& orbs_H = ent_H.orbs;
  for (auto & orb_H: orbs_H) {
    Face f(n_act);
    for (auto & pt : orb_H) {
      f[pt] = 1;
    }
    StabChain<Telt,Tidx_label> sub = Kernel_Stabilizer_OnSets<Telt,Tidx_label,Tint>(ent_G.g, f);
    Tint size_sub = Order<Telt,Tidx_label,Tint>(sub);
    if (size_sub != size_G && size_sub != size_H) {
      return sub;
    }
  }
  return {};
}


// The combined use of the 4 methods for finding intermediate subgroups.
// All of those methods can be iterated since they all can potentially
// give more intermediate subgroups when reapplied.
template <typename Telt, typename Tidx_label, typename Tint>
std::optional<StabChain<Telt, Tidx_label>> Kernel_AscendingChain_All(AscendingEntry<Telt,Tidx_label,Tint> const& ent_H,
                                                                     AscendingEntry<Telt,Tidx_label,Tint> const& ent_G) {
#ifdef TIMINGS_ASCENDING_CHAINS_COSETS
  MicrosecondTime_perm time;
#endif
  //
  // First the alternating method as it is the cheapest one.
  //
  std::optional<StabChain<Telt, Tidx_label>> opt2 =
    Kernel_AscendingChain_Alt(ent_H, ent_G);
  if (opt2) {
#ifdef TIMINGS_ASCENDING_CHAINS_COSETS
    std::cerr << "|ACC: Kernel_AscendingChain_Alt, success|=" << time << "\n";
#endif
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
    std::cerr << "ACC: Finding an intermediate subgroup by Alt method\n";
#endif
    return *opt2;
  }
#ifdef TIMINGS_ASCENDING_CHAINS_COSETS
  std::cerr << "|ACC: Kernel_AscendingChain_Alt, failure|=" << time << "\n";
#endif
  //
  // Next the Block method as it is not that expensive, especially when it fails.
  //
  std::optional<StabChain<Telt, Tidx_label>> opt3 =
    Kernel_AscendingChain_Block(ent_H, ent_G);
  if (opt3) {
#ifdef TIMINGS_ASCENDING_CHAINS_COSETS
    std::cerr << "|ACC: Kernel_AscendingChain_Block, success|=" << time << "\n";
#endif
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
    std::cerr << "ACC: Finding an intermediate subgroup by the Block method\n";
#endif
    return *opt3;
  }
#ifdef TIMINGS_ASCENDING_CHAINS_COSETS
  std::cerr << "|ACC: Kernel_AscendingChain_Block, failure|=" << time << "\n";
#endif
  //
  // Then the subset method as it is very powerful but relatively expensive.
  //
  std::optional<StabChain<Telt, Tidx_label>> opt1 =
    Kernel_AscendingChain_Subset(ent_H, ent_G);
  if (opt1) {
#ifdef TIMINGS_ASCENDING_CHAINS_COSETS
    std::cerr << "|ACC: Kernel_AscendingChain_Subset, success|=" << time << "\n";
#endif
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
    std::cerr << "ACC: Finding an intermediate subgroup by the Subset method\n";
#endif
    return *opt1;
  }
#ifdef TIMINGS_ASCENDING_CHAINS_COSETS
  std::cerr << "|ACC: Kernel_AscendingChain_Subset, failure|=" << time << "\n";
#endif
  //
  // Latest is the gens method which works rarely and is relatively expensive
  //
  std::optional<StabChain<Telt, Tidx_label>> opt4 =
    Kernel_AscendingChain_Gens(ent_H, ent_G);
  if (opt4) {
#ifdef TIMINGS_ASCENDING_CHAINS_COSETS
    std::cerr << "|ACC: Kernel_AscendingChain_Gens, success|=" << time << "\n";
#endif
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
    std::cerr << "ACC: Finding an intermediate subgroup by Gens method\n";
#endif
    return *opt4;
  }
#ifdef TIMINGS_ASCENDING_CHAINS_COSETS
  std::cerr << "|ACC: Kernel_AscendingChain_Gens, failure|=" << time << "\n";
#endif
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
  std::cerr << "ACC: Failing to find an intermediate subgroup |H|=" << ent_H.ord << " |G|=" << ent_G.ord << "\n";
#endif
  //
  //
  //
  return {};
}


  /*
    Computes an ascending chain of subgroups from H to G.
    It is an heuristic program and there is no guarantee that we will get an optimal one.
    Order of function call (H, G) vs (G, H) is reverted compared to GAP.
                      --- What GAP is doing ---
    In GAP the following is done:
    ---Compute an ascending chain from the trivial subgroup and if one can restricy
    it to one for the pair of subgroup, do it. Otherwise apply AscendingChainOp.
    ---For the case of Cyclic Group C_n in the symmetric group S_n, we find in GAP
    some strange patterns. But they are better than what
    ---Apply normalizing steps, that seems completely natural.
    ---For H < G Compute Centralizer(G, Centre(H)). Could work. Require to implement
       Centralizer and Centre.
    ---Can we find index 2 subgroups easily?
       Probably not. But since we work with symmetric group, we can find easily if
       the permutations are all alternating or not. But then, computing G \cap A_n
       could be problematic.
    ---IntermediateSubgroup function that uses
       ---Maximal subgroups (technique excluded)
       ---Block decomposition on the action on the right cosets.

                      -------- Tactics --------
    Apply the following tactics:
    ---Orbit set stabilizers and see what we obtain. Basic strategy but should be powerful.
    ---Block decompositions strategies.
    ---GeneratorsOfGroup of course could help us get some intermediate.
    ---Computing the alternating subgroups of the action on orbits.
    Objections:
    ---If the index is prime then nothing more can be done. That is obvious. But could we
    have more sophisticated criterion?
    ---We have to factor out the cost of computing. Maybe sometimes it is not worth to
    overly refine the decomposition. So, we may have something similar to the
    ---The contruction Centralizer(G, Centre(H)) seems to be pretty efficient in finding
    intermediate subgroups. It does not return a normalizing group, it just returns something
    that it invariant under conjugation.
    ---By contrast the Normalizer idea does not seem to be very powerful.
    ---The alternating subgroup technique does not appear to ever succeed.
   */
template <typename Telt, typename Tidx_label, typename Tint>
std::vector<StabChain<Telt, Tidx_label>> Kernel_AscendingChainPair(StabChain<Telt, Tidx_label> const &H,
                                                                   StabChain<Telt, Tidx_label> const &G) {
#ifdef TIMINGS_ASCENDING_CHAINS_COSETS
  MicrosecondTime_perm time_total;
#endif
  AscendingEntry<Telt,Tidx_label,Tint> ent_H = get_ascending_entry<Telt,Tidx_label,Tint>(H);
  AscendingEntry<Telt,Tidx_label,Tint> ent_G = get_ascending_entry<Telt,Tidx_label,Tint>(G);
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
  if (!Kernel_IsSubgroup(ent_G.g, ent_H.g)) {
    std::cerr << "ACC: Kernel_AscendingChainPair, H should be a subgroup of G\n";
  }
#endif

#ifdef DEBUG_ASCENDING_CHAINS_COSETS
  std::cerr << "ACC: Kernel_AscendingChainPair, ent_H.ord=" << ent_H.ord << " ent_G.ord=" << ent_G.ord << "\n";
#endif
  if (ent_H.ord == ent_G.ord) {
    return {H};
  }
  std::vector<AscendingEntry<Telt,Tidx_label,Tint>> l_grp{ent_H, ent_G};
  size_t pos = 0;
  while(true) {
    auto get_intermediate=[&]() -> std::optional<StabChain<Telt, Tidx_label>> {
      Tint index = l_grp[pos+1].ord / l_grp[pos].ord;
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
      std::cerr << "ACC: Kernel_AscendingChainPair, l_grp[pos+1].ord=" << l_grp[pos+1].ord << " l_grp[pos].ord=" << l_grp[pos].ord << "\n";
      std::cerr << "ACC: Kernel_AscendingChainPair, Before IsPrime_loc, pos=" << pos << " index=" << index << "\n";
#endif
      bool is_prime = IsPrime_loc(index);
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
      std::cerr << "ACC: index=" << index << " is_prime=" << is_prime << "\n";
#endif
      if (is_prime) { // Cannot improve when the index is prime
        return {};
      }
      return Kernel_AscendingChain_All<Telt,Tidx_label,Tint>(l_grp[pos], l_grp[pos+1]);
    };
    std::optional<StabChain<Telt, Tidx_label>> opt = get_intermediate();
    if (opt) {
      auto iter = l_grp.begin();
      std::advance(iter, pos + 1);
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
      std::cerr << "ACC: Before get_ascending_entry\n";
#endif
      AscendingEntry<Telt,Tidx_label,Tint> ent = get_ascending_entry<Telt,Tidx_label,Tint>(*opt);
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
      std::cerr << "ACC: After get_ascending_entry\n";
#endif
      l_grp.insert(iter, ent);
    } else { // No method work, going to the next one.
      pos++;
      if (pos + 1 == l_grp.size()) {
        break;
      }
    }
  }
  std::vector<StabChain<Telt, Tidx_label>> chain;
  for (auto & ent : l_grp) {
    chain.push_back(ent.g);
  }
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
  for (size_t i_grp=1; i_grp<chain.size(); i_grp++) {
    if (!Kernel_IsSubgroup(chain[i_grp], chain[i_grp-1])) {
      std::cerr << "ACC: Error in Kernel_AscendingChainPair at i_grp=" << i_grp << "\n";
      throw PermutalibException{1};
    }
  }
#endif
#ifdef TIMINGS_ASCENDING_CHAINS_COSETS
  std::cerr << "|ACC: Kernel_AscendingChainPair|=" << time_total << "\n";
#endif
  return chain;
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

template <typename Telt, typename Tidx_label, typename Tint>
void KernelCheckLeftCosets(StabChain<Telt, Tidx_label> const &G,
                           StabChain<Telt, Tidx_label> const &H,
                           std::vector<Telt> const& ListLeftTransversal) {
  std::vector<Telt> l_elt_g = get_all_elements(G);
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
  std::cerr << "ACC: l_elt_g\n";
#endif
  std::vector<Telt> l_elt_h = get_all_elements(H);
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
  std::cerr << "ACC: l_elt_h\n";
#endif
  std::unordered_set<Telt> set_elt;
  size_t ProdSize = l_elt_h.size() * ListLeftTransversal.size();
  if (ProdSize != l_elt_g.size()) {
    std::cerr << "ACC: |l_elt_g|=" << l_elt_g.size() << "\n";
    std::cerr << "ACC: |l_elt_h|=" << l_elt_h.size() << "\n";
    std::cerr << "ACC: |ListLeftTransversal|=" << ListLeftTransversal.size() << "\n";
    std::cerr << "ACC: ProdSize=" << ProdSize << "\n";
    std::cerr << "ACC: Discrepancy at the order level\n";
    throw PermutalibException{1};
  }
  std::vector<std::unordered_set<Telt>> l_cos;
  for (auto & eElt : ListLeftTransversal) {
    std::unordered_set<Telt> set;
    for (auto & e_h : l_elt_h) {
      Telt eProd = eElt * e_h;
      if (set_elt.count(eProd) == 1) {
        std::cerr << "ACC: The element eProd is already present\n";
        throw PermutalibException{1};
      }
      set_elt.insert(eProd);
      set.insert(eProd);
    }
    l_cos.push_back(set);
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
    std::cerr << "ACC: Now |l_cos|=" << l_cos.size() << "\n";
#endif
  }
  for (size_t i_cos=0; i_cos<l_cos.size(); i_cos++) {
    for (size_t j_cos=i_cos+1; j_cos<l_cos.size(); j_cos++) {
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
      std::cerr << "ACC: Now i_cos=" << i_cos << " j_cos=" << j_cos << "\n";
#endif
      size_t the_int = 0;
      for (auto & val : l_cos[i_cos]) {
        if (l_cos[j_cos].count(val) == 1) {
          the_int += 1;
        }
      }
      if (the_int > 0) {
        std::cerr << "ACC: Intersection between i_cos=" << i_cos << " j_cos=" << j_cos << " has size " << the_int << "\n";
        throw PermutalibException{1};
      }
    }
  }
}

template <typename Telt, typename Tidx_label, typename Tint>
void KernelCheckRightCosets(StabChain<Telt, Tidx_label> const &G,
                            StabChain<Telt, Tidx_label> const &H,
                            std::vector<Telt> const& ListRightTransversal) {
  std::vector<Telt> l_elt_g = get_all_elements(G);
  std::vector<Telt> l_elt_h = get_all_elements(H);
  std::unordered_set<Telt> set_elt;
  size_t ProdSize = l_elt_h.size() * ListRightTransversal.size();
  if (ProdSize != l_elt_g.size()) {
    std::cerr << "ACC: |l_elt_g|=" << l_elt_g.size() << "\n";
    std::cerr << "ACC: |l_elt_h|=" << l_elt_h.size() << "\n";
    std::cerr << "ACC: |ListRightTransversal|=" << ListRightTransversal.size() << "\n";
    std::cerr << "ACC: ProdSize=" << ProdSize << "\n";
    std::cerr << "ACC: Discrepancy at the order level\n";
    throw PermutalibException{1};
  }
  std::vector<std::unordered_set<Telt>> l_cos;
  for (auto & eCos : ListRightTransversal) {
    std::unordered_set<Telt> set;
    for (auto & e_h : l_elt_h) {
      Telt eProd = e_h * eCos;
      if (set_elt.count(eProd) == 1) {
        std::cerr << "ACC: The element eProd is already present\n";
        throw PermutalibException{1};
      }
      set_elt.insert(eProd);
      set.insert(eProd);
    }
    l_cos.push_back(set);
  }
  for (size_t i_cos=0; i_cos<l_cos.size(); i_cos++) {
    for (size_t j_cos=i_cos+1; j_cos<l_cos.size(); j_cos++) {
      size_t the_int = 0;
      for (auto & val : l_cos[i_cos]) {
        if (l_cos[j_cos].count(val) == 1) {
          the_int += 1;
        }
      }
      if (the_int > 0) {
        std::cerr << "ACC: Intersection between i_cos=" << i_cos << " j_cos=" << j_cos << " has size " << the_int << "\n";
        throw PermutalibException{1};
      }
    }
  }
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
struct RightCosetIterator {
private:
  std::vector<std::vector<Telt>> ll_cos;
  std::vector<size_t> l_size;
  std::vector<size_t> l_pos;
  size_t n_level;
  bool is_end;
  Telt result;
  void compute_position() {
#ifdef DEBUG_ASCENDING_CHAINS_COSETS_DISABLE
    std::cerr << "ACC: compute_position, start\n";
    std::cerr << "ACC: l_pos/l_size =";
    for (size_t i_level=0; i_level<n_level; i_level++) {
      std::cerr << " (" << l_pos[i_level] << "|" << l_size[i_level] << ")";
    }
    std::cerr << "\n";
#endif
    result = ll_cos[0][l_pos[0]];
    for (size_t i_level=1; i_level<n_level; i_level++) {
      result *= ll_cos[i_level][l_pos[i_level]];
    }
#ifdef DEBUG_ASCENDING_CHAINS_COSETS_DISABLE
    std::cerr << "ACC: compute_position, end\n";
#endif
  }
  void single_increase() {
    for (size_t i_level=0; i_level<n_level; i_level++) {
      if (l_pos[i_level] < l_size[i_level] - 1) {
        for (size_t j_level=0; j_level<i_level; j_level++) {
          l_pos[j_level] = 0;
        }
        l_pos[i_level] += 1;
        compute_position();
        return;
      }
    }
    is_end = true;
  }
public:
  RightCosetIterator(StabChain<Telt,Tidx_label> const& H, StabChain<Telt,Tidx_label> const& G) {
#ifdef TIMINGS_ASCENDING_CHAINS_COSETS
  MicrosecondTime_perm time;
#endif
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
    std::cerr << "ACC: RightCosetIterator, begin constructor\n";
#endif
    std::vector<StabChain<Telt,Tidx_label>> chain = Kernel_AscendingChainPair<Telt,Tidx_label,Tint>(H, G);
#ifdef TIMINGS_ASCENDING_CHAINS_COSETS
    std::cerr << "|ACC: RightCosetIterator, chain|=" << time << "\n";
#endif
    n_level = chain.size() - 1;
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
    std::cerr << "ACC: n_level=" << n_level << "\n";
    for (size_t i_level=0; i_level<=n_level; i_level++) {
      std::cerr << "ACC: i_level=" << i_level << " ord=" << Order<Telt,Tidx_label,Tint>(chain[i_level]) << "\n";
    }
#endif
    result = G->comm->identity;
#ifdef TIMINGS_ASCENDING_CHAINS_COSETS
    std::cerr << "|ACC: RightCosetIterator, result|=" << time << "\n";
#endif
    for (size_t i_level=0; i_level<n_level; i_level++) {
      std::vector<Telt> l_cos = Kernel_RightTransversal_Direct<Telt,Tidx_label,Tint>(chain[i_level + 1], chain[i_level]);
#ifdef DEBUG_ASCENDING_CHAINS_COSETS_DISABLE
      std::cerr << "ACC: i_level=" << i_level << " |l_cos|=" << l_cos.size()
                << " ord1=" << Order<Telt,Tidx_label,Tint>(chain[i_level])
                << " ord2=" << Order<Telt,Tidx_label,Tint>(chain[i_level + 1])
                << "\n";
      KernelCheckRightCosets<Telt,Tidx_label,Tint>(chain[i_level + 1], chain[i_level], l_cos);
#endif
      result *= l_cos[0];
      ll_cos.push_back(l_cos);
      l_size.push_back(l_cos.size());
      l_pos.push_back(0);
    }
#ifdef TIMINGS_ASCENDING_CHAINS_COSETS
    std::cerr << "|ACC: RightCosetIterator, ll_cos / l_size / l_pos|=" << time << "\n";
#endif
    is_end = false;
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
    std::cerr << "ACC: RightCosetIterator, exit\n";
#endif
  }
  RightCosetIterator() {
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
    std::cerr << "ACC: RightCosetIterator, end constructor\n";
#endif
    is_end = true;
  }
  Telt const& operator*() const {
    return result;
  }
  bool operator==(const RightCosetIterator<Telt,Tidx_label,Tint>& rci) const {
    if (is_end == rci.is_end) {
      return true;
    }
    if (is_end || rci.is_end) {
      return false;
    }
    if (n_level != rci.n_level) {
      return false;
    }
    for (size_t i_level=0; i_level<n_level; i_level++) {
      if (l_pos[i_level] != rci.l_pos[i_level]) {
        return false;
      }
    }
    return true;
  }
  bool operator!=(const RightCosetIterator<Telt,Tidx_label,Tint>& rci) const {
    if (is_end == rci.is_end) {
      return false;
    }
    if (is_end || rci.is_end) {
      return true;
    }
    if (n_level != rci.n_level) {
      return true;
    }
    for (size_t i_level=0; i_level<n_level; i_level++) {
      if (l_pos[i_level] != rci.l_pos[i_level]) {
        return true;
      }
    }
    return false;
  }
  Telt operator++() {
    single_increase();
    return result;
  }
  Telt operator++(int) {
    Telt tmp = result;
    single_increase();
    return tmp;
  }
};

template <typename Telt, typename Tidx_label, typename Tint>
struct KernelRightCosets {
private:
  StabChain<Telt,Tidx_label> H;
  StabChain<Telt,Tidx_label> G;
public:
  using iterator = RightCosetIterator<Telt,Tidx_label,Tint>;
  using const_iterator = RightCosetIterator<Telt,Tidx_label,Tint>;
  KernelRightCosets(StabChain<Telt,Tidx_label> const& _H, StabChain<Telt,Tidx_label> const& _G) : H(_H), G(_G) {
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
    std::cerr << "ACC: KernelRightCosets, constructor\n";
#endif
  }
  RightCosetIterator<Telt,Tidx_label,Tint> begin() const {
    return RightCosetIterator<Telt,Tidx_label,Tint>(H, G);
  }
  RightCosetIterator<Telt,Tidx_label,Tint> end() const {
    return RightCosetIterator<Telt,Tidx_label,Tint>();
  }
};

template <typename Telt, typename Tidx_label, typename Tint>
std::vector<Telt> enumerate_right_cosets(StabChain<Telt,Tidx_label> const& H, StabChain<Telt,Tidx_label> const& G) {
  std::vector<Telt> l_cos;
  KernelRightCosets<Telt, Tidx_label, Tint> rc(H, G);
  for (auto & eCos: rc) {
    l_cos.push_back(eCos);
  }
  return l_cos;
}

template <typename Telt, typename Tidx_label, typename Tint>
struct LeftCosetIterator {
private:
  std::vector<std::vector<Telt>> ll_cos;
  std::vector<size_t> l_size;
  std::vector<size_t> l_pos;
  size_t n_level;
  bool is_end;
  Telt result;
  void compute_position() {
#ifdef DEBUG_ASCENDING_CHAINS_COSETS_DISABLE
    std::cerr << "ACC: compute_position, start\n";
    std::cerr << "ACC: l_pos/l_size =";
    for (size_t i_level=0; i_level<n_level; i_level++) {
      std::cerr << " (" << l_pos[i_level] << "|" << l_size[i_level] << ")";
    }
    std::cerr << "\n";
#endif
    result = ll_cos[0][l_pos[0]];
    for (size_t i_level=1; i_level<n_level; i_level++) {
      result = ll_cos[i_level][l_pos[i_level]] * result;
    }
#ifdef DEBUG_ASCENDING_CHAINS_COSETS_DISABLE
    std::cerr << "ACC: compute_position, end\n";
#endif
  }
  void single_increase() {
    for (size_t i_level=0; i_level<n_level; i_level++) {
      if (l_pos[i_level] < l_size[i_level] - 1) {
        for (size_t j_level=0; j_level<i_level; j_level++) {
          l_pos[j_level] = 0;
        }
        l_pos[i_level] += 1;
        compute_position();
        return;
      }
    }
    is_end = true;
  }
public:
  LeftCosetIterator(StabChain<Telt,Tidx_label> const& H, StabChain<Telt,Tidx_label> const& G) {
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
    std::cerr << "ACC: RightCosetIterator, begin constructor\n";
#endif
    std::vector<StabChain<Telt,Tidx_label>> chain = Kernel_AscendingChainPair<Telt,Tidx_label,Tint>(H, G);
    n_level = chain.size() - 1;
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
    for (size_t i_level=0; i_level<=n_level; i_level++) {
      std::cerr << "ACC: i_level=" << i_level << " ord=" << Order<Telt,Tidx_label,Tint>(chain[i_level]) << "\n";
    }
#endif
    result = G->comm->identity;
    for (size_t i_level=0; i_level<n_level; i_level++) {
      std::vector<Telt> l_cos = Kernel_LeftTransversal_Direct<Telt,Tidx_label,Tint>(chain[i_level + 1], chain[i_level]);
#ifdef DEBUG_ASCENDING_CHAINS_COSETS_DISABLE
      std::cerr << "ACC: i_level=" << i_level << " |l_cos|=" << l_cos.size()
                << " ord1=" << Order<Telt,Tidx_label,Tint>(chain[i_level])
                << " ord2=" << Order<Telt,Tidx_label,Tint>(chain[i_level + 1])
                << "\n";
      KernelCheckLeftCosets<Telt,Tidx_label,Tint>(chain[i_level + 1], chain[i_level], l_cos);
#endif
      result = l_cos[0] * result;
      ll_cos.push_back(l_cos);
      l_size.push_back(l_cos.size());
      l_pos.push_back(0);
    }
    is_end = false;
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
    std::cerr << "ACC: LeftCosetIterator, exit\n";
#endif
  }
  LeftCosetIterator() {
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
    std::cerr << "ACC: LeftCosetIterator, end constructor\n";
#endif
    is_end = true;
  }
  Telt const& operator*() const {
    return result;
  }
  bool operator==(const LeftCosetIterator<Telt,Tidx_label,Tint>& rci) const {
    if (is_end == rci.is_end) {
      return true;
    }
    if (is_end || rci.is_end) {
      return false;
    }
    if (n_level != rci.n_level) {
      return false;
    }
    for (size_t i_level=0; i_level<n_level; i_level++) {
      if (l_pos[i_level] != rci.l_pos[i_level]) {
        return false;
      }
    }
    return true;
  }
  bool operator!=(const LeftCosetIterator<Telt,Tidx_label,Tint>& rci) const {
    if (is_end == rci.is_end) {
      return false;
    }
    if (is_end || rci.is_end) {
      return true;
    }
    if (n_level != rci.n_level) {
      return true;
    }
    for (size_t i_level=0; i_level<n_level; i_level++) {
      if (l_pos[i_level] != rci.l_pos[i_level]) {
        return true;
      }
    }
    return false;
  }
  Telt operator++() {
    single_increase();
    return result;
  }
  Telt operator++(int) {
    Telt tmp = result;
    single_increase();
    return tmp;
  }
};

template <typename Telt, typename Tidx_label, typename Tint>
struct KernelLeftCosets {
private:
  StabChain<Telt,Tidx_label> H;
  StabChain<Telt,Tidx_label> G;
public:
  using iterator = LeftCosetIterator<Telt,Tidx_label,Tint>;
  using const_iterator = LeftCosetIterator<Telt,Tidx_label,Tint>;
  KernelLeftCosets(StabChain<Telt,Tidx_label> const& _H, StabChain<Telt,Tidx_label> const& _G) : H(_H), G(_G) {
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
    std::cerr << "ACC: KernelLeftCosets, constructor\n";
#endif
  }
  LeftCosetIterator<Telt,Tidx_label,Tint> begin() const {
    return LeftCosetIterator<Telt,Tidx_label,Tint>(H, G);
  }
  LeftCosetIterator<Telt,Tidx_label,Tint> end() const {
    return LeftCosetIterator<Telt,Tidx_label,Tint>();
  }
};

template <typename Telt, typename Tidx_label, typename Tint>
std::vector<Telt> enumerate_left_cosets(StabChain<Telt,Tidx_label> const& H, StabChain<Telt,Tidx_label> const& G) {
  std::vector<Telt> l_cos;
  KernelLeftCosets<Telt, Tidx_label, Tint> rc(H, G);
  for (auto & eCos: rc) {
    l_cos.push_back(eCos);
  }
  return l_cos;
}

/*
  Computes the information for G as the union of U g_i V.
  The decomposition is from U = H0 \subset H1 \subset ... Hm = G.
  For that we need to do decompositions of the kind
  H(i+1) = union H(i) g_k
 */
template<typename Telt, typename Tidx_label>
struct DoubleCosetSplitEntry {
  StabChain<Telt,Tidx_label> grp;
  std::vector<Telt> l_cos;
  std::unordered_map<Telt, size_t> map;
  bool is_normal;
};

template<typename Telt>
struct KernelDccEntry {
  Telt cos;
  std::vector<Telt> stab_gens;
};


/*
  This is the central function for computing the double cosets.
  Input:
  ---The DCSE is the description of the level (group, cosets, map)
  ---The KernelDccEntry is the double coset entry.
     de.cos is the coset
     de.stab_gens is a generator set for a group of permutation x
        such that G de.cos x = G de.cos
  ---compute_stabs is whether to compute the stabilizers.
  Output:
  ---The vector of KernelDccEntry is returned into output.
           ------------
 */
template<typename Telt, typename Tidx_label, typename Tint>
std::vector<KernelDccEntry<Telt>> span_double_cosets(DoubleCosetSplitEntry<Telt,Tidx_label> const& dcse, KernelDccEntry<Telt> const& de, bool const& compute_stabs, Telt const& id) {
#ifdef DEBUG_SPAN_DOUBLE_COSETS
  std::cerr << "ACC: ---------------- span_double_cosets |dcse.grp|=" << Order<Telt,Tidx_label,Tint>(dcse.grp) << " ----------------\n";
#endif
  auto f_can=[&](Telt const& u) -> Telt {
    return MinimalElementCosetStabChain(dcse.grp, u);
  };
  std::vector<Telt> stab_gens_std;
  Telt cos_inv = Inverse(de.cos);
  for (auto &eGen : de.stab_gens) {
    Telt gen_std = de.cos * eGen * cos_inv;
    stab_gens_std.push_back(gen_std);
  }
  std::vector<std::vector<size_t>> list_perm;
#ifdef DEBUG_SPAN_DOUBLE_COSETS
  size_t iGen = 0;
#endif
  for (auto &eGen : stab_gens_std) {
    std::vector<size_t> perm;
    for (auto & eCos : dcse.l_cos) {
      Telt prod = eCos * eGen;
      Telt prod_can = f_can(prod);
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
      if (dcse.map.count(prod_can) == 0) {
        std::cerr << "ACC: span_double_cosets, missing entry in list_perm construction\n";
        throw PermutalibException{1};
      }
#endif
      size_t pos = dcse.map.at(prod_can);
      perm.push_back(pos);
    }
#ifdef DEBUG_SPAN_DOUBLE_COSETS
    std::cerr << "ACC: iGen=" << iGen << "/" << stab_gens_std.size() << " perm=[";
    for (size_t i=0; i<perm.size(); i++) {
      if (i>0)
        std::cerr << ",";
      std::cerr << perm[i];
    }
    std::cerr << "]\n";
    iGen += 1;
#endif
    list_perm.emplace_back(std::move(perm));
  }
  size_t n_gen = list_perm.size();
#ifdef DEBUG_SPAN_DOUBLE_COSETS
  std::cerr << "ACC: span_double_cosets |list_perm|=" << list_perm.size() << "\n";
#endif
  size_t n_cos = dcse.l_cos.size();
  std::vector<KernelDccEntry<Telt>> dcc_entries;
  if (!compute_stabs) {
    // No need to compute the stabilizers here.
    Face f_done(n_cos);
#ifdef DEBUG_SPAN_DOUBLE_COSETS
    std::cerr << "ACC: compute_stabs=false n_cos=" << n_cos << "\n";
#endif
    for (size_t i=0; i<n_cos; i++) {
#ifdef DEBUG_SPAN_DOUBLE_COSETS
      std::cerr << "ACC: i=" << i << "/" << n_cos << "\n";
#endif
      if (f_done[i] == 0) {
        Telt new_cos = dcse.l_cos[i] * de.cos;
        Telt new_cos_can = f_can(new_cos);
        KernelDccEntry<Telt> new_de{new_cos_can,{}};
        dcc_entries.push_back(new_de);
        std::vector<size_t> l_idx;
        auto f_insert=[&](size_t const& pos) -> void {
          l_idx.push_back(pos);
          f_done[pos] = 1;
        };
        size_t start = 0;
        f_insert(i);
        while(true) {
          size_t len = l_idx.size();
#ifdef DEBUG_SPAN_DOUBLE_COSETS
          std::cerr << "ACC: compute_stabs=false start=" << start << " len=" << len << "\n";
#endif
          for (auto & perm : list_perm) {
            for (size_t u=start; u<len; u++) {
              size_t img = perm[l_idx[u]];
              if (f_done[img] == 0) {
                f_insert(img);
              }
            }
          }
          start = len;
          if (start == l_idx.size()) {
            break;
          }
        }
      }
    }
  } else {
    // We go to the next step, so we need the stabilizers
    // Not sure what to do for in the normal case.
    Face f_done(n_cos);
#ifdef DEBUG_SPAN_DOUBLE_COSETS
    std::cerr << "ACC: n_cos=" << n_cos << "\n";
#endif
    for (size_t i=0; i<n_cos; i++) {
#ifdef DEBUG_SPAN_DOUBLE_COSETS
      std::cerr << "ACC: i=" << i << "/" << n_cos << "\n";
#endif
      if (f_done[i] == 0) {
        Telt new_cos = dcse.l_cos[i] * de.cos;
        Telt new_cos_can = f_can(new_cos);
        std::unordered_set<Telt> set_gens;
        auto f_insert_gen=[&](Telt const& eGen) -> void {
          if (!eGen.isIdentity()) {
#ifdef DEBUG_SPAN_DOUBLE_COSETS
            Telt imgElt = new_cos_can * eGen;
            Telt imgElt_can = f_can(imgElt);
            if (imgElt_can != new_cos_can) {
              std::cerr << "ACC: The element new_cos_can is not preserved\n";
              throw PermutalibException{1};
            }
#endif
            set_gens.insert(eGen);
          }
        };
        std::vector<std::pair<size_t,Telt>> l_idx;
        std::unordered_map<size_t, size_t> map;
        auto f_insert=[&](std::pair<size_t,Telt> const& ent) -> void {
          size_t len = l_idx.size();
          l_idx.push_back(ent);
          f_done[ent.first] = 1;
          map[ent.first] = len;
        };
        size_t start = 0;
        f_insert(std::pair<size_t,Telt>{i, id});
        while(true) {
          size_t len = l_idx.size();
#ifdef DEBUG_SPAN_DOUBLE_COSETS
          std::cerr << "ACC: compute_stabs=true start=" << start << " len=" << len << "\n";
#endif
          for (size_t i_gen=0; i_gen<n_gen; i_gen++) {
            std::vector<size_t> const& perm = list_perm[i_gen];
            Telt const& eGen = de.stab_gens[i_gen];
            for (size_t u=start; u<len; u++) {
              size_t img = perm[l_idx[u].first];
              Telt imgGen = l_idx[u].second * eGen;
              if (f_done[img] == 0) {
                f_insert(std::pair<size_t,Telt>{img, imgGen});
              } else {
                size_t pos = map[img];
                Telt newStabElt = l_idx[pos].second * Inverse(imgGen);
                f_insert_gen(newStabElt);
              }
            }
          }
          start = len;
          if (start == l_idx.size()) {
            break;
          }
        }
        std::vector<Telt> vect_gens(set_gens.begin(), set_gens.end());
        auto get_reduced_vect_gens=[&]() -> std::vector<Telt> {
          if (vect_gens.size() > 2) {
            StabChainOptions<Tint, Telt> options = GetStandardOptions<Tint, Telt>(id);
            StabChain<Telt,Tidx_label> g = StabChainOp_listgen<Telt, Tidx_label, Tint>(vect_gens, options);
#ifdef DEBUG_SPAN_DOUBLE_COSETS
            StabChain<Telt,Tidx_label> g_de_stabgens = StabChainOp_listgen<Telt, Tidx_label, Tint>(de.stab_gens, options);
            Tint ord_g_de_sg = Order<Telt,Tidx_label,Tint>(g_de_stabgens);
            Tint ord_g = Order<Telt,Tidx_label,Tint>(g);
            Tint ord_l_cos = l_idx.size();
            std::cerr << "ACC: span_double_cosets |vect_gens|=" << vect_gens.size() << " |g|=" << ord_g << " |l_cos|=" << ord_l_cos << " |de.stab_gens|=" << ord_g_de_sg << "\n";
            if (ord_g * ord_l_cos != ord_g_de_sg) {
              std::cerr << "ACC: incoherence of order : ord_g_de_sg=" << ord_g_de_sg << "\n";
              throw PermutalibException{1};
            }
#endif
            return Kernel_SmallGeneratingSet<Telt,Tidx_label,Tint>(g);
          } else {
            return vect_gens;
          }
        };
        std::vector<Telt> v_gens = get_reduced_vect_gens();
#ifdef DEBUG_SPAN_DOUBLE_COSETS
        std::cerr << "ACC: |v_gens|=" << v_gens.size() << "\n";
#endif
        KernelDccEntry<Telt> new_de{std::move(new_cos_can), std::move(v_gens)};
        dcc_entries.emplace_back(std::move(new_de));
      }
    }
  }
#ifdef DEBUG_SPAN_DOUBLE_COSETS
  std::cerr << "ACC: Returning |dcc_entries|=" << dcc_entries.size() << "\n";
#endif
  return dcc_entries;
}




template<typename Telt, typename Tidx_label, typename Tint>
struct InnerDoubleCosetComputer {
private:
  std::vector<DoubleCosetSplitEntry<Telt,Tidx_label>> levels;
  size_t n_level;
  Telt id;
public:
  InnerDoubleCosetComputer(StabChain<Telt,Tidx_label> const& G, StabChain<Telt,Tidx_label> const& U) {
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
    std::cerr << "ACC: InnerDoubleCosetComputer, constructor\n";
#endif
    std::vector<StabChain<Telt,Tidx_label>> chain = Kernel_AscendingChainPair<Telt,Tidx_label,Tint>(U, G);
    n_level = chain.size() - 1;
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
    std::cerr << "ACC: InnerDoubleCosetComputer, n_level=" << n_level << "\n";
#endif
    id = U->comm->identity;
    for (size_t i_level=0; i_level<n_level; i_level++) {
      std::vector<Telt> l_cos = Kernel_RightTransversal_Direct<Telt,Tidx_label,Tint>(chain[i_level + 1], chain[i_level]);
      std::unordered_map<Telt, size_t> map;
      for (size_t u=0; u<l_cos.size(); u++) {
        map[l_cos[u]] = u;
      }
      bool is_normal = Kernel_IsNormalSubgroup(chain[i_level + 1], chain[i_level]);
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
      std::cerr << "ACC: i_level=" << i_level
                << " ord1=" << Order<Telt,Tidx_label,Tint>(chain[i_level])
                << " ord2=" << Order<Telt,Tidx_label,Tint>(chain[i_level + 1])
                << " |l_cos|=" << l_cos.size() << " is_normal=" << is_normal << "\n";
#endif
      DoubleCosetSplitEntry<Telt,Tidx_label> level{chain[i_level], l_cos, map, is_normal};
      levels.push_back(level);
    }
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
    std::cerr << "ACC: ---------------------------------------------------------\n";
#endif
  }
  std::vector<KernelDccEntry<Telt>> double_cosets_kernel(StabChain<Telt,Tidx_label> const& V, bool const& do_last) const {
    std::vector<Telt> small_gens = Kernel_SmallGeneratingSet<Telt,Tidx_label,Tint>(V);
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
    std::cerr << "ACC: double_cosets |V|=" << Order<Telt,Tidx_label,Tint>(V) << " |small_gens|=" << small_gens.size() << "\n";
#endif
    KernelDccEntry<Telt> de{id, small_gens};
    std::vector<KernelDccEntry<Telt>> l_de{de};
    for (size_t i_level=0; i_level<n_level; i_level++) {
      size_t j_level = n_level - 1 - i_level;
      DoubleCosetSplitEntry<Telt,Tidx_label> const& dcse = levels[j_level];
      bool compute_stabs = true;
      if (j_level == 0 && !do_last) {
        compute_stabs = false;
      }
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
      std::cerr << "ACC: i_level=" << i_level << " compute_stabs=" << compute_stabs << " |l_de|=" << l_de.size() << "\n";
#endif
      std::vector<KernelDccEntry<Telt>> new_l_de;
      for (auto & de: l_de) {
        std::vector<KernelDccEntry<Telt>> elist = span_double_cosets<Telt,Tidx_label,Tint>(dcse, de, compute_stabs, id);
        new_l_de.insert(new_l_de.end(), elist.begin(), elist.end());
      }
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
      std::cerr << "ACC: |new_l_de|=" << new_l_de.size() << "\n";
#endif
      l_de = new_l_de;
    }
    return l_de;
  }
  std::vector<KernelDccEntry<Telt>> double_cosets_and_stabilizers(StabChain<Telt,Tidx_label> const& V) const {
    return double_cosets_kernel(V, true);
  }
  std::vector<Telt> double_cosets(StabChain<Telt,Tidx_label> const& V) const {
    std::vector<KernelDccEntry<Telt>> l_de = double_cosets_kernel(V, false);
    std::vector<Telt> l_cos;
    for (auto & de: l_de) {
      l_cos.push_back(de.cos);
    }
    return l_cos;
  }
};

  /*
    Checking the list of double cosets.
    * We check that the pairwise intersection is empty.
    * We check that the sum of the sizes of the double cosets is equal to the size of the full group.
    The check is very expensive and with large groups will run forever.
    But it is also very elementary so bugs revealed there will be in the code, not in the checking
    code.
   */
template<typename Telt, typename Tidx_label>
void ExhaustiveCheck_DoubleCosets(StabChain<Telt,Tidx_label> const& G, StabChain<Telt,Tidx_label> const& U, StabChain<Telt,Tidx_label> const& V, std::vector<Telt> const& list_dcc) {
  std::vector<Telt> l_elt_g = get_all_elements(G);
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
  std::cerr << "ACC: We have l_elt_g\n";
#endif
  std::vector<Telt> l_elt_u = get_all_elements(U);
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
  std::cerr << "ACC: We have l_elt_u\n";
#endif
  std::vector<Telt> l_elt_v = get_all_elements(V);
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
  std::cerr << "ACC: We have l_elt_v\n";
#endif
  std::vector<std::unordered_set<Telt>> listfull_dcc;
  size_t n_elt_dcc = 0;
  for (auto & dcc: list_dcc) {
    std::unordered_set<Telt> set;
    for (auto & e_u: l_elt_u) {
      for (auto & e_v: l_elt_v) {
        Telt full_elt = e_u * dcc * e_v;
        set.insert(full_elt);
      }
    }
    n_elt_dcc += set.size();
    listfull_dcc.push_back(set);
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
    std::cerr << "ACC: |listfull_dcc|=" << listfull_dcc.size() << " set=" << set.size() << "\n";
#endif
  }
  for (size_t i_dcc=0; i_dcc<list_dcc.size(); i_dcc++) {
    for (size_t j_dcc=i_dcc+1; j_dcc<list_dcc.size(); j_dcc++) {
#ifdef DEBUG_ASCENDING_CHAINS_COSETS
      std::cerr << "ACC: i_dcc=" << i_dcc << " j_dcc=" << j_dcc << "\n";
#endif
      size_t the_int = 0;
      for (auto & val : listfull_dcc[i_dcc]) {
        if (listfull_dcc[j_dcc].count(val)) {
          the_int += 1;
        }
      }
      if (the_int > 0) {
        std::cerr << "ACC: Non-trivial insersection between i_dcc=" << i_dcc << " and j_dcc=" << j_dcc << " the_int=" << the_int << "\n";
        throw PermutalibException{1};
      }
    }
  }
  if (n_elt_dcc != l_elt_g.size()) {
    std::cerr << "ACC: n_elt_dcc=" << n_elt_dcc << " |l_elt_g|=" << l_elt_g.size() << "\n";
    std::cerr << "ACC: The double cosets do not cover the full group\n";
    throw PermutalibException{1};
  }
}

  /*
    Check that the sum of the sizes is equal to the size of the group.
    It uses the computation of intersection. So, it is exposed to
    possible errors in the computation of group Intersection.
   */
template<typename Telt, typename Tidx_label, typename Tint>
void FastCheckSizes_DoubleCosets(StabChain<Telt,Tidx_label> const& G, StabChain<Telt,Tidx_label> const& U, StabChain<Telt,Tidx_label> const& V, std::vector<Telt> const& list_dcc) {
  Telt id = V->comm->identity;
  StabChainOptions<Tint, Telt> options = GetStandardOptions<Tint, Telt>(id);
  std::vector<Telt> LGen_V = Kernel_SmallGeneratingSet<Telt,Tidx_label,Tint>(V);
  Tint sum_sizes = 0;
  Tint size_G = Order<Telt, Tidx_label, Tint>(G);
  Tint size_U = Order<Telt, Tidx_label, Tint>(U);
  Tint size_V = Order<Telt, Tidx_label, Tint>(V);
  for (auto & dcc : list_dcc) {
    Telt dcc_inv = Inverse(dcc);
    std::vector<Telt> NewLGen;
    for (auto & eGen: LGen_V) {
      Telt NewGen = dcc * eGen * dcc_inv;
      NewLGen.push_back(NewGen);
    }
    StabChain<Telt,Tidx_label> ConjV = StabChainOp_listgen<Telt, Tidx_label, Tint>(NewLGen, options);
    StabChain<Telt,Tidx_label> eInt = Kernel_Intersection<Telt, Tidx_label, Tint>(U, ConjV);
    Tint size_Int = Order<Telt, Tidx_label, Tint>(eInt);
    Tint n_cos = size_V / size_Int;
    Tint double_cos_size = size_U * n_cos;
    sum_sizes += double_cos_size;
  }
  if (size_G != sum_sizes) {
    std::cerr << "ACC: size_G=" << size_G << "\n";
    std::cerr << "ACC: sum_sizes=" << sum_sizes << "\n";
    std::cerr << "ACC: The double cosets do not cover the full group\n";
    throw PermutalibException{1};
  }
}

  /*
    Check that the pairwise intersection of the double cosets
    is empty.
   */
template<typename Telt, typename Tidx_label, typename Tint>
void FastCheckIntersection_DoubleCosets(StabChain<Telt,Tidx_label> const& G, StabChain<Telt,Tidx_label> const& U, StabChain<Telt,Tidx_label> const& V, std::vector<Telt> const& list_dcc) {
  Telt id = V->comm->identity;
  StabChainOptions<Tint, Telt> options = GetStandardOptions<Tint, Telt>(id);
  std::vector<Telt> LGen_V = Kernel_SmallGeneratingSet<Telt,Tidx_label,Tint>(V);
  Tint size_G = Order<Telt, Tidx_label, Tint>(G);
  Tint size_U = Order<Telt, Tidx_label, Tint>(U);
  Tint size_V = Order<Telt, Tidx_label, Tint>(V);
  std::vector<Telt> l_elt_u = get_all_elements(U);
  std::vector<std::unordered_set<Telt>> l_elts_dcc;
  for (auto & dcc: list_dcc) {
    Telt dcc_inv = Inverse(dcc);
    std::vector<Telt> NewLGen;
    for (auto & eGen: LGen_V) {
      Telt NewGen = dcc * eGen * dcc_inv;
      NewLGen.push_back(NewGen);
    }
    StabChain<Telt,Tidx_label> ConjV = StabChainOp_listgen<Telt, Tidx_label, Tint>(NewLGen, options);
    StabChain<Telt,Tidx_label> eInt = Kernel_Intersection<Telt, Tidx_label, Tint>(U, ConjV);
    //
    std::vector<Telt> l_cos = enumerate_right_cosets<Telt,Tidx_label,Tint>(eInt, ConjV);
    std::unordered_set<Telt> set;
    for (auto & eU : l_elt_u) {
      for (auto & eCos : l_cos) {
        Telt prod = eU * eCos;
        set.insert(prod);
      }
    }
    l_elts_dcc.push_back(set);
  }
  size_t n_dcc = list_dcc.size();
  for (size_t i_dcc=0; i_dcc<n_dcc; i_dcc++) {
    for (size_t j_dcc=i_dcc+1; j_dcc<n_dcc; j_dcc++) {
      for (auto & eX : l_elts_dcc[i_dcc]) {
        if (l_elts_dcc[j_dcc].count(eX) == 1) {
          std::cerr << "ACC: The intersection is not empty i_dcc=" << i_dcc << " j_dcc=" << j_dcc << "\n";
          throw PermutalibException{1};
        }
      }
    }
  }
}


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
                # NC is safe (initializing as TrivialSubgroup(G)
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
          Add(nstab, st);

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
