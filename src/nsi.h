// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GAP_NSI_H_
#define SRC_GAP_NSI_H_

#include "NumberTheory.h"
#include "StabChain.h"
#include "stbcbckt.h"
#include <algorithm>
#include <limits>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

/*
#
# V 1.1 - Bug fixed
#
# By Steve Linton
#
# Bad documentation by Chris Jefferson
#
# Finds the minimal image of a Set set under a group G.
#
# Usage: NewSmallestImage(G, set, stab, x -> x);
#
# (ignore the last argument!)
#
# Where stab should be a subgroup of
# Stabilizer(G,set);
#
# If in doubt, the best way to invoke this algorithm
# is:
# NewSmallestImage(G, set, Stabilizer(G, set, OnSets), x -> x);
#
# Returns a pair [image, stabilizer], where stabilizer is a subgroup of
Stabilizer(G, set), possibly larger than the one given into the function.
#
# Note: The return type of this is NOT a set, but provides the pointwise mapping
of the input set. # This means the permutation which actually provides the
smallest image can be found cheaply as follows:
#
# res := NewSmallestImage(G, set, Group(()), x -> x);
# perm := RepresentativeAction(G, set, res[1], OnTuples);



#
# Search node data:
#
#  selected -- indices in set of points being mapped to minimal image at this
node #  image    -- sequence-wise image of set under element represented by this
node #  substab --  Stab_K(selected) sequence stabilizer #  children --  nodes
corresponding to extensions of selected #  parent -- node corr to all but last
element of selected.] #  childno #  next -- across row #  prev -- likewise #
deleted
#

# At each level

# Find the next pt of minimum image and all the corresponding nodes

# If at any time in this process we notice a potential prune, we have a
# generator of K -- construct it, close K with it and delete all tree
# branches that are now non-canonical -- propagate new K down through
# tree. Continue with surviving nodes at current level.

*/

namespace permutalib {

template <typename T> void Remove(std::vector<T> &eV, int const &pos) {
  eV.erase(eV.begin() + pos);
}

template <typename Telt> struct ResultCanonicalization {
  std::vector<int> set;
  Telt g;
};

template <typename Telt, typename T>
Telt PermListList(std::vector<T> const &list1, std::vector<T> const &list2) {
  using Tidx = typename Telt::Tidx;
  std::unordered_map<T, int> eMap;
  size_t siz = list2.size();
  for (size_t i = 0; i < siz; i++) {
    eMap[list1[i]] = i;
  }
  //
  std::vector<Tidx> eList(siz);
  for (size_t i = 0; i < siz; i++)
    eList[i] = eMap[list2[i]];
  return Telt(std::move(eList));
}

template <typename Telt, typename Tidx_label, typename Tint>
StabChain<Telt, Tidx_label>
Action(StabChain<Telt, Tidx_label> const &S,
       std::vector<typename Telt::Tidx> const &set) {
  using Tidx = typename Telt::Tidx;
  Tidx n = S->comm->n;
  Tidx siz = Tidx(set.size());
  std::vector<Tidx> map_idx(n, 0);
  for (Tidx i = 0; i < siz; i++)
    map_idx[set[i]] = i;
  //
  std::vector<Telt> LGen;
  for (auto &eGen : Kernel_GeneratorsOfGroup(S)) {
    std::vector<Tidx> eList(siz);
    for (Tidx i = 0; i < siz; i++) {
      Tidx ePt = set[i];
      Tidx fPt = OnPoints(ePt, eGen);
      Tidx j = map_idx[fPt];
      eList[i] = j;
    }
    Telt ePerm = Telt(std::move(eList));
    LGen.emplace_back(std::move(ePerm));
  }
  return FCT_Group<Telt, Tidx_label, Tint>(LGen, siz);
}

template <typename Telt, typename Tidx>
std::vector<Tidx> OnTuples(std::vector<Tidx> const &V, Telt const &g) {
  Tidx len = V.size();
  std::vector<Tidx> retV(len);
  for (Tidx i = 0; i < len; i++)
    retV[i] = PowAct(V[i], g);
  return retV;
}

template <typename Telt, typename Tidx>
void OnTuples_inplace(std::vector<Tidx> &V, Telt const &g) {
  size_t len = V.size();
  for (size_t i = 0; i < len; i++)
    V[i] = PowAct(V[i], g);
}

template <typename T> std::vector<T> Set(std::vector<T> const &V) {
  std::vector<T> retV = V;
  std::sort(retV.begin(), retV.end());
  return retV;
}

/*
  Modification done:
  --- skip_fnuc eliminated as it is the identity in the case that interest us.
 */
template <typename Telt, typename Tidx_label, typename Tint>
std::pair<std::vector<typename Telt::Tidx>, StabChain<Telt, Tidx_label>>
NewCanonicImage(StabChain<Telt, Tidx_label> const &g,
                std::vector<typename Telt::Tidx> const &set,
                StabChain<Telt, Tidx_label> const &k_group) {
  using Tidx = typename Telt::Tidx;
#ifdef DEBUG_NSI
  std::cerr << "CPP NewCanonicImage : beginning\n";
#endif
  Tidx infinity = std::numeric_limits<Tidx>::max();
  Tidx initial_upb = infinity;
  Tidx max_val_type = std::numeric_limits<Tidx>::max();

  auto calculateBestOrbit = [&](std::vector<Tidx> const &orbmins,
                                std::vector<Tidx> const &orbitCounts,
                                std::vector<Tidx> const &orbsizes) -> Tidx {
#ifdef DEBUG_NSI
    std::cerr << "CPP Beginning of calculateBestOrbit\n";
#endif
    struct typeCnt {
      Tidx orbSize;
      Tidx orbCount;
    };
    // The comparison goes
    // log(a.orbCount) / a.orbSize < log(b.orbCount) / b.orbSize
    // which is equivalent to
    // a.orbCount ^ b.orbSize    <    b.orbCount ^ a.orbSize
    auto comparisonLower = [](typeCnt const &a, typeCnt const &b) -> bool {
#ifdef DEBUG_NSI
      std::cerr << "CPP compLower a=" << a.orbCount << " / " << a.orbSize
                << " b=" << b.orbCount << " / " << b.orbSize << "\n";
#endif
      if (a.orbSize == 1) {
        if (b.orbSize == 1) {
          return a.orbCount < b.orbCount;
        } else {
          return true;
        }
      } else {
        if (b.orbSize == 1)
          return false;
        Tint pow1 = ComputePower(Tint(a.orbCount), b.orbSize);
        //        std::cerr << "pow1=" << pow1 << " a.orbCount=" <<
        //        Tidx(a.orbCount) << " b.orbSize=" << size_t(b.orbSize) <<
        //        "\n";
        Tint pow2 = ComputePower(Tint(b.orbCount), a.orbSize);
        //        std::cerr << "pow2=" << pow2 << " b.orbCount=" <<
        //        Tidx(b.orbCount) << " a.orbSize=" << size_t(a.orbSize) <<
        //        "\n";
#ifdef DEBUG_NSI
        std::cerr << "CPP pow1=" << pow1 << " pow2=" << pow2 << "\n";
#endif
        return pow1 < pow2;
      }
    };
    auto comparisonEqual = [](typeCnt const &a, typeCnt const &b) -> bool {
#ifdef DEBUG_NSI
      std::cerr << "CPP compEqual a=" << a.orbCount << " / " << a.orbSize
                << " b=" << b.orbCount << " / " << b.orbSize << "\n";
#endif
      if (a.orbSize == 1) {
        if (b.orbSize == 1) {
          return a.orbCount == b.orbCount;
        } else {
          return false;
        }
      } else {
        if (b.orbSize == 1)
          return false;
        Tint pow1 = ComputePower(Tint(a.orbCount), b.orbSize);
        //        std::cerr << "pow1=" << pow1 << " a.orbCount=" <<
        //        Tidx(a.orbCount) << " b.orbSize=" << size_t(b.orbSize) <<
        //        "\n";
        Tint pow2 = ComputePower(Tint(b.orbCount), a.orbSize);
        //        std::cerr << "pow2=" << pow2 << " b.orbCount=" <<
        //        Tidx(b.orbCount) << " a.orbSize=" << size_t(a.orbSize) <<
        //        "\n";
        return pow1 == pow2;
      }
    };
    auto selector = [&](Tidx const &jdx) -> typeCnt {
      return {orbsizes[jdx], orbitCounts[jdx]};
    };
    Tidx index = 0;
    typeCnt result_0 = selector(0);
    Tidx result_1 = orbmins[0];
#ifdef DEBUG_NSI
    std::cerr << "CPP result_0=" << result_0.orbCount << " / "
              << result_0.orbSize << "   result_1=" << int(result_1 + 1)
              << "\n";
#endif
    for (Tidx i = 1; i < Tidx(orbmins.size()); i++) {
      typeCnt ret_0 = selector(i);
      Tidx ret_1 = orbmins[i];
      bool lower = false;
#ifdef DEBUG_NSI
      std::cerr << "CPP i=" << int(i + 1) << " ret_0=" << ret_0.orbCount
                << " / " << ret_0.orbSize << "   ret_1=" << int(ret_1 + 1)
                << "\n";
#endif
      if (comparisonLower(ret_0, result_0)) {
        lower = true;
      } else {
        if (comparisonEqual(ret_0, result_0)) {
          if (ret_1 < result_1)
            lower = true;
        }
      }
#ifdef DEBUG_NSI
      std::cerr << "CPP lower=" << lower << "\n";
#endif
      //
      if ((orbitCounts[index] == 0) || (lower && orbitCounts[i] != 0)) {
        index = i;
        result_0 = ret_0;
        result_1 = ret_1;
      }
    }
    return index;
  };

  struct Node {
    std::vector<Tidx> selected;
    std::vector<Tidx> image;
    StabChain<Telt, Tidx_label> substab;
    bool deleted;
    std::shared_ptr<Node> next;
    std::shared_ptr<Node> prev;
    std::shared_ptr<Node> parent;
    // children
    Tidx childno;
    bool IsBoundChildren;
    std::vector<std::shared_ptr<Node>> children;
    std::vector<Tidx> validkids;
  };
  using NodePtr = std::shared_ptr<Node>;
  std::vector<NodePtr> ListPtr;

  Tidx n = set[set.size() - 1] + 1;
  Tidx n_largest = LargestMovedPoint(StrongGeneratorsStabChain(g));
  if (n_largest > n)
    n = n_largest;
#ifdef DEBUG_NSI
  std::cerr << "CPP n=" << int(n) << "\n";
  std::cerr << "DEBUG MATCH set=" << GapStringIntVector(set) << "\n";
#endif
  StabChain<Telt, Tidx_label> s = CopyStabChain(g);
  StabChain<Telt, Tidx_label> l_group = Action<Telt, Tidx_label, Tint>(k_group, set);
  Tidx m = Tidx(set.size());
  Node root_v;
  root_v.image = set;
  root_v.substab = l_group;
  root_v.deleted = false;
  root_v.next = nullptr;
  root_v.prev = nullptr;
  root_v.parent = nullptr;
  // unset values
  //  root_v.selectedbaselength = max_val_type;
  root_v.IsBoundChildren = false;
  NodePtr root = std::make_shared<Node>(root_v);
  // no need to put root in the list of nodes to be deleted as the setting of
  // all to nullptr eventually kills it.
  //  ListPtr.push_back(root);

  // Node exploration functions
  auto leftmost_node = [&](Tidx const &depth) -> NodePtr {
#ifdef DEBUG_NSI
    std::cerr << "CPP Beginning of leftmost_node\n";
#endif
    NodePtr n = root;
    while (Tidx(n->selected.size()) < depth)
      n = n->children[0];
    return n;
  };
  auto next_node = [&](NodePtr const &node) -> NodePtr {
#ifdef DEBUG_NSI
    std::cerr << "CPP Beginning of next_node\n";
#endif
    NodePtr n = node;
    while (true) {
      n = n->next;
      if (n == nullptr || !n->deleted)
        break;
    }
    return n;
  };
  // Delete a node, and recursively deleting all it's children.
  std::function<void(NodePtr &)> delete_node = [&](NodePtr &node) -> void {
#ifdef DEBUG_NSI
    std::cerr << "CPP Beginning of delete_node\n";
#endif
    if (node->deleted) {
      return;
    }
    if (node->prev != nullptr) {
      node->prev->next = node->next;
    }
    if (node->next != nullptr) {
      node->next->prev = node->prev;
    }
    node->deleted = true;
    if (node->parent != nullptr) {
      Remove(node->parent->children, node->childno);
      if (node->parent->children.size() == 0) {
        delete_node(node->parent);
      } else {
        for (Tidx i = node->childno; i < Tidx(node->parent->children.size());
             i++) {
          node->parent->children[i]->childno = i;
        }
      }
    }
    if (node->IsBoundChildren) {
      for (auto &enode : node->children)
        delete_node(enode);
    }
  };

  // Given a group 'gp' and a set 'set', find orbit representatives
  // of 'set' in 'gp' simply.
  std::vector<Tidx> q_sor;
  q_sor.reserve(n);
  Face b_sor(n);
  auto simpleOrbitReps =
      [&](StabChain<Telt, Tidx_label> const &gp,
          std::vector<Tidx> const &set) -> std::vector<Tidx> {
#ifdef DEBUG_NSI
    std::cerr << "CPP Beginning of simpleOrbitReps\n";
#endif
    Tidx m = Tidx(set.size());
    Tidx n = set[m - 1] + 1;
    for (Tidx i = 0; i < n; i++)
      b_sor[i] = 0;
#ifdef DEBUG_NSI
    std::cerr << "DEBUG n=" << int(n) << "\n";
#endif
    for (auto &eVal : set) {
#ifdef DEBUG_NSI
      std::cerr << "DEBUG eVal=" << int(eVal) << " n=" << int(n) << "\n";
#endif
      b_sor[eVal] = 1;
    }
    boost::dynamic_bitset<>::size_type seed = b_sor.find_first();
    Tidx seed_tidx = Tidx(seed);
    std::vector<Tidx> reps;
    std::vector<Telt> gens = Kernel_GeneratorsOfGroup(gp);
    while (seed != boost::dynamic_bitset<>::npos && seed_tidx < n) {
#ifdef DEBUG_NSI
      std::cerr << "DEBUG seed=" << int(seed) << " n=" << int(n) << "\n";
#endif
      b_sor[seed] = 0;
      q_sor.clear();
      q_sor.push_back(seed_tidx);
      reps.push_back(seed_tidx);
      size_t pos = 0;
      while (true) {
        size_t idx, qsiz = q_sor.size();
        for (idx = pos; idx < qsiz; idx++) {
          Tidx pt = q_sor[idx];
          for (auto &gen : gens) {
            Tidx im = PowAct(pt, gen);
            if (b_sor[im] == 1) {
#ifdef DEBUG_NSI
              std::cerr << "DEBUG im=" << int(im) << " n=" << int(n) << "\n";
#endif
              b_sor[im] = 0;
              q_sor.push_back(im);
            }
          }
        }
        if (idx == q_sor.size())
          break;
        pos = idx;
      }
      seed = b_sor.find_next(seed);
      seed_tidx = Tidx(seed);
    }
    return reps;
  };
  auto DifferenceVect_local =
      [](Tidx const &m, std::vector<Tidx> const &LIdx) -> std::vector<Tidx> {
    Face b(m);
    for (auto &eVal : LIdx)
      b[eVal] = 1;
    size_t len_ret = size_t(m) - LIdx.size();
    size_t pos = 0;
    std::vector<Tidx> cands(len_ret);
    for (Tidx i = 0; i < m; i++)
      if (b[i] == 0) {
        cands[pos] = i;
        pos++;
      }
    return cands;
  };
  // We need to break all the cycles in order to the memory free to happen
  // correctly. We use a hack in order to get that behavior: A vector of all the
  // nodes, then set the pointer to zero and so all cycles eliminated. Maybe we
  // could do better, but the hack should be adequate.
  auto free_all_nodes = [&]() -> void {
    //    std::cerr << "|ListPtr|=" << ListPtr.size() << "\n";
    for (auto &e_node : ListPtr) {
      e_node->prev = nullptr;
      e_node->next = nullptr;
      e_node->parent = nullptr;
    }
  };
  if (set.size() == 0) {
    free_all_nodes();
    return {{}, g};
  }
  Tidx depth;
  std::vector<Tidx> orbmins;
  orbmins.reserve(n);
  std::vector<Tidx> orbsizes;
  orbsizes.reserve(n);
  std::vector<Tidx> q;
  q.reserve(n);
  std::vector<Tidx> orbnums(n, max_val_type);
  for (depth = 0; depth < m; depth++) {
#ifdef DEBUG_NSI
    std::cerr << "CPP depth=" << int(depth + 1) << "\n";
#endif
    std::vector<Telt> gens = GetListGenerators(s);
    for (Tidx i = 0; i < n; i++)
      orbnums[i] = max_val_type;
    orbmins.clear();
    orbsizes.clear();
    Tidx upb = initial_upb;
    // Make orbit of x, updating orbnums, orbmins and orbsizes as approriate.
    auto make_orbit = [&](Tidx const &x) -> Tidx {
#ifdef DEBUG_NSI
      if (orbnums[x] == max_val_type) {
        std::cerr << "CPP Beginning of make_orbit x=" << int(x + 1)
                  << " orbnums[x]=" << int(-1) << "\n";
      } else {
        std::cerr << "CPP Beginning of make_orbit x=" << int(x + 1)
                  << " orbnums[x]=" << int(orbnums[x] + 1) << "\n";
      }
#endif
      if (orbnums[x] != max_val_type) {
        return orbnums[x];
      }
#ifdef DEBUG_NSI
      std::cerr << "CPP After check\n";
#endif
      q.clear();
      q.push_back(x);
      Tidx rep = x;
      Tidx num = Tidx(orbmins.size());
      orbnums[x] = num;
#ifdef DEBUG_NSI
      std::cerr << "CPP 1 : x=" << int(x + 1)
                << " Assign orbnums[x]=" << int(orbnums[x] + 1) << "\n";
#endif
      size_t pos = 0;
      while (true) {
        size_t idx, qsiz = q.size();
        for (idx = pos; idx < qsiz; idx++) {
          Tidx pt = q[idx];
          for (auto &gen : gens) {
            Tidx img = PowAct(pt, gen);
            if (orbnums[img] == max_val_type) {
              orbnums[img] = num;
#ifdef DEBUG_NSI
              std::cerr << "CPP 2 : img=" << int(img + 1)
                        << " Assign orbnums[img]=" << int(orbnums[img] + 1)
                        << "\n";
#endif
              q.push_back(img);
              if (img < rep)
                rep = img;
            }
          }
        }
        if (idx == q.size())
          break;
        pos = idx;
      }
#ifdef DEBUG_NSI
      std::cerr << "CPP make_orbit rep=" << int(rep + 1) << " |q|=" << q.size()
                << " num=" << int(num + 1) << "\n";
#endif
      orbmins.push_back(rep);
      orbsizes.push_back(Tidx(q.size()));
      return num;
    };
    /*
      # At this point, all bottom nodes are blue
      # first pass creates appropriate set of virtual red nodes
    */

    std::vector<Tidx> minOrbitMset = {infinity};
    NodePtr node = leftmost_node(depth);
    while (node != nullptr) {
#ifdef DEBUG_NSI
      std::cerr << "CPP m=" << m
                << " node.selected=" << GapStringIntVector(node->selected)
                << "\n";
#endif
      std::vector<Tidx> cands = DifferenceVect_local(m, node->selected);
#ifdef DEBUG_NSI
      std::cerr << "CPP 1 : cands=" << GapStringIntVector(cands) << "\n";
#endif

      std::vector<Tidx> orbitMset;
      for (auto &y : cands) {
        Tidx x = node->image[y];
        Tidx num = make_orbit(x);
#ifdef DEBUG_NSI
        std::cerr << "CPP x=" << int(x + 1) << " num=" << int(num + 1) << "\n";
#endif
        orbitMset.push_back(orbmins[num]);
      }
#ifdef DEBUG_NSI
      std::cerr << "CPP bef orbitMset=" << GapStringIntVector(orbitMset)
                << "\n";
#endif
      std::sort(orbitMset.begin(), orbitMset.end());
#ifdef DEBUG_NSI
      std::cerr << "CPP aft orbitMset=" << GapStringIntVector(orbitMset)
                << "\n";
#endif
      if (orbitMset < minOrbitMset) {
#ifdef DEBUG_NSI
        std::cerr << "CPP orbitMset comparison case 1\n";
#endif
        minOrbitMset = orbitMset;
        NodePtr node2 = node->prev;
        while (node2 != nullptr) {
          delete_node(node2);
          node2 = node2->prev;
        }
      } else {
        if (orbitMset > minOrbitMset) {
#ifdef DEBUG_NSI
          std::cerr << "CPP orbitMset comparison case 2\n";
#endif
          delete_node(node);
        }
      }
      node = next_node(node);
    }

    std::vector<Tidx> globalOrbitCounts(orbmins.size(), 0);
    node = leftmost_node(depth);
    while (node != nullptr) {
      std::vector<Tidx> cands = DifferenceVect_local(m, node->selected);
#ifdef DEBUG_NSI
      std::cerr << "CPP 2 : cands=" << GapStringIntVector(cands) << "\n";
#endif
      if (cands.size() > 1 && !IsTrivial(node->substab)) {
        cands = simpleOrbitReps(node->substab, cands);
      }
      /*
        # These index the children of node that will
        # not be immediately deleted under rule C
      */
      for (auto &y : cands) {
        Tidx x = node->image[y];
#ifdef DEBUG_NSI
        std::cerr << "CPP y=" << int(y + 1) << " x=" << int(x + 1) << "\n";
#endif
        Tidx num = make_orbit(x);
        Tidx siz = Tidx(globalOrbitCounts.size());
        if (num < siz) {
          globalOrbitCounts[num]++;
        } else {
          for (Tidx u = siz; u <= num; u++)
            globalOrbitCounts.push_back(0);
          globalOrbitCounts[num] = 1;
        }
#ifdef DEBUG_NSI
        std::cerr << "CPP globalOrbitCounts : num=" << int(num + 1)
                  << " cnt=" << globalOrbitCounts[num] << "\n";
#endif
      }
      node = next_node(node);
    }
#ifdef DEBUG_NSI
    std::cerr << "CPP globalOrbitCounts=" << GapStringTVector(globalOrbitCounts)
              << "\n";
#endif
    Tidx globalBestOrbit =
        calculateBestOrbit(orbmins, globalOrbitCounts, orbsizes);
    upb = orbmins[globalBestOrbit];
#ifdef DEBUG_NSI
    std::cerr << "CPP globalBestOrbit=" << int(globalBestOrbit + 1)
              << " upb=" << int(upb + 1) << "\n";
#endif

    node = leftmost_node(depth);
    while (node != nullptr) {
      std::vector<Tidx> cands = DifferenceVect_local(m, node->selected);
#ifdef DEBUG_NSI
      std::cerr << "CPP 3 : cands=" << GapStringIntVector(cands) << "\n";
#endif
      if (cands.size() > 1 && !IsTrivial(node->substab)) {
        cands = simpleOrbitReps(node->substab, cands);
      }
      /*
        # These index the children of node that will
        # not be immediately deleted under rule C
      */
      node->validkids.clear();
      for (auto &y : cands) {
        Tidx x = node->image[y];
        Tidx num = orbnums[x];
        if (num == max_val_type) {
          /*
            # Need a new orbit. Also require the smallest point
            # as the rep.
            #
            #
            # If there is no prospect of the new orbit being
            # better than the current best then go on to the next candidate
          */
          num = make_orbit(x);
          Tidx rep = orbmins[num];
          if (rep < upb) {
            upb = rep;
            NodePtr node2 = node->prev;
            while (node2 != nullptr) {
              delete_node(node2);
              node2 = node2->prev;
            }
            node->validkids.clear();
            node->validkids.push_back(y);
#ifdef DEBUG_NSI
            std::cerr << "CPP validkids set to {y} with y=" << int(y + 1)
                      << "\n";
#endif
          }
        } else {
          Tidx rep = orbmins[num];
#ifdef DEBUG_NSI
          std::cerr << "CPP before insertion rep=" << int(rep + 1)
                    << " num=" << int(num + 1) << " upb=" << int(upb + 1)
                    << "\n";
#endif
          if (rep == upb) {
            node->validkids.push_back(y);
#ifdef DEBUG_NSI
            std::cerr << "CPP validkids inserting y=" << int(y + 1) << "\n";
#endif
          }
        }
      }
      if (node->validkids.size() == 0) {
        delete_node(node);
      }
      node = next_node(node);
    }
    /*
      # Second pass. Actually make all the red nodes and turn them blue
    */
#ifdef DEBUG_NSI
    std::cerr << "CPP Before ChangeStabChain\n";
#endif
    ChangeStabChain(s, {upb}, false);
#ifdef DEBUG_NSI
    std::cerr << "CPP After ChangeStabChain\n";
#endif
    bool do_continue = false;
    if (s->orbit.size() == 1) {
      /*
        # In this case nothing much can happen. Each surviving node will have
        exactly one child # and none of the imsets will change # so we mutate
        the nodes in-place
      */
      node = leftmost_node(depth);
      while (node != nullptr) {
        //        if (node->selectedbaselength == max_val_type)
        //          node->selectedbaselength = Tidx(node->selected.size());
        node->selected.push_back(node->validkids[0]);
#ifdef DEBUG_NSI
        std::cerr << "CPP Now node.selected="
                  << GapStringIntVector(node->selected) << "\n";
#endif
        node = next_node(node);
      }
      s = s->stabilizer;
      if (Tidx(leftmost_node(depth + 1)->selected.size()) != m) {
        do_continue = true;
      }
    }

    if (!do_continue) {
      node = leftmost_node(depth);
      NodePtr prevnode = nullptr;
      while (node != nullptr) {
        node->IsBoundChildren = true;
        node->children.clear();
#ifdef DEBUG_NSI
        std::cerr << "CPP node.validkids="
                  << GapStringIntVector(node->validkids) << "\n";
#endif
        for (auto &x : node->validkids) {
          Node newnode_v;
          newnode_v.selected = node->selected;
          newnode_v.selected.push_back(x);
          //          PrintStabChain(node->substab);
#ifdef DEBUG_NSI
          std::cerr << "DEBUG Before Stabilize_OnPoints x=" << int(x + 1)
                    << "\n";
#endif
          newnode_v.substab =
              Kernel_Stabilizer_OnPoints<Telt, Tidx_label, Tint>(node->substab,
                                                                 x);
#ifdef DEBUG_NSI
          std::cerr << "DEBUG After Stabilize_OnPoints\n";
#endif
          newnode_v.parent = node;
          newnode_v.childno = Tidx(node->children.size());
          newnode_v.next = nullptr;
          newnode_v.prev = prevnode;
          newnode_v.deleted = false;
          newnode_v.IsBoundChildren = false;
#ifdef DEBUG_NSI
          std::cerr << "CPP newnode.selected="
                    << GapStringIntVector(newnode_v.selected) << "\n";
#endif
          NodePtr newnode = std::make_shared<Node>(newnode_v);
          ListPtr.push_back(newnode);
          if (prevnode != nullptr) {
            prevnode->next = newnode;
          }
          prevnode = newnode;
          node->children.push_back(newnode);

          std::vector<Tidx> image = node->image;
          if (image[x] != upb) {
            while (true) {
              const Telt &g = s->comm->labels[s->transversal[image[x]]];
              OnTuples_inplace(image, g);
              if (image[x] == upb)
                break;
            }
            newnode->image = image;
          } else {
            newnode->image = image;
          }
        }
        node = next_node(node);
      }

#ifdef DEBUG_NSI
      std::cerr << "CPP Before s:=s.stabilizer operation\n";
#endif
      s = s->stabilizer;
      if (Tidx(leftmost_node(depth + 1)->selected.size()) == m) {
        break;
      }
    }
  }
  free_all_nodes();
  NodePtr node = leftmost_node(depth + 1);
  return {node->image, node->substab};
}

template <typename Telt, typename Tidx_label, typename Tint>
std::optional<std::pair<std::vector<typename Telt::Tidx>,size_t>>
NewCanonicImageInitialTriv(StabChain<Telt, Tidx_label> const &g,
                           std::vector<typename Telt::Tidx> const &set,
                           size_t const& max_size) {
  using Tidx = typename Telt::Tidx;
#ifdef DEBUG_NSI
  std::cerr << "CPP NewCanonicImage : beginning\n";
#endif
  Tidx infinity = std::numeric_limits<Tidx>::max();
  Tidx initial_upb = infinity;
  Tidx max_val_type = std::numeric_limits<Tidx>::max();

  auto calculateBestOrbit = [&](std::vector<Tidx> const &orbmins,
                                std::vector<Tidx> const &orbitCounts,
                                std::vector<Tidx> const &orbsizes) -> Tidx {
#ifdef DEBUG_NSI
    std::cerr << "CPP Beginning of calculateBestOrbit\n";
#endif
    struct typeCnt {
      Tidx orbSize;
      Tidx orbCount;
    };
    // The comparison goes
    // log(a.orbCount) / a.orbSize < log(b.orbCount) / b.orbSize
    // which is equivalent to
    // a.orbCount ^ b.orbSize    <    b.orbCount ^ a.orbSize
    auto comparisonLower = [](typeCnt const &a, typeCnt const &b) -> bool {
#ifdef DEBUG_NSI
      std::cerr << "CPP compLower a=" << a.orbCount << " / " << a.orbSize
                << " b=" << b.orbCount << " / " << b.orbSize << "\n";
#endif
      if (a.orbSize == 1) {
        if (b.orbSize == 1) {
          return a.orbCount < b.orbCount;
        } else {
          return true;
        }
      } else {
        if (b.orbSize == 1)
          return false;
        Tint pow1 = ComputePower(Tint(a.orbCount), b.orbSize);
        //        std::cerr << "pow1=" << pow1 << " a.orbCount=" <<
        //        Tidx(a.orbCount) << " b.orbSize=" << size_t(b.orbSize) <<
        //        "\n";
        Tint pow2 = ComputePower(Tint(b.orbCount), a.orbSize);
        //        std::cerr << "pow2=" << pow2 << " b.orbCount=" <<
        //        Tidx(b.orbCount) << " a.orbSize=" << size_t(a.orbSize) <<
        //        "\n";
#ifdef DEBUG_NSI
        std::cerr << "CPP pow1=" << pow1 << " pow2=" << pow2 << "\n";
#endif
        return pow1 < pow2;
      }
    };
    auto comparisonEqual = [](typeCnt const &a, typeCnt const &b) -> bool {
#ifdef DEBUG_NSI
      std::cerr << "CPP compEqual a=" << a.orbCount << " / " << a.orbSize
                << " b=" << b.orbCount << " / " << b.orbSize << "\n";
#endif
      if (a.orbSize == 1) {
        if (b.orbSize == 1) {
          return a.orbCount == b.orbCount;
        } else {
          return false;
        }
      } else {
        if (b.orbSize == 1)
          return false;
        Tint pow1 = ComputePower(Tint(a.orbCount), b.orbSize);
        //        std::cerr << "pow1=" << pow1 << " a.orbCount=" <<
        //        Tidx(a.orbCount) << " b.orbSize=" << size_t(b.orbSize) <<
        //        "\n";
        Tint pow2 = ComputePower(Tint(b.orbCount), a.orbSize);
        //        std::cerr << "pow2=" << pow2 << " b.orbCount=" <<
        //        Tidx(b.orbCount) << " a.orbSize=" << size_t(a.orbSize) <<
        //        "\n";
        return pow1 == pow2;
      }
    };
    auto selector = [&](Tidx const &jdx) -> typeCnt {
      return {orbsizes[jdx], orbitCounts[jdx]};
    };
    Tidx index = 0;
    typeCnt result_0 = selector(0);
    Tidx result_1 = orbmins[0];
#ifdef DEBUG_NSI
    std::cerr << "CPP result_0=" << result_0.orbCount << " / "
              << result_0.orbSize << "   result_1=" << int(result_1 + 1)
              << "\n";
#endif
    for (Tidx i = 1; i < Tidx(orbmins.size()); i++) {
      typeCnt ret_0 = selector(i);
      Tidx ret_1 = orbmins[i];
      bool lower = false;
#ifdef DEBUG_NSI
      std::cerr << "CPP i=" << int(i + 1) << " ret_0=" << ret_0.orbCount
                << " / " << ret_0.orbSize << "   ret_1=" << int(ret_1 + 1)
                << "\n";
#endif
      if (comparisonLower(ret_0, result_0)) {
        lower = true;
      } else {
        if (comparisonEqual(ret_0, result_0)) {
          if (ret_1 < result_1)
            lower = true;
        }
      }
#ifdef DEBUG_NSI
      std::cerr << "CPP lower=" << lower << "\n";
#endif
      //
      if ((orbitCounts[index] == 0) || (lower && orbitCounts[i] != 0)) {
        index = i;
        result_0 = ret_0;
        result_1 = ret_1;
      }
    }
    return index;
  };

  struct Node {
    std::vector<Tidx> selected;
    std::vector<Tidx> image;
    bool deleted;
    std::shared_ptr<Node> next;
    std::shared_ptr<Node> prev;
    std::shared_ptr<Node> parent;
    // children
    Tidx childno;
    bool IsBoundChildren;
    std::vector<std::shared_ptr<Node>> children;
    std::vector<Tidx> validkids;
  };
  using NodePtr = std::shared_ptr<Node>;
  std::vector<NodePtr> ListPtr;

  Tidx n = set[set.size() - 1] + 1;
  Tidx n_largest = LargestMovedPoint(StrongGeneratorsStabChain(g));
  if (n_largest > n)
    n = n_largest;
#ifdef DEBUG_NSI
  std::cerr << "CPP n=" << int(n) << "\n";
  std::cerr << "DEBUG MATCH set=" << GapStringIntVector(set) << "\n";
#endif
  StabChain<Telt, Tidx_label> s = CopyStabChain(g);
  Tidx m = Tidx(set.size());
  Node root_v;
  root_v.image = set;
  root_v.deleted = false;
  root_v.next = nullptr;
  root_v.prev = nullptr;
  root_v.parent = nullptr;
  // unset values
  //  root_v.selectedbaselength = max_val_type;
  root_v.IsBoundChildren = false;
  NodePtr root = std::make_shared<Node>(root_v);
  // no need to put root in the list of nodes to be deleted as the setting of
  // all to nullptr eventually kills it.
  //  ListPtr.push_back(root);

  // Node exploration functions
  auto leftmost_node = [&](Tidx const &depth) -> NodePtr {
#ifdef DEBUG_NSI
    std::cerr << "CPP Beginning of leftmost_node\n";
#endif
    NodePtr n = root;
    while (Tidx(n->selected.size()) < depth)
      n = n->children[0];
    return n;
  };
  auto next_node = [&](NodePtr const &node) -> NodePtr {
#ifdef DEBUG_NSI
    std::cerr << "CPP Beginning of next_node\n";
#endif
    NodePtr n = node;
    while (true) {
      n = n->next;
      if (n == nullptr || !n->deleted)
        break;
    }
    return n;
  };
  // Delete a node, and recursively deleting all it's children.
  std::function<void(NodePtr &)> delete_node = [&](NodePtr &node) -> void {
#ifdef DEBUG_NSI
    std::cerr << "CPP Beginning of delete_node\n";
#endif
    if (node->deleted) {
      return;
    }
    if (node->prev != nullptr) {
      node->prev->next = node->next;
    }
    if (node->next != nullptr) {
      node->next->prev = node->prev;
    }
    node->deleted = true;
    if (node->parent != nullptr) {
      Remove(node->parent->children, node->childno);
      if (node->parent->children.size() == 0) {
        delete_node(node->parent);
      } else {
        for (Tidx i = node->childno; i < Tidx(node->parent->children.size());
             i++) {
          node->parent->children[i]->childno = i;
        }
      }
    }
    if (node->IsBoundChildren) {
      for (auto &enode : node->children)
        delete_node(enode);
    }
  };

  // Given a group 'gp' and a set 'set', find orbit representatives
  // of 'set' in 'gp' simply.
  std::vector<Tidx> q_sor;
  q_sor.reserve(n);

  auto DifferenceVect_local =
      [](Tidx const &m, std::vector<Tidx> const &LIdx) -> std::vector<Tidx> {
    Face b(m);
    for (auto &eVal : LIdx)
      b[eVal] = 1;
    size_t len_ret = size_t(m) - LIdx.size();
    size_t pos = 0;
    std::vector<Tidx> cands(len_ret);
    for (Tidx i = 0; i < m; i++)
      if (b[i] == 0) {
        cands[pos] = i;
        pos++;
      }
    return cands;
  };
  // We need to break all the cycles in order to the memory free to happen
  // correctly. We use a hack in order to get that behavior: A vector of all the
  // nodes, then set the pointer to zero and so all cycles eliminated. Maybe we
  // could do better, but the hack should be adequate.
  auto free_all_nodes = [&]() -> void {
    //    std::cerr << "|ListPtr|=" << ListPtr.size() << "\n";
    for (auto &e_node : ListPtr) {
      e_node->prev = nullptr;
      e_node->next = nullptr;
      e_node->parent = nullptr;
    }
  };
  if (set.size() == 0) {
    free_all_nodes();
    return {};
  }
  Tidx depth;
  std::vector<Tidx> orbmins;
  orbmins.reserve(n);
  std::vector<Tidx> orbsizes;
  orbsizes.reserve(n);
  std::vector<Tidx> q;
  q.reserve(n);
  std::vector<Tidx> orbnums(n, max_val_type);
  for (depth = 0; depth < m; depth++) {
#ifdef DEBUG_NSI
    std::cerr << "CPP depth=" << int(depth + 1) << "\n";
#endif
    std::vector<Telt> gens = GetListGenerators(s);
    for (Tidx i = 0; i < n; i++)
      orbnums[i] = max_val_type;
    orbmins.clear();
    orbsizes.clear();
    Tidx upb = initial_upb;
    // Make orbit of x, updating orbnums, orbmins and orbsizes as approriate.
    auto make_orbit = [&](Tidx const &x) -> Tidx {
#ifdef DEBUG_NSI
      if (orbnums[x] == max_val_type) {
        std::cerr << "CPP Beginning of make_orbit x=" << int(x + 1)
                  << " orbnums[x]=" << int(-1) << "\n";
      } else {
        std::cerr << "CPP Beginning of make_orbit x=" << int(x + 1)
                  << " orbnums[x]=" << int(orbnums[x] + 1) << "\n";
      }
#endif
      if (orbnums[x] != max_val_type) {
        return orbnums[x];
      }
#ifdef DEBUG_NSI
      std::cerr << "CPP After check\n";
#endif
      q.clear();
      q.push_back(x);
      Tidx rep = x;
      Tidx num = Tidx(orbmins.size());
      orbnums[x] = num;
#ifdef DEBUG_NSI
      std::cerr << "CPP 1 : x=" << int(x + 1)
                << " Assign orbnums[x]=" << int(orbnums[x] + 1) << "\n";
#endif
      size_t pos = 0;
      while (true) {
        size_t idx, qsiz = q.size();
        for (idx = pos; idx < qsiz; idx++) {
          Tidx pt = q[idx];
          for (auto &gen : gens) {
            Tidx img = PowAct(pt, gen);
            if (orbnums[img] == max_val_type) {
              orbnums[img] = num;
#ifdef DEBUG_NSI
              std::cerr << "CPP 2 : img=" << int(img + 1)
                        << " Assign orbnums[img]=" << int(orbnums[img] + 1)
                        << "\n";
#endif
              q.push_back(img);
              if (img < rep)
                rep = img;
            }
          }
        }
        if (idx == q.size())
          break;
        pos = idx;
      }
#ifdef DEBUG_NSI
      std::cerr << "CPP make_orbit rep=" << int(rep + 1) << " |q|=" << q.size()
                << " num=" << int(num + 1) << "\n";
#endif
      orbmins.push_back(rep);
      orbsizes.push_back(Tidx(q.size()));
      return num;
    };
    /*
      # At this point, all bottom nodes are blue
      # first pass creates appropriate set of virtual red nodes
    */

    std::vector<Tidx> minOrbitMset = {infinity};
    NodePtr node = leftmost_node(depth);
    while (node != nullptr) {
#ifdef DEBUG_NSI
      std::cerr << "CPP m=" << m
                << " node.selected=" << GapStringIntVector(node->selected)
                << "\n";
#endif
      std::vector<Tidx> cands = DifferenceVect_local(m, node->selected);
#ifdef DEBUG_NSI
      std::cerr << "CPP 1 : cands=" << GapStringIntVector(cands) << "\n";
#endif

      std::vector<Tidx> orbitMset;
      for (auto &y : cands) {
        Tidx x = node->image[y];
        Tidx num = make_orbit(x);
#ifdef DEBUG_NSI
        std::cerr << "CPP x=" << int(x + 1) << " num=" << int(num + 1) << "\n";
#endif
        orbitMset.push_back(orbmins[num]);
      }
#ifdef DEBUG_NSI
      std::cerr << "CPP bef orbitMset=" << GapStringIntVector(orbitMset)
                << "\n";
#endif
      std::sort(orbitMset.begin(), orbitMset.end());
#ifdef DEBUG_NSI
      std::cerr << "CPP aft orbitMset=" << GapStringIntVector(orbitMset)
                << "\n";
#endif
      if (orbitMset < minOrbitMset) {
#ifdef DEBUG_NSI
        std::cerr << "CPP orbitMset comparison case 1\n";
#endif
        minOrbitMset = orbitMset;
        NodePtr node2 = node->prev;
        while (node2 != nullptr) {
          delete_node(node2);
          node2 = node2->prev;
        }
      } else {
        if (orbitMset > minOrbitMset) {
#ifdef DEBUG_NSI
          std::cerr << "CPP orbitMset comparison case 2\n";
#endif
          delete_node(node);
        }
      }
      node = next_node(node);
    }

    std::vector<Tidx> globalOrbitCounts(orbmins.size(), 0);
    node = leftmost_node(depth);
    while (node != nullptr) {
      std::vector<Tidx> cands = DifferenceVect_local(m, node->selected);
#ifdef DEBUG_NSI
      std::cerr << "CPP 2 : cands=" << GapStringIntVector(cands) << "\n";
#endif
      /*
        # These index the children of node that will
        # not be immediately deleted under rule C
      */
      for (auto &y : cands) {
        Tidx x = node->image[y];
#ifdef DEBUG_NSI
        std::cerr << "CPP y=" << int(y + 1) << " x=" << int(x + 1) << "\n";
#endif
        Tidx num = make_orbit(x);
        Tidx siz = Tidx(globalOrbitCounts.size());
        if (num < siz) {
          globalOrbitCounts[num]++;
        } else {
          for (Tidx u = siz; u <= num; u++)
            globalOrbitCounts.push_back(0);
          globalOrbitCounts[num] = 1;
        }
#ifdef DEBUG_NSI
        std::cerr << "CPP globalOrbitCounts : num=" << int(num + 1)
                  << " cnt=" << globalOrbitCounts[num] << "\n";
#endif
      }
      node = next_node(node);
    }
#ifdef DEBUG_NSI
    std::cerr << "CPP globalOrbitCounts=" << GapStringTVector(globalOrbitCounts)
              << "\n";
#endif
    Tidx globalBestOrbit =
        calculateBestOrbit(orbmins, globalOrbitCounts, orbsizes);
    upb = orbmins[globalBestOrbit];
#ifdef DEBUG_NSI
    std::cerr << "CPP globalBestOrbit=" << int(globalBestOrbit + 1)
              << " upb=" << int(upb + 1) << "\n";
#endif

    node = leftmost_node(depth);
    while (node != nullptr) {
      std::vector<Tidx> cands = DifferenceVect_local(m, node->selected);
#ifdef DEBUG_NSI
      std::cerr << "CPP 3 : cands=" << GapStringIntVector(cands) << "\n";
#endif
      /*
        # These index the children of node that will
        # not be immediately deleted under rule C
      */
      node->validkids.clear();
      for (auto &y : cands) {
        Tidx x = node->image[y];
        Tidx num = orbnums[x];
        if (num == max_val_type) {
          /*
            # Need a new orbit. Also require the smallest point
            # as the rep.
            #
            #
            # If there is no prospect of the new orbit being
            # better than the current best then go on to the next candidate
          */
          num = make_orbit(x);
          Tidx rep = orbmins[num];
          if (rep < upb) {
            upb = rep;
            NodePtr node2 = node->prev;
            while (node2 != nullptr) {
              delete_node(node2);
              node2 = node2->prev;
            }
            node->validkids.clear();
            node->validkids.push_back(y);
#ifdef DEBUG_NSI
            std::cerr << "CPP validkids set to {y} with y=" << int(y + 1)
                      << "\n";
#endif
          }
        } else {
          Tidx rep = orbmins[num];
#ifdef DEBUG_NSI
          std::cerr << "CPP before insertion rep=" << int(rep + 1)
                    << " num=" << int(num + 1) << " upb=" << int(upb + 1)
                    << "\n";
#endif
          if (rep == upb) {
            node->validkids.push_back(y);
#ifdef DEBUG_NSI
            std::cerr << "CPP validkids inserting y=" << int(y + 1) << "\n";
#endif
          }
        }
      }
      if (node->validkids.size() == 0) {
        delete_node(node);
      }
      node = next_node(node);
    }
    /*
      # Second pass. Actually make all the red nodes and turn them blue
    */
#ifdef DEBUG_NSI
    std::cerr << "CPP Before ChangeStabChain\n";
#endif
    ChangeStabChain(s, {upb}, false);
#ifdef DEBUG_NSI
    std::cerr << "CPP After ChangeStabChain\n";
#endif
    bool do_continue = false;
    if (s->orbit.size() == 1) {
      /*
        # In this case nothing much can happen. Each surviving node will have
        exactly one child # and none of the imsets will change # so we mutate
        the nodes in-place
      */
      node = leftmost_node(depth);
      while (node != nullptr) {
        //        if (node->selectedbaselength == max_val_type)
        //          node->selectedbaselength = Tidx(node->selected.size());
        node->selected.push_back(node->validkids[0]);
#ifdef DEBUG_NSI
        std::cerr << "CPP Now node.selected="
                  << GapStringIntVector(node->selected) << "\n";
#endif
        node = next_node(node);
      }
      s = s->stabilizer;
      if (Tidx(leftmost_node(depth + 1)->selected.size()) != m) {
        do_continue = true;
      }
    }

    if (!do_continue) {
      node = leftmost_node(depth);
      NodePtr prevnode = nullptr;
      while (node != nullptr) {
        node->IsBoundChildren = true;
        node->children.clear();
#ifdef DEBUG_NSI
        std::cerr << "CPP node.validkids="
                  << GapStringIntVector(node->validkids) << "\n";
#endif
        for (auto &x : node->validkids) {
          Node newnode_v;
          newnode_v.selected = node->selected;
          newnode_v.selected.push_back(x);
#ifdef DEBUG_NSI
          std::cerr << "DEBUG Before Stabilize_OnPoints x=" << int(x + 1)
                    << "\n";
#endif
#ifdef DEBUG_NSI
          std::cerr << "DEBUG After Stabilize_OnPoints\n";
#endif
          newnode_v.parent = node;
          newnode_v.childno = Tidx(node->children.size());
          newnode_v.next = nullptr;
          newnode_v.prev = prevnode;
          newnode_v.deleted = false;
          newnode_v.IsBoundChildren = false;
#ifdef DEBUG_NSI
          std::cerr << "CPP newnode.selected="
                    << GapStringIntVector(newnode_v.selected) << "\n";
#endif
          NodePtr newnode = std::make_shared<Node>(newnode_v);
          ListPtr.push_back(newnode);
          if (ListPtr.size() == max_size) {
            free_all_nodes();
            return {};
          }
          if (prevnode != nullptr) {
            prevnode->next = newnode;
          }
          prevnode = newnode;
          node->children.push_back(newnode);

          std::vector<Tidx> image = node->image;
          if (image[x] != upb) {
            while (true) {
              const Telt &g = s->comm->labels[s->transversal[image[x]]];
              OnTuples_inplace(image, g);
              if (image[x] == upb)
                break;
            }
            newnode->image = image;
          } else {
            newnode->image = image;
          }
        }
        node = next_node(node);
      }

#ifdef DEBUG_NSI
      std::cerr << "CPP Before s:=s.stabilizer operation\n";
#endif
      s = s->stabilizer;
      if (Tidx(leftmost_node(depth + 1)->selected.size()) == m) {
        break;
      }
    }
  }
  size_t n_node = ListPtr.size();
  free_all_nodes();
  NodePtr node = leftmost_node(depth + 1);
  std::pair<std::vector<typename Telt::Tidx>,size_t> pair{node->image,n_node};
  return pair;
}


template<typename Tidx>
std::pair<std::vector<Tidx>,bool> get_input_info(Face const& set) {
  size_t siz = set.size();
  Face ret(siz);
  size_t cnt = set.count();
  if (2 * cnt <= siz) {
    std::vector<Tidx> set_i(cnt);
    size_t pos = 0;
    boost::dynamic_bitset<>::size_type aRow = set.find_first();
    while (aRow != boost::dynamic_bitset<>::npos) {
      set_i[pos] = Tidx(aRow);
      pos++;
      aRow = set.find_next(aRow);
    }
    return {std::move(set_i), true};
  } else {
    std::vector<Tidx> set_i(siz - cnt);
    Tidx siz_i = Tidx(siz);
    size_t pos = 0;
    for (Tidx i = 0; i < siz_i; i++) {
      int val = set[i];
      if (val == 0) {
        set_i[pos] = i;
        pos++;
      }
    }
    return {std::move(set_i), false};
  }
}

template<typename Tidx>
Face convert_vector_out(size_t siz, std::vector<Tidx> const& out_v, bool const& cnt_bool) {
  Face ret(siz);
  if (cnt_bool) {
    for (auto &eVal : out_v) {
      ret[eVal] = 1;
    }
  } else {
    for (size_t i = 0; i < siz; i++)
      ret[i] = 1;
    for (auto &eVal : out_v)
      ret[eVal] = 0;
  }
  return ret;
}



template <typename Telt, typename Tidx_label, typename Tint, typename F>
std::pair<Face,StabChain<Telt,Tidx_label>> Kernel_GeneralCanonicalImagePair(StabChain<Telt, Tidx_label> const &g,
                                                                            Face const &set, F f) {
  using Tidx = typename Telt::Tidx;
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
  if (g->comm->n != Tidx(set.size())) {
    std::cerr << "The set set should have the same size as the number of "
                 "elements on which g acts\n";
    std::cerr << "g->comm->n=" << int(g->comm->n) << " |set|=" << set.size()
              << "\n";
    throw PermutalibException{1};
  }
#endif
  if (set.count() == 0 || set.count() == set.size())
    return {set, g};
  StabChain<Telt, Tidx_label> k_group =
    Kernel_Stabilizer_OnSets<Telt, Tidx_label, Tint>(g, set);
  size_t siz = set.size();
  std::pair<std::vector<Tidx>,bool> pair = get_input_info<Tidx>(set);
  std::pair<std::vector<Tidx>,StabChain<Telt,Tidx_label>> pairCan =
    NewCanonicImage<Telt, Tidx_label, Tint>(g, pair.first, k_group);
  Face ret = convert_vector_out(siz, pairCan.first, pair.second);
  return f(k_group, ret, pairCan);
}

template <typename Telt, typename Tidx_label, typename Tint>
std::pair<Face,size_t> Kernel_GeneralCanonicalInitialTriv(StabChain<Telt, Tidx_label> const &g, Face const &set, size_t const& max_size) {
  using Tidx = typename Telt::Tidx;
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
  if (g->comm->n != Tidx(set.size())) {
    std::cerr << "The set set should have the same size as the number of "
                 "elements on which g acts\n";
    std::cerr << "g->comm->n=" << int(g->comm->n) << " |set|=" << set.size()
              << "\n";
    throw PermutalibException{1};
  }
#endif
  if (set.count() == 0 || set.count() == set.size())
    return {set, 0};
  size_t siz = set.size();
  std::pair<std::vector<Tidx>,bool> pair = get_input_info<Tidx>(set);
  std::optional<std::pair<std::vector<Tidx>,size_t>> opt_can =
    NewCanonicImageInitialTriv<Telt, Tidx_label, Tint>(g, pair.first, max_size);
  if (opt_can) {
    std::pair<std::vector<Tidx>,size_t> const& pair_can = *opt_can;
    Face ret = convert_vector_out(siz, pair_can.first, pair.second);
    std::pair<Face,size_t> pair_ret{std::move(ret), pair_can.second};
    return pair_ret;
  }
#ifdef PERMUTALIB_TRACK_METHOD
  Tint order = Order<Telt, Tidx_label, Tint>(g);
  std::cerr << "Method change in Kernel_GeneralCanonicalInitialTriv |set|=" << set.size() << "/" << set.count() << " |g|=" << order << "\n";
#endif
  // The happy path failed, now using the computational strategy
  StabChain<Telt, Tidx_label> k_group =
    Kernel_Stabilizer_OnSets<Telt, Tidx_label, Tint>(g, set);
  std::pair<std::vector<Tidx>,StabChain<Telt,Tidx_label>> pairCan =
    NewCanonicImage<Telt, Tidx_label, Tint>(g, pair.first, k_group);
  Face ret = convert_vector_out(siz, pairCan.first, pair.second);
  std::pair<Face,size_t> pair_ret{std::move(ret), max_size};
  return pair_ret;
}

template <typename Telt, typename Tidx_label, typename Tint>
std::pair<Face,StabChain<Telt,Tidx_label>> CanonicalImage_SubgroupStabilizer(StabChain<Telt, Tidx_label> const &g,
                                                                             Face const &set) {
  using Tidx = typename Telt::Tidx;
  auto f=[&]([[maybe_unused]] StabChain<Telt, Tidx_label> const& k, Face const& ret, std::pair<std::vector<Tidx>,StabChain<Telt,Tidx_label>> const& pair) -> std::pair<Face,StabChain<Telt,Tidx_label>> {
    return {ret, pair.second};
  };
  return Kernel_GeneralCanonicalImagePair<Telt,Tidx_label,Tint>(g, set, f);
}

template <typename Telt, typename Tidx_label, typename Tint>
std::pair<Face,StabChain<Telt,Tidx_label>> CanonicalImage_ConjugateStabilizer(StabChain<Telt, Tidx_label> const &g,
                                                                             Face const &set) {
  using Tidx = typename Telt::Tidx;
  auto f=[&](StabChain<Telt, Tidx_label> const& k, Face const& ret, [[maybe_unused]] std::pair<std::vector<Tidx>,StabChain<Telt,Tidx_label>> const& pair) -> std::pair<Face,StabChain<Telt,Tidx_label>> {
    return {ret, k};
  };
  return Kernel_GeneralCanonicalImagePair<Telt,Tidx_label,Tint>(g, set, f);
}

template <typename Telt, typename Tidx_label, typename Tint>
Face Kernel_CanonicalImage(StabChain<Telt, Tidx_label> const &g,
                           Face const &set) {
  return CanonicalImage_SubgroupStabilizer<Telt,Tidx_label,Tint>(g, set).first;
}

// clang-format off
}  // namespace permutalib
#endif  // SRC_GAP_NSI_H_
// clang-format on
