#ifndef DEFINE_PERMUTALIB_NSI_H
#define DEFINE_PERMUTALIB_NSI_H

#include "StabChain.h"
#include "stbcbckt.h"
#include <unordered_map>


//#define DEBUG_NSI

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
# Returns a pair [image, stabilizer], where stabilizer is a subgroup of Stabilizer(G, set), possibly larger than the one given into the function.
#
# Note: The return type of this is NOT a set, but provides the pointwise mapping of the input set.
# This means the permutation which actually provides the smallest image can be found cheaply as follows:
#
# res := NewSmallestImage(G, set, Group(()), x -> x);
# perm := RepresentativeAction(G, set, res[1], OnTuples);



#
# Search node data:
#
#  selected -- indices in set of points being mapped to minimal image at this node
#  image    -- sequence-wise image of set under element represented by this node
#  substab --  Stab_K(selected) sequence stabilizer
#  children --  nodes corresponding to extensions of selected
#  parent -- node corr to all but last element of selected.]
#  childno
#  next -- across row
#  prev -- likewise
#  deleted
#

# At each level

# Find the next pt of minimum image and all the corresponding nodes

# If at any time in this process we notice a potential prune, we have a
# generator of K -- construct it, close K with it and delete all tree
# branches that are now non-canonical -- propagate new K down through
# tree. Continue with surviving nodes at current level.

*/


namespace permutalib {

template<typename T>
void Remove(std::vector<T> & eV, int const& pos)
{
  eV.erase(eV.begin() + pos);
}


template<typename Telt>
struct ResultCanonicalization {
  std::vector<int> set;
  Telt g;
};


template<typename Telt,typename T>
Telt PermListList(std::vector<T> const& list1, std::vector<T> const& list2)
{
  std::unordered_map<T,int> eMap;
  int siz = list2.size();
  for (int i=0; i<siz; i++) {
    eMap[list1[i]] = i;
  }
  //
  std::vector<int> eList;
  for (int i=0; i<siz; i++) {
    int val = eMap[list2[i]];
    eList[i] = val;
  }
  return Telt(eList);
}




template<typename Telt, typename Tint>
StabChain<Telt> Action(StabChain<Telt> const& S, std::vector<typename Telt::Tidx> const& set)
{
  using Tidx = typename Telt::Tidx;
  int n = S->comm->n;
  Tidx siz = set.size();
  std::vector<Tidx> map_idx(n, 0);
  for (Tidx i=0; i<siz; i++)
    map_idx[set[i]] = i;
  //
  std::vector<Telt> LGen;
  for (auto & eGen : Kernel_GeneratorsOfGroup(S)) {
    std::vector<Tidx> eList(siz);
    for (Tidx i=0; i<siz; i++) {
      Tidx ePt = set[i];
      Tidx fPt = OnPoints(ePt, eGen);
      Tidx j = map_idx[fPt];
      eList[i] = j;
    }
    Telt ePerm = Telt(eList);
    LGen.emplace_back(ePerm);
  }
  return FCT_Group<Telt,Tint>(LGen, siz);
}

template<typename Telt, typename Tidx>
std::vector<Tidx> OnTuples(std::vector<Tidx> const& V, Telt const& g)
{
  Tidx len = V.size();
  std::vector<Tidx> retV(len);
  for (Tidx i=0; i<len; i++)
    retV[i] = PowAct(V[i], g);
  return retV;
}


template<typename T>
std::vector<T> Set(std::vector<T> const& V)
{
  std::vector<T> retV = V;
  std::sort(retV.begin(), retV.end());
  return retV;
}

/*
  Modification done:
  --- skip_fnuc eliminated as it is the identity in the case that interest us.
 */
template<typename Telt, typename Tint>
std::vector<typename Telt::Tidx> NewCanonicImage(StabChain<Telt> const& g, std::vector<typename Telt::Tidx> const& set, StabChain<Telt> const& k)
{
  using Tidx = typename Telt::Tidx;
#ifdef DEBUG_NSI
  std::cerr << "CPP NewCanonicImage : beginning\n";
#endif
  Tidx infinity = std::numeric_limits<Tidx>::max();
  Tidx initial_upb = infinity;
  Tidx max_val_type = std::numeric_limits<Tidx>::max();


  auto calculateBestOrbit=[&](PreAllocatedVector<Tidx> const& orbmins, std::vector<Tidx> const& orbitCounts, PreAllocatedVector<Tidx> const& orbsizes) -> Tidx {
#ifdef DEBUG_NSI
    std::cerr << "CPP Beginning of calculateBestOrbit\n";
#endif
    struct typeCnt {
      Tidx orbSize;
      Tidx orbCount;
    };
    // The comparison goes
    // log(a.orbCount) / a.orbSize < log(b.orbCount) / b.orbSize
    // whish is equivalent to
    // a.orbCount ^ b.orbSize    <    b.orbCount ^ a.orbSize
    auto comparisonLower=[](typeCnt const& a, typeCnt const& b) -> bool {
#ifdef DEBUG_NSI
      std::cerr << "CPP compLower a=" << a.orbCount << " / " << a.orbSize << " b=" << b.orbCount << " / " << b.orbSize << "\n";
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
        Tint pow1 = 1;
        Tint eV = a.orbCount;
        for (Tidx idx=0; idx<b.orbSize; idx++)
          pow1 *= eV;
        //
        Tint pow2 = 1;
        eV = b.orbCount;
        for (Tidx idx=0; idx<a.orbSize; idx++)
          pow2 *= eV;
        //
#ifdef DEBUG_NSI
        std::cerr << "CPP pow1=" << pow1 << " pow2=" << pow2 << "\n";
#endif
        return pow1 < pow2;
      }
    };
    auto comparisonEqual=[](typeCnt const& a, typeCnt const& b) -> bool {
#ifdef DEBUG_NSI
      std::cerr << "CPP compEqual a=" << a.orbCount << " / " << a.orbSize << " b=" << b.orbCount << " / " << b.orbSize << "\n";
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
        Tint pow1 = 1;
        Tint eV = a.orbCount;
        for (Tidx idx=0; idx<b.orbSize; idx++)
          pow1 *= eV;
        //
        Tint pow2 = 1;
        eV = b.orbCount;
        for (Tidx idx=0; idx<a.orbSize; idx++)
          pow2 *= eV;
        //
        return pow1 == pow2;
      }
    };
    auto selector=[&](Tidx const& jdx) -> typeCnt {
      return {orbsizes[jdx], orbitCounts[jdx]};
      /*
      if (orbsizes[i] == 1) {
        return - std::pow(2.0, 32.0) + orbitCounts[i];
      } else {
        return (log(double(orbitCounts[i]))/log(2.0)) / double(orbsizes[i]);
      }*/
    };
    Tidx index = 0;
    typeCnt result_0 = selector(0);
    Tidx result_1 = orbmins[0];
#ifdef DEBUG_NSI
    std::cerr << "CPP result_0=" << result_0.orbCount << " / " << result_0.orbSize << "   result_1=" << int(result_1+1) << "\n";
#endif
    for (size_t i=1; i<orbmins.size(); i++) {
      typeCnt ret_0 = selector(i);
      Tidx ret_1 = orbmins[i];
      bool lower=false;
#ifdef DEBUG_NSI
      std::cerr << "CPP i=" << int(i+1) << " ret_0=" << ret_0.orbCount << " / " << ret_0.orbSize << "   ret_1=" << int(ret_1+1) << "\n";
#endif
      if (comparisonLower(ret_0, result_0)) {
        lower=true;
      } else {
        if (comparisonEqual(ret_0, result_0)) {
          if (ret_1 < result_1)
            lower=true;
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
    StabChain<Telt> substab;
    bool deleted;
    std::shared_ptr<Node> next;
    std::shared_ptr<Node> prev;
    std::shared_ptr<Node> parent;
    Tidx selectedbaselength;
    // children
    Tidx childno;
    bool IsBoundChildren;
    std::vector<std::shared_ptr<Node>> children;
    std::vector<Tidx> validkids;
  };
  using NodePtr = std::shared_ptr<Node>;
  std::vector<NodePtr> ListPtr;

  Tidx n = set[set.size()-1] + 1;
  Tidx n_largest = LargestMovedPoint(StrongGeneratorsStabChain(g));
  if (n_largest > n)
    n = n_largest;
#ifdef DEBUG_NSI
  std::cerr << "CPP n=" << int(n) << "\n";
  std::cerr << "DEBUG MATCH set=" << GapStringIntVector(set) << "\n";
#endif
  StabChain<Telt> s = CopyStabChain(g);
  StabChain<Telt> l = Action<Telt,Tint>(k, set);
  Tidx m = set.size();
  Node root_v;
  root_v.image = set;
  root_v.substab = l;
  root_v.deleted = false;
  root_v.next = nullptr;
  root_v.prev = nullptr;
  root_v.parent = nullptr;
  // unset values
  root_v.selectedbaselength = max_val_type;
  root_v.IsBoundChildren = false;
  NodePtr root = std::make_shared<Node>(root_v);
  // no need to put root in the list of nodes to be deleted as the setting of all to nullptr eventually kills it.
  //  ListPtr.push_back(root);

  // Node exploration functions
  auto leftmost_node =[&](Tidx const& depth) -> NodePtr {
#ifdef DEBUG_NSI
    std::cerr << "CPP Beginning of leftmost_node\n";
#endif
    NodePtr n = root;
    while (Tidx(n->selected.size()) < depth)
      n = n->children[0];
    return n;
  };
  auto next_node =[&](NodePtr const& node) -> NodePtr {
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
  //Delete a node, and recursively deleting all it's children.
  std::function<void(NodePtr &)> delete_node=[&](NodePtr & node) -> void {
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
        for (Tidx i=node->childno; i<Tidx(node->parent->children.size()); i++) {
          node->parent->children[i]->childno = i;
        }
      }
    }
    if (node->IsBoundChildren) {
      for (auto & enode : node->children)
        delete_node(enode);
    }
  };
  /*
  auto delete_nodes=[&](std::vector<NodePtr> & nodes) {
    for (auto & e_node : nodes)
      delete_node(e_node);
  };
  */
  
  // Filter nodes by stabilizer group,
  // Updates the stabilizer group of the node,
  /*
  std::function<void(NodePtr &)> clean_subtree =[&](NodePtr & node) -> void {
#ifdef DEBUG_NSI
    std::cerr << "CPP Beginning of clean_subtree\n";
#endif
    if (!node->IsBoundChildren)
      return;
    std::vector<NodePtr> bad;

    Face seen(m);
    Tidx x;
    for (auto & c : node->children) {
      if (c->selectedbaselength != max_val_type) {
        x = c->selected[c->selectedbaselength];
      } else {
        x = c->selected[c->selected.size() - 1];
      }
      if (seen[x] == 1) {
        bad.push_back(c);
      } else {
        std::vector<Tidx> q = {x};
        std::vector<Telt> gens = Kernel_GeneratorsOfGroup(node->substab);
        size_t olen = 1;
        size_t pos = 0;
        seen[x] = 1;
        while (true) {
          size_t idx;
          for (idx=pos; idx<olen; idx++) {
            Tidx pt = q[idx];
            for (auto & gen : gens) {
              Tidx im = PowAct(pt,gen);
              if (seen[im] == 0) {
                seen[im] = 1;
                q.emplace_back(im);
                olen++;
              }
            }
          }
          if (idx == olen)
            break;
          pos = idx;
        }
        Tint quot = Order<Telt,Tint>(node->substab) / Order<Telt,Tint>(c->substab);
        if (Tint(olen) < quot) {
          c->substab = Kernel_Stabilizer_OnPoints<Telt,Tint>(node->substab, x);
          clean_subtree(c);
        }
      }
    }
    delete_nodes(bad);
  };
  */
  //Add a new stabilizer element, mapping node1 to node2, and then call
  // clean_subtree to remove any new subtrees.
  // unused in this specific code.
  /*
  auto handle_new_stabilizer_element=[&](NodePtr & node1, NodePtr & node2) -> void {
    // so node1 and node2 represnet group elements that map set to the same
    // place in two different ways
    Telt perm1 = PermListList<Telt,int>(node1->image, node2->image);
    StabChain<Telt> l_copy = CopyStabChain(l);
    ClosureGroup<Telt,Tint>(l_copy,perm1);
    root->substab = l_copy;
    clean_subtree(root);
  };
  */

  // Given a group 'gp' and a set 'set', find orbit representatives
  // of 'set' in 'gp' simply.
  PreAllocatedVector<Tidx> q_sor(n);
  Face b_sor(n);
  auto simpleOrbitReps=[&](StabChain<Telt> const& gp, std::vector<Tidx> const& set) -> std::vector<Tidx> {
#ifdef DEBUG_NSI
    std::cerr << "CPP Beginning of simpleOrbitReps\n";
#endif
    Tidx m = set.size();
    Tidx n = set[m-1] + 1;
    for (Tidx i=0; i<n; i++)
      b_sor[i] = 0;
#ifdef DEBUG_NSI
    std::cerr << "DEBUG n=" << int(n) << "\n";
#endif
    for (auto & eVal : set) {
#ifdef DEBUG_NSI
      std::cerr << "DEBUG eVal=" << int(eVal) << " n=" << int(n) << "\n";
#endif
      b_sor[eVal] = 1;
    }
    boost::dynamic_bitset<>::size_type seed=b_sor.find_first();
    Tidx seed_tidx = Tidx(seed);
    std::vector<Tidx> reps;
    std::vector<Telt> gens = Kernel_GeneratorsOfGroup(gp);
    while (seed != boost::dynamic_bitset<>::npos && seed_tidx < n) {
#ifdef DEBUG_NSI
      std::cerr << "DEBUG seed=" << int(seed) << " n=" << int(n) << "\n";
#endif
      b_sor[seed]=0;
      q_sor.clear();
      q_sor.push_back(seed_tidx);
      reps.push_back(seed_tidx);
      size_t pos=0;
      while(true) {
        size_t idx, qsiz=q_sor.size();
        for (idx=pos; idx<qsiz; idx++) {
          Tidx pt = q_sor[idx];
          for (auto & gen : gens) {
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
  auto DifferenceVect_local=[](Tidx const& m, std::vector<Tidx> const& LIdx) -> std::vector<Tidx> {
     Face b(m);
     for (auto & eVal : LIdx)
       b[eVal] = 1;
     std::vector<Tidx> cands;
     for (Tidx i=0; i<m; i++)
       if (b[i] == 0)
         cands.push_back(i);
     return cands;
  };
  // We need to break all the cycles in order to the memory free to happen correctly.
  // We use a hack in order to get that behavior: A vector of all the nodes, then
  // set the pointer to zero and so all cycles eliminated.
  // Maybe we could do better, but the hack should be adequate.
  auto free_all_nodes=[&]() -> void {
                        //    std::cerr << "|ListPtr|=" << ListPtr.size() << "\n";
    for (auto & e_node : ListPtr) {
      e_node->prev=nullptr;
      e_node->next=nullptr;
      e_node->parent=nullptr;
    }
  };
  if (set.size() == 0) {
    free_all_nodes();
    return {};
  }
  Tidx depth;
  PreAllocatedVector<Tidx> orbmins(n);
  PreAllocatedVector<Tidx> orbsizes(n);
  PreAllocatedVector<Tidx> q(n);
  std::vector<Tidx> orbnums(n,max_val_type);
  for (depth=0; depth<m; depth++) {
#ifdef DEBUG_NSI
    std::cerr << "CPP depth=" << int(depth+1) << "\n";
#endif
    std::vector<Telt> gens = GetListGenerators(s);
    for (Tidx i=0; i<n; i++)
      orbnums[i] = max_val_type;
    orbmins.clear();
    orbsizes.clear();
    Tidx upb = initial_upb;
    // Make orbit of x, updating orbnums, orbmins and orbsizes as approriate.
    auto make_orbit=[&](Tidx const& x) -> Tidx {
#ifdef DEBUG_NSI
      if (orbnums[x] == max_val_type) {
        std::cerr << "CPP Beginning of make_orbit x=" << int(x+1) << " orbnums[x]=" << int(-1) << "\n";
      } else {
        std::cerr << "CPP Beginning of make_orbit x=" << int(x+1) << " orbnums[x]=" << int(orbnums[x]+1) << "\n";
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
      Tidx num = orbmins.size();
      orbnums[x] = num;
#ifdef DEBUG_NSI
      std::cerr << "CPP 1 : x=" << int(x+1) << " Assign orbnums[x]=" << int(orbnums[x] + 1) << "\n";
#endif
      size_t pos=0;
      while (true) {
        size_t idx, qsiz = q.size();
        for (idx=pos; idx<qsiz; idx++) {
          Tidx pt = q[idx];
          for (auto& gen : gens) {
            Tidx img = PowAct(pt, gen);
            if (orbnums[img] == max_val_type) {
              orbnums[img] = num;
#ifdef DEBUG_NSI
              std::cerr << "CPP 2 : img=" << int(img+1) << " Assign orbnums[img]=" << int(orbnums[img] + 1) << "\n";
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
      std::cerr << "CPP make_orbit rep=" << int(rep+1) << " |q|=" << q.size() << " num=" << int(num+1) << "\n";
#endif
      orbmins.push_back(rep);
      orbsizes.push_back(q.size());
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
      std::cerr << "CPP m=" << m << " node.selected=" << GapStringIntVector(node->selected) << "\n";
#endif
      std::vector<Tidx> cands = DifferenceVect_local(m, node->selected);
#ifdef DEBUG_NSI
      std::cerr << "CPP 1 : cands=" << GapStringIntVector(cands) << "\n";
#endif

      std::vector<Tidx> orbitMset;
      for (auto & y : cands) {
        Tidx x = node->image[y];
        Tidx num = make_orbit(x);
#ifdef DEBUG_NSI
        std::cerr << "CPP x=" << int(x+1) << " num=" << int(num+1) << "\n";
#endif
        orbitMset.push_back(orbmins[num]);
      }
#ifdef DEBUG_NSI
      std::cerr << "CPP bef orbitMset=" << GapStringIntVector(orbitMset) << "\n";
#endif
      std::sort(orbitMset.begin(), orbitMset.end());
#ifdef DEBUG_NSI
      std::cerr << "CPP aft orbitMset=" << GapStringIntVector(orbitMset) << "\n";
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
      for (auto & y : cands) {
        Tidx x = node->image[y];
#ifdef DEBUG_NSI
        std::cerr << "CPP y=" << int(y+1) << " x=" << int(x+1) << "\n";
#endif
        Tidx num = make_orbit(x);
        Tidx siz = globalOrbitCounts.size();
        if (num < siz) {
          globalOrbitCounts[num]++;
        } else {
          for (Tidx u=siz; u<=num; u++)
            globalOrbitCounts.push_back(0);
          globalOrbitCounts[num] = 1;
        }
#ifdef DEBUG_NSI
        std::cerr << "CPP globalOrbitCounts : num=" << int(num+1) << " cnt=" << globalOrbitCounts[num] << "\n";
#endif
      }
      node = next_node(node);
    }
#ifdef DEBUG_NSI
    std::cerr << "CPP globalOrbitCounts=" << GapStringTVector(globalOrbitCounts) << "\n";
#endif
    Tidx globalBestOrbit = calculateBestOrbit(orbmins, globalOrbitCounts, orbsizes);
    upb = orbmins[globalBestOrbit];
#ifdef DEBUG_NSI
    std::cerr << "CPP globalBestOrbit=" << int(globalBestOrbit+1) << " upb=" << int(upb+1) << "\n";
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
      for (auto & y : cands) {
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
            std::cerr << "CPP validkids set to {y} with y=" << int(y+1) << "\n";
#endif
          }
        } else {
          Tidx rep = orbmins[num];
#ifdef DEBUG_NSI
          std::cerr << "CPP before insertion rep=" << int(rep+1) << " num=" << int(num+1) << " upb=" << int(upb+1) << "\n";
#endif
          if (rep == upb) {
            node->validkids.push_back(y);
#ifdef DEBUG_NSI
            std::cerr << "CPP validkids inserting y=" << int(y+1) << "\n";
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
        # In this case nothing much can happen. Each surviving node will have exactly one child
        # and none of the imsets will change
        # so we mutate the nodes in-place
      */
      node = leftmost_node(depth);
      while (node != nullptr) {
        if (node->selectedbaselength == max_val_type)
          node->selectedbaselength = node->selected.size();
        node->selected.push_back(node->validkids[0]);
#ifdef DEBUG_NSI
        std::cerr << "CPP Now node.selected=" << GapStringIntVector(node->selected) << "\n";
#endif
        node = next_node(node);
      }
      s = s->stabilizer;
      if (Tidx(leftmost_node(depth+1)->selected.size()) != m) {
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
        std::cerr << "CPP node.validkids=" << GapStringIntVector(node->validkids) << "\n";
#endif
        for (auto & x : node->validkids) {
          Node newnode_v;
          std::vector<Tidx> selected = node->selected;
          selected.push_back(x);
          newnode_v.selected = selected;
          //          PrintStabChain(node->substab);
#ifdef DEBUG_NSI
          std::cerr << "DEBUG Before Stabilize_OnPoints x=" << int(x+1) << "\n";
#endif
          newnode_v.substab = Kernel_Stabilizer_OnPoints<Telt,Tint>(node->substab, x);
#ifdef DEBUG_NSI
          std::cerr << "DEBUG After Stabilize_OnPoints\n";
#endif
          newnode_v.parent = node;
          newnode_v.childno = node->children.size();
          newnode_v.next = nullptr;
          newnode_v.prev = prevnode;
          newnode_v.deleted = false;
          newnode_v.IsBoundChildren = false;
#ifdef DEBUG_NSI
          std::cerr << "CPP newnode.selected=" << GapStringIntVector(selected) << "\n";
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
              Telt g = s->comm->labels[s->transversal[image[x]]];
              image = OnTuples(image, g);
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
      if (Tidx(leftmost_node(depth+1)->selected.size()) == m) {
        break;
      }
    }
  }
  free_all_nodes();
  return leftmost_node(depth+1)->image;
}



template<typename Telt, typename Tint>
Face Kernel_CanonicalImage(StabChain<Telt> const& g, Face const& set)
{
  using Tidx = typename Telt::Tidx;
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
  if (Tidx(g->comm->n) != Tidx(set.size())) {
    std::cerr << "The set set should have the same size as the number of elements on which g acts\n";
    throw PermutalibException{1};
  }
#endif
  if (set.count() == 0 || set.count() == set.size())
    return set;
  std::vector<Tidx> set_i;
  Face ret(set.size());
  if (2 * set.count() <= set.size()) {
    boost::dynamic_bitset<>::size_type aRow = set.find_first();
    while (aRow != boost::dynamic_bitset<>::npos) {
      set_i.push_back(aRow);
      aRow = set.find_next(aRow);
    }
    StabChain<Telt> k = Kernel_Stabilizer_OnSets<Telt,Tint>(g, set);
    std::vector<Tidx> eSetCan = NewCanonicImage<Telt,Tint>(g, set_i, k);
#ifdef DEBUG_NSI
    std::cerr << "CPP eSetCan=" << GapStringIntVector(eSetCan) << "\n";
#endif
    for (auto & eVal : eSetCan) {
      ret[eVal] = 1;
    }
  } else {
    Face setC(set.size());
    for (size_t i=0; i<set.size(); i++)
      setC[i] = 1 - set[i];
    boost::dynamic_bitset<>::size_type aRow = setC.find_first();
    while (aRow != boost::dynamic_bitset<>::npos) {
      set_i.push_back(aRow);
      aRow = setC.find_next(aRow);
    }
    StabChain<Telt> k = Kernel_Stabilizer_OnSets<Telt,Tint>(g, setC);
    std::vector<Tidx> eSetCan = NewCanonicImage<Telt,Tint>(g, set_i, k);
    for (size_t i=0; i<set.size(); i++)
      ret[i] = 1;
    for (auto & eVal : eSetCan)
      ret[eVal] = 0;
  }
  return ret;
}



}
#endif
