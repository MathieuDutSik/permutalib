#ifndef INCLUDE_NEW_SMALLEST_IMAGE
#define INCLUDE_NEW_SMALLEST_IMAGE

#include "StabChain.h"
#include "stbcbckt.h"

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
  return PermList(eList);
}


template<typename Telt, typename Tint>
StabChain<Telt> Action(StabChain<Telt> const& S, std::vector<int> const& set)
{
  int n = S->comm->n;
  int siz = set.size();
  std::vector<int> map_idx(n, -1);
  for (int i=0; i<siz; i++)
    map_idx[set[i]] = i;
  //
  std::vector<Telt> LGen;
  for (auto & eGen : GeneratorsOfGroup(S)) {
    std::vector<int> eList(siz);
    for (int i=0; i<siz; i++) {
      int ePt = set[i];
      int fPt = OnPoints(ePt, eGen);
      int j = map_idx[fPt];
      eList[i] = j;
    }
    Telt ePerm = PermList(eList);
    LGen.push_back(ePerm);
  }
  return Group<Telt,Tint>(LGen, siz);
}

template<typename Telt>
std::vector<int> OnTuples(std::vector<int> const& V, Telt const& g)
{
  std::vector<int> retV(V.size());
  for (size_t i=0; i<V.size(); i++) {
    retV[i] = PowAct(V[i], g);
  }
  return retV;
}



/*
  Modification done:
  --- skip_fnuc eliminated as it is the identity in the case that interest us.
 */
template<typename Telt, typename Tint>
ResultCanonicalization<Telt> NewSmallestImage(StabChain<Telt> const& g, std::vector<int> const& set, StabChain<Telt> const& k)
{

  int infinity = 10 + g->comm->n;
  int initial_upb = infinity;


  auto calculateBestOrbit=[&](std::vector<int> const& orbmins, std::vector<int> const& orbitCounts, std::vector<int> const& orbsizes) -> int {
    auto selector=[&](int const& i) -> double {
      if (orbsizes[i] == 1) {
        return - std::pow(2.0, 32.0) + orbitCounts[i];
      } else {
        return (log(double(orbitCounts[i]))/log(2.0)) / double(orbsizes[i]);
      }
    };
    int index = 0;
    double result_0 = selector(0);
    int result_1 = orbmins[0];
    for (size_t i=1; i<orbmins.size(); i++) {
      double ret_0 = selector(1);
      int ret_1 = orbmins[i];
      bool lower=false;
      if (ret_0 < result_0) {
        lower=true;
      } else {
        if (ret_0 == result_0) {
          if (ret_1 < result_1)
            lower=true;
        }
      }
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
    std::vector<int> selected;
    std::vector<int> image;
    std::vector<int> imset;
    StabChain<Telt> substab;
    bool deleted;
    std::shared_ptr<Node> next;
    std::shared_ptr<Node> prev;
    std::shared_ptr<Node> parent;
    int selectedbaselength;
    // children
    int childno;
    bool IsBoundChildren;
    std::vector<std::shared_ptr<Node>> children;
    std::vector<int> validkids;
  };
  using NodePtr = std::shared_ptr<Node>;

  int n = std::max(LargestMovedPoint(g), set[set.size()-1]);
  StabChain<Telt> s = CopyStabChain(g);
  StabChain<Telt> l = Action<Telt,Tint>(k, set);
  int m = set.size();
  Node root_v;
  root_v.image = set;
  root_v.imset = set;
  root_v.substab = l;
  root_v.deleted = false;
  root_v.next = nullptr;
  root_v.prev = nullptr;
  root_v.parent = nullptr;
  // unset values
  root_v.selectedbaselength = -1;
  root_v.IsBoundChildren = false;
  NodePtr root = std::make_shared<Node>(root_v);

  // Node exploration functions
  auto leftmost_node =[&](int const& depth) -> NodePtr {
    NodePtr n = root;
    while (int(n->selected.size()) < depth - 1) {
      n = n->children[0];
    }
    return n;
  };
  auto next_node =[&](NodePtr const& node) -> NodePtr {
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
        for (int i=node->childno; i<int(node->parent->children.size()); i++) {
          node->parent->children[i]->childno = i;
        }
      }
    }
    if (node->IsBoundChildren) {
      for (auto & enode : node->children)
        delete_node(enode);
    }
  };
  auto delete_nodes=[&](std::vector<NodePtr> & nodes) {
    for (auto & e_node : nodes)
      delete_node(e_node);
  };

  // Filter nodes by stabilizer group,
  // Updates the stabilizer group of the node,
  std::function<void(NodePtr &)> clean_subtree =[&](NodePtr & node) -> void {
    if (!node->IsBoundChildren)
      return;
    std::vector<NodePtr> bad;

    std::vector<int> range= ClosedInterval(0,m);
    Face seen = BlistList(range,{});
    int x;
    for (auto & c : node->children) {
      if (c->selectedbaselength != -1) {
        x = c->selected[c->selectedbaselength];
      } else {
        x = c->selected[c->selected.size() - 1];
      }
      if (seen[x] == 1) {
        bad.push_back(c);
      } else {
        std::vector<int> q = {x};
        std::vector<Telt> gens = GeneratorsOfGroup(node->substab);
        size_t olen = 1;
        size_t pos = 0;
        seen[x] = 1;
        while (true) {
          size_t idx;
          for (idx=pos; idx<olen; idx++) {
            int pt = q[idx];
            for (auto & gen : gens) {
              int im = PowAct(pt,gen);
              if (seen[im] == 0) {
                seen[im] = 1;
                q.push_back(im);
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
          c->substab = Stabilizer_OnPoints<Telt,Tint>(node->substab, x);
          clean_subtree(c);
        }
      }
    }
    delete_nodes(bad);
  };

  //Add a new stabilizer element, mapping node1 to node2, and then call
  // clean_subtree to remove any new subtrees.
  auto handle_new_stabilizer_element=[&](NodePtr & node1, NodePtr & node2) -> void {
    // so node1 and node2 represnet group elements that map set to the same
    // place in two different ways
    Telt perm1 = PermListList<Telt,int>(node1->image, node2->image);
    StabChain<Telt> l_copy = CopyStabChain(l);
    ClosureGroup<Telt,Tint>(l_copy,perm1);
    root->substab = l_copy;
    clean_subtree(root);
  };

  // Given a group 'gp' and a set 'set', find orbit representatives
  // of 'set' in 'gp' simply.
  auto simpleOrbitReps=[&](StabChain<Telt> const& gp, std::vector<int> const& set) -> std::vector<int> {
    int m = set.size();
    int n = set[m-1];
    Face b = BlistList(ClosedInterval(0,n+1), set); // Check the n+1 here
    int seed = set[0];
    std::vector<int> reps;
    std::vector<Telt> gens = GeneratorsOfGroup(gp);
    while (seed != -1 && seed <= n) {
      b[seed]=0;
      std::vector<int> q = {seed};
      reps.push_back(seed);
      size_t pos=0;
      while(true) {
        size_t idx, qsiz=q.size();
        for (idx=pos; idx<qsiz; idx++) {
          int pt = q[idx];
          for (auto & gen : gens) {
            int im = PowAct(pt, gen);
            if (b[im] == 1) {
              b[im] = 0;
              q.push_back(im);
            }
          }
        }
        if (idx == q.size())
          break;
        pos = idx;
      }
      seed = b.find_next(seed);
    }
    return reps;
  };


  if (set.size() == 0) {
    return {{}, {}};
  }
  size_t depth;
  for (depth=0; depth<m; depth++) {
    std::vector<Telt> gens = s->generators;
    std::vector<int> orbnums(n,-1);
    std::vector<int> orbmins;
    std::vector<int> orbsizes;
    int upb = initial_upb;
    // Make orbit of x, updating orbnums, orbmins and orbsizes as approriate.
    auto make_orbit=[&](int const& x) -> int {
      if (orbnums[x] != -1) {
        return orbnums[x];
      }
      std::vector<int> q = {x};
      int rep = x;
      int num = orbmins.size() + 1;
      orbnums[x] = num;
      size_t pos=0;
      while (true) {
        size_t idx, qsiz = q.size();
        for (idx=pos; idx<qsiz; idx++) {
          int pt = q[idx];
          for (auto& gen : gens) {
            int img = PowAct(pt, gen);
            if (orbnums[img] == -1) {
              orbnums[img] = num;
              q.push_back(img);
              if (img < rep) {
                rep = img;
              }
            }
          }
        }
        if (idx == q.size())
          break;
        pos = idx;
      }
      orbmins.push_back(rep);
      orbsizes.push_back(q.size());
      return num;
    };
    /*
      # At this point, all bottom nodes are blue
      # first pass creates appropriate set of virtual red nodes
    */


    std::vector<int> minOrbitMset = {infinity};
    NodePtr node = leftmost_node(depth);
    while (node != nullptr) {
      std::vector<int> cands = DifferenceVect(ClosedInterval(0,m), node->selected);

      std::vector<int> orbitMset;
      for (auto & y : cands) {
        int x = node->image[y];
        int num = make_orbit(x);
        orbitMset.push_back(orbmins[num]);
      }
      std::sort(orbitMset.begin(), orbitMset.end());
      if (IsBound(bestOrbitMset)) {
        if (orbitMset != bestOrbitMset) {
          delete_node(node);
        }
      } else {
        if (orbitMset < minOrbitMset) {
          minOrbitMset = orbitMset;
          NodePtr node2 = node->prev;
          while (node2 != nullptr) {
            delete_node(node2);
            node2 = node2->prev;
          }
        } else {
          if (orbitMset > minOrbitMset) {
            delete_node(node);
          }
        }
      }
      node = next_node(node);
    }

    std::vector<int> globalOrbitCounts(orbmins.size(), 0);
    node = leftmost_node(depth);
    while (node != nullptr) {
      cands = DifferenceVect(ClosedInterval(0, m), node->selected);
      if (cands.size() > 1 && !IsTrivial(node->substab)) {
        cands = simpleOrbitReps(node->substab, cands);
      }
      /*
        # These index the children of node that will
        # not be immediately deleted under rule C
      */
      for (auto & y : cands) {
        int x = node->image[y];
        int num = make_orbit(x);
        globalOrbitCounts[num]++;
      }
      node = next_node(node);
    }
    int globalBestOrbit = calculateBestOrbit(orbmins, globalOrbitCounts, orbsizes);
    upb = orbmins[globalBestOrbit];


    node = leftmost_node(depth);
    while (node != nullptr) {
      cands = DifferenceVect(ClosedInterval(0,m), node->selected);
      if (cands.size() > 1 && !IsTrivial(node->substab)) {
        cands = simpleOrbitReps(node->substab, cands);
      }
      /*
        # These index the children of node that will
        # not be immediately deleted under rule C
      */
      node->validkids.clear();
      for (auto & y : cands) {
        int x = node->image[y];
        int num = orbnums[x];
        if (num == -1) {
          /*
            # Need a new orbit. Also require the smallest point
            # as the rep.
            #
            #
            # If there is no prospect of the new orbit being
            # better than the current best then go on to the next candidate
          */
          num = make_orbit(x);
          int rep = orbmins[num];
          if (rep < upb) {
            upb = rep;
            NodePtr node2 = node->prev;
            while (node2 != nullptr) {
              delete_node(node2);
              node2 = node2->prev;
            }
            node->validkids.clear();
            node->validkids.push_back(y);
          }
        } else {
          int rep = orbmins[num];
          if (rep == upb) {
            node->validkids.push_back(y);
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
    ChangeStabChain(s, {upb}, false);
    bool do_continue = false;
    if (s->orbit.size() == 1) {
      /*
        # In this case nothing much can happen. Each surviving node will have exactly one child
        # and none of the imsets will change
        # so we mutate the nodes in-place
      */
      node = leftmost_node(depth);
      while (node != nullptr) {
        if (!IsBound(node->selectedbaselength)) {
          node->selectedbaselength = node->selected.size();
        }
        node->selected.push_back(node->validkids[0]);
        node = next_node(node);
      }
      s = s->stabilizer;
      if (int(leftmost_node(depth+1)->selected.size()) != m) {
        do_continue = true;
      }
    }

    if (!do_continue) {
      node = leftmost_node(depth);
      NodePtr prevnode = nullptr;
      int nodect = 0;
      while (node != nullptr) {
        node->IsBoundChildren = false;
        node->children.clear();
        for (auto & x : node->validkids) {
          Node newnode_v;
          std::vector<int> selected = node->selected;
          selected.push_back(x);
          newnode_v.selected = selected;
          newnode_v.substab = Stabilizer_OnPoints<Telt,Tint>(node->substab, x);
          newnode_v.parent = node;
          newnode_v.childno = node->children.size();
          newnode_v.next = nullptr;
          newnode_v.prev = prevnode;
          newnode_v.deleted = false;
          NodePtr newnode = std::make_shared<Node>(newnode_v);
          nodect = nodect + 1;
          if (prevnode != nullptr) {
            prevnode->next = newnode;
          }
          prevnode = newnode;
          node->children.push_back(newnode);

          std::vector<int> image = node->image;
          if (image[x] != upb) {
            while (true) {
              image = OnTuples(image, s->transversal[image[x]]);
              if (image[x] == upb)
                break;
            }
            newnode->image = image;
            newnode->imset = Set(image);
          } else {
            newnode->image = image;
            newnode->imset = node->imset;
          }
        }
        node = next_node(node);
      }

      s = s->stabilizer;
      if (int(leftmost_node(depth+1)->selected.size()) == m) {
        break;
      }
    }
  }
  return {leftmost_node(depth+1)->image, l};
}



template<typename Telt, typename Tint>
Face CanonicalImage(StabChain<Telt> const& g, Face const& set)
{
  std::vector<int> set_i;
  boost::dynamic_bitset<>::size_type aRow = set.find_first();
  while (aRow != boost::dynamic_bitset<>::npos) {
    set_i.push_back(aRow);
    aRow = set.find_next(aRow);
  }
  StabChain<Telt> k = Stabilizer_OnSets<Telt,Tint>(g, set);
  Face ret(set.size());
  for (auto & eVal : NewSmallestImage<Telt,Tint>(g, set, k).set) {
    ret[eVal] = 1;
  }
  return ret;
}



}
#endif
