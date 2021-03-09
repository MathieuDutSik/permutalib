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


_IMAGES_NSI_HASH_LIMIT :=100;

_IMAGES_RATIO := function(selector)
    return function(orbmins, orbitCounts, orbsizes)
        local index, result, i, ret;
        index := 1;
        result := [selector(1, orbmins, orbitCounts, orbsizes), orbmins[1]];
        for i in [2..Length(orbmins)] do
            ret := [selector(i, orbmins, orbitCounts, orbsizes), orbmins[i]];
            if (orbitCounts[index] = 0) or (ret < result and orbitCounts[i] <> 0) then
                index := i;
                result := ret;
            fi;
        od;
        return index;
    end;
end;

template<typename T>
void Remove(std::vector<T> & eV, int const& pos)
{
  std::erase(eV.begin() + pos);
}



_IMAGES_RARE_RATIO_ORBIT_FIX := _IMAGES_RATIO(
    function(i, orbmins, orbitCounts, orbsizes)
        if(orbsizes[i]) = 1 then return Float(-(2^32)+orbitCounts[i]); fi;
        return (Log2(Float(orbitCounts[i])))/orbsizes[i];
    end
);



/*
  Modification done:
  --- skip_fnuc eliminated as it is the identity in the case that interest us.
 */
template<typename Telt, typename Tint>
_NewSmallestImage := function(g, set, k, early_exit, config_option)
    local   leftmost_node,  next_node,  delete_node,  delete_nodes,
            clean_subtree,  handle_new_stabilizer_element,
            simpleOrbitReps,  make_orbit,  n,  s,  l,  m,  hash,
            root,  depth,  gens,  orbnums,  orbmins,
            orbsizes,  upb,  imsets,  imsetnodes,  node,  cands,  y,
            x,  num,  rep,  node2,  prevnode,  nodect,  changed,
            newnode,  image,  dict,  seen,  he,  bestim,  bestnode,
            imset,  p, config,
            globalOrbitCounts, globalBestOrbit, minOrbitMset, orbitMset,
            savedArgs, countOrbDict, bestOrbitMset;

    if IsString(config_option) then
        config_option := ValueGlobal(config_option);
    fi;

    savedArgs := rec(perminv := ());

    config := rec(skipNewOrbit := ReturnFalse);
    config.initial_upb := infinity;
    config.getQuality := pt -> orbmins[pt];
    config.calculateBestOrbit := _IMAGES_RARE_RATIO_ORBIT_FIX;

    struct Node {
      std::vector<int> selected;
      std::vector<int> image;
      std::vector<int> imset;
      int substab;
      bool deleted;
      std::shared_ptr<Node> next;
      std::shared_ptr<Node> prev;
      std::shared_ptr<Node> parent;
      int selectedbaselength;
      // children
      int childno;
      bool IsBoundChildren;
      std::vector<std::shared_ptr<Node>> children;
    };
    using NodePtr = std::shared_ptr<Node>;

    n := Maximum(LargestMovedPoint(g), Maximum(set));
    s := StabChainMutable(g);
    l := Action(k,set);
    m := Length(set);
    root := rec(selected := [],
                image := set,
                imset := Immutable(Set(set)),
                substab := l,
                deleted := false,
                next := fail,
                prev := fail,
                parent := fail);

    // Node exploration functions
    auto leftmost_node =[&](int const& depth) -> NodePtr {
        NodePtr n = root;
        while (n.selected.size() < depth - 1) {
            n = n->children[0];
        }
        return n;
    };
    auto next_node =[&](NodePtr const& node) -> NodePtr {
      NodePtr n = node;
      while (true) {
        n = n->next;
        if (n == nullptr || !n.deleted)
          break;
      }
      return n;
    };
    //Delete a node, and recursively deleting all it's children.
    std::function<void(NodePtr)> delete_node=[&](NodePtr & node) -> void {
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
        Remove(node.parent.children, node.childno);
        if (node->parent->children.size() == 0) {
          delete_node(node->parent);
        } else {
          for (int i=node->childno; i<node->parent->children.size(); i++) {
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
    std::function<void(NodePtr)> clean_subtree =[&](NodePtr & node) -> void {
      if (!node-IsBoundChildren)
        return;
      std::vectror<NodePtr> bad;

      std::vector<int> range= ClosedInterval(0,m);
      Face seen = BlistList(range,{});
      int x;
      for (auto & c : node.children) {
        if (c->selectedbaselength != -1) {
          x = c.selected[c.selectedbaselength];
        } else {
          x = c.selected[c->selected.size() - 1];
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
              int p = q[idx];
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
          mpz_class quot = Size<mpz_class>(node->substab) / Size<mpz_class>(c->substab);
          if (mpz_class(olen) < quot) {
            c->substab = Stabilize(node->substab,x);
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
      Telt perm1 = PermListList(node1.image, node2.image);
      l = ClosureGroup(l,perm1);
      root->substab = l;
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
        seed = PositionVect(b,true,seed);
      }
      return reps;
    };

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

    if (set.size() == 0) {
    }

    if set = [] then
      return [ [], k^(savedArgs.perminv)];
    fi;
    for depth in [1..m] do
        gens := s.generators;
        orbnums := ListWithIdenticalEntries(n,-1);
        orbmins := [];
        orbsizes := [];
        upb := config.initial_upb;
        imsets := [];
        imsetnodes := [];
        /*
        # At this point, all bottom nodes are blue
        # first pass creates appropriate set of virtual red nodes
        */


        minOrbitMset := [infinity];
        node := leftmost_node(depth);
        while node <> fail do
            cands := Difference([1..m],node.selected);

            orbitMset := [];
            for y in cands do
                x := node.image[y];
                num := make_orbit(x);
                Add(orbitMset, orbmins[num]);
            od;
            Sort(orbitMset);
            if IsBound(bestOrbitMset) then
                if orbitMset <> bestOrbitMset then
                    delete_node(node);
                fi;
            else
                if orbitMset < minOrbitMset then
                    minOrbitMset := orbitMset;
                    node2 := node.prev;
                    while node2 <> fail do
                        delete_node(node2);
                        node2 := node2.prev;
                    od;
                elif orbitMset > minOrbitMset then
                    delete_node(node);
                fi;
            fi;
            node := next_node(node);
        od;

        globalOrbitCounts := ListWithIdenticalEntries(Length(orbmins), 0) ;
        node := leftmost_node(depth);
        while node <> fail do
            cands := Difference([1..m],node.selected);
            if Length(cands) > 1 and not IsTrivial(node.substab) then
                cands := simpleOrbitReps(node.substab,cands);
            fi;
            /*
              # These index the children of node that will
              # not be immediately deleted under rule C
            */
            for y in cands do
                x := node.image[y];
                num := make_orbit(x);
                if IsBound(globalOrbitCounts[num]) then
                    globalOrbitCounts[num] := globalOrbitCounts[num] + 1;
                else
                    globalOrbitCounts[num] := 1;
                fi;
            od;
            node := next_node(node);
        od;
        globalBestOrbit := config.calculateBestOrbit(orbmins, globalOrbitCounts, orbsizes);
        upb := orbmins[globalBestOrbit];


        node := leftmost_node(depth);
        while node <> fail do

            cands := Difference([1..m],node.selected);
            if Length(cands) > 1 and not IsTrivial(node.substab) then
                cands := simpleOrbitReps(node.substab,cands);
            fi;
            /*
            # These index the children of node that will
            # not be immediately deleted under rule C
            */
            node.validkids := [];
            for y in cands do
                x := node.image[y];

                num := orbnums[x];
                if num = -1 then
                  /*
                    # Need a new orbit. Also require the smallest point
                    # as the rep.
                    #
                    #
                    # If there is no prospect of the new orbit being
                    # better than the current best then go on to the next candidate
                  */
                    if config.skipNewOrbit() then
                        continue;
                    fi;
                    num := make_orbit(x);
                    rep := config.getQuality(num);
                    if rep < upb then
                             // CAJ - Support bailing out early when a smaller
                             // set is found
                        if early_exit[1] and rep < early_exit[2][depth] then
                            return [MinImage.Smaller, l^(savedArgs.perminv)];
                        fi;
                             // END of bailing out early
                        upb := rep;
                        node2 := node.prev;
                        while node2 <> fail do
                            delete_node(node2);
                            node2 := node2.prev;
                        od;
                        node.validkids := [y];
                    fi;
                else
                    rep := config.getQuality(num);
                    if rep = upb then
                        Add(node.validkids,y);
                    fi;
                fi;
            od;
            if node.validkids = [] then
                delete_node(node);
            fi;
            node := next_node(node);
        od;
        // CAJ - Support bailing out early when a larger set is found
        if early_exit[1] and upb > early_exit[2][depth] then
            return [MinImage.Larger, l^(savedArgs.perminv)];
        fi;
        /*
        # Second pass. Actually make all the red nodes and turn them blue
        */
        ChangeStabChain(s, [upb], false);
        if Length(s.orbit) = 1 then
           /*
             # In this case nothing much can happen. Each surviving node will have exactly one child
             # and none of the imsets will change
             # so we mutate the nodes in-place
           */
            node := leftmost_node(depth);
            while node <> fail do
                if not IsBound(node.selectedbaselength) then
                    node.selectedbaselength := Length(node.selected);
                fi;
                Assert(1, Length(node.validkids)=1);
                Add(node.selected, node.validkids[1]);
                node := next_node(node);
            od;
            s := s.stabilizer;
            if Size(leftmost_node(depth+1).selected) = m then
                break;
            fi;

            continue; // to the next depth
        fi;
        node := leftmost_node(depth);
        prevnode := fail;
        nodect := 0;
        changed := false;
        while node <> fail do
            node.children := [];
            for x in node.validkids do
                newnode := rec( selected := Concatenation(node.selected,[x]),
                                substab := Stabilizer(node.substab,x),
                                parent := node,
                                childno := Length(node.children)+1,
                                next := fail,
                                prev := prevnode,
                                deleted := false);
                nodect := nodect+1;
                if prevnode <> fail then
                    prevnode.next := newnode;
                fi;
                prevnode := newnode;
                Add(node.children,newnode);

                image := node.image;
                if image[x] <> upb then
                    repeat
                        image := OnTuples(image, s.transversal[image[x]]);
                    until image[x] = upb;
                    newnode.image := image;
                    newnode.imset := Set(image);
                    MakeImmutable(newnode.imset);
                    changed := true;
                else
                    newnode.image := image;
                    newnode.imset := node.imset;
                fi;
            od;
            node := next_node(node);
        od;

        s := s.stabilizer;
        if Size(leftmost_node(depth+1).selected) = m then
            break;
        fi;
    od;
    return [OnTuples(leftmost_node(depth+1).image,savedArgs.perminv),l^savedArgs.perminv];
end;
