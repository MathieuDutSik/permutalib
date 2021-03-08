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

_IMAGES_Get_Hash := function(m)
    local jenkins_hash;
    if IsBoundGlobal("JENKINS_HASH") then
        jenkins_hash := ValueGlobal("JENKINS_HASH");
         return s->jenkins_hash(s,GAPInfo.BytesPerVariable*m+GAPInfo.BytesPerVariable);
     else
       return s->HashKeyBag(s,57,0,GAPInfo.BytesPerVariable*m+GAPInfo.BytesPerVariable);
    fi;
end;

/*
  GAP dictionaries don't (currently) provide a way of getting the values
  stored in them, so here we cache them seperately
*/
_countingDict := function(dictexample)
    local data;
    data := rec(
        d := NewDictionary(dictexample, true),
        l := []
    );

    return rec(
        add := function(list)
            local val;
            val := LookupDictionary(data.d, list);
            if val = fail then
                val := 0;
                Add(data.l, list);
            fi;
            val := val + 1;
            AddDictionary(data.d, list, val);
        end,

        findElement := function(comp)
            local smallval, smalllist, val, i;
            smalllist := data.l[1];
            smallval := LookupDictionary(data.d, smalllist);
            for i in data.l do
                val := LookupDictionary(data.d, i);
                if comp(val, smallval) or (val = smallval and i < smalllist) then
                    smallval := val;
                    smalllist := i;
                fi;
            od;
            return smalllist;
        end,

        dump := function() return data; end
        );
end;

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

_IMAGES_RARE_RATIO_ORBIT := _IMAGES_RATIO(
    function(i, orbmins, orbitCounts, orbsizes)
        return (Log2(Float(orbitCounts[i])))/orbsizes[i];
    end
);

_IMAGES_COMMON_RATIO_ORBIT := _IMAGES_RATIO(
    function(i, orbmins, orbitCounts, orbsizes)
        return -(Log2(Float(orbitCounts[i])))/orbsizes[i];
    end
);

_IMAGES_RARE_RATIO_ORBIT_FIX := _IMAGES_RATIO(
    function(i, orbmins, orbitCounts, orbsizes)
        if(orbsizes[i]) = 1 then return Float(-(2^32)+orbitCounts[i]); fi;
        return (Log2(Float(orbitCounts[i])))/orbsizes[i];
    end
);

_IMAGES_COMMON_RATIO_ORBIT_FIX := _IMAGES_RATIO(
    function(i, orbmins, orbitCounts, orbsizes)
        if(orbsizes[i]) = 1 then return Float(-(2^32)-orbitCounts[i]); fi;
        return -(Log2(Float(orbitCounts[i])))/orbsizes[i];
    end
);

_IMAGES_RARE_ORBIT := _IMAGES_RATIO(
    function(i, orbmins, orbitCounts, orbsizes)
        return orbitCounts[i];
    end
);

_IMAGES_COMMON_ORBIT := _IMAGES_RATIO(
    function(i, orbmins, orbitCounts, orbsizes)
        return -orbitCounts[i];
    end
);


_NewSmallestImage := function(g,set,k,skip_func, early_exit, disableStabilizerCheck_in, config_option)
    local   leftmost_node,  next_node,  delete_node,  delete_nodes,
            clean_subtree,  handle_new_stabilizer_element,
            simpleOrbitReps,  make_orbit,  n,  s,  l,  m,  hash,
            lastupb,  root,  depth,  gens,  orbnums,  orbmins,
            orbsizes,  upb,  imsets,  imsetnodes,  node,  cands,  y,
            x,  num,  rep,  node2,  prevnode,  nodect,  changed,
            newnode,  image,  dict,  seen,  he,  bestim,  bestnode,
            imset,  p,
            config,
            globalOrbitCounts, globalBestOrbit, minOrbitMset, orbitMset,
            savedArgs,
            countOrbDict,
            bestOrbitMset
            ;

    if IsString(config_option) then
        config_option := ValueGlobal(config_option);
    fi;

    if config_option.branch = "static" then
            savedArgs := rec( config_option := config_option, g := g, k := k, set := set );
        if config_option.order = "MinOrbit" then
            savedArgs.perm := MinOrbitPerm(g);
        elif config_option.order = "MaxOrbit" then
            savedArgs.perm := MaxOrbitPerm(g);
        else
            ErrorNoReturn("Invalid 'order' when branch = 'static' in CanonicalImage");
        fi;
        savedArgs.perminv := savedArgs.perm^-1;
        g := g^savedArgs.perm;
        k := k^savedArgs.perm;
        set := OnTuples(set, savedArgs.perm);
        config_option := rec(branch := "minimum");
    else
        savedArgs := rec(perminv := ());
    fi;

    if config_option.branch = "minimum" then
        config := rec(
                   skipNewOrbit := function() return (upb <= lastupb + 1); end,
                   getQuality := pt -> orbmins[pt],
                   getBasePoint := IdFunc,
                   initial_lastupb := 0,
                   initial_upb := infinity,
                   countRareOrbits := false,
                   tryImproveStabilizer := true,
                   preFilterByOrbMset := false,
               );
    elif config_option.branch = "dynamic" then
        config := rec(skipNewOrbit := ReturnFalse,
                      preFilterByOrbMset := false);
        if config_option.order in ["MinOrbit", "MaxOrbit", "SingleMaxOrbit"] then
            config.getBasePoint := pt->pt[2];
            config.initial_lastupb := [-infinity, -infinity];
            config.initial_upb := [infinity, infinity];
            config.countRareOrbits := false;
            config.tryImproveStabilizer := true;

            if config_option.order = "MinOrbit" then
                config.getQuality := pt -> [orbsizes[pt], orbmins[pt]];
            elif config_option.order = "MaxOrbit" then
                config.getQuality := pt -> [-orbsizes[pt], orbmins[pt]];
            elif config_option.order = "SingleMaxOrbit" then
                config.getQuality := function(pt)
                                    if orbsizes[pt] = 1 then
                                        return [-(2^64), orbmins[pt]];
                                    else
                                        return [-orbsizes[pt], orbmins[pt]];
                                    fi;
                                 end;
            else
                ErrorNoReturn("?");
            fi;
        elif config_option.order in ["RareOrbit", "CommonOrbit", "RareRatioOrbit", "CommonRatioOrbit",
                                     "RareRatioOrbitFix", "CommonRatioOrbitFix"] then
            config.getBasePoint := IdFunc;
            config.initial_lastupb := 0;
            config.initial_upb := infinity;
            config.countRareOrbits := true;
            config.tryImproveStabilizer := false;
            config.getQuality := pt -> orbmins[pt];
            if config_option.order = "RareOrbit" then
                config.calculateBestOrbit := _IMAGES_RARE_ORBIT;
            elif config_option.order = "CommonOrbit" then
                config.calculateBestOrbit := _IMAGES_COMMON_ORBIT;
            elif config_option.order = "RareRatioOrbit" then
                config.calculateBestOrbit := _IMAGES_RARE_RATIO_ORBIT;
            elif config_option.order = "RareRatioOrbitFix" then
                config.calculateBestOrbit := _IMAGES_RARE_RATIO_ORBIT_FIX;
            elif config_option.order = "CommonRatioOrbit" then
                config.calculateBestOrbit := _IMAGES_COMMON_RATIO_ORBIT;
            elif config_option.order = "CommonRatioOrbitFix" then
                config.calculateBestOrbit := _IMAGES_COMMON_RATIO_ORBIT_FIX;
            else
                ErrorNoReturn("?");
            fi;
        else
            ErrorNoReturn("Invalid ordering: ", config_option.order);
        fi;

        if IsBound(config_option.orbfilt) then
            if config_option.orbfilt = "Min" then
              // This space intensionally blank
            elif config_option.orbfilt = "Rare" then
                config.findBestOrbMset := function(x,y) return x < y; end;
            elif config_option.orbfilt = "Common" then
                config.findBestOrbMset := function(x,y) return x > y; end;
            else
                Error("Invalid 'orbfilt' option");
            fi;
            config.preFilterByOrbMset := true;
        fi;
    else
        ErrorNoReturn("'branch' must be minimum, static or dynamic");
    fi;

    if disableStabilizerCheck_in = true then
        config.tryImproveStabilizer := false;
    fi;
    struct Node {
      std::vector<int> selected;
      std::vector<int> image;
      std::vector<int> imset;
      int substab;
      bool deleted;
      std::shared_ptr<Node> next;
      std::shared_ptr<Node> prev;
      std::shared_ptr<Node> parent;
      // children
      bool IsBoundChildreM
      std::vector<std::shared_ptr<Node>> children;
    };
    using NodePtr = std::shared_ptr<Node>;

    // Node exploration functions
    auto leftmost_node =[&](int const& depth) -> NodePtr {
        local   n,  i;
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
    delete_node := function(node)
        local   i;
        if node.deleted then
            return;
        fi;
        if node.prev <> fail then
            node.prev.next := node.next;
        fi;
        if node.next <> fail then
            node.next.prev := node.prev;
        fi;
        node.deleted := true;
        if node.parent <> fail then
            Remove(node.parent.children, node.childno);
            if Length(node.parent.children) = 0 then
                delete_node(node.parent);
            else
                for i in [node.childno..Length(node.parent.children)] do
                    node.parent.children[i].childno := i;
                od;
            fi;
        fi;
        if IsBound(node.children) then
            delete_nodes(ShallowCopy(node.children));
        fi;
    end;
    delete_nodes := function(nodes)
        local   node;
        for node in nodes do
            delete_node(node);
        od;
    end;

    # Filter nodes by stabilizer group,
    # Updates the stabilizer group of the node,
    clean_subtree := function(node)
        local   bad,  seen,  c,  x,  q,  gens,  olen,  pt,  gen,  im;
        if not IsBound(node.children) then
            return;
        fi;
        bad := [];

        seen := BlistList([1..m],[]);
        for c in node.children do
            if IsBound(c.selectedbaselength) then
                x := c.selected[c.selectedbaselength];
            else
                x := c.selected[Length(c.selected)];
            fi;
            if seen[x] then
                Add(bad,c);
            else
                q := [x];
                gens := GeneratorsOfGroup(node.substab);
                olen := 1;
                seen[x] := true;
                for pt in q do
                    for gen in gens do
                        im := pt^gen;
                        if not seen[im] then
                            seen[im] := true;
                            Add(q,im);
                            olen := olen+1;
                        fi;
                    od;
                od;
                if olen < Size(node.substab)/Size(c.substab) then
                    c.substab := Stabilizer(node.substab,x);
                    clean_subtree(c);
                fi;
            fi;
        od;
        delete_nodes(bad);
    end;

    //Add a new stabilizer element, mapping node1 to node2, and then call
    // clean_subtree to remove any new subtrees.
    handle_new_stabilizer_element := function(node1,node2)
        local   perm1,  i;
        # so node1 and node2 represnet group elements that map set to the same
        # place in two different ways
        perm1 := PermListList(node1.image, node2.image);
        Assert(1, not perm1 in l);
        l := ClosureGroup(l,perm1);
        root.substab := l;
        clean_subtree(root);
    end;

    // Given a group 'gp' and a set 'set', find orbit representatives
    // of 'set' in 'gp' simply.
    simpleOrbitReps := function(gp,set)
        local   m,  n,  b,  seed,  reps,  gens,  q,  pt,  gen,  im;
        m := Length(set);
        n := set[m];
        b := BlistList([1..n],set);
        seed := set[1];
        reps := [];
        gens := GeneratorsOfGroup(gp);
        while seed <> fail and seed <= n do
            b[seed] := false;
            q := [seed];
            Add(reps,seed);
            for pt in q do
                for gen in gens do
                    im := pt^gen;
                    if b[im] then
                        b[im] := false;
                        Add(q,im);
                    fi;
                od;
            od;
            seed := Position(b,true,seed);
        od;
        return reps;
    end;

    // Make orbit of x, updating orbnums, orbmins and orbsizes as approriate.
    make_orbit := function(x)
        local   q,  rep,  num,  pt,  gen,  img;
        if orbnums[x] <> -1 then
            return orbnums[x];
        fi;
        q := [x];
        rep := x;
        num := Length(orbmins)+1;
        orbnums[x] := num;
        for pt in q do
            for gen in gens do
                img := pt^gen;
                if orbnums[img] = -1 then
                    orbnums[img] := num;
                    Add(q,img);
                    if img < rep then
                        rep := img;
                    fi;
                fi;
            od;
        od;
        Add(orbmins,rep);
        Add(orbsizes,Length(q));
        return num;
    end;

    if set = [] then
      return [ [], k^(savedArgs.perminv)];
    fi;

    n := Maximum(LargestMovedPoint(g), Maximum(set));
    s := StabChainMutable(g);
    l := Action(k,set);
    m := Length(set);
    hash := _IMAGES_Get_Hash(m);
    lastupb := config.initial_lastupb;
    root := rec(selected := [],
                image := set,
                imset := Immutable(Set(set)),
                substab := l,
                deleted := false,
                next := fail,
                prev := fail,
                parent := fail);
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

        if IsBound(config.findBestOrbMset) then
            countOrbDict := _countingDict([1,2,3]);
            node := leftmost_node(depth);
            while node <> fail do
                cands := Difference([1..m],skip_func(node.selected));

                orbitMset := [];
                for y in cands do
                    x := node.image[y];
                    num := make_orbit(x);
                    Add(orbitMset, orbmins[num]);
                od;
                Sort(orbitMset);
                countOrbDict.add(orbitMset);
                node := next_node(node);
            od;

            bestOrbitMset := countOrbDict.findElement(config.findBestOrbMset);
            Unbind(countOrbDict); # Free memory
        fi;

        if config.preFilterByOrbMset then
            minOrbitMset := [infinity];
            node := leftmost_node(depth);
            while node <> fail do
                cands := Difference([1..m],skip_func(node.selected));

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
        fi;

        if config.countRareOrbits then
            globalOrbitCounts := ListWithIdenticalEntries(Length(orbmins), 0) ;
            node := leftmost_node(depth);
            while node <> fail do
                cands := Difference([1..m],skip_func(node.selected));
                if Length(cands) > 1 and not IsTrivial(node.substab) then
                    cands := simpleOrbitReps(node.substab,cands);
                fi;
                #
                # These index the children of node that will
                # not be immediately deleted under rule C
                #
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
        fi;


        node := leftmost_node(depth);
        while node <> fail do

            cands := Difference([1..m],skip_func(node.selected));
            if Length(cands) > 1 and not IsTrivial(node.substab) then
                cands := simpleOrbitReps(node.substab,cands);
            fi;
            #
            # These index the children of node that will
            # not be immediately deleted under rule C
            #
            node.validkids := [];
            for y in cands do
                x := node.image[y];

                num := orbnums[x];
                if num = -1 then
                    #
                    # Need a new orbit. Also require the smallest point
                    # as the rep.
                    #
                    #
                    # If there is no prospect of the new orbit being
                    # better than the current best then go on to the next candidate
                    #
                    if config.skipNewOrbit() then
                        continue;
                    fi;
                    num := make_orbit(x);
                    rep := config.getQuality(num);
                    if rep < upb then
                        ### CAJ - Support bailing out early when a smaller
                        # set is found
                        if early_exit[1] and rep < early_exit[2][depth] then
                            return [MinImage.Smaller, l^(savedArgs.perminv)];
                        fi;
                        ### END of bailing out early
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
        ### CAJ - Support bailing out early when a larger set is found
        if early_exit[1] and upb > early_exit[2][depth] then
            return [MinImage.Larger, l^(savedArgs.perminv)];
        fi;
        /*
        # Second pass. Actually make all the red nodes and turn them blue
        */
        lastupb := upb;
        ChangeStabChain(s,[config.getBasePoint(upb)],false);
        if Length(s.orbit) = 1 then
            #
            # In this case nothing much can happen. Each surviving node will have exactly one child
            # and none of the imsets will change
            # so we mutate the nodes in-place
            #
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
            if Size(skip_func(leftmost_node(depth+1).selected)) = m then
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
                if image[x] <> config.getBasePoint(upb) then
                    repeat
                        image := OnTuples(image, s.transversal[image[x]]);
                    until image[x] = config.getBasePoint(upb);
                    newnode.image := image;
                    newnode.imset := Set(image);
                    MakeImmutable(newnode.imset);
                    changed := true;
                else
                    newnode.image := image;
                    newnode.imset := node.imset;
                fi;
#                Print("Made a node ",newnode.selected, " ",newnode.image,"\n");
            od;
            node := next_node(node);
        od;
        /*
        # Third pass detect stabilizer elements
        */

        if  changed and config.tryImproveStabilizer then
            node := leftmost_node(depth+1);
            if nodect > _IMAGES_NSI_HASH_LIMIT then
                dict := SparseHashTable(hash);
                seen := [];
                while node <> fail do
                    he := GetHashEntry(dict,node.imset);
                    if  fail <> he then
                        handle_new_stabilizer_element(node, he);
                    else
                        AddHashEntry(dict, node.imset, node);
                    fi;
                    node := next_node(node);
                od;
                s := s.stabilizer;
                if Length(s.generators) = 0 then
                    node := leftmost_node(depth+1);
                    bestim := node.imset;
                    bestnode := node;
                    node := next_node(node);
                    while node <> fail do
                        if node.imset < bestim then
                            bestim := node.imset;
                            bestnode := node;
                        fi;
                        node := next_node(node);
                    od;
                    return [OnTuples(bestnode.image,savedArgs.perminv),l^savedArgs.perminv];
                fi;
            else
                while node <> fail do
                    imset := node.imset;
                    p := PositionSorted(imsets, imset);
                    if p <= Length(imsets) and imsets[p] = imset then
                        handle_new_stabilizer_element(node, imsetnodes[p]);
                    else
                        Add(imsets,imset,p);
                        Add(imsetnodes,node,p);
                    fi;
                    node := next_node(node);
                od;
                s := s.stabilizer;
                if Length(s.generators) = 0 then
                    return [OnTuples(imsetnodes[1].image,savedArgs.perminv),l^savedArgs.perminv];
                fi;
            fi;
        else
            s := s.stabilizer;
        fi;
        if Size(skip_func(leftmost_node(depth+1).selected)) = m then
            break;
        fi;
    od;
    return [OnTuples(leftmost_node(depth+1).image,savedArgs.perminv),l^savedArgs.perminv];
end;
