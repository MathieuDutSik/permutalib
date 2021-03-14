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

_IMAGES_NSI_HASH_LIMIT :=100;


_IMAGES_TIME_CLASSES := [];

_IMAGES_DeclareTimeClass := function(name)
    BindGlobal(name, Length(_IMAGES_TIME_CLASSES)+1);
    Add(_IMAGES_TIME_CLASSES,MakeImmutable(name));
end;


_IMAGES_DeclareTimeClass("pass1");
_IMAGES_DeclareTimeClass("pass2");
_IMAGES_DeclareTimeClass("pass3");
_IMAGES_DeclareTimeClass("shortcut");
_IMAGES_DeclareTimeClass("changeStabChain");
_IMAGES_DeclareTimeClass("orbit");
_IMAGES_DeclareTimeClass("skippedorbit");
_IMAGES_DeclareTimeClass("getcands");
_IMAGES_DeclareTimeClass("improve");
_IMAGES_DeclareTimeClass("check1");
_IMAGES_DeclareTimeClass("check2");
_IMAGES_DeclareTimeClass("prune");
_IMAGES_DeclareTimeClass("ShallowNode");
_IMAGES_DeclareTimeClass("DeepNode");
_IMAGES_DeclareTimeClass("FilterOrbCount");

_IMAGES_TIME_CLASSES := MakeImmutable(_IMAGES_TIME_CLASSES);

_IMAGES_nsi_stats := ListWithIdenticalEntries(Length(_IMAGES_TIME_CLASSES),0);

_IMAGES_DO_TIMING := true;
if IsBound( MakeThreadLocal ) then
    MakeThreadLocal("_IMAGES_DO_TIMING");
fi;

if _IMAGES_DO_TIMING then
    _IMAGES_StartTimer := function(cat)
        _IMAGES_nsi_stats[cat] := _IMAGES_nsi_stats[cat] - Runtime();
    end;

    _IMAGES_StopTimer := function(cat)
        _IMAGES_nsi_stats[cat] := _IMAGES_nsi_stats[cat] + Runtime();
    end;

    _IMAGES_IncCount := function(cat)
        _IMAGES_nsi_stats[cat] := _IMAGES_nsi_stats[cat] + 1;
    end;

    _IMAGES_ResetStats := function()
        _IMAGES_nsi_stats := ListWithIdenticalEntries(Length(_IMAGES_TIME_CLASSES),0);
    end;

    _IMAGES_ResetStats();
    if IsBound( MakeThreadLocal ) then
        MakeThreadLocal("_IMAGES_nsi_stats");
    fi;

    _IMAGES_GetStats := function()
        local   r,  c;
        r := rec();
        for c in _IMAGES_TIME_CLASSES do
            r.(c) := _IMAGES_nsi_stats[ValueGlobal(c)];
        od;
        return r;
    end;

else
    _IMAGES_StartTimer := function(cat)
        return;
    end;

    _IMAGES_StopTimer := function(cat)
        return;
    end;

    _IMAGES_IncCount := function(cat)
        return;
    end;

    _IMAGES_ResetStats := function()
        return;
    end;

    _IMAGES_GetStats := function()
        return fail;
    end;

fi;

_IMAGES_Get_Hash := function(m)
    local jenkins_hash;
    if IsBoundGlobal("JENKINS_HASH") then
        jenkins_hash := ValueGlobal("JENKINS_HASH");
         return s->jenkins_hash(s,GAPInfo.BytesPerVariable*m+GAPInfo.BytesPerVariable);
     else
       return s->HashKeyBag(s,57,0,GAPInfo.BytesPerVariable*m+GAPInfo.BytesPerVariable);
    fi;
end;


# GAP dictionaries don't (currently) provide a way of getting the values
# stored in them, so here we cache them seperately
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


if not IsBound(InfoNSI) then
    DeclareInfoClass("InfoNSI");
fi;

_IMAGES_RATIO := function(selector)
    return function(orbmins, orbitCounts, orbsizes)
        local index, result, i, ret;
        Print("GAP Beginning of calculateBestOrbit\n");
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
    Print("GAP NewSmallestImage : beginning\n");
    if IsString(config_option) then
        config_option := ValueGlobal(config_option);
    fi;
#    Print("GAP config_option.branch=", config_option.branch, "\n");
#    Print("GAP config_option=", config_option, "\n");
#    Print("GAP early_exit=", early_exit, "\n");

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
                # This space intensionally blank
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
#    Print("GAP 1 : peminv=", savedArgs.perminv, "\n");

    if disableStabilizerCheck_in = true then
        config.tryImproveStabilizer := false;
    fi;

    ## Node exploration functions
    leftmost_node := function(depth)
        local   n,  i;
        Print("GAP Beginning of leftmost_node\n");
        n := root;
        while Length(n.selected) < depth -1 do
            n := n.children[1];
        od;
        return n;
    end;

    next_node := function(node)
        local   n;
        Print("GAP Beginning of next_node\n");
        n := node;
        repeat
            n := n.next;
        until n = fail or not n.deleted;
        return n;
    end;

    # Delete a node, and recursively deleting all it's children.
    delete_node := function(node)
        local   i;
        Print("GAP Beginning of delete_node\n");
        if node.deleted then
            return;
        fi;
        Info(InfoNSI,3,"Deleting ",node.selected);
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
        Print("GAP Beginning of clean_subtree\n");
        Info(InfoNSI,3,"Cleaning at ",node.selected);
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
                Info(InfoNSI, 5, "Removing ", c, " because ", x);
                Add(bad,c);
            else
                q := [x];
                gens := GeneratorsOfGroup(node.substab);
                olen := 1;
                Info(InfoNSI, 5, "Keeping ", c, " because ", x);
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

    # Add a new stabilizer element, mapping node1 to node2, and then call
    # clean_subtree to remove any new subtrees.
    handle_new_stabilizer_element := function(node1,node2)
        local   perm1,  i;
        # so node1 and node2 represnet group elements that map set to the same
        # place in two different ways
        perm1 := PermListList(node1.image, node2.image);
        Info(InfoNSI, 2, "Can map ",node1.image, " to ", node2.image, " : ", perm1);
        Assert(1, not perm1 in l);
        l := ClosureGroup(l,perm1);
        Info(InfoNSI,2,"Found new stabilizer element. Stab now ",Size(l));
        root.substab := l;
        clean_subtree(root);
    end;

    # Given a group 'gp' and a set 'set', find orbit representatives
    # of 'set' in 'gp' simply.
    simpleOrbitReps := function(gp,set)
        local   m,  n,  b,  seed,  reps,  gens,  q,  pt,  gen,  im;
        Print("GAP Beginning of simpleOrbitReps\n");
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

    # Make orbit of x, updating orbnums, orbmins and orbsizes as approriate.
    make_orbit := function(x)
        local   q,  rep,  num,  pt,  gen,  img;
        Print("GAP Beginning of make_orbit x=", x, " orbnums[x]=", orbnums[x], "\n");
        if orbnums[x] <> -1 then
            return orbnums[x];
        fi;
        Print("GAP After check\n");
        q := [x];
        rep := x;
        num := Length(orbmins)+1;
        orbnums[x] := num;
        Print("GAP 1 : x=", x, " Assign orbnums[x]=", orbnums[x], "\n");
        for pt in q do
            for gen in gens do
                img := pt^gen;
                if orbnums[img] = -1 then
                    orbnums[img] := num;
                    Print("GAP 2 : img=", img, " Assign orbnums[img]=", orbnums[img], "\n");
                    Add(q,img);
                    if img < rep then
                        rep := img;
                    fi;
                fi;
            od;
        od;
        Print("GAP make_orbit rep=", rep, " |q|=", Length(q), " num=", num, "\n");
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
#    Print("GAP NSI k=", k, " set=", set, " l=", l, "\n");
#    Print("GAP |k|=", Order(k), "\n");

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
        Print("CPP depth=", depth, "\n");
        Info(InfoNSI, 3, "Stabilizer is :", s.generators);
        gens := s.generators;
        orbnums := ListWithIdenticalEntries(n,-1);
        orbmins := [];
        orbsizes := [];
        upb := config.initial_upb;
        imsets := [];
        imsetnodes := [];
        #
        # At this point, all bottom nodes are blue
        # first pass creates appropriate set of virtual red nodes
        #
        _IMAGES_StartTimer(pass1);

        if IsBound(config.findBestOrbMset) then
            countOrbDict := _countingDict([1,2,3]);
            node := leftmost_node(depth);
            while node <> fail do
                _IMAGES_StartTimer(getcands);
                cands := Difference([1..m],skip_func(node.selected));

                _IMAGES_StopTimer(getcands);
                orbitMset := [];
                for y in cands do
                    _IMAGES_IncCount(check1);
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
                _IMAGES_StartTimer(getcands);
                Print("GAP m=", m, " node.selected=", skip_func(node.selected), "\n");
                cands := Difference([1..m],skip_func(node.selected));
                Print("GAP 1 : cands=", List(cands, x->x), "\n");

                _IMAGES_StopTimer(getcands);
                orbitMset := [];
                for y in cands do
                    _IMAGES_IncCount(check1);
                    x := node.image[y];
                    num := make_orbit(x);
                    Print("GAP x=", x, " num=", num, "\n");
                    Add(orbitMset, orbmins[num]);
                od;
                Print("GAP bef orbitMset=", orbitMset, "\n");
                Sort(orbitMset);
                Print("GAP aft orbitMset=", orbitMset, "\n");
                Info(InfoNSI, 5, "Considering: ", orbitMset, "::",node.selected);
                if IsBound(bestOrbitMset) then
                    if orbitMset <> bestOrbitMset then
                        delete_node(node);
                    fi;
                else
                    if orbitMset < minOrbitMset then
                        Info(InfoNSI, 4, "New min: ", orbitMset);
                        Print("GAP orbitMset comparison case 1\n");
                        minOrbitMset := orbitMset;
                        node2 := node.prev;
                        while node2 <> fail do
                            Info(InfoNSI, 4, "Clean up old big set");
                            _IMAGES_IncCount(FilterOrbCount);
                            delete_node(node2);
                            node2 := node2.prev;
                        od;
                    elif orbitMset > minOrbitMset then
                        Info(InfoNSI, 4, "Too big!");
                        Print("GAP orbitMset comparison case 2\n");
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
                _IMAGES_StartTimer(getcands);
                cands := Difference([1..m],skip_func(node.selected));
                Print("GAP 2 : cands=", List(cands, x->x), "\n");
                if Length(cands) > 1 and not IsTrivial(node.substab) then
                    cands := simpleOrbitReps(node.substab,cands);
                fi;
                #
                # These index the children of node that will
                # not be immediately deleted under rule C
                #
                _IMAGES_StopTimer(getcands);
                for y in cands do
                    _IMAGES_IncCount(check1);
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

            _IMAGES_StartTimer(getcands);
            cands := Difference([1..m],skip_func(node.selected));
            Print("GAP 3 : cands=", List(cands, x->x), "\n");
            if Length(cands) > 1 and not IsTrivial(node.substab) then
                cands := simpleOrbitReps(node.substab,cands);
            fi;
            #
            # These index the children of node that will
            # not be immediately deleted under rule C
            #
            _IMAGES_StopTimer(getcands);
            node.validkids := [];
            for y in cands do
                _IMAGES_IncCount(check1);
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
                        _IMAGES_IncCount(skippedorbit);
                        continue;
                    fi;
                    _IMAGES_StartTimer(orbit);
                    num := make_orbit(x);
                    _IMAGES_StopTimer(orbit);
                    rep := config.getQuality(num);
                    if rep < upb then
                        _IMAGES_StartTimer(improve);
                        ### CAJ - Support bailing out early when a smaller
                        # set is found
                        if early_exit[1] and rep < early_exit[2][depth] then
#                            Print("GAP 2 : peminv=", savedArgs.perminv, "\n");
                            return [MinImage.Smaller, l^(savedArgs.perminv)];
                        fi;
                        ### END of bailing out early
                        upb := rep;
                        node2 := node.prev;
                        while node2 <> fail do
                            delete_node(node2);
                            node2 := node2.prev;
                        od;
                        _IMAGES_IncCount(ShallowNode);
                        node.validkids := [y];
                        Info(InfoNSI,3,"Best down to ",upb);
                        _IMAGES_StopTimer(improve);
                    fi;
                else
                    _IMAGES_IncCount(check2);
                    rep := config.getQuality(num);
                    if rep = upb then
                        _IMAGES_IncCount(ShallowNode);
                        Add(node.validkids,y);
                    fi;
                fi;
            od;
            if node.validkids = [] then
                _IMAGES_StartTimer(prune);
                delete_node(node);
                _IMAGES_StopTimer(prune);
            fi;
            node := next_node(node);
        od;
        ### CAJ - Support bailing out early when a larger set is found
        if early_exit[1] and upb > early_exit[2][depth] then
#            Print("GAP 3 : peminv=", savedArgs.perminv, "\n");
            return [MinImage.Larger, l^(savedArgs.perminv)];
        fi;
        ###
        Info(InfoNSI,2,"Layer ",depth," pass 1 complete. Best is ",upb);
        _IMAGES_StopTimer(pass1);
        #
        # Second pass. Actually make all the red nodes and turn them blue
        #
        lastupb := upb;
        Info(InfoNSI, 2, "Branch on ", upb);
        _IMAGES_StartTimer(changeStabChain);
        Print("GAP Before ChangeStabChain\n");
        ChangeStabChain(s,[config.getBasePoint(upb)],false);
        Print("GAP After ChangeStabChain\n");
        _IMAGES_StopTimer(changeStabChain);
        if Length(s.orbit) = 1 then
            #
            # In this case nothing much can happen. Each surviving node will have exactly one child
            # and none of the imsets will change
            # so we mutate the nodes in-place
            #
            _IMAGES_StartTimer(shortcut);
            node := leftmost_node(depth);
            while node <> fail do
                if not IsBound(node.selectedbaselength) then
                    node.selectedbaselength := Length(node.selected);
                fi;
                Assert(1, Length(node.validkids)=1);
                Add(node.selected, node.validkids[1]);
                Print("GAP Now node.selected=", node.selected, "\n");
                node := next_node(node);
            od;
            Info(InfoNSI,2,"Nothing can happen, short-cutting");
            s := s.stabilizer;
            _IMAGES_StopTimer(shortcut);
            if Size(skip_func(leftmost_node(depth+1).selected)) = m then
                Info(InfoNSI,2,"Skip would skip all remaining points");
                break;
            fi;

            continue; # to the next depth
        fi;
        _IMAGES_StartTimer(pass2);
        node := leftmost_node(depth);
        prevnode := fail;
        nodect := 0;
        changed := false;
        while node <> fail do
            node.children := [];
            Print("GAP node.validkids=", node.validkids, "\n");
            for x in node.validkids do
                _IMAGES_IncCount(DeepNode);
                newnode := rec( selected := Concatenation(node.selected,[x]),
                                substab := Stabilizer(node.substab,x),
                                parent := node,
                                childno := Length(node.children)+1,
                                next := fail,
                                prev := prevnode,
                                deleted := false);
                Print("GAP newnode.selected=", newnode.selected, "\n");
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
        _IMAGES_StopTimer(pass2);
        Info(InfoNSI,2,"Layer ",depth," pass 2 complete. ",nodect," new nodes");
        #
        # Third pass detect stabilizer elements
        #

        _IMAGES_StartTimer(pass3);
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
#                    if hash(node.imset) in seen then
#                        Error("");
#                    fi;
#                    AddSet(seen, hash(node.imset));
                    fi;
                    node := next_node(node);
                od;
                Info(InfoNSI,2,"Layer ",depth," pass 3 complete. Used hash table");
                s := s.stabilizer;
                if Length(s.generators) = 0 then
                    Info(InfoNSI,2,"Run out of group, return best image");
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
                    _IMAGES_StopTimer(pass3);
#                    Print("GAP 4 : peminv=", savedArgs.perminv, "\n");
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
                Info(InfoNSI,2,"Layer ",depth," pass 3 complete. ",Length(imsets)," images");
                s := s.stabilizer;
                if Length(s.generators) = 0 then
                    Info(InfoNSI,2,"Run out of group, return best image");
                    _IMAGES_StopTimer(pass3);
#                    Print("GAP 5 : peminv=", savedArgs.perminv, "\n");
                    return [OnTuples(imsetnodes[1].image,savedArgs.perminv),l^savedArgs.perminv];
                fi;
            fi;
        else
            Print("GAP Before s:=s.stabilizer operation\n");
            s := s.stabilizer;
        fi;
        _IMAGES_StopTimer(pass3);
        if Size(skip_func(leftmost_node(depth+1).selected)) = m then
            Info(InfoNSI,2,"Skip would skip all remaining points");
            break;
        fi;
    od;
#    Print("GAP 6 : peminv=", savedArgs.perminv, "\n");
    return [OnTuples(leftmost_node(depth+1).image,savedArgs.perminv),l^savedArgs.perminv];
end;
