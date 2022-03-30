#############################################################################
##
#W  stbc.gi                     GAP library                    Heiko Theißen
#W                                                               Ákos Seress
##
##
#Y  Copyright (C)  1996,  Lehrstuhl D für Mathematik,  RWTH Aachen,  Germany
#Y  (C) 1998 School Math and Comp. Sci., University of St Andrews, Scotland
#Y  Copyright (C) 2002 The GAP Group
##

GetStringExpressionOfStabChain:=function(eRec)
  local eStab, strRet;
  if IsInt(eRec) then
    return Concatenation("integer_", String(eRec));
  fi;
  if IsRecord(eRec) then
    eStab:=eRec;
    strRet:="record_";
    while(true)
    do
      if IsBound(eStab.orbit) then
        strRet:=Concatenation(strRet, "orbit_", String(eStab.orbit));
      fi;
      if IsBound(eStab.transversal) then
        strRet:=Concatenation(strRet, "_transversal_", String(eStab.transversal));
      fi;
      if IsBound(eStab.stabilizer) then
        eStab:=eStab.stabilizer;
      else
        break;
      fi;
    od;
    return strRet;
  fi;
  Print("Need to write further code for GetStringExpressionOfStabChain for handling another case\n");
  Print(NullMat(5));
end;


SortVector:=function(eV)
    local eVV;
    eVV:=eV;
    Sort(eVV);
    return eVV;
end;




GetListStabCommPartition:=function(ListS)
    local len, Status, ListStr, i, LVal, j, labels_i, labels_j;
    len:=Length(ListS);
    Status:=ListWithIdenticalEntries(len,0);
    ListStr:=[];
    for i in [1..len]
    do
        if IsBound(ListS[i].labels) then
            labels_i:=ListS[i].labels;
        else
            labels_i:="unset";
        fi;
        if Status[i]=0 then
            LVal:=[];
            for j in [1..len]
            do
                if IsBound(ListS[j].labels) then
                    labels_j:=ListS[j].labels;
                else
                    labels_j:="unset";
                fi;
                if IsIdenticalObj(labels_i, labels_j) or i=j then
                    Add(LVal, j);
                    Status[j]:=1;
                fi;
            od;
            Add(ListStr, LVal);
        fi;
    od;
    return ListStr;
end;


PrintStabChain:=function(eRec)
  local eStab, iLev;
  #
  Print("GAP Reference Partition=", GetListStabCommPartition(ListStabChain(eRec)), "\n");
  eStab:=eRec;
  iLev:=0;
  while(true)
  do
    Print("GAP iLev=", iLev, "\n");
    if IsBound(eStab.orbit) then
      Print("GAP   orbit=", eStab.orbit, "\n");
    else
      Print("GAP   orbit=[  ]\n");
    fi;
    if IsBound(eStab.transversal) then
      Print("GAP   transversal=", eStab.transversal, "\n");
    else
      Print("GAP   transversal=[  ]\n");
    fi;
    Print("XXX ELIMINATE begin\n");
    if IsBound(eStab.cycles) then
        Print("GAP   cycles=", eStab.cycles, "\n");
    else
        Print("GAP   No cycles\n");
    fi;
    Print("XXX ELIMINATE end\n");
    if IsBound(eStab.stabilizer) then
      eStab:=eStab.stabilizer;
    else
      break;
    fi;
    iLev:=iLev+1;
  od;
end;


PrintListStabCommPartition:=function(mesg, ListS)
    local ListStr;
    ListStr:=GetListStabCommPartition(ListS);
    Print(mesg, " ListStabCommPartition=", ListStr, "\n");
end;



GetCompleteListLabels:=function(ListS)
    return List(ListS, x->x.labels);
end;





PrintTopOrbit:=function(eRec)
  if IsBound(eRec.orbit) then
    return eRec.orbit;
  else
    return "[  ]";
  fi;
end;




PrintStabChainTransversals:=function(S)
  local Swork, iLevel;
  Swork := S;
  iLevel:=0;
  while(true)
  do
    if IsBound(Swork.transversal) then
      Print("GAP i=", iLevel, " transversal=", Swork.transversal, "\n");
    else
      Print("GAP i=", iLevel, " transversal=[  ]\n");
    fi;
    iLevel:=iLevel+1;
    if IsBound(Swork.stabilizer) then
      Swork:=Swork.stabilizer;
    else
      break;
    fi;
  od;
end;



PrintStabChainOrbits:=function(S)
  local Swork, iLevel;
  Swork := S;
  iLevel:=0;
  while(true)
  do
    if IsBound(Swork.orbit) then
      Print("GAP i=", iLevel, " orbit=", Swork.orbit, "\n");
    else
      Print("GAP i=", iLevel, " orbit=[  ]\n");
    fi;
    iLevel:=iLevel+1;
    if IsBound(Swork.stabilizer) then
      Swork:=Swork.stabilizer;
    else
      break;
    fi;
  od;
end;




#############################################################################
##
#F  StabChain( <G>, <options> ) . . . . . . . . . . . . make stabilizer chain
##
InstallGlobalFunction( StabChain, function( arg )
    if Length( arg ) = 1  then
        return StabChainImmutable( arg[ 1 ] );
    else
        return Immutable( StabChainOp( arg[ 1 ], arg[ 2 ] ) );
    fi;
end );

InstallMethod( StabChainImmutable,"use StabChainMutable",
  true, [ IsObject ], 0, StabChainMutable );

InstallMethod( StabChainMutable,"call StabChainOp", true, [ IsGroup ], 0,
    function(G)
    Print("GAP Call to StabChainOp (call StabChainOp)\n");
    return StabChainOp( G, rec(  ) );
end );

InstallOtherMethod( StabChainOp,"with base", true, [ IsPermGroup,
        IsList and IsCyclotomicCollection ], 0,
    function( G, base )
    Print("GAP Call to StabChainOp (with base)\n");
    return StabChainOp( G, rec( base := base ) );
end );

InstallOtherMethod( StabChainOp,"empty base", true,
  [ IsPermGroup, IsList and IsEmpty ], 0,
    function( G, base )
    Print("GAP Call to StabChainOp (empty base)\n");
    return StabChainOp( G, rec( base := base ) );
end );

InstallMethod( StabChainOp,"trivial group",
  [ IsPermGroup and IsTrivial, IsRecord ],
    function( G, options )
    local   S,  T,  pnt;
    Print("GAP Call to StabChainOp (trivial group)\n");

#    Print(NullMat(67));
    
    S := EmptyStabChain( [  ], One( G ) );
    if     IsBound( options.base )
       and (        IsBound( options.reduced )
                and not options.reduced
             or     not IsBound( options.reduced )
                and not DefaultStabChainOptions.reduced )  then
        T := S;
        for pnt  in options.base  do
            InsertTrivialStabilizer( T, pnt );
            T := T.stabilizer;
        od;
    fi;
    return S;
end );

InstallMethod( StabChainOp,"group and option",
  [ IsPermGroup, IsRecord ],
    function( G, options )
    local   S,  T,  degree,  pcgs, UseNonPortedMethods;

    UseNonPortedMethods:=false;

    Print("GAP Call to StabChainOp (group and option)\n");
    # If a stabilizer chain <S> is already known, modify it.
    if HasStabChainMutable( G )  then
        Print("GAP Begin of HaseStabChainMutable=true section\n");
        S := StructuralCopy( StabChainMutable( G ) );
        if IsBound( options.base )  then
            if not IsBound( options.reduced )  then
                options.reduced := DefaultStabChainOptions.reduced;
            fi;
            if not ChangeStabChain( S, options.base, options.reduced )
               then
                return false;
            fi;
        elif IsBound( options.reduced )  and  options.reduced  then
            Print("GAP Before ReduceStabChain\n");
            ReduceStabChain( S );
        fi;
        Print("GAP End of HaseStabChainMutable=true section\n");
    # Otherwise construct a new GAP object <S>.
    else
        Print("GAP Case HaseStabChainMutable=false\n");
        CopyOptionsDefaults( G, options );

        # For solvable groups, use the pcgs algorithm.
        pcgs := [  ];
        if  UseNonPortedMethods and   options.tryPcgs and (not IsBound(options.base))
           and (# the group is know to be solvable
	     (HasIsSolvableGroup(G) and IsSolvableGroup(G))
		# or the degree is small and the group is not known to be
		# insolvable
	        or (Length(MovedPoints(G))<=100 and not
		    (HasIsSolvableGroup(G) and not IsSolvableGroup(G))
		    )) then
            S := EmptyStabChain( [  ], One( G ) );
            if IsBound( options.base )  then  S.base := options.base;
                                        else  S.base := [  ];          fi;
	    if HasPcgs(G) and IsBound(Pcgs(G)!.stabChain) then
	      # is there already a pcgs with a stabchain?
	      # the translation to a record is  necessary to be able to copy
	      # the stab chain.
	      pcgs:=rec(stabChain:=CopyStabChain(Pcgs(G)!.stabChain));
	    else
	      Print("GAP Before call to TryPcgsPermGroup\n");
	      pcgs := TryPcgsPermGroup( [ G, GroupStabChain( G, S, true ) ],
			      # get the series elementary abelian -- its much
			      # better
                            false, false, true );
	    fi;
        fi;
        if IsPcgs( pcgs ) and UseNonPortedMethods then
	  options.random := 1000;
	  S := pcgs!.stabChain;

	  if not HasPcgs(G) then
	    # remember the pcgs
	    SetPcgs(G,pcgs);
	    SetPcgsElementaryAbelianSeries(G,pcgs);
	    S := CopyStabChain(S); # keep the pcgs' pristine stabchain
	    if IsBound(options.base) then
              Print("Before ChangeStabChain, step 2\n");
	      ChangeStabChain( S, options.base, options.reduced );
	    fi;
	  fi;
	elif UseNonPortedMethods and IsRecord(pcgs) then
	  S:=pcgs.stabChain;
	  if IsBound(options.base) then
            Print("Before ChangeStabChain, step 3\n");
	    ChangeStabChain( S, options.base, options.reduced );
	  fi;
        else
            degree := LargestMovedPoint( G );
            Print("DEBUG degree=", degree, "\n");
            if degree > 100  then

                # random Schreier-Sims
		Print("GAP Before call to StabChainRandomPermGroup\n");
                S := StabChainRandomPermGroup(
                         ShallowCopy( GeneratorsOfGroup( G ) ), One( G ),
                         options );

            else

                # ordinary Schreier Sims
                S := EmptyStabChain( [  ], One( G ) );
                Unbind( S.generators );
                if not IsTrivial( G )  then
                    if not IsBound( options.base )  then
                        options.base := [  ];
                    fi;
                    Print("GAP options=", options, "\n");
                    S.cycles := [  ];
                    Print("GAP Before call to StabChainStrong\n");
                    StabChainStrong( S, GeneratorsOfGroup( G ), options );
		    Print("GAP After main call to StabChainStrong\n");
		    PrintStabChain(S);
                    T := S;
                    while IsBound( T.stabilizer )  do
                        T.generators := T.labels{ T.genlabels };
                        Unbind( T.cycles );
                        T := T.stabilizer;
                    od;
                    T.generators := T.labels{ T.genlabels };
                    Unbind( T.cycles );
		else
		  S.generators:=[];
                fi;

            fi; # random / deterministic
        fi;

        # Now extend <S>, if desired.
        if not options.reduced  and  IsBound( options.base )  then
            ExtendStabChain( S, options.base );
        fi;

    fi;

    # if the parent is random, this group should be also
    # at base change or strong gens constr, may be no info about random
    if IsBound( options.random ) then
        if IsBound( StabChainOptions( Parent( G ) ).random )  then
            options.random := Minimum( StabChainOptions( Parent( G ) ).random,
                                      options.random );
        fi;
        StabChainOptions( G ).random := options.random;
    fi;
    SetStabChainMutable( G, S );
    return S;
end );

#############################################################################
##
#F  TrimStabChain( <C>,<n> )
##
##
InstallGlobalFunction(TrimStabChain,function( C,n )
local i;
  # typically all permutations in a stabilizer chain are just links to the
  # `labels' component. Thus reducing here will make them all small.
  for i in C.labels do
    if IsInternalRep(i) then
      TRIM_PERM(i,n);
    fi;
  od;
end);

#############################################################################
##
#F  CopyStabChain( <C> )  . . . . . . . . . . . . . . . . . . . copy function
##
##  This function produces a memory-disjoint copy of a stabilizer chain, with
##  `labels'  components   possibly   shared by  several  levels,   but  with
##  superfluous labels  removed. An entry  in  `labels' is superfluous  if it
##  does not occur among  the  `genlabels' or `translabels'   on any  of  the
##  levels which share that `labels' component.
##
##  This is useful for  stabiliser sub-chains that  have been obtained as the
##  (iterated) `stabilizer' component of a bigger chain.
##
InstallGlobalFunction(CopyStabChain,function( C1 )
    local   C,Xlabels,  S,  len,  xlab,  need,  poss,  i;

    # To begin with, make a deep copy.
    C := StructuralCopy( C1 );

    # First pass: Collect the necessary genlabels.
    Xlabels := [  ];
    S := C;
    while IsBound( S.stabilizer )  do
        len := Length( S.labels );
        if len = 0  or  IsPerm( S.labels[ len ] )  then
            Add( S.labels, [ 1 ] );  len := len + 1;
            Add( Xlabels, S.labels );
        fi;
        UniteSet( S.labels[ len ], S.genlabels );
        UniteSet( S.labels[ len ], S.translabels );
        S := S.stabilizer;
    od;

    # Second pass: Find the new positions of the labels.
    for xlab  in Xlabels  do
        need := xlab[ Length( xlab ) ];

        # If all labels are needed, change nothing.
        if Length( need ) = Length( xlab ) - 1  then
            Unbind( xlab[ Length( xlab ) ] );

        else
            poss := [  ];
            for i  in [ 1 .. Length( need ) ]  do
                poss[ need[ i ] ] := i;
            od;
            Add( xlab, poss );
        fi;
    od;

    # Third pass: Update the genlabels and translabels.
    S := C;
    while IsBound( S.stabilizer )  do
        len := Length( S.labels );
        if len <> 0  and  not IsPerm( S.labels[ len ] )  then
            poss := S.labels[ len ];
            S.genlabels := poss{ S.genlabels };
            S.translabels{ S.orbit } := poss{ S.translabels{ S.orbit } };
        fi;
        S := S.stabilizer;
    od;

    # Fourth pass: Update the labels.
    for xlab  in Xlabels  do
        len := Length( xlab );
        if len <> 0  and  not IsPerm( xlab[ len ] )  then
            need := xlab[ Length( xlab ) - 1 ];
            xlab{ [ 1 .. Length( need ) ] } := xlab{ need };
            for i  in [ Length( need ) + 1 .. Length( xlab ) ]  do
                Unbind( xlab[ i ] );
            od;
        fi;
    od;

    return C;
end);

#############################################################################
##
#M  StabChainOptions( <G> ) . . . . . . . . . . . . . for a permutation group
##
InstallMethod( StabChainOptions, true, [ IsPermGroup ], 0,
    G -> rec(  ) );

#############################################################################
##
#V  DefaultStabChainOptions . . . . . .  options record for stabilizer chains
##
InstallValue( DefaultStabChainOptions,rec( reduced := true,
                                 random := 1000,
                                tryPcgs := true ));

#############################################################################
##
#F  CopyOptionsDefaults( <G>, <options> ) . . . . . . . copy options defaults
##
InstallGlobalFunction(CopyOptionsDefaults,function( G, options )
    local   P,  name;

    # See whether we know a base for <G>.
    if not IsBound( options.knownBase )  then
	if HasStabChainMutable(G) then
	  options.knownBase := BaseStabChain(StabChainMutable(G));
        elif   HasBaseOfGroup( G )  then
	  options.knownBase := BaseOfGroup( G );
        else
	  P := Parent( G );
	  while     not HasBaseOfGroup( P )
		and not IsIdenticalObj( P, Parent( P ) )  do
	    P := Parent( P );
	  od;
	  if HasStabChainMutable(P) then
	    options.knownBase := BaseStabChain(StabChainMutable(P));
	  elif HasBaseOfGroup( P )  then
	    options.knownBase := BaseOfGroup( P );
	  fi;
        fi;
    fi;

    # See whether we know the exact size.
    if not IsBound( options.size )  then
        if HasSize( G )  then
            options.size := Size( G );
        elif IsBound( StabChainOptions( G ).size )  then
            options.size := StabChainOptions( G ).size;
        fi;
    fi;

    # Copy the default values.
    for name  in RecNames( DefaultStabChainOptions )  do
        if not IsBound( options.( name ) )  then
            if IsBound( StabChainOptions( G ).( name ) )  then
                options.( name ) := StabChainOptions( G ).( name );
            else
                options.( name ) := DefaultStabChainOptions.( name );
            fi;
        fi;
    od;

    # In the case of random construction, see whether  we know an upper limit
    # for the size.
    if IsBound( options.size ) then
	options.limit := options.size;
    elif not IsBound( options.limit )  then
	if IsBound( StabChainOptions( G ).limit )  then
	    options.limit := StabChainOptions( G ).limit;
	else
	    P := Parent( G );
	    while     not HasSize( P )
		  and not IsBound( StabChainOptions( P ).limit )
		  and not IsIdenticalObj( P, Parent( P ) )  do
		P := Parent( P );
	    od;
	    if HasSize( P )  then
		options.limit := Size( P );
	    elif IsBound( StabChainOptions( P ).limit )  then
		options.limit := StabChainOptions( P ).limit;
	    fi;
	fi;
    fi;

end);

#############################################################################
##
#F  StabChainBaseStrongGenerators( <base>, <sgs>[, <one>] )
##
InstallGlobalFunction(StabChainBaseStrongGenerators,function(arg)
local   base,sgs,one,S,  T,  pnt;

    Print("DEBUG begin of StabChainBaseStrongGenerators\n");
    base:=arg[1];
    sgs:=arg[2];
    if Length(arg)=3 then
      one:=arg[3];
    else
      one:= One(arg[2][1]);
    fi;
    Print("DEBUG one=", one, " sgs=", sgs, "\n");
    S := EmptyStabChain( [  ], one );
    T := S;
    for pnt  in base  do
        InsertTrivialStabilizer( T, pnt );
        AddGeneratorsExtendSchreierTree( T, sgs );
        sgs := Filtered( sgs, g -> pnt ^ g = pnt );
        T := T.stabilizer;
    od;
    Print("DEBUG exiting of StabChainBaseStrongGenerators\n");
    return S;
end);

#############################################################################
##
#F  GroupStabChain( <arg> ) . . . . . . make (sub)group from stabilizer chain
##
InstallGlobalFunction(GroupStabChain,function( arg )
local   S,  G,  P,L;

    if Length( arg ) = 1  then
        S := arg[ 1 ];
	if not IsBound(S.generators) then
	  G := GroupByGenerators( [], S.identity );
	else
	  G := GroupByGenerators( S.generators, S.identity );
	fi;
    else
        P := arg[ 1 ];
        S := arg[ 2 ];
	if not IsBound(S.generators) then
	  L := [];
	else
	  L := S.generators;
	fi;
        if Length( arg ) = 3  and  arg[ 3 ] = true  then
            G := SubgroupNC( P, L );
        else
            G := Subgroup( P, L );
        fi;
    fi;
    SetStabChainMutable( G, S );

    return G;
end);

#############################################################################
##
#F  DepthSchreierTrees( <S> ) . . . . . . . . . . . . depth of Schreier trees
##
InstallGlobalFunction( DepthSchreierTrees, function( S )
    local   depths,  gens,  dep,  pnt,  sum,  i;

    depths := "";
    gens := [  ];
    while IsBound( S.stabilizer )  do
        UniteSet( gens, S.labels{ S.genlabels } );
        dep := [  ];  dep[ S.orbit[ 1 ] ] := -1;
        for pnt  in S.orbit  do
            dep[ pnt ] := dep[ pnt ^ S.transversal[ pnt ] ] + 1;
        od;
        sum := 0;
        for i  in dep  do
            sum := sum + i;
        od;
        i := sum / Length( S.orbit );
        Append( depths, Concatenation( String( Int( i ) ), "." ) );
        i := Int( 100 * ( i - Int( i ) ) );
        if i < 10  then
            Append( depths, "0" );
        fi;
        Append( depths, Concatenation( String( i ), "-",
                String( Maximum( Compacted( dep ) ) ), " " ) );
        S := S.stabilizer;
    od;
    Append( depths, Concatenation( String( Length( gens ) ), " gens" ) );
    return depths;
end );

#############################################################################
##
#F  AddGeneratorsExtendSchreierTree( <S>, <new> ) . . . . . .  add generators
##
##  This function may be called with a generatorless <S>.
##
InstallGlobalFunction( AddGeneratorsExtendSchreierTree, function( S, new )
    local   gen,  pos,  # new generator and its position in <S>.labels
            old,  ald,  # genlabels before extension
            len,        # initial length of the orbit of <S>
            img,        # image during orbit algorithm
            debug_fct,  # whether to debug or not the program
            i,  j;      # loop variable

    debug_fct:=true;
    if debug_fct then
      Print("GAP AGEST : Beginning of AddGeneratorsExtendSchreierTree\n");
#      Print("GAP AGEST 1: genlabels=", S.genlabels, "\n");
    fi;
    # Put in the new labels.
    old := BlistList( [ 1 .. Length( S.labels ) ], S.genlabels );
    old[ 1 ] := true;
    ald := StructuralCopy( old );
    if debug_fct then
        Print("GAP AGEST newgens=", new, "\n");
        Print("DEBUG AGEST 1: old=", old, "\n");
        Print("DEBUG AGEST 1: ald=", ald, "\n");
        Print("DEBUG AGEST labels=", S.labels, "\n");
        Print("DEBUG AGEST 2: genlabels=", S.genlabels, "\n");
    fi;
    for gen  in new  do
        pos := Position( S.labels, gen );
        if pos = fail  then
            Add( S.labels, gen );
            Add( old, false );
            Add( ald, true );
            if debug_fct then
              Print("GAP AGEST  genlabels insert 1:\n");
#              Print("GAP AGEST  genlabels insert 1: pos=", Length(S.labels), "\n");
#              Print("GAP AGEST  genlabels insert X: pos=", Length(S.labels), "\n");
            fi;
            Add( S.genlabels, Length( S.labels ) );
        elif not ald[ pos ]  then
            if debug_fct then
              Print("GAP AGEST  genlabels insert 2:\n");
#              Print("GAP AGEST  genlabels insert 2: pos=", pos, "\n");
#              Print("GAP AGEST  genlabels insert X: pos=", pos, "\n");
            fi;
            Add( S.genlabels, pos );
        fi;
        if     IsBound( S.generators )
           and pos <> 1 and not gen in S.generators  then
            Add( S.generators, gen );
        fi;
    od;
    if debug_fct then
        Print("DEBUG AGEST 2: old=", old, "\n");
        Print("DEBUG AGEST 2: ald=", ald, "\n");
    fi;

    # Extend the orbit and the transversal with the new labels.
    len := Length( S.orbit );
    if debug_fct then
        Print("GAP AGEST len=", len, "\n");
        Print("GAP AGEST S->orbit=", S.orbit, "\n");
    fi;
    i := 1;
    if debug_fct then
#        Print("XXX ELIMINATE begin\n");
    fi;
    if IsBound( S.cycles )  then
        if debug_fct then
          Print("GAP AGEST Cycles S.cycles=", S.cycles, "\n");
        fi;
        while i <= Length( S.orbit )  do
            if debug_fct then
              Print("GAP   AGEST i=", i, " |cycles|=", Length(S.cycles), "\n");
            fi;
            for j  in S.genlabels  do

                # Use new labels for old points, all labels for new points.
                if i > len  or  not old[ j ]  then
                    img := S.orbit[ i ] / S.labels[ j ];
                    if debug_fct then
                      Print("GAP     AGEST img=", img, " g=", S.labels[j], "\n");
                    fi;
                    if IsBound( S.translabels[ img ] )  then
                        if debug_fct then
                            if i > Length(S.cycles) then
                                Print("DEBUG CycleLength i=", i, " |S.cycles|=", Length(S.cycles), "\n");
                            fi;
                            Print("GAP       |S->cycles|=", Length(S.cycles), " i=", i, "\n");
                        fi;
                        S.cycles[ i ] := true;
                        if debug_fct then
                            Print("GAP       AGEST assign true\n");
                        fi;
                    else
                        S.translabels[ img ] := j;
                        S.transversal[ img ] := S.labels[ j ];
                        if debug_fct then
                            Print("GAP       AGEST S.transversal[img]=", S.transversal[img], "\n");
                        fi;
                        Add( S.orbit, img );
                        if debug_fct then
                            Print("GAP       AGEST insert img\n");
                        fi;
                        Add( S.cycles, false );
                    fi;
                fi;

            od;
            i := i + 1;
        od;

    else
        if debug_fct then
          Print("GAP AGEST No Cycles\n");
        fi;
        while i <= Length( S.orbit )  do
#            if debug_fct then
#              Print("GAP AGEST i=", i, "\n");
#            fi;
            for j  in S.genlabels  do

                # Use new labels for old points, all labels for new points.
                if i > len  or  not old[ j ]  then
                    img := S.orbit[ i ] / S.labels[ j ];
                    if not IsBound( S.translabels[ img ] )  then
                        S.translabels[ img ] := j;
                        S.transversal[ img ] := S.labels[ j ];
                        Add( S.orbit, img );
                    fi;
                fi;

            od;
            i := i + 1;
        od;
    fi;
    if debug_fct then
#        Print("XXX ELIMINATE end\n");
    fi;
end );

#############################################################################
##
#F  ChooseNextBasePoint( <S>, <base>, <newgens> ) . . . . . . . . . . . local
##
InstallGlobalFunction( ChooseNextBasePoint, function( S, base, newgens )
    local   i,  pnt,  bpt,  pos;

    Print("GAP base = ", base, "\n");
    i := 1;
    while     i <= Length( base )
          and ForAll( newgens, gen -> base[ i ] ^ gen = base[ i ] )  do
        i := i + 1;
    od;
    if i <= Length( base )  then
        pnt := base[ i ];
    elif IsPermCollection( newgens )  then
        pnt := SmallestMovedPoint( newgens );
    else
        pnt := 1;
    fi;

    # If <pnt> is before  the base point <bpt>  of  <S>, insert a  new level.
    # `Before' means (1) <pnt> before <bpt> in <base> or (2) <pnt> in <base>,
    # <bpt> not in <base> or (3) <pnt> less than <bpt> both not in <base>.
    if IsBound( S.orbit )  then
        bpt := S.orbit[ 1 ];
        pos := Position( base, bpt );
    else
        bpt := infinity;
        pos := fail;
    fi;
    Print("GAP pnt=", pnt, " bpt=", bpt, " pos=", pos, "\n");
    if    pos <> fail  and  i < pos              # (1)
       or pos =  fail  and  i <= Length( base )  # (2)
       or pos =  fail  and  pnt < bpt  then      # (3)
        Print("GAP matching test\n");
#        Print("ChooseNextBasePoint: Before InsertTrivialStabilizer, S=\n");
#        MyPrintStabChain(S);
        InsertTrivialStabilizer( S, pnt );
#        Print("ChooseNextBasePoint:  After InsertTrivialStabilizer, S=\n");
#        MyPrintStabChain(S);
        if IsBound( S.stabilizer.cycles )  then
            S.cycles := [ false ];
            Print("GAP   Initializing cycles\n");
        elif IsBound( S.stabilizer.relativeOrders )  then
            Unbind( S.stabilizer.relativeOrders );
            Unbind( S.stabilizer.base           );
        fi;
    fi;
    Print("GAP Exiting ChooseNextBasePoint\n");
end );

#############################################################################
##
#F  StabChainStrong( <S>, <newgens>, <options> )  . . Schreier-Sims algorithm
##
InstallGlobalFunction( StabChainStrong, function( S, newgens, options )
    local   gen,        # one generator from <newgens>
            pnt,        # next base point to use
            len,        # length of orbit of <S>
            pnts,       # points to use for Schreier generators
            p,          # point in orbit of <S>
            rep, r, rr, # representative of <p>
            gen1, old,  # numbers of labels to be used for Schreier gens
            g,          # one of these labels
            sch,        # Schreier generator for '<S>.stabilizer'
            img,  i,  j;# loop variables
    Print("GAP Begin StabChainStrong : newgens=", newgens, "\n");
    Print("GAP StabChainStrong 1: genlabels=", S.genlabels, "\n");
    # It is possible to prescribe a new operation domain for each level.
    if IsPermOnEnumerator( S.identity )  then
        newgens := List( newgens, gen -> PermOnEnumerator
                         ( Enumerator( S.identity ), gen ) );
    fi;

    # Determine the next base point.
    if IsBound( options.nextBasePoint )  then
        if not IsBound( S.orbit )  then
            pnt := options.nextBasePoint( S );
            InsertTrivialStabilizer( S, pnt );
        fi;
    else
        ChooseNextBasePoint( S, options.base, newgens );
    fi;
    Print("GAP StabChainStrong 2: genlabels=", S.genlabels, "\n");

    # Add the new generators to <S>.
    pnt := S.orbit[ 1 ];
    len := Length( S.orbit );
    old := Length( S.genlabels );
    Print("GAP Before AddGeneratorsExtendSchreierTree\n");
    AddGeneratorsExtendSchreierTree( S, newgens );

    # If a new generator fixes the base point, put it into the stabilizer.
    Print("GAP newgens=", newgens, "\n");
    for gen  in newgens  do
        Print("GAP eGen=", gen, " eGen=", gen, "\n");
        if gen <> S.identity  and  pnt ^ gen = pnt  then
            Print("CPP   1: Calling StabChainStrong with eGen=", gen, "\n");
            StabChainStrong( S.stabilizer, [ gen ], options );
        fi;
    od;

    # Compute the Schreier generators (seems to work better backwards).
    if IsBound( S.cycles )  then
        pnts := ListBlist( [ 1 .. Length( S.orbit ) ], S.cycles );
    else
        pnts := [ 1 .. Length( S.orbit ) ];
    fi;
    Print("GAP pnts = ", pnts, "\n");
    Print("GAP Usecycle=", IsBound(S.cycles), "\n");
    if IsBound(S.cycles) then
      Print("GAP cycles=", S.cycles, "\n");
    fi;
    gen1 := 1;
    Print("GAP StabChainStrong O=", S.orbit, "\n");
    for i  in Reversed( pnts )  do
        p := S.orbit[ i ];
	Print("GAP StabChainStrong i=", i, " p=", p, "\n");
        if IsBound( options.knownBase )  then
            rep := InverseRepresentativeWord( S, p );
        else
            rep := InverseRepresentative( S, p );
        fi;

        # Take only new generators for old, all generators for new points.
        if i <= len  then
            gen1 := old + 1;
        fi;
	Print("GAP StabChainStrong gen1=", gen1, " rep=", rep, "\n");
        for j  in [ gen1 .. Length( S.genlabels ) ]  do
          g := S.labels[ S.genlabels[ j ] ];
          Print("GAP StabChainStrong   j=", j, " g=", g, "\n");

          # Avoid computing Schreier generators that will be trivial.
          if S.translabels[ p / g ] <> S.genlabels[ j ]  then
            # If a base is known, use it to test the Schreier generator.
            if IsBound( options.knownBase )  then
                if not MembershipTestKnownBase( S,
                           options.knownBase, [ rep, [ g ] ] )  then

                    # If this is the  first Schreier generator for this orbit
                    # point, multiply the representative.
                    if IsList( rep )  then
                        r := S.identity;
                        for rr  in rep  do
                            r := LeftQuotient( rr, r );
                        od;
                        rep := r;
                    fi;

                    sch := rep / g;
                    img := pnt ^ sch;
                    while img <> pnt  do
                        sch := sch * S.transversal[ img ];
                        img := img ^ S.transversal[ img ];
                    od;
                    StabChainStrong( S.stabilizer, [ sch ], options );
                fi;

            # If no  base is known, construct the  Schreier generator and put
            # it in the chain if it is non-trivial.
            else
                sch := SiftedPermutation( S, ( g * rep ) ^ -1 );
		Print("GAP sch=", sch, " g=", g, " rep=", rep, "\n");
                if sch <> S.identity  then
                    StabChainStrong( S.stabilizer, [ sch ], options );
                fi;
            fi;

          fi;
        od;
    od;
    Print("GAP exiting StabChainStrong\n");
end );


#############################################################################
##
#F  StabChainForcePoint( <S>, <pnt> ) . . . . . . . .  force <pnt> into orbit
##
InstallGlobalFunction( StabChainForcePoint, function( S, pnt )
    Print("GAP Beginning of StabChainForcePoint pnt=", pnt, "\n");
    PrintStabChain(S);
    # Do nothing if <pnt> is already in the orbit of <S>.
    if    not IsBound( S.translabels )
       or not IsBound( S.translabels[ pnt ] )  then
        Print("GAP Matching the first test\n");

        # If all generators of <S> fix <pnt>, insert a trivial stabilizer.
        if IsFixedStabilizer( S, pnt )  then
            Print("GAP Matching the second test\n");
            InsertTrivialStabilizer( S, pnt );

        # Get  <pnt> in   the orbit   of   the stabilizer and  swap   the two
        # stabilizers.  If this  is  unsuccessful, the   stabilizer  chain is
        # incorrect, return `false' then.
        elif    not StabChainForcePoint( S.stabilizer, pnt )
             or not StabChainSwap( S )  then
            Print("GAP StabChainForcePoint, return false\n");
            return false;
        fi;

    fi;
    Print("GAP StabChainForcePoint, return true\n");
    PrintStabChain(S);
    return true;
end );

#############################################################################
##
#F  StabChainSwap( <S> )  . . . . . . . . . . . . . . .  swap two base points
##
InstallGlobalFunction( StabChainSwap, function( S )
    local   a,  b,      # basepoints that are to be switched
            T,  Tstab,  # copy of $S$ with $Tstab$ becomes $S_b$
            len,        # length of $Tstab.orbit$ to be reached
            pnt,        # one point from $a^S$ not yet in $a^{T_b}$
            ind,        # index of <pnt> in $S.orbit$
            img,        # image $b^{Rep(S,pnt)^-}$
            gen,        # new generator of $T_b$
            i;          # loop variable
    Print("GAP Beginning of StabChainSwap\n");
    PrintStabChain(S);
    # get the two basepoints $a$ and $b$ that we have to switch
    a := S.orbit[ 1 ];
    b := S.stabilizer.orbit[ 1 ];
    Print("GAP StabChainSwap a=", a, " b=", b, "\n");

    # set $T = S$ and compute $b^T$ and a transversal $T/T_b$
    T := EmptyStabChain( S.labels, S.identity, b );
    Unbind( T.generators );
    Print("GAP LGens=", S.generators, "\n");
    AddGeneratorsExtendSchreierTree( T, S.generators );
    Print("GAP StabChainSwap : after first AGEST\n");

    # initialize $Tstab$, which will become $T_b$
    Tstab := EmptyStabChain( [  ], S.identity, a );
    Unbind( Tstab.generators );
    Print("GAP StabChainSwap : before second AGEST gens=", S.stabilizer.stabilizer.generators, "\n");
    AddGeneratorsExtendSchreierTree( Tstab,
            S.stabilizer.stabilizer.generators );
    Print("GAP StabChainSwap : after second AGEST\n");
    Print("GAP After AddGeneratorsExtendSchreierTree 1 : Tstab=\n");
    PrintStabChain(Tstab);

    # in the end $|b^T||a^{T_b}| = [T:T_{ab}] = [S:S_{ab}] = |a^S||b^{S_a}|$
    ind := 1;
    len := Length( S.orbit ) * Length( S.stabilizer.orbit )
           / Length( T.orbit );
    Print("GAP StabChainSwap |Tstab->orbit|=", Length(Tstab.orbit), " len=", len, "\n");
    while Length( Tstab.orbit ) < len  do
        Print("GAP beginning of loop\n");

        # choose a point $pnt \in a^S \ a^{T_b}$ with representative $s$
        repeat
            ind := ind + 1;

            # If <ind> exceeds the length of  <S>.orbit, <S> was an incorrect
            # stabilizer chain.
            if ind > Length( S.orbit )  then
                return false;
            fi;
            Print("GAP |orbit|=", Length(S.orbit), " ind=", ind, "\n");

            pnt := S.orbit[ ind ];
            Print("GAP ind=", ind, " pnt=", pnt, "\n");
        until not IsBound( Tstab.translabels[ pnt ] );
        Print("GAP ind=", ind, " pnt=", pnt, "\n");

        # find out what $s^-$ does with $b$ (without computing $s$!)
        img := b;
        i := pnt;
        while i <> a  do
            img := img ^ S.transversal[ i ];
            i   := i   ^ S.transversal[ i ];
        od;
        Print("GAP i=", i, " img=", img, "\n");

        # if $b^{s^-}} \in b^{S_a}$ with representative $r \in S_a$
        if IsBound( S.stabilizer.translabels[ img ] )  then
            Print("GAP Determining gen, step 1\n");

            # with $gen = s^- r^-$ we have
            # $b^gen = {b^{s^-}}^{r^-} = img^{r-} = b$, so $gen \in S_b$
            # and $pnt^gen = {pnt^{s^-}}^{r^-} = a^{r-} = a$, so $gen$ is new
            gen := S.identity;
            while pnt ^ gen <> a  do
                Print("GAP pnt^gen=", pnt^gen, "\n");
                gen := gen * S.transversal[ pnt ^ gen ];
            od;
            Print("GAP Determining gen, step 2 gen=", gen, "\n");
            while b ^ gen <> b  do
                Print("GAP b^gen=", b^gen, "\n");
                gen := gen * S.stabilizer.transversal[ b ^ gen ];
            od;
            Print("GAP Determining gen, step 3 gen=", gen, "\n");
            AddGeneratorsExtendSchreierTree( Tstab, [ gen ] );
            Print("GAP After AddGeneratorsExtendSchreierTree 2 : Tstab=\n");
            PrintStabChain(Tstab);

        fi;

    od;
    Print("GAP After while loop\n");
    Print("GAP T=\n");
    PrintStabChain(T);
    Print("GAP Tstab=\n");
    PrintStabChain(Tstab);

    # copy everything back into the stabchain
    S.labels      := T.labels;
    S.genlabels   := T.genlabels;
    S.orbit       := T.orbit;
    S.translabels := T.translabels;
    S.transversal := T.transversal;
    PrintStabChain(S);
    Print("GAP StabChainSwap 1: |orbit|=", Length(Tstab.orbit), "\n");
    if Length( Tstab.orbit ) = 1  then
        S.stabilizer := S.stabilizer.stabilizer;
        Print("GAP StabChainSwap 2:\n");
    else
        S.stabilizer.labels      := Tstab.labels;
        S.stabilizer.genlabels   := Tstab.genlabels;
	if not IsBound(Tstab.generators) then
	  Tstab.generators:=Tstab.labels{Tstab.genlabels};
	fi;
        S.stabilizer.generators  := Tstab.generators;
        S.stabilizer.orbit       := Tstab.orbit;
        S.stabilizer.translabels := Tstab.translabels;
        S.stabilizer.transversal := Tstab.transversal;
        Print("GAP StabChainSwap 3:\n");
    fi;
    PrintStabChain(S);
    return true;
end );


#############################################################################
##
#F  LabsLims( <lab>, <hom>, <labs>, <lims> )  . . . .  help for next function
##
InstallGlobalFunction( LabsLims, function( lab, hom, labs, lims )
    local   pos;

    pos := Position( labs, lab );
    if pos = fail  then
        AddSet( labs, lab );
        pos := Position( labs, lab );
        if IsFunction( hom )  then
            Add(lims, hom(lab), pos);
        else
            Add(lims, lab ^ hom, pos);
        fi;
    fi;
    return lims[ pos ];
end );

#############################################################################
##
#F  ConjugateStabChain( <arg> ) . . . . . . . . .  conjugate stabilizer chain
##
InstallGlobalFunction( ConjugateStabChain, function( arg )
    local   S,  T,  hom,  map,  cond,           # arguments
            newlevs,                            # new labels lists
            len,                                # number of labels in <S>
            labels,  labpos,  orbit,  edges,    # conjugated components
            labs,  lims,                        # list of all labels/images
            img,  pnt, nb_assign,               # image of label and point
            pos,  L,  l,  i;                    # loop variables

    # Get the arguments.
    S := arg[ 1 ];  T := arg[ 2 ];  hom := arg[ 3 ];  map := arg[ 4 ];
    if Length( arg ) > 4  then  cond := arg[ 5 ];
                          else  cond := S -> IsBound( S.stabilizer );  fi;
    newlevs := [  ];

    # Prepare common  lists for the labels and  their images at the different
    # levels.
    labs := [ S.identity ];
    lims := [ T.identity ];
#    Print(NullMat(5));


    # ND MDS: With the existing code what is done in the code is mapping the
    # ND MDS: labels from the one in the S.labels to the one in return.
    # ND MDS: The last S.labels is not mapped for some strange reason.
    # ND MDS: Maybe because it is not actually used in the computation.

    # Loop over the stabilizer chain.
    nb_assign:=0;
    while cond( S )  do
        len := Length( S.labels );
        Print("DEBUG ConjugateStabChain, Begin loop\n");

        # If this is a  new  labels component, map  the  labels and mark  the
        # component so that it can be recognized at deeper levels.
        if len = 0  or  IsPerm( S.labels[ len ] )  then
            Print("DEBUG ConjugateStabChain, matching test. len=", len, "\n");
            if IsPerm( hom )  then
                # ND MDS. This is the case that applies to us when calling with permutation mapping
                labels := OnTuples( S.labels, hom );
                Print("DEBUG ConjugateStabChain, 1 : labels=", labels, "\n");
                labpos := [ 1 .. len ];
            else
                if IsIdenticalObj( S, T )  then  labels := [ T.identity ];
                                        else  labels := T.labels;        fi;
                Print("DEBUG ConjugateStabChain, 2 : labels=", labels, "\n");
                labpos := ListWithIdenticalEntries( len, 0 );
                labpos[ 1 ] := 1;
            fi;
            # ND MDS: The genlabels are actually not used.
            Add( S.labels, rec( labels := labels,
                                labpos := labpos,
                             genlabels := Set( S.genlabels ) ) );
            Add( newlevs, S.labels );

        # The current labels component is not  new, so take the mapped labels
        # from the   mark  that was  inserted  when the  component  was first
        # encountered.
        else
            labels := S.labels[ len ].labels;
            Print("DEBUG ConjugateStabChain, 3 : labels=", labels, "\n");
#            Print(NullMat(5));

            labpos := S.labels[ len ].labpos;
            UniteSet( S.labels[ len ].genlabels, S.genlabels );
        fi;

        # Map the orbit and edges.
        edges := [  ];
        if IsPerm( map )  and  IsPerm( hom )  then
            # This is the case that matters to us.
            orbit := OnTuples( S.orbit, map );
            edges{ orbit } := S.translabels{ S.orbit };
        else
            orbit := [  ];
            for pnt  in S.orbit  do
                if   IsFunction( map )  then  img := map( pnt );
                elif IsList    ( map )  then  img := map[ pnt ];
                else                          img := pnt ^ map;   fi;
                if not IsBound( edges[ img ] )  then
                    Add( orbit, img );
                    pos := labpos[ S.translabels[ pnt ] ];

                    # We can  already map  the  labels that appear   as edges
                    # because we know that their images  will be distinct and
                    # non-trivial.
                    if pos = 0  then
                        # ND MDS: this case should not happen for us because value 0 occurs only when in the non-perm case
                        Add( labels, LabsLims( S.transversal[ pnt ], hom,
                                labs, lims ) );
                        Print("DEBUG ConjugateStabChain, 4 : labels=", labels, "\n");
                        pos := Length( labels );
                        labpos[ S.translabels[ pnt ] ] := pos;
                    fi;

                    edges[ img ] := pos;
                fi;
            od;
            if not IsPerm( hom )  then
                T.labpos := labpos;
            fi;
        fi;

        # Build a level of <T>  (`genlabels' will be completed when  `labpos'
        # is complete).
        Print("DEBUG ConjugateStabChain, 5 : labels=", labels, "\n");
        nb_assign:=nb_assign+1;
        T.labels      := labels;
        T.genlabels   := S.genlabels;
        T.orbit       := orbit;
        T.translabels := edges;
        T.transversal := [  ];
        T.transversal{ orbit } := labels{ edges{ orbit } };
#        Print("ConjugateStabChain 1: Now T.transversal{orbit}=", T.transversal{orbit}, "\n");

        # Step down to the stabilizer.
        S := S.stabilizer;
        if not IsBound( T.stabilizer )  then
            T.stabilizer := EmptyStabChain( T.labels, T.identity );
        fi;
        T := T.stabilizer;

    od;
    Print("DEBUG nb_assign=", nb_assign, "\n");

    # For   the distinct labels  components  of  the original  chain, map the
    # labels that  did  not appear     as edges  and   remove the   auxiliary
    # components.
    for L  in newlevs  do
        l := L[ Length( L ) ];
        i := Position( l.labpos, 0 );
        # ND MDS: In the perm case, we do not have the case of a 0 in labpos
        # ND MDS: So, this does not apply here.
        while i <> fail  do
            if i in l.genlabels  then
                img := LabsLims( L[ i ], hom, labs, lims );
                pos := Position( l.labels, img );
                if pos = fail  then
                    Add( l.labels, img );
                    pos := Length( l.labels );
                fi;
                l.labpos[ i ] := pos;
            fi;
            i := Position( l.labpos, 0, i );
        od;
        Unbind( L[ Length( L ) ] );
    od;

    # Now that all labels have been mapped, complete  the `genlabels' and put
    # in `generators'.
    # ND MDS: The code below obviously does not apply to the IsPerm(hom) case.
    if not IsPerm( hom )  then
        L := arg[ 2 ];
        while IsBound( L.labpos )  do
            L.genlabels := Set( L.labpos{ L.genlabels } );
            RemoveSet( L.genlabels, 0 );
            RemoveSet( L.genlabels, 1 );
            Unbind( L.labpos );
            L := L.stabilizer;
        od;
    fi;
    # ND MDS: The code below does not apply as we do not use the generators entry.
    L := arg[ 2 ];
    while IsBound( L.stabilizer )  do
        L.generators := L.labels{ L.genlabels };
        L := L.stabilizer;
    od;

    # Return the mapped stabilizer from the first level  where <cond> was not
    # satisfied (i.e., the ``end'' of the original chain).
    return T;

end );


GetStabilizerDepth:=function(S)
  local T, dep;
  T := S;
#  Print("GetStabilizerDepth\n");
  dep:=1;
  while(true)
  do
    if IsBound(T.stabilizer)=false then
      break;
    fi;
    dep:=dep+1;
    T := T.stabilizer;
  od;
  return dep;
end;


#############################################################################
##
#F  ChangeStabChain(<G>,<base>[,<reduced>])  change/extend a stabilizer chain
##
##  reduced = -1    : extend stabilizer chain
##  reduced = false : change stabilizer chain, do not reduce it
##  reduced = true  : change stabilizer chain, reduce it
##
InstallGlobalFunction(ChangeStabChain,function( arg )
local   G,  base,  reduced,
        cnj,  S,  newBase,  old,  new,  i, dep1, dep2, 
        strG_orig, strG_current, strS_current, strG_final, strS, idx, KeyUpdating;
    # Get the arguments.
    G := arg[ 1 ];
    strG_orig:=GetStringExpressionOfStabChain(G);
    strG_current:=strG_orig;
    Print("GAP Beginning ChangeStabChain, GetStabilizerDepth = ", GetStabilizerDepth(G), "\n");
    base := arg[ 2 ];
    if Length( arg ) > 2  then
        reduced := arg[ 3 ];
    else
        reduced := true;
    fi;

    cnj := G.identity;
    S := G;
    strS_current:=GetStringExpressionOfStabChain(S);
    Print("GAP we have strS_current\n");
    idx:=0;
    KeyUpdating:=function(str)
      local strGloc, strSloc, DoPrint;
      strGloc:=GetStringExpressionOfStabChain(G);
      strSloc:=GetStringExpressionOfStabChain(S);
      idx:=idx+1;
      DoPrint:=false;
      if DoPrint then
        Print("GAP KU At ", idx, " of ", str, " dep(G)/dep(S)=", GetStabilizerDepth(G), "/", GetStabilizerDepth(S), "\n");
        Print("GAP KU At step ", idx, " of ", str, "\n");
        Print("GAP sgs(G) = ", StrongGeneratorsStabChain(G), "\n");
        Print("GAP sgs(S) = ", StrongGeneratorsStabChain(S), "\n");
        if strG_current=strGloc then
          Print("GAP   KU At step ", idx, " of ", str, " no change of G\n");
        else
          Print("GAP   KU At step ", idx, " of ", str, " CHANGE of G\n");
          strG_current:=strGloc;
        fi;
        if strS_current=strSloc then
          Print("GAP   KU At step ", idx, " of ", str, " no change of S\n");
        else
          Print("GAP   KU At step ", idx, " of ", str, " CHANGE of S\n");
          strS_current:=strSloc;
        fi;
      fi;
      Print("GAP Gptr at str=", str, "\n");
      PrintStabChain(G);
      Print("GAP Sptr at str=", str, "\n");
      PrintStabChain(S);
    end;
    KeyUpdating("Begin ChangeStabChain");
    newBase := [  ];
    i := 1;
    Print("GAP ChangeStabChain base = ", base, "\n");
    Print("GAP ChangeStabChain 1 orbit=", PrintTopOrbit(G), "\n");
    while IsBound( S.stabilizer )  or  i <= Length( base )  do
        KeyUpdating("Before BasePoint");
        old := BasePoint( S );
        Print("GAP ChangeStabChain old=", old, " i=", i, " |base|=", Length(base), "\n");
        KeyUpdating("After BasePoint");
        # Cut off unwanted trivial stabilizers at the end.
        if     Length( S.genlabels ) = 0
           and ( reduced = true  or  i > Length( base ) )  then
            Print("GAP Before RemoveStabChain\n");
            KeyUpdating("Before RemoveStabChain");
            RemoveStabChain( S );
            KeyUpdating("After RemoveStabChain");
            i := Length( base ) + 1;

        # Determine the new base point for this level.
        elif i <= Length( base )  then
            new := base[ i ] / cnj;
            Print("GAP newpnt=", new, "\n");
            i := i + 1;

            # Stabilizer chain extension.
            if reduced = -1  then
                AddSet( newBase, new );
                if new <> old  then
                    if IsFixedStabilizer( S, new )  then
                        InsertTrivialStabilizer( S, new );
                        KeyUpdating("After InsertTrivialStabilizer");
                    else
                        Error("<base> must be an extension of base of <G>");
                    fi;
                fi;
                S := S.stabilizer;
                KeyUpdating("After S:=S.stabilizer 1");

            # Base change. Return `false' if <S> turns out to be incorrect.
            elif reduced = false  or  not IsFixedStabilizer( S, new )  then
                if IsBound( S.stabilizer )  then
                    KeyUpdating("Before StabChainForcePoint");
                    if not StabChainForcePoint( S, new )  then
                        KeyUpdating("After StabChainForcePoint, return false");
                        return false;
                    fi;
                    KeyUpdating("After StabChainForcePoint, return true");
                    Print("GAP 1: cnj=", cnj, "\n");
                    cnj := LeftQuotient( InverseRepresentative( S, new ),
                                   cnj );
                    Print("GAP 2: cnj=", cnj, "\n");
                else
                    InsertTrivialStabilizer( S, new );
                    KeyUpdating("After InsertTrivialStabilizer");
                fi;
                AddSet( newBase, S.orbit[ 1 ] );
                S := S.stabilizer;
                KeyUpdating("After S:=S.stabilizer 2");
            fi;

        # Delete unwanted trivial  stabilizers  (e.g., double  points  in the
        # base).
        elif    old in newBase
             or reduced = true  and  Length( S.orbit ) = 1  then
            Print("GAP Stabilizer shift in ChangeStabChain\n");
            S.labels     := S.stabilizer.labels;
            S.genlabels  := S.stabilizer.genlabels;
            S.generators := S.stabilizer.generators;
            KeyUpdating("MissCase 1");
            if IsBound( S.stabilizer.orbit )  then
                S.orbit       := S.stabilizer.orbit;
                S.translabels := S.stabilizer.translabels;
                S.transversal := S.stabilizer.transversal;
            else
                Unbind( S.orbit );
                Unbind( S.translabels );
                Unbind( S.transversal );
            fi;
            KeyUpdating("MissCase 2");
            if IsBound( S.stabilizer.stabilizer )  then
                S.stabilizer := S.stabilizer.stabilizer;
            else
                Unbind( S.stabilizer );
            fi;
            KeyUpdating("MissCase 3");

        # Simply move down the stabilizer chain (to look for double points).
        else
            AddSet( newBase, old );
            S := S.stabilizer;
            KeyUpdating("After S:=S.stabilizer 3");
        fi;

    od;
    Print("GAP LEAVE i=", i, " |base|=", Length(base), "\n");
    KeyUpdating("After the loop");
    strG_final:=GetStringExpressionOfStabChain(G);
    strS:=GetStringExpressionOfStabChain(S);

    # Conjugate to move all the points to the beginning of their orbit.
    Print("GAP Before ConjugateStabChain cnj=", cnj, "\n");
    if cnj <> S.identity  then
        ConjugateStabChain( G, G, cnj, cnj );
    fi;
    KeyUpdating("After ConjugateStabChain");
    Print("GAP Leaving ChangeStabChain\n");
    return true;
end);

#############################################################################
##
#F  ExtendStabChain( <S>, <base> )  . . . . . . . . extend a stabilizer chain
##
InstallGlobalFunction(ExtendStabChain,function( S, base )
    ChangeStabChain( S, base, -1 );
end);


#############################################################################
##
#F  ReduceStabChain( <S> )  . . . . . . . . . . . . reduce a stabilizer chain
##
InstallGlobalFunction(ReduceStabChain,function( S )
    ChangeStabChain( S, [  ], true );
end);

#############################################################################
##
#F  EmptyStabChain( <labels>,<id>[,<limgs>,<idimg>][,<pnt>] ) . .  stab chain
##
InstallGlobalFunction(EmptyStabChain,function( arg )
local   S;

    S := rec(  labels := arg[ 1 ],
            genlabels := [  ],
           generators := [  ],
             identity := arg[ 2 ] );
    if Length( S.labels ) = 0  or  S.labels[ 1 ] <> S.identity  then
        Add( S.labels, S.identity, 1);
    fi;
    if Length( arg ) >= 4  then
        S.labelimages := arg[ 3 ];
        S.genimages   := [  ];
        S.idimage     := arg[ 4 ];
        if Length( S.labelimages ) = 0  then
            Add( S.labelimages, S.idimage );
        fi;
    fi;
    if Length( arg ) mod 2 = 1  then
        InitializeSchreierTree( S, arg[ Length( arg ) ] );
    fi;
    return S;
end);

#############################################################################
##
#F  InitializeSchreierTree( <S>, <pnt> )  . . . .  initialize a Schreier tree
##
InstallGlobalFunction( InitializeSchreierTree, function( S, pnt )
    S.orbit       := [ pnt ];
    S.translabels := [  ];  S.translabels[ pnt ] := 1;
    S.transversal := [  ];  S.transversal[ pnt ] := S.identity;
    if IsBound( S.idimage )  then
        S.transimages := [  ];  S.transimages[ pnt ] := S.idimage;
    fi;
end );

#############################################################################
##
#F  InsertTrivialStabilizer( <S>, <pnt> ) . .  add redundant base point <pnt>
##
InstallGlobalFunction( InsertTrivialStabilizer, function( S, pnt )
    Print("GAP begin InsertTrivialStabilizer\n");
    S.stabilizer := ShallowCopy( S );
    S.genlabels  := ShallowCopy( S.stabilizer.genlabels );
    if IsBound( S.generators )  then
        S.generators := ShallowCopy( S.stabilizer.generators );
        if IsBound( S.idimage )  then
            S.genimages := ShallowCopy( S.stabilizer.genimages );
        fi;
    fi;
    InitializeSchreierTree( S, pnt );
end );

#############################################################################
##
#F  RemoveStabChain( <S> )  . . . . . . . .  cut off rest of stabilizer chain
##
InstallGlobalFunction(RemoveStabChain,function( S )
local  name;

    for name  in RecNames( S )  do
        if name <> "identity"  and  name <> "idimage"  then
            Unbind( S.( name ) );
        fi;
    od;
    S.labels     := [ S.identity ];
    S.genlabels  := [  ];
    S.generators := [  ];
end);

#############################################################################
##
#F  BasePoint( <S> )  . . . . . . . . . . . . . . . . . . . base point of <S>
##
InstallGlobalFunction( BasePoint, function( S )
    if IsBound( S.orbit )  then  return S.orbit[ 1 ];
                           else  return false;         fi;
end );

#############################################################################
##
#F  IsInBasicOrbit( <S>, <pnt> )  . . . . . . . . .  is <pnt> in basic orbit?
##
InstallGlobalFunction( IsInBasicOrbit, function( S, pnt )
    return     IsBound( S.translabels )
           and IsBound( S.translabels[ pnt ] );
end );

#############################################################################
##
#F  IsFixedStabilizer( <S>, <pnt> ) . . . . . . . . .  is <pnt> fixed by <S>?
##
InstallGlobalFunction( IsFixedStabilizer, function( S, pnt )
    return ForAll( S.generators, gen -> pnt ^ gen = pnt );
end );

#############################################################################
##
#F  InverseRepresentative( <S>, <pnt> ) . .  perm mapping <pnt> to base point
##
InstallGlobalFunction( InverseRepresentative, function( S, pnt )
    local   bpt,  rep, te, fct_debug;

    bpt := S.orbit[ 1 ];
    rep := S.identity;
    fct_debug:=false;
    if fct_debug then
      Print("GAP INVREP bpt=", bpt, " pnt=", pnt, "\n");
    fi;
    while pnt <> bpt  do
	te:=S.transversal[pnt];
        if fct_debug then
          Print("GAP INVREP te=", te, "\n");
        fi;
	pnt:=pnt^te;
        if fct_debug then
          Print("GAP INVREP   pnt=", pnt, "\n");
        fi;
        rep := rep * te;
        if fct_debug then
          Print("GAP INVREP   rep=", rep, "\n");
        fi;
    od;
    if fct_debug then
      Print("GAP INVREP return rep=", rep, "\n");
    fi;
    return rep;
end );

#############################################################################
##
#F  QuickInverseRepresentative( <S>, <pnt> )  . . . . . . . same, but quicker
##
InstallGlobalFunction( QuickInverseRepresentative, function( S, pnt )
local   bpt,  rep,  lab,  pow;

    bpt := S.orbit[ 1 ];
    rep := S.identity;
    lab := S.translabels[ pnt ];
    pow := 1;
    while pnt <> bpt  do
        pnt := pnt ^ S.transversal[ pnt ];
        if S.translabels[ pnt ] = lab  then
            pow := pow + 1;
        else
            rep := rep * S.labels[ lab ] ^ pow;
            lab := S.translabels[ pnt ];
            pow := 1;
        fi;
    od;
    return rep;
end );

#############################################################################
##
#F  InverseRepresentativeWord( <S>, <pnt> ) . . . . . . . inverse rep as word
##
InstallGlobalFunction( InverseRepresentativeWord, function( S, pnt )
local   word,  bpt;

    word := [  ];
    bpt := S.orbit[ 1 ];
    while pnt <> bpt  do
        Add( word,   S.transversal[ pnt ] );
        pnt := pnt ^ S.transversal[ pnt ];
    od;
    return word;
end );

#############################################################################
##
#F  SiftedPermutation( <S>, <g> ) . . . . . . . . . . . .  sifted permutation
##
##  This function may be called with a generatorless <S>.
##
InstallGlobalFunction(SiftedPermutation,function( S, g )
local   bpt,  img;

    # The  following     condition   tests     `IsBound(S.stabilizer)',   not
    # `IsEmpty(S.genlabels)'. This is necessary because  the function may  be
    # called with an inconsistent chain from `NormalClosure'.
    while     IsBound( S.stabilizer )
          and g <> S.identity  do
        bpt := S.orbit[ 1 ];
        img := bpt ^ g;
        if IsBound( S.transversal[ img ] )  then
            while img <> bpt  do
                g := g * S.transversal[ img ];
                img := bpt ^ g;
            od;
            S := S.stabilizer;
        else
            return g;
        fi;
    od;
    return g;
end);

#############################################################################
##
#F  MinimalElementCosetStabChain( <S>, <g> )  . . .  minimal element of coset
##
##  This function may be called with a generatorless <S>.
##
InstallGlobalFunction(MinimalElementCosetStabChain,function( S, g )
local   p,i,a,bp,pp;

    while not IsEmpty( S.genlabels )  do

	if IsPlistRep(S.orbit) and IsPosInt(S.orbit[1])
	  and IsInternalRep(g) then
	  p:=SMALLEST_IMG_TUP_PERM(S.orbit,g);
	else
	  p:=infinity;
	  for i in S.orbit do
	    a:=i^g;
	    if a<p then
	      p:=a;
	    fi;
	  od;
	fi;

	bp:=S.orbit[1];
	pp:=p/g;
        while bp<>pp  do
            g:=LeftQuotient(S.transversal[pp],g);
	    pp:=p/g;
        od;
#        while S.orbit[ 1 ] ^ g <> p  do
#            g := LeftQuotient( S.transversal[ p / g ], g );
#        od;
        S := S.stabilizer;
    od;
    return g;
end);

#############################################################################
##
#M  MembershipTestKnownBase( <S>, <knownBase>, <word> ) . . . with known base
##
##  This function may be called with a generatorless <S>.
##
InstallMethod( MembershipTestKnownBase, "stabchain, base, word",true,
  [ IsRecord, IsList and IsCyclotomicCollection, IsList ], 0,
    function( S, knownBase, word )
    local   base,  g,  i,  j,  bpt;

    base := Concatenation( BaseStabChain( S ), knownBase );
    for g  in word  do
        if IsPerm( g )  then
            base := OnTuples( base, g );
        else
            for i  in Reversed( [ 1 .. Length( g ) ] )  do
                for j  in [ 1 .. Length( base ) ]  do
                    base[ j ] := base[ j ] / g[ i ];
                od;
            od;
        fi;
    od;
    while     not IsEmpty( S.genlabels )
          and IsBound( S.translabels[ base[ 1 ] ] )  do
        bpt := S.orbit[ 1 ];
        while base[ 1 ] <> bpt  do
            base := OnTuples( base, S.transversal[ base[ 1 ] ] );
        od;
        base := base{ [ 2 .. Length( base ) ] };
        S := S.stabilizer;
    od;
    return base = knownBase;
end );

InstallOtherMethod( MembershipTestKnownBase, true, [ IsRecord,
        IsGroup, IsPerm ], 0,
    function( S, G, t )
    return SiftedPermutation( S, t ) = S.identity;
end );

# the base `BaseOfGroup' does not need to confirm to the stabilizer chain
# `S'. therefore the following method is invalid. AH
#InstallOtherMethod( MembershipTestKnownBase, true, [ IsRecord,
#        IsGroup and HasBaseOfGroup, IsPerm ], 0,
#    function( S, G, t )
#    Error("this method may not be used!");
#    return MembershipTestKnownBase( S, BaseOfGroup( G ), [ t ] );
#end );

#############################################################################
##
#F  BaseStabChain( <S> )  . . . . . . . . . . . . . . . . . . . . . . .  base
##
##  This function may be called with a generatorless <S>.
##
InstallGlobalFunction(BaseStabChain,function( S )
    local   base;

    base := [  ];
    while IsBound( S.stabilizer )  do
        Add( base, S.orbit[ 1 ] );
        S := S.stabilizer;
    od;
    return base;
end);

#############################################################################
##
#F  SizeStabChain( <S> )  . . . . . . . . . . . . . . . . . . . . . . .  size
##
##  This function may be called with a generatorless <S>.
##
InstallGlobalFunction(SizeStabChain,function( S )
    local   size;

    size := 1;
    while not IsEmpty( S.genlabels )  do
        size := size * Length( S.orbit );
        S := S.stabilizer;
    od;
    return size;
end);

#############################################################################
##
#F  StrongGeneratorsStabChain( <S> )  . . . . . . . . . . . strong generators
##
InstallGlobalFunction(StrongGeneratorsStabChain,function( S )
    local   sgs;

    sgs := [  ];
    while not IsEmpty( S.generators )  do
        UniteSet( sgs, S.generators );
        S := S.stabilizer;
    od;
    return sgs;
end);

#############################################################################
##
#F  IndicesStabChain( <S> ) . . . . . . . . . . . . . . . . . . . . . indices
##
##  This function may be called with a generatorless <S>.
##
InstallGlobalFunction(IndicesStabChain,function( S )
    local   ind;

    ind := [  ];
    while IsBound( S.stabilizer )  do
        Add( ind, Length( S.orbit ) );
        S := S.stabilizer;
    od;
    return ind;
end);

#############################################################################
##
#F  ListStabChain( <S> )  . . . . . . . . . . . . .  stabilizer chain as list
##
##  This function may be called with a generatorless <S>.
##
InstallGlobalFunction(ListStabChain,function( S )
    local   list;

    list := [  ];
    while IsBound( S.stabilizer )  do
        Add( list, S );
        S := S.stabilizer;
    od;
    Add( list, S );
    return list;
end);

#############################################################################
##
#F  OrbitStabChain( <S>, <pnt> )  . . . . . . . . . orbit of stabilizer chain
##
InstallGlobalFunction(OrbitStabChain,function( S, pnt )
    if IsBound( S.edges )  and  IsBound( S.edges[ pnt ] )  then
        return ShallowCopy( S.orbit );
    else
        return OrbitPerms( S.generators, pnt );
    fi;
end);

#############################################################################
##
#M  MinimalStabChain( <G> ) . . . . . . . . . . . .  minimal stabilizer chain
##
InstallMethod( MinimalStabChain, "Perm", true, [ IsPermGroup] , 0,
    function( G )
    local TheGRP;
    Print("GAP Before call to MinimalStabChain\n");
    TheGRP:=StabChainOp( G, rec( base := [ 1 .. LargestMovedPoint( G ) ] ) );
    Print("GAP After call to MinimalStabChain\n");
    return TheGRP;
end );


#############################################################################
##
#F  SCMinSmaGens(<G>,<S>,<emptyset>,<identity element>,<flag>)
##
##  This function computes a stabilizer chain for a minimal base image and
##  a smallest generating set wrt. this base for a permutation
##  group.
##
##  <G> must be a permutation group and <S> a mutable stabilizer chain for
##  <G> that defines a base <bas>. Let <mbas> the smallest image (OnTuples)
##  of <G>. Then this operation changes <S> to a stabilizer chain wrt.
##  <mbas>.
##  The arguments <emptyset> and <identity element> are needed
##  only for the recursion.
##
##  The function returns a record whose component `gens' is a list whose
##  first element is the smallest element wrt. <bas>. (i.e. an element which
##  maps <bas> to <mbas>. If <flag> is `true', `gens' is  the smallest
##  generating set wrt. <bas>. (If <flag> is `false' this will not be
##  computed.)
InstallGlobalFunction(SCMinSmaGens,function (G,S,bas,pre,flag)
local   Sgens,      # smallest generating system of <S>, result
        gens,       # smallest generating system of <S>, result
	span,	    # <gens>
	stb,	    # Stab_span(bas{1..n-1})
	min,        # minimum in orbit
	nbas,	    # bas+[min]
	rep,        # representative mapping minimal
	gen,        # one generator in <gens>
	orb,        # basic orbit of <S>
	pnt,        # one point in <orb>
	T;          # stabilizer in <S>

  Sgens:=S.generators;
  # handle the anchor case
  if Length(Sgens) = 0  then
      return rec(gens:=[pre],span:=SubgroupNC(Parent(G),[pre]));
  fi;

  # the new ``base'' point is the point to which the current level base
  # point is mapped under pre.
  pnt:=S.orbit[1];

  # find a representative that moves the point as small as possible
  rep:=S.identity;
  min:=Minimum(S.orbit);
  nbas:=Concatenation(bas,[min]);
  while pnt<>min do
    gen:=S.transversal[min];
    rep:=LeftQuotient(gen,rep);
    min:=min^gen;
  od;

  # now we want orbit and stabilizer wrt. this point.
  ConjugateStabChain(S,S,rep,rep);

  # this element will have to be pre-multiplied to all generators
  # generated on this or lower level
  pre:=pre*rep;

  # recursive call to change base below and compute smallest mapping
  # generators there
  gens:=SCMinSmaGens(G,S.stabilizer,nbas,pre,flag);

  # do we want to compute the minimal generating set?
  if flag=false then
    return rec(gens:=[pre]);
  fi;

  span:=gens.span;
  gens:=gens.gens;
  pre:=gens[1]; # the smallest generators is the premul. element

  # get the sorted orbit (the basepoint will be the first point)
  orb := Set( S.orbit );

  # compute the stabilizer in `span' of the first base points (that's the
  # group we're extending at this level)
  stb:=Stabilizer(span,bas,OnTuples);

  # this stabilizer will cover already some points
  SubtractSet( orb, Orbit( stb, S.orbit[1]));

  # handle the points in the orbit
  while Length(orb) > 0  do

    # take the smallest remaining point (coset) and get one representative
    # for it
    pnt := orb[1];
    gen := S.identity;
    while S.orbit[1] ^ gen <> pnt  do
	gen := LeftQuotient( S.transversal[ pnt / gen ], gen );
    od;

    # now change gen by elements in the lower stabilizers  to
    # find the minimal element in its coset.
    T := S.stabilizer;
    while Length(T.generators) <> 0  do
      pnt := Minimum( OnTuples( T.orbit, gen ) );
      while T.orbit[1] ^ gen <> pnt  do
	gen := LeftQuotient( T.transversal[ pnt / gen ], gen );
      od;
      T := T.stabilizer;
    od;

    # pre-multiply with the element mapping the base to the smallest base
    gen:=pre*gen;

    # add this generator to the generators list
    Add( gens, gen );
    #NC is safe -- always use parent(G)
    span:=ClosureSubgroupNC(span,gen);
    stb:=Stabilizer(span,bas,OnTuples);

    # test which cosets we can now cover: reduce orbit
    SubtractSet( orb, Orbit( stb, S.orbit[1]));

  od;

  # return the smallest generating system
  return rec(gens:=gens,span:=span);
end);


#############################################################################
##
#F  LargestElementStabChain(<S>,<id>)
##
InstallGlobalFunction(LargestElementStabChain,function(S,rep)
local   min,    # minimum in orbit
	pnt,    # one point in <orb>
	i,	# loop
	val,	# point image
	gen;    # gen. in transversal

  # handle the anchor case
  if Length(S.generators) = 0  then
      return rep;
  fi;

  # the new ``base'' point is the point to which the current level base
  # point is mapped under pre.
  pnt:=S.orbit[1];

  # find a representative that moves the point as large as possible
  min:=0;
  val:=0;
  for i in S.orbit do
    if i^rep>val then
      min:=i;
      val:=i^rep;
    fi;
  od;

  while pnt<>min do
    gen:=S.transversal[min];
    rep:=LeftQuotient(gen,rep);
    min:=min^gen;
  od;

  # recursive call to change base below and compute smallest mapping
  # generators there
  return LargestElementStabChain(S.stabilizer,rep);

end);

#############################################################################
##
#F  ElementsStabChain(<S>)
##
InstallGlobalFunction(ElementsStabChain,function ( S )
    local   elms,               # element list, result
            stb,                # elements of the stabilizer
            pnt,                # point in the orbit of <S>
            rep;                # inverse representative for that point

    # if <S> is trivial then it is easy
    if Length(S.generators) = 0  then
        elms := [ S.identity ];

    # otherwise
    else

        # start with the empty list
        elms := [];

        # compute the elements of the stabilizer
        stb := ElementsStabChain( S.stabilizer );

        # loop over all points in the orbit
        for pnt  in S.orbit  do

           # add the corresponding coset to the set of elements
           rep := S.identity;
           while S.orbit[1] ^ rep <> pnt  do
                rep := LeftQuotient( S.transversal[pnt/rep], rep );
           od;
           UniteSet( elms, stb * rep );

        od;

   fi;

   # return the result
   return elms;
end);

InstallMethod( ViewObj,"stabilizer chain records", true,
  [ IsRecord ], 0,
function(r)
local sz;
  if not (IsBound(r.stabilizer) and IsBound(r.generators) and
          IsBound(r.orbit) and IsBound(r.identity) and
          IsBound(r.transversal)) then
    TryNextMethod();
  fi;
  sz:= SizeStabChain(r);

  Print("<stabilizer chain record, Base ",BaseStabChain(r),
        ", Orbit length ",Length(r.orbit),", Size: ",sz,">");
end);

#############################################################################
##
#E
