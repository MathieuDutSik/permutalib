#ifndef INCLUDE_PERMUTALIB_ASCENDING_CHAIN_H
#define INCLUDE_PERMUTALIB_ASCENDING_CHAIN_H


#############################################################################
##
##  IntermediateGroup(<G>,<U>)  . . . . . . . . . subgroup of G containing U
##
##  This routine tries to find a subgroup E of G, such that G>E>U. If U is
##  maximal, it returns fail. This is done by using the maximal subgroups machinery or
##  finding minimal blocks for
##  the operation of G on the Right Cosets of U.
##
InstallGlobalFunction( IntermediateGroup, function(G,U)
local o,b,img,G1,c,m,mt,hardlimit,gens,t,k,intersize;

  if U=G then
    return fail;
  fi;

  intersize:=Size(G);
  m:=ValueOption("intersize");
  if IsInt(m) and m<=intersize then
    return fail; # avoid infinite recursion
  fi;

  # use maximals, use `Try` as we call with limiting options
  IsNaturalAlternatingGroup(G);
  IsNaturalSymmetricGroup(G);
  m:=TryMaximalSubgroupClassReps(G:cheap,intersize:=intersize,nolattice);
  if m<>fail and Length(m)>0 then

    m:=Filtered(m,x->Size(x) mod Size(U)=0 and Size(x)>Size(U));
    SortBy(m,x->Size(G)/Size(x));
    
    gens:=SmallGeneratingSet(U);
    for c in m do
      if Index(G,c)<50000 then
        t:=RightTransversal(G,c:noascendingchain); # conjugates
        for k in t do
          if ForAll(gens,x->k*x/k in c) then
            Info(InfoCoset,2,"Found Size ",Size(c));
            # U is contained in c^k
            return c^k;
          fi;
        od;
      else
        t:=DoConjugateInto(G,c,U,true:intersize:=intersize,onlyone:=true);
        if t<>fail and t<>[] then 
          Info(InfoCoset,2,"Found Size ",Size(c));
          return c^(Inverse(t));
        fi;
      fi;
    od;

    Info(InfoCoset,2,"Found no intermediate subgroup ",Size(G)," ",Size(U));
    return fail;
  fi;

  # old code -- obsolete

  c:=ValueOption("refineChainActionLimit");
  if IsInt(c) then
    hardlimit:=c;
  else
    hardlimit:=1000000;
  fi;

  if Index(G,U)>hardlimit then return fail;fi;

  if IsPermGroup(G) and Length(GeneratorsOfGroup(G))>3 then
    G1:=Group(SmallGeneratingSet(G));
    if HasSize(G) then
      SetSize(G1,Size(G));
    fi;
    G:=G1;
  fi;
  o:=ActionHomomorphism(G,RightTransversal(G,U:noascendingchain),
    OnRight,"surjective");
  img:=Range(o);
  b:=Blocks(img,MovedPoints(img));
  if Length(b)=1 then
    return fail;
  else
    b:=StabilizerOfBlockNC(img,First(b,i->1 in i));
    b:=PreImage(o,b);
    return b;
  fi;
end );


#############################################################################
##
#F  RefinedChain(<G>,<c>) . . . . . . . . . . . . . . . .  refine chain links
##
InstallGlobalFunction(RefinedChain,function(G,cc)
local bound,a,b,c,cnt,r,i,j,bb,normalStep,gens,hardlimit,cheap,olda;
  bound:=(10*LogInt(Size(G),10)+1)*Maximum(Factors(Size(G)));
  bound:=Minimum(bound,20000);
  cheap:=ValueOption("cheap")=true;
  c:=ValueOption("refineIndex");
  if IsInt(c) then
    bound:=c;
  fi;

  c:=[];  
  for i in [2..Length(cc)] do  
    Add(c,cc[i-1]);
    if Index(cc[i],cc[i-1]) > bound then
      a:=AsSubgroup(Parent(cc[i]),cc[i-1]);
      olda:=TrivialSubgroup(a);
      while Index(cc[i],a)>bound and Size(a)>Size(olda) do
        olda:=a;
        # try extension via normalizer
        b:=Normalizer(cc[i],a);
        if Size(b)>Size(a) then
         # extension by normalizer surely is a normal step
          normalStep:=true;
          bb:=b;
        else
          bb:=cc[i];
          normalStep:=false;
          b:=Centralizer(cc[i],Centre(a));
        fi;
        if Size(b)=Size(a) or Index(b,a)>bound then
          cnt:=8+2^(LogInt(Index(bb,a),9));
          #if cheap then cnt:=Minimum(cnt,50);fi;
          cnt:=Minimum(cnt,40); # as we have better intermediate
          repeat
            if cnt<20 and not cheap then
              # if random failed: do hard work
              b:=IntermediateGroup(bb,a);
              if b=fail then
                b:=bb;
              fi;
              cnt:=0;
            else
            # larger indices may take more tests...
              Info(InfoCoset,5,"Random");
              repeat
                r:=Random(bb);
              until not(r in a);
              if normalStep then
                # NC is safe
                b:=ClosureSubgroupNC(a,r);
              else
                # self normalizing subgroup: thus every element not in <a>
                     # will surely map one generator out
                j:=0;
                gens:=GeneratorsOfGroup(a);
                repeat
                  j:=j+1;
                until not(gens[j]^r in a);
                r:=gens[j]^r;

                # NC is safe
                b:=ClosureSubgroupNC(a,r);
              fi;
              if Size(b)<Size(bb) then
                Info(InfoCoset,1,"improvement found ",Size(bb)/Size(b));
                bb:=b;
              fi;
              cnt:=cnt-1;
            fi;
          until Index(bb,a)<=bound or cnt<1;
        fi;
        if Index(b,a)>bound and Length(c)>1 then
          bb:=IntermediateGroup(b,c[Length(c)-1]);
          if bb<>fail and Size(bb)>Size(c[Length(c)]) then
            c:=Concatenation(c{[1..Length(c)-1]},[bb],Filtered(cc,x->Size(x)>=Size(b)));
            return RefinedChain(G,c);
          fi;
        fi;

        a:=b;
        if a<>cc[i] then #not upper level
          Add(c,a);
        fi;

      od;
    fi;
  od;
  Add(c,cc[Length(cc)]);
  a:=c[Length(c)];
  for i in [Length(c)-1,Length(c)-2..1] do
    #enforce parent relations
    if not HasParent(c[i]) then
      SetParent(c[i],a);
      a:=c[i];
    else
      a:=AsSubgroup(a,c[i]);
      c[i]:=a;
    fi;
  od;
  return c;
end);




#############################################################################
##
#M  AscendingChainOp(<G>,<pnt>) . . . approximation of
##
InstallMethod( AscendingChainOp, "PermGroup", IsIdenticalObj,
  [IsPermGroup,IsPermGroup],0,
function(G,U)
local s,c,mp,o,i,step,a;
  s:=G;
  c:=[G];
  repeat
    mp:=MovedPoints(s);
    o:=ShallowCopy(OrbitsDomain(s,mp));
    Sort(o,function(a,b) return Length(a)<Length(b);end);
    i:=1;
    step:=false;
    while i<=Length(o) and step=false do
      if not IsTransitive(U,o[i]) then
	Info(InfoCoset,2,"AC: orbit");
	o:=ShallowCopy(OrbitsDomain(U,o[i]));
	Sort(o,function(a,b) return Length(a)<Length(b);end);
	# union of same length -- smaller index
	a:=Union(Filtered(o,x->Length(x)=Length(o[1])));
	if Length(a)=Sum(o,Length) then
	  a:=Set(o[1]);
	fi;
	s:=Stabilizer(s,a,OnSets);
	step:=true;
      elif Index(G,U)>NrMovedPoints(U) 
	  and IsPrimitive(s,o[i]) and not IsPrimitive(U,o[i]) then
	Info(InfoCoset,2,"AC: blocks");
	s:=Stabilizer(s,Set(List(MaximalBlocks(U,o[i]),Set)),
                      OnSetsDisjointSets);
	step:=true;
      else
	i:=i+1;
      fi;
    od;
    if step then
      Add(c,s);
    fi;
  until step=false or Index(s,U)=1; # we could not refine better
  if Index(s,U)>1 then
    Add(c,U);
  fi;
  Info(InfoCoset,2,"Indices",List([1..Length(c)-1],i->Index(c[i],c[i+1])));
  return RefinedChain(G,Reversed(c));
end);



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
    Info(InfoCoset,1,"DoubleCosetFlip");
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
  #c:=AscendingChain(G,a:refineChainActionLimit:=Index(G,a));
  Print(NullMat(5));
  
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
      Info(InfoCoset,1,"flip to get better chain");
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
          #Print(maxidx(c),obj,Length(Orbit(G,obj,act))," ",
          #          Length(Orbit(a,obj,act)),"\n");
          G1:=Stabilizer(G,obj,act);
          if Index(G,G1)<maxidx(c) then
            a1:=Stabilizer(a,obj,act);
          else
            a1:=G;
          fi;
        fi;
        Info(InfoCoset,4,"attempt up step ",obj," index:",Size(a)/Size(a1));
        if Index(G,G1)<maxidx(c) and Index(a,a1)<=uplimit and (
          maxidx(c)>avoidlimit or Size(a1)>Size(c[1])) then
          c1:=AscendingChain(G1,a1:refineIndex:=avoidlimit,
                                   refineChainActionLimit:=actlimit,
                                   indoublecoset);
          if maxidx(c1)<maxidx(c) then
            c:=Concatenation(c1,[G]);
            cano:=true;
            Info(InfoCoset,1,"improved chain with up step ",obj,
            " index:",Size(a)/Size(a1)," maxidx=",maxidx(c));
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
            Info(InfoCoset,3,"Try maximal of Indices ",Index(G,i),":",
              Index(a,p));
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
        #Print("j=",Size(j),"\n");
        t:=AscendingChain(G,j:refineIndex:=avoidlimit,
                              refineChainActionLimit:=actlimit,indoublecoset);
        Info(InfoCoset,4,"maxidx ",Index(a,j)," yields ",maxidx(t),": ",
          List(t,Size));
        if maxidx(t)<maxidx(c) and (maxidx(c)>badlimit or
          # only increase up-step if index gets better by extra index
          (maxidx(c)>maxidx(t)*Size(c[1])/Size(t[1])) ) then
          c:=t;
          cano:=true;
          Info(InfoCoset,1,"improved chain with up step index:",
                Size(a)/Size(j));
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
  if ValueOption("noupfuse")=true then cano:=false;fi;

  Info(InfoCoset,1,"Chosen series is ",List(c,Size));
  #if ValueOption("indoublecoset")<>true then Error("GNASH");fi;

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
    if (step>1 or cano) and Length(t)<badlimit and IsPermGroup(G) and 
      not normal then
      # in this case, we can beneficially compute the action once and then use
      # homomorphism methods to obtain the permutation image
      Info(InfoCoset,2,"using perm action on step ",step,": ",Length(t));
      hom:=Subgroup(G,SmallGeneratingSet(a1));
      hom:=ActionHomomorphism(hom,t,OnRight,"surjective");
    else
      hom:=fail;
    fi;
    homs[step]:=hom;
  od;

  omi:=[];
  omiz:=[];

  for step in [1..Length(c)-1] do
    a1:=c[Length(c)-step+1];
    a2:=c[Length(c)-step];
    normal:=IsNormal(a1,a2);
    indx:=Index(a1,a2);
    if normal then
      Info(InfoCoset,1,"Normal Step :",indx,": ",Length(r)," double cosets");
    else
      Info(InfoCoset,1,"Step :",indx,": ",Length(r)," double cosets");
    fi;
    

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
      Info(InfoCoset,4,"number ",nu,", |stab|=",Size(lst));
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
        cnt:=cnt-1;

        # compute orbit and stabilizers for the next step
        # own Orbitalgorithm and stabilizer computation

        #while blist[posi] do posi:=posi+1;od;
        posi:=Position(blist,false,posi);
        ps:=posi;
        blist[ps]:=true;
        bsz:=bsz-1;
        e:=t[ps];
        mop:=1;
        mo:=ps;

        rep := [ One(b) ];
        st := TrivialSubgroup(lst);

        o:=[ps];
        if cano or compst then
          oi:=[];
          oi[ps]:=1; # reverse index
        fi;
        orbcnt:=orbcnt+1;

        i:=1;
        while i<=Length(o) 
          # will not grab if nonreg,. orbiut and stabilizer not computed,
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
                #NC is safe (initializing as TrivialSubgroup(G)
                st := ClosureSubgroupNC(st,rep[i]*lstgens[j]/rep[oi[ps]]);
              fi;
            else
              # new image
              blist[ps]:=true;
              bsz:=bsz-1;
              Add(o,ps);
              if cano or compst then
                Add(rep,rep[i]*lstgens[j]);
                if cano and ps<mo then
                  mo:=ps;
                  mop:=Length(rep);
                fi;
                oi[ps]:=Length(o);
              fi;
            fi;
          od;
          i:=i+1;
        od;
        Info(InfoCoset,5,"|o|=",Length(o));

        ep:=e*rep[mop]*p;
        Add(nr,ep);

        if compst then
          st:=st^rep[mop];
          Add(nstab,st);
        fi;

        if cano and step=1 and not normal then 
          Add(omi,mo);
          Add(omiz,Length(o));
        fi;

        siz:=sifa*Length(o); #order

        if unten then
          if flip then
            Add(dcs,[ep^(-1),siz]);
          else
            Add(dcs,[ep,siz]);
          fi;
        fi;

      od;
      Info(InfoCoset,4,"Get ",orbcnt," orbits");

      if normal then
        # in the normal case, we can obtain the other orbits easily via
        # the orbit theorem (same stabilizer)
        rt:=RightTransversal(lst,st);
        Assert(1,Length(rt)=Length(o));

        while bsz>0 do
          ps:=Position(blist,false);
          e:=t[ps];
          blist[ps]:=true;

          ep:=e*p;
          mo:=ep;
          mop:=ps;
          # tick off the orbit
          for i in rt do
            #ps:=PositionCanonical(t,e*p*i/p);
            j:=ep*i/p;
            ps:=PositionCanonical(t,ep*i/p);
            if cano then
              if ps<mop then
                mop:=ps;
                mo:=j;
              fi;
            fi;
            blist[ps]:=true;
          od;
          bsz:=bsz-Length(rt);

          Add(nr,mo);
          Add(nstab,st);

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
    Info(InfoCoset,3,Length(r)," double cosets so far.");
  od;

  if cano then
    # do the final up step

    IsSSortedList(omi);

    # canonization fct
    canrep:=function(x)
    local stb, p, pinv, t, hom,ps, mop, mo, o, oi, rep, st, lstgens, lstgensop,
          i, img, step, j,calcs;
      stb:=b;
      p:=One(G);
      for step in [1..Length(c)-1] do
        calcs:=step<Length(c)-1;
        pinv:=p^-1;
        t:=tra[step];
        hom:=homs[step];
        # orbit-stabilizer algorithm
        ps:=PositionCanonical(t,x);
        mop:=1;
        mo:=ps;
        o:=[ps];
        oi:=[];
        oi[ps]:=1;
        rep:=[One(stb)];
        st:=TrivialSubgroup(b);

        lstgens:=GeneratorsOfGroup(stb);
        if Length(lstgens)>4 and
          Length(lstgens)/(AbelianRank(stb)+1)*2>5 then
          lstgens:=SmallGeneratingSet(stb);
        fi;
        lstgensop:=List(lstgens,i->i^pinv); # conjugate generators: operation

        if hom<>fail then
          lstgensop:=List(lstgensop,i->Image(hom,i));
        fi;
        i:=1;
        while i<=Length(o) do
          for j in [1..Length(lstgensop)] do
            if hom=fail then
              img:=t[o[i]]*lstgensop[j];
              ps:=PositionCanonical(t,img);
            else
              ps:=o[i]^lstgensop[j];
            fi;
            if IsBound(oi[ps]) then
              # known image

              # if there is only one orbit on the top step, we know the
              # stabilizer!
              if calcs then
                #NC is safe (initializing as TrivialSubgroup(G)
                st := ClosureSubgroupNC(st,rep[i]*lstgens[j]/rep[oi[ps]]);
                if Size(st)*Length(o)=Size(b) then i:=Length(o);fi;
              fi;
              #fi;
            else
              Add(o,ps);
              Add(rep,rep[i]*lstgens[j]);
              if ps<mo then
                mo:=ps;
                mop:=Length(rep);
                if step=1 and mo in omi then
                  #Print("found\n");
                  if Size(st)*omiz[Position(omi,mo)]=Size(stb) then
                    # we have the minimum and the right stabilizer: break
                    #Print("|Orbit|=",Length(o),
                    #" of ",omiz[Position(omi,mo)]," min=",mo,"\n");
                    i:=Length(o);
                  fi;
                fi;
              fi;
              oi[ps]:=Length(o);
              if Size(st)*Length(o)=Size(b) then i:=Length(o);fi;
            fi;
          od;
          i:=i+1;
        od;
        
        if calcs then
          stb:=st^(rep[mop]);
        fi;
        #if HasSmallGeneratingSet(st) then
        #  SetSmallGeneratingSet(stb,List(SmallGeneratingSet(st),x->x^rep[mop]));
        #fi;

        #else
        #  stb:=omis;
        #fi;
        x:=x*(rep[mop]^pinv)/t[mo];
        p:=t[mo]*p;
        #Print("step ",step," |Orbit|=",Length(o),"nmin=",mo,"\n");

        #if ForAny(GeneratorsOfGroup(stb),
        #     i->not x*p*i/p in t!.subgroup) then
        #     Error("RRR");
        #fi;

      od;
      return p;
    end;

    # now fuse orbits under the left action of a
    indx:=Index(a,a2);
    Info(InfoCoset,2,"fusion index ",indx);
    #t:=Filtered(RightTransversal(a,a2),x->not x in a2);
    t:=RightTransversal(a,a2);
    sifa:=Size(a2)*Size(b);

    # cluster according to A-double coset sizes and C lengths
    #sizes:=List(r,x->Size(a)*Size(b)/Size(Intersection(b,a^x)));
    hom:=ActionHomomorphism(a,t,OnRight,"surjective");
    sizes:=[];
    for i in [1..Length(r)] do
      lr:=Intersection(a,b^(r[i]^-1));
      # size of double coset and 
      Add(sizes,[Size(a)*Size(b)/Size(lr),
                 Length(OrbitsDomain(Image(hom,lr),[1..Length(t)],OnPoints))]);
    od;
    ps:=ShallowCopy(sizes);
    sizes:=Set(sizes); # sizes corresponding to clusters
    cluster:=List(sizes,s->Filtered([1..Length(r)],x->ps[x]=s));

    # now process per cluster
    for i in [1..Length(sizes)] do
      sel:=cluster[i];
      lr:=r{sel};
      lstabs:=stabs{sel};
      SortParallel(lr,lstabs); # quick find
      IsSSortedList(lr);
      ssizes:=List(lstabs,x->sifa/Size(x));
      num:=Sum(ssizes)/sizes[i][1]; # number of double cosets to be created
      if num>1 and sizes[i][1]/Size(a)<=10*Index(a,a2)^2 then
        # fuse orbits together
        lr:=List(lr,x->CanonicalRightCosetElement(a,x));
        o:=DCFuseSubgroupOrbits(G,b,lr,function(r,g)
            return CanonicalRightCosetElement(a,r*g);
          end,1000,num);
        for j in o do
          # record double coset
          if flip then
            Add(dcs,[lr[j[1]]^(-1),sizes[i][1]]);
          else
            Add(dcs,[lr[j[1]],sizes[i][1]]);
          fi;
          Info(InfoCoset,2,"orbit fusion ",Length(dcs)," orblen=",Length(j));
        od;
        lr:=[];lstabs:=[];
      else
        while num>1 do
          # take first representative as rep for double coset
          #stab:=Intersection(b,a^lr[1]);

          # check how does its double coset a*lr[1]*b split up into a2-DC's
          o:=OrbitsDomain(Image(hom,Intersection(a,b^(lr[1]^-1))),
                [1..Length(t)],OnPoints);

          # identify which of the a2-cosets they are they are (so we can
          # remove them)
          o:=List(o,x->Position(lr,canrep(t[x[1]]*lr[1])));

          # record double coset
          if flip then
            Add(dcs,[lr[1]^(-1),sizes[i][1]]);
          else
            Add(dcs,[lr[1],sizes[i][1]]);
          fi;
          sel:=Difference([1..Length(lr)],o);
          lr:=lr{sel};lstabs:=lstabs{sel};
          Info(InfoCoset,2,"new fusion ",Length(dcs)," orblen=",Length(o),
              " remainder ",Length(lr));

          num:=num-1;
        od;

        # remainder must be a single double coset
        if flip then
          Add(dcs,[lr[1]^(-1),sizes[i][1]]);
        else
          Add(dcs,[lr[1],sizes[i][1]]);
        fi;
        Info(InfoCoset,2,"final fusion ",Length(dcs)," orblen=",Length(lr),
            " remainder ",0);

      fi;

    od;
  fi;

  if AssertionLevel()>2 then
    # test
    bsz:=Size(G);
    t:=[];
    if flip then
      # flip back
      c:=a;
      a:=b;
      b:=c;
    fi;
    for i in dcs do
      bsz:=bsz-i[2];
      if AssertionLevel()>0 then
        r:=CanonicalRightCosetElement(a,i[1]);
        if ForAny(t,j->r in RepresentativesContainedRightCosets(j)) then
          Error("duplicate!");
        fi;
      fi;
      r:=DoubleCoset(a,i[1],b);
      if AssertionLevel()>0 and Size(r)<>i[2] then
        Error("size error!");
      fi;
      Add(t,r);
    od;
    if bsz<>0 then
      Error("number");
    fi;
  fi;

  return dcs;
end);






#endif
