Local_RandomSubset:=function(eSet, k)
    local i, sSet, V, h;
    sSet:=[];
    V:=ListWithIdenticalEntries(Length(eSet), 1);
    for i in [1..k]
    do
        while(true)
        do
            h:=Random([1..Length(eSet)]);
            if V[h]=1 then
                V[h]:=0;
                Add(sSet, eSet[h]);
                break;
            fi;
        od;
    od;
    return Set(sSet);
end;



GetListCandidateGroups:=function()
    local ListGroup, n, ePow, siz, eNB, i;
    ListGroup:=[];
    for n in [3..8]
    do
        Add(ListGroup, SymmetricGroup(n));
        Add(ListGroup, AlternatingGroup(n));
    od;
    for n in [9, 10, 11, 12, 21, 22, 23, 24]
    do
        Add(ListGroup, MathieuGroup(n));
    od;
    for ePow in [4..50]
    do
        if Length(Set(FactorsInt(ePow)))=1 then
            Add(ListGroup, PSL(2,ePow));
        fi;
    od;
    for siz in [4..15]
    do
        eNB:=NrTransitiveGroups(siz);
        for i in [1..eNB]
        do
            Add(ListGroup, TransitiveGroup(siz, i));
        od;
    od;
    return ListGroup;
end;


TestSpecificGroupSet_Stabilizer:=function(nbMov, eGRP, eSet)
    local eDir, FileName, output, LGen, eGen, iMov, eImg, pos, eVal, eBinary, FileErr, FileRes, eCommand, eStab1, eStab2, test;
    Print("Treating one pair Group/Set eSet=", eSet, "\n");
    if Maximum(eSet) > nbMov then
        Error("The eSet is too large");
    fi;
    eDir:="/tmp/DebugStabOnSets_datarun/";
    eCommand:=Concatenation("mkdir -p ", eDir);
    Exec(eCommand);
    #
    FileName:=Concatenation(eDir, "Input");
    RemoveFileIfExist(FileName);
    output:=OutputTextFile(FileName, true);
    LGen:=GeneratorsOfGroup(eGRP);
    AppendTo(output, Length(LGen), " ", nbMov, "\n");
    for eGen in LGen
    do
        for iMov in [1..nbMov]
        do
            eImg:=OnPoints(iMov, eGen);
            AppendTo(output, " ", eImg-1);
        od;
        AppendTo(output, "\n");
    od;
    #
    for iMov in [1..nbMov]
    do
        pos:=Position(eSet, iMov);
        if pos=fail then
            eVal:=0;
        else
            eVal:=1;
        fi;
        AppendTo(output, " ", eVal);
    od;
    AppendTo(output, "\n");
    CloseStream(output);
    #
    eBinary:="/home/mathieu/GITall/GIT/permutalib/src_gap/GapStabilizerOnSet";
    FileErr:=Concatenation(eDir, "CppError");
    FileRes:=Concatenation(eDir, "GapOutput");
    eCommand:=Concatenation(eBinary, " ", FileName, " 2> ", FileErr, " ", FileRes);
    Exec(eCommand);
    #
    eStab1:=Stabilizer(eGRP, eSet, OnSets);
    eStab2:=ReadAsFunction(FileRes)();
    test:=eStab1=eStab2;
    if test=false then
        Error("Found some error. Please debug");
    fi;
    RemoveFileIfExist(FileName);
    RemoveFileIfExist(FileErr);
    RemoveFileIfExist(FileRes);
end;






TestSpecificGroupSet_Equivalence:=function(nbMov, eGRP, eSet, fSet)
    local eDir, FileName, output, LGen, eGen, iMov, eImg, pos, eVal, eBinary, FileErr, FileRes, eCommand, eTest1, eTest2;
    Print("Treating one pair Group/Set eSet=", eSet, " fSet=", fSet, "\n");
    if Maximum(eSet) > nbMov then
        Error("The eSet is too large");
    fi;
    eDir:="/tmp/DebugStabOnSets_datarun/";
    eCommand:=Concatenation("mkdir -p ", eDir);
    Exec(eCommand);
    #
    FileName:=Concatenation(eDir, "Input");
    RemoveFileIfExist(FileName);
    output:=OutputTextFile(FileName, true);
    LGen:=GeneratorsOfGroup(eGRP);
    AppendTo(output, Length(LGen), " ", nbMov, "\n");
    for eGen in LGen
    do
        for iMov in [1..nbMov]
        do
            eImg:=OnPoints(iMov, eGen);
            AppendTo(output, " ", eImg-1);
        od;
        AppendTo(output, "\n");
    od;
    #
    for iMov in [1..nbMov]
    do
        pos:=Position(eSet, iMov);
        if pos=fail then
            eVal:=0;
        else
            eVal:=1;
        fi;
        AppendTo(output, " ", eVal);
    od;
    AppendTo(output, "\n");
    #
    for iMov in [1..nbMov]
    do
        pos:=Position(fSet, iMov);
        if pos=fail then
            eVal:=0;
        else
            eVal:=1;
        fi;
        AppendTo(output, " ", eVal);
    od;
    AppendTo(output, "\n");
    CloseStream(output);
    #
    eBinary:="/home/mathieu/GITall/GIT/permutalib/src_gap/GapRepresentativeActionOnSet";
    FileErr:=Concatenation(eDir, "CppError");
    FileRes:=Concatenation(eDir, "GapOutput");
    eCommand:=Concatenation(eBinary, " ", FileName, " 2> ", FileErr, " ", FileRes);
    Exec(eCommand);
    #
    eTest1:=RepresentativeAction(eGRP, eSet, fSet, OnSets);
    eTest2:=ReadAsFunction(FileRes)();
    if eTest1=fail then
        if eTest2<>fail then
            Error("Found some error in TestSpecificGroupSet_Equivalence, case 1\n");
        fi;
    else
        if fSet<>OnSets(eSet, eTest2) then
            Error("Found some error in TestSpecificGroupSet_Equivalence, case 2\n");
        fi;
    fi;
    RemoveFileIfExist(FileName);
    RemoveFileIfExist(FileErr);
    RemoveFileIfExist(FileRes);
end;



TestSpecificGroup:=function(method, nbMov, eGRP)
    local iMov, sizSet, i, eSet, fSet, eElt;
    Print("ListGens(eGRP)=", GeneratorsOfGroup(eGRP), "\n");
    for iMov in [1..5]
    do
        if nbMov<4 then
            sizSet:=2;
        else
            sizSet:=Random([2..nbMov-2]);
        fi;
        if method="stabilizer" then
            for i in [1..5]
            do
                eSet:=Local_RandomSubset([1..nbMov], sizSet);
                TestSpecificGroupSet_Stabilizer(nbMov, eGRP, eSet);
            od;
        fi;
        if method="equivalence" then
            for i in [1..5]
            do
                eSet:=Local_RandomSubset([1..nbMov], sizSet);
                fSet:=Local_RandomSubset([1..nbMov], sizSet);
                TestSpecificGroupSet_Equivalence(nbMov, eGRP, eSet, fSet);
            od;
            for i in [1..5]
            do
                eSet:=Local_RandomSubset([1..nbMov], sizSet);
                eElt:=Random(eGRP);
                fSet:=OnSets(eSet, eElt);
                TestSpecificGroupSet_Equivalence(nbMov, eGRP, eSet, fSet);
            od;
        fi;
    od;
end;



TestAllGroups:=function(method)
    local ListGroups, nbGroups, eGRP, nMov, iGRP;
    ListGroups:=GetListCandidateGroups();
#    ListGroups:=Filtered(ListGroups, x->IsSolvable(x)=false);
    nbGroups:=Length(ListGroups);
    Print("|ListGroups|=", nbGroups, "\n");

    for iGRP in [1..nbGroups]
    do
        eGRP:=ListGroups[iGRP];
        nMov:=LargestMovedPoint(eGRP);
        if Position([1,2], iGRP) = fail then
            Print("    iGRP=", iGRP, " / ", nbGroups, " nMov=", nMov, "\n");
            TestSpecificGroup(method, nMov, eGRP);
        fi;
    od;
end;


#TestAllGroups("stabilizer");
TestAllGroups("equivalence");
