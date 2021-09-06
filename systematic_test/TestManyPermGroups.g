MyRemoveFileIfExist:=function(FileName)
    if IsExistingFile(FileName) then
        RemoveFile(FileName);
    fi;
end;




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



MergeTransitiveGroup:=function(GRP1, GRP2)
    local LMoved1, LMoved2, n1, n2, LGens, eGen, eList;
    LMoved1:=MovedPoints(GRP1);
    LMoved2:=MovedPoints(GRP2);
    n1:=Maximum(LMoved1);
    n2:=Maximum(LMoved2);
    LGens:=[];
    for eGen in GeneratorsOfGroup(GRP1)
    do
        Add(LGens, eGen);
    od;
    for eGen in GeneratorsOfGroup(GRP2)
    do
        eList:=Concatenation([1..n1], List([1..n2], x->n1 + OnPoints(x, eGen)));
        Add(LGens, PermList(eList));
    od;
    return Group(LGens);
end;




GetListCandidateGroups:=function()
    local ListGroup, n, ePow, siz, eNB, i, j;
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
    # Some non-transitive cases
    for i in [2..6]
    do
        for j in [2..6]
        do
            Add(ListGroup, MergeTransitiveGroup(SymmetricGroup(i), SymmetricGroup(j)));
            Add(ListGroup, MergeTransitiveGroup(SymmetricGroup([3..i+2]), SymmetricGroup(j)));
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
    MyRemoveFileIfExist(FileName);
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
    eBinary:="../src_gap/GapStabilizerOnSet";
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
    MyRemoveFileIfExist(FileName);
    MyRemoveFileIfExist(FileErr);
    MyRemoveFileIfExist(FileRes);
end;





TestSpecificGroupSet_Canonical:=function(nbMov, eGRP, eSet)
    local eDir, FileName, output, LGen, eGen, iMov, eImg, pos, eVal, eBinary, FileErr, FileRes, eCommand, eSet1, eSet2, test;
    Print("Treating one pair Group/Set eSet=", eSet, "\n");
    if Maximum(eSet) > nbMov then
        Error("The eSet is too large");
    fi;
    eDir:="/tmp/DebugStabOnSets_datarun/";
    eCommand:=Concatenation("mkdir -p ", eDir);
    Exec(eCommand);
    #
    FileName:=Concatenation(eDir, "Input");
    MyRemoveFileIfExist(FileName);
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
    eBinary:="../src_gap/GapCanonicalImage";
    FileErr:=Concatenation(eDir, "CppError");
    FileRes:=Concatenation(eDir, "GapOutput");
    eCommand:=Concatenation(eBinary, " ", FileName, " 2> ", FileErr, " ", FileRes);
    Exec(eCommand);
    #
    eSet1:=CanonicalImage(eGRP, eSet, OnSets);
    eSet2:=ReadAsFunction(FileRes)();
    test:=eSet1=eSet2;
    if test=false then
        Error("Found some error. Please debug");
    fi;
    MyRemoveFileIfExist(FileName);
    MyRemoveFileIfExist(FileErr);
    MyRemoveFileIfExist(FileRes);
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
    MyRemoveFileIfExist(FileName);
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
    eBinary:="../src_gap/GapRepresentativeActionOnSet";
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
    MyRemoveFileIfExist(FileName);
    MyRemoveFileIfExist(FileErr);
    MyRemoveFileIfExist(FileRes);
end;




TestPropertiesGroup:=function(nbMov, eGRP)
    local eDir, FileName, output, LGen, eGen, iMov, eImg, pos, eVal, eBinary, FileErr, FileRes, eCommand, Result1, Result2, test;
    Print("Treating a group for its properties\n");
    eDir:="/tmp/DebugGroupProperties_datarun/";
    eCommand:=Concatenation("mkdir -p ", eDir);
    Exec(eCommand);
    #
    FileName:=Concatenation(eDir, "Input");
    MyRemoveFileIfExist(FileName);
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
    CloseStream(output);
    #
    eBinary:="../src_gap/GroupProperties";
    FileErr:=Concatenation(eDir, "CppError");
    FileRes:=Concatenation(eDir, "GapOutput");
    eCommand:=Concatenation(eBinary, " ", FileName, " ", FileRes, " 2> ", FileErr);
#    Print("eCommand=", eCommand, "\n");
    Exec(eCommand);
    #
    Result1:=rec(IsPrimitive:=IsPrimitive(eGRP), IsTransitive:=IsTransitive(eGRP), IsCommutative:=IsCommutative(eGRP), IsCyclic:=IsCyclic(eGRP));
    Result2:=ReadAsFunction(FileRes)();
    test:=Result1 = Result2;
    if test=false then
        Error("Found some error. Please debug");
    fi;
    MyRemoveFileIfExist(FileName);
    MyRemoveFileIfExist(FileErr);
    MyRemoveFileIfExist(FileRes);
end;




TestSmallGeneratingSet:=function(nbMov, eGRP)
    local eDir, FileName, output, LGen, eGen, iMov, eImg, pos, eVal, eBinary, FileErr, FileRes, eCommand, Result1, Result2, test;
    Print("Checking SmallGeneratingSet feature\n");
    eDir:="/tmp/DebugSmallGeneratingSet_datarun/";
    eCommand:=Concatenation("mkdir -p ", eDir);
    Exec(eCommand);
    #
    FileName:=Concatenation(eDir, "Input");
    MyRemoveFileIfExist(FileName);
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
    CloseStream(output);
    #
    eBinary:="../src_gap/GapSmallGeneratingSet";
    FileErr:=Concatenation(eDir, "CppError");
    FileRes:=Concatenation(eDir, "GapOutput");
    eCommand:=Concatenation(eBinary, " ", FileName, " ", FileRes, " 2> ", FileErr);
#    Print("eCommand=", eCommand, "\n");
    Exec(eCommand);
    #
    LGen:=ReadAsFunction(FileRes)();
    if Group(LGen) <> eGRP then
        Error("Found some error. Please debug");
    fi;
    MyRemoveFileIfExist(FileName);
    MyRemoveFileIfExist(FileErr);
    MyRemoveFileIfExist(FileRes);
end;






TestDerivedSubgroup:=function(nbMov, eGRP)
    local eDir, FileName, output, LGen, eGen, iMov, eImg, pos, eVal, eBinary, FileErr, FileRes, eCommand, Result1, Result2, test, eGRP_der;
    Print("Checking DerivedSubgroup feature\n");
    eDir:="/tmp/DebugDerivedSubgroup_datarun/";
    eCommand:=Concatenation("mkdir -p ", eDir);
    Exec(eCommand);
    #
    FileName:=Concatenation(eDir, "Input");
    MyRemoveFileIfExist(FileName);
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
    CloseStream(output);
    #
    eBinary:="../src_gap/GapDerivedSubgroup";
    FileErr:=Concatenation(eDir, "CppError");
    FileRes:=Concatenation(eDir, "GapOutput");
    eCommand:=Concatenation(eBinary, " ", FileName, " ", FileRes, " 2> ", FileErr);
    Exec(eCommand);
    #
    eGRP_der:=ReadAsFunction(FileRes)();
    if eGRP_der <> DerivedSubgroup(eGRP) then
        Error("Found some error. Please debug");
    fi;
    MyRemoveFileIfExist(FileName);
    MyRemoveFileIfExist(FileErr);
    MyRemoveFileIfExist(FileRes);
end;



TestCentreSubgroup:=function(nbMov, eGRP)
    local eDir, FileName, output, LGen, eGen, iMov, eImg, pos, eVal, eBinary, FileErr, FileRes, eCommand, Result1, Result2, test, eGRP_cent;
    Print("Checking CentreSubgroup feature\n");
    eDir:="/tmp/DebugCentreSubgroup_datarun/";
    eCommand:=Concatenation("mkdir -p ", eDir);
    Exec(eCommand);
    #
    FileName:=Concatenation(eDir, "Input");
    MyRemoveFileIfExist(FileName);
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
    CloseStream(output);
    #
    eBinary:="../src_gap/GapCentreSubgroup";
    FileErr:=Concatenation(eDir, "CppError");
    FileRes:=Concatenation(eDir, "GapOutput");
    eCommand:=Concatenation(eBinary, " ", FileName, " ", FileRes, " 2> ", FileErr);
    Exec(eCommand);
    #
    eGRP_cent:=ReadAsFunction(FileRes)();
    if eGRP_cent <> Centre(eGRP) then
        Error("Found some error. Please debug");
    fi;
    MyRemoveFileIfExist(FileName);
    MyRemoveFileIfExist(FileErr);
    MyRemoveFileIfExist(FileRes);
end;



TestCentralizerElt:=function(nbMov, eGRP, g)
    local eDir, FileName, output, LGen, eGen, iMov, eImg, pos, eVal, eBinary, FileErr, FileRes, eCommand, Result1, Result2, test, PrtElement, eCentrElt;
    Print("Checking CentralizerElt feature\n");
    eDir:="/tmp/DebugCentralizerElt_datarun/";
    eCommand:=Concatenation("mkdir -p ", eDir);
    Exec(eCommand);
    #
    FileName:=Concatenation(eDir, "Input");
    MyRemoveFileIfExist(FileName);
    output:=OutputTextFile(FileName, true);
    LGen:=GeneratorsOfGroup(eGRP);
    AppendTo(output, Length(LGen), " ", nbMov, "\n");
    PrtElement:=function(u)
        for iMov in [1..nbMov]
        do
            eImg:=OnPoints(iMov, u);
            AppendTo(output, " ", eImg-1);
        od;
        AppendTo(output, "\n");
    end;
    for eGen in LGen
    do
        PrtElement(eGen);
    od;
    PrtElement(g);
    CloseStream(output);
    #
    eBinary:="../src_gap/GapCentralizerElt";
    FileErr:=Concatenation(eDir, "CppError");
    FileRes:=Concatenation(eDir, "GapOutput");
    eCommand:=Concatenation(eBinary, " ", FileName, " ", FileRes, " 2> ", FileErr);
    Exec(eCommand);
    #
    eCentrElt:=ReadAsFunction(FileRes)();
    if eCentrElt <> Centralizer(eGRP, g) then
        Error("Found some error. Please debug");
    fi;
    MyRemoveFileIfExist(FileName);
    MyRemoveFileIfExist(FileErr);
    MyRemoveFileIfExist(FileRes);
end;





TestSpecificGroup:=function(method, size_opt, nbMov, eGRP)
    local iMov, sizSet, i, eSet, fSet, eElt, iter, g, idx1, idx2;
    Print("ListGens(eGRP)=", GeneratorsOfGroup(eGRP), "\n");
    for iMov in [1..5]
    do
        if size_opt=1 then
            sizSet:=1;
        fi;
        if size_opt=2 then
            if nbMov<4 then
                sizSet:=2;
            else
                sizSet:=Random([2..nbMov-2]);
            fi;
        fi;
        if size_opt=3 then # We want the size to be at most nMov / 2
            if nbMov<4 then
                return;
            fi;
            sizSet:=Random([2..Int(nbMov/2)]);
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
        if method="canonical" then
            for i in [1..5]
            do
                eSet:=Local_RandomSubset([1..nbMov], sizSet);
                TestSpecificGroupSet_Canonical(nbMov, eGRP, eSet);
            od;
        fi;
    od;
    if method="properties" then
        TestPropertiesGroup(nbMov, eGRP);
    fi;
    if method="smallgeneratingset" then
        TestSmallGeneratingSet(nbMov, eGRP);
    fi;
    if method="derivedsubgroup" then
        TestDerivedSubgroup(nbMov, eGRP);
    fi;
    if method="centresubgroup" then
        TestCentreSubgroup(nbMov, eGRP);
    fi;
    if method="centralizerelt" then
        for iter in [1..50]
        do
            g:=Random(eGRP);
            TestCentralizerElt(nbMov, eGRP, g);
            TestCentralizerElt(2*nbMov, eGRP, g);
            #
            g:=Random(SymmetricGroup(2*nbMov));
            TestCentralizerElt(2*nbMov, eGRP, g);
            g:=Random(SymmetricGroup(nbMov));
            TestCentralizerElt(nbMov, eGRP, g);
            #
            idx1:=Random([1..nbMov]);
            idx2:=Random([1..nbMov]);
            if idx1<>idx2 then
                g:=(idx1,idx2);
                TestCentralizerElt(nbMov, eGRP, g);
            fi;
        od;
    fi;
end;



TestAllGroups:=function(method, size_opt)
    local ListGroups, nbGroups, eGRP, nMov, iGRP;
    ListGroups:=GetListCandidateGroups();
#    ListGroups:=Filtered(ListGroups, x->IsSolvable(x)=false);
    nbGroups:=Length(ListGroups);
    Print("|ListGroups|=", nbGroups, "\n");

    for iGRP in [1..nbGroups]
    do
        eGRP:=ListGroups[iGRP];
        nMov:=LargestMovedPoint(eGRP);
        Print("    iGRP=", iGRP, " / ", nbGroups, " nMov=", nMov, "\n");
        TestSpecificGroup(method, size_opt, nMov, eGRP);
    od;
end;


WriteAllGroupsInFile:=function(eFile)
    local ListGroups, output, eGRP, LGen, nbMov, eGen, iMov, eImg;
    ListGroups:=GetListCandidateGroups();
    if IsExistingFile(eFile) then
        RemoveFile(eFile);
    fi;
    output:=OutputTextFile(eFile, true);
    AppendTo(output, Length(ListGroups), "\n");
    for eGRP in ListGroups
    do
        LGen:=GeneratorsOfGroup(eGRP);
        nbMov:=LargestMovedPoint(eGRP);
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
    od;
    CloseStream(output);
end;


TestAllGroups("centralizerelt", 1);
#TestAllGroups("derivedsubgroup", 1);
#TestAllGroups("centresubgroup", 1);
#TestAllGroups("smallgeneratingset", 1);
#TestAllGroups("properties", 1);
#TestAllGroups("stabilizer", 1);
#TestAllGroups("equivalence", 2);
#TestAllGroups("canonical", 3);

#WriteAllGroupsInFile("AllFileTest");

