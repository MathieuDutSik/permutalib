GetListCandidateGroups:=function()
    local ListGroup;
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
    for ePow in [1..50]
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


TestSpecificGroupSet:=function(nbMov, eGRP, eSet)
    local output;
    FileName:="/tmp/Input";
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
    FileRes:="/tmp/GapOutput";
    eCommand:=Concatenation(eBinary, " ", FileName, " ", FileRes);
    Exec(eCommand);
    #
    eStab1:=Stabilizer(eGRP, eSet, OnSets);
    eStab2:=ReadAsFunction(FileRes)();
    test:=eStab1=eStab2;
    if test=false then
        Error("Found some error. Please debug");
    fi;
end;



TestSpecificGroup:=function(nbMov, eGRP)
    for iMov in [1..5]
    do
        sizSet:=Random([2..nMov]);
        for i in [1..5]
        do
            eSet:=Local_RandomSubset([1..nbMov], sizSet);
            TestSpecificGroupSet(nbMov, eGRP, eSet);
        od;
    od;
end;



TestAllGroups:=function()
    local ListGroups;
    ListGroup:=GetListCandidateGroups();
    for eGRP in ListGroup
    do
        nMov:=LargestMovedPoint(eGRP);
        TestSpecificGroup(nbMov, eGRP);
    od;
end;
