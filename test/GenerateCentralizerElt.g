RequirePackage("gapcommon");
Local_RemoveFileIfExist:=function(FileName)
  if IsExistingFile(FileName)=true then
    RemoveFile(FileName);
  fi;
end;



GenerateCentralizerEltExample:=function(FileName, GRP, eElt)
    local LGen, eGen, nbMov, output, iMov, eImg, pos, eVal, eCentrElt;
    nbMov:=LargestMovedPoint(GRP);
    #
    Local_RemoveFileIfExist(FileName);
    output:=OutputTextFile(FileName, true);
    LGen:=GeneratorsOfGroup(GRP);
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
        eImg:=OnPoints(iMov, eElt);
        AppendTo(output, " ", eImg-1);
    od;
    AppendTo(output, "\n");
    CloseStream(output);
    #
    Print("Before Centralizer\n");
    eCentrElt:=Centralizer(GRP, eElt);
    Print("|eCentrElt|=", Order(eCentrElt), "\n");
end;


# We need to consider cases where the group
# is not solvable. Because otherwise GAP uses
# a different algorithm.



  eFile:="ExampleCentralizerElt";


# inadequate since there are special method for symmetric group
#GenerateCentralizerEltExample(eFile, SymmetricGroup(5), (1,5)(2,4)); 

GenerateCentralizerEltExample(eFile, Group([(1,2,3,4,5), (3,4,5)]), (1,5)(2,4));
Print("eFile=", eFile, "\n");
