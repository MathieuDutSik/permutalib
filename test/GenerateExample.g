RequirePackage("gapcommon");
Local_RemoveFileIfExist:=function(FileName)
  if IsExistingFile(FileName)=true then
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



CreateExampleOnSetCase:=function(FileName, GRP, eSet)
  local LGen, eGen, nbMov, output, iMov, eImg, pos, eVal, eStab;
  nbMov:=Maximum(LargestMovedPoint(GRP), Maximum(eSet));
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
  Print("|eG|=", Order(GRP), "\n");
  Print("Before Stabilizer_OnSets\n");
  eStab:=Stabilizer(GRP, eSet, OnSets);
  Print("|eStab|=", Order(eStab), "\n");
end;


# The case of Sym4 is simpler. However it is also solvable
# and in that case GAP uses the PCGS algorithms that we do not
# want to implement.
#
# Also, we need to put Group([(1,2,3,4,5),(1,2)])
#       instead of SymmetricGroup(5)
# because this triggers a different algo for the Symmetric group
#


DoSym5:=true;
#DoSym5:=false;
if DoSym5 then
  eSet:=Local_RandomSubset([1..5], 3);
  eFile:="ExampleGRP_Set";
#  CreateExampleOnSetCase(eFile, Group([(1,2,3,4,5,6,7,8,9),(1,2)]), eSet);
#  CreateExampleOnSetCase(eFile, Group([(1,2,3,4,5),(4,5)]), eSet);
#  CreateExampleOnSetCase(eFile, Group([(1,2,3,4,5),(3,4,5)]), eSet);
  CreateExampleOnSetCase(eFile, Group([ (1,4,5,9,3)(2,8,10,7,6)(12,15,16,20,14)(13,19,21,18,17), (1,21,5,12,20)(2,16,3,4,17)(6,18,7,19,15)(8,13,9,14,11) ]), [ 2, 4, 5, 9, 11, 12, 13, 14, 18 ]);
  Print("eFile=", eFile, "\n");
fi;
