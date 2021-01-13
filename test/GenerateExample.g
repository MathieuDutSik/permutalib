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



CreateExampleOnSetCase:=function(FileName, GRP, sizSet)
  local LGen, SetMovedPt, eGen, nbMov, output, iMov, eImg, eSet, pos, eVal, eStab;
#  PrintStabChain(GRP);
  LGen:=GeneratorsOfGroup(GRP);
  SetMovedPt:=[];
  for eGen in LGen
  do
    SetMovedPt:=Union(SetMovedPt, MovedPoints(eGen));
  od;
  nbMov:=Length(SetMovedPt);
  if SetMovedPt<>[1..nbMov] then
    Print("Some assumption need to be rethough. The set of moving points is a closed interval [1..N]\n");
    Error("Please correct");
  fi;
  #
  Local_RemoveFileIfExist(FileName);
  output:=OutputTextFile(FileName, true);
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
  eSet:=Local_RandomSubset([1..nbMov], sizSet);
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


#DoMathieu24:=true;
DoMathieu24:=false;
if DoMathieu24 then
  eSize:=12;
  eFile:=Concatenation("ExampleM24_", String(eSize));
  CreateExampleOnSetCase(eFile, MathieuGroup(24), eSize);
fi;

DoMathieu12:=false;
if DoMathieu12 then
  eSize:=6;
  eFile:=Concatenation("ExampleM12_", String(eSize));
  CreateExampleOnSetCase(eFile, MathieuGroup(12), eSize);
fi;

DoSym6:=false;
if DoSym6 then
  for eSize in [1..3]
  do
    eFile:=Concatenation("ExampleSym6_", String(eSize));
    CreateExampleOnSetCase(eFile, SymmetricGroup(6), eSize);
  od;
  #
fi;

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
  eSize:=3;
  eFile:="ExampleGRP_Set";
#  CreateExampleOnSetCase(eFile, Group([(1,2,3,4,5,6,7,8,9),(1,2)]), eSize);
#  CreateExampleOnSetCase(eFile, Group([(1,2,3,4,5),(4,5)]), eSize);
  CreateExampleOnSetCase(eFile, Group([(1,2,3,4,5),(3,4,5)]), eSize);
  Print("eFile=", eFile, "\n");
fi;
