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
  PrintStabChain(GRP);
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

  eStab:=Stabilizer(GRP, eSet, OnSets);
  Print("|eStab|=", Order(eStab), "\n");
end;


DoMathieu:=false;
if DoMathieu then
  for eSize in [1..12]
  do
    eFile:=Concatenation("ExampleM24_", String(eSize));
    CreateExampleOnSetCase(eFile, MathieuGroup(24), eSize);
  od;
  #
  for eSize in [1..6]
  do
    eFile:=Concatenation("ExampleM12_", String(eSize));
    CreateExampleOnSetCase(eFile, MathieuGroup(12), eSize);
  od;
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

DoSym4:=true;
if DoSym4 then
  for eSize in [2..2]
  do
    eFile:=Concatenation("ExampleSym4_", String(eSize));
    CreateExampleOnSetCase(eFile, SymmetricGroup(4), eSize);
  od;
fi;
