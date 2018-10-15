CreateExampleOnSetCase:=function(FileName, GRP, sizSet)
  local LGen, SetMovedPt, eGen, nbMov, output, iMov, eImg, eSet, pos, eVal, eStab;
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
  RemoveFileIfExist(FileName);
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
  eSet:=RandomSubset([1..nbMov], sizSet);
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



CreateExampleOnSetCase("ExampleM24", MathieuGroup(24), 12);


