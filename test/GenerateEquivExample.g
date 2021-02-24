RequirePackage("gapcommon");
Local_RemoveFileIfExist:=function(FileName)
  if IsExistingFile(FileName)=true then
    RemoveFile(FileName);
  fi;
end;


CreateExampleOnSetCase:=function(FileName, GRP, eSet, fSet)
  local LGen, eGen, nbMov, output, iMov, eImg, pos, eVal, eEquiv;
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
  Print("GAP |eG|=", Order(GRP), "\n");
  Print("GAP Before RepresentativeAction_OnSets\n");
  eEquiv:=RepresentativeAction(GRP, eSet, fSet, OnSets);
  Print("eEquiv=", eEquiv, "\n");
end;


# We need to consider cases where the group
# is not solvable. Because otherwise GAP uses
# a different algorithm.



  eFile:="ExampleEquivGRP_Set";
  CreateExampleOnSetCase(eFile, Group([(1,2,3,4,5),(3,4,5)]), [1, 5], [2, 5]);
Print("eFile=", eFile, "\n");
