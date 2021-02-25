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
  Print("GAP eEquiv=", eEquiv, "\n");
end;


# We need to consider cases where the group
# is not solvable. Because otherwise GAP uses
# a different algorithm.



  eFile:="ExampleEquivGRP_Set";
#  CreateExampleOnSetCase(eFile, Group([(1,2,3,4,5),(3,4,5)]), [1, 5], [2, 5]);
#  CreateExampleOnSetCase(eFile, Group( [ ( 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23), ( 3,17,10, 7, 9)( 4,13,14,19, 5)( 8,18,11,12,23)(15,20,22,21,16) ] ), [ 5, 6, 7, 9, 11, 12, 14, 15, 18, 22, 23 ], [ 5, 6, 8, 11, 13, 14, 18, 19, 20, 21, 23 ]);
#  CreateExampleOnSetCase(eFile, Group( [ ( 1, 9, 6, 7, 5)( 2,10, 3, 8, 4), ( 1,10, 7, 8)( 2, 9, 4, 6) ] ), [ 1, 2, 5, 6, 9 ], [ 3, 6, 8, 9, 10 ]);
  CreateExampleOnSetCase(eFile, Group( [ ( 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23), ( 3,17,10, 7, 9)( 4,13,14,19, 5)( 8,18,11,12,23)(15,20,22,21,16),
                                         ( 1,24)( 2,23)( 3,12)( 4,16)( 5,18)( 6,10)( 7,20)( 8,14)( 9,21)(11,17)(13,22)(15,19) ] ),
                         [ 1, 4, 5, 6, 7, 8, 9, 12, 14, 15, 16, 21, 24 ], [ 1, 3, 4, 5, 7, 8, 10, 12, 13, 18, 19, 22, 24 ]);
  
  
  
Print("eFile=", eFile, "\n");
