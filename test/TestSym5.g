RequirePackage("gapcommon");

LGen:=[(1,2,3,4,5),(1,2)];
#LGen:=GeneratorsOfGroup(MathieuGroup(24));
#LGen:=[(1,2,3,4)];


GRP:=Group(LGen);

TheOrd:=Order(GRP);
Print("TheOrd=", TheOrd, "\n");
