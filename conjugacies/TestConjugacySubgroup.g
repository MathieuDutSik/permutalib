GRP:=SymmetricGroup(7);

CJ:=ConjugacyClassesSubgroups(GRP);
ListOrd:=List(CJ, x->Order(Representative(x)));
pos:=Position(ListOrd, 16);

eRepr:=Representative(CJ[pos]);
fRepr:=ConjugateGroup(eRepr, Random(GRP));

g:=RepresentativeAction(GRP, eRepr, fRepr);
