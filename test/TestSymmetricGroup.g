
n:=5;

ePerm1:=PermList(Concatenation([2..n], [1]));
ePerm2:=(1,2);

eG:=Group([ePerm1, ePerm2]);
wS:=MinimalStabChain(eG);

#u:=Random(eG);

