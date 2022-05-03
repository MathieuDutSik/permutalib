// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include <iostream>
#include <libgap-api.h>

int main(int argc, char **argv) {
  std::cerr << "argc=" << argc << "\n";
  GAP_Initialize(argc, argv, 0, 0, 1);

  std::string comm = "G:=Group([ (1,2), (2,3) ]);;  f1:=[1,3];; f2:=[1,4];;  "
                     "RepresentativeAction(G, f1, f2, OnSets)<>fail;";
  Int ok = GAP_Enter();
  std::cerr << "ok=" << ok << "\n";
  Obj res = GAP_EvalString(comm.c_str());
  std::cerr << "res=" << res << "\n";
  Int rc = GAP_LenList(res);
  std::cerr << "rc=" << rc << "\n";
  for (Int i = 1; i <= rc; i++) {
    Obj ires = GAP_ElmList(res, i);
    if (GAP_ElmList(ires, 1) == GAP_True) {
      Char *buffer = GAP_CSTR_STRING(GAP_ElmList(ires, 5));
      if (buffer)
        printf("%s\n", buffer);
    }
  }
}
