// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GAP_FACE_BASIC_H_
#define SRC_GAP_FACE_BASIC_H_

#include "exception.h"
#include <functional>
#include <iostream>
#include <vector>

#include <bitset>
#include <boost/dynamic_bitset.hpp>

namespace permutalib {

using Face = boost::dynamic_bitset<>;

template <typename Tidx> std::vector<Tidx> FaceToVector(Face const &eSet) {
  size_t nbVert = eSet.count();
  std::vector<Tidx> eList(nbVert);
  boost::dynamic_bitset<>::size_type aRow = eSet.find_first();
  for (size_t i = 0; i < nbVert; i++) {
    eList[i] = Tidx(aRow);
    aRow = eSet.find_next(aRow);
  }
  return eList;
}

template<typename Telt>
Face FaceAct(Face const& f, Telt const& u) {
  size_t len = f.size();
  Face fret(len);
  boost::dynamic_bitset<>::size_type aRow = f.find_first();
  while (aRow != boost::dynamic_bitset<>::npos) {
    size_t pos = u.at(aRow);
    fret[pos] = 1;
    aRow = f.find_next(aRow);
  }
  return fret;
}

bool TestEqual(Face const& f1, Face const& f2) {
  size_t len = f1.size();
  if (len != f2.size())
    return false;
  for (size_t i=0; i<len; i++)
    if (f1[i] != f2[i])
      return false;
  return true;
}

Face ReadFace(std::istream &is) {
  if (!is.good()) {
    std::cerr << "ReadFace operation failed because stream is not valid\n";
    throw PermutalibException{1};
  }
  size_t len;
  int eVal;
  is >> len;
  Face eFace(len);
  for (size_t i = 0; i < len; i++) {
    is >> eVal;
    eFace[i] = eVal;
  }
  return eFace;
}

// We require x and y to be of the same size
bool operator<(Face const &x, Face const &y) {
  size_t len = x.size();
  for (size_t i = 0; i < len; i++) {
    if (x[i] == 0 && y[i] == 1)
      return true;
    if (x[i] == 1 && y[i] == 0)
      return false;
  }
  return false;
}

// clang-format off
}  // namespace permutalib
#endif  // SRC_GAP_FACE_BASIC_H_
// clang-format on
