#ifndef DEFINE_PERMUTALIB_FACE_BASIC_H
#define DEFINE_PERMUTALIB_FACE_BASIC_H

#include <iostream>
#include "exception.h"
#include <functional>

#include <bitset>
#include <boost/dynamic_bitset.hpp>

namespace permutalib {


typedef boost::dynamic_bitset<> Face;

std::vector<int> FaceToVector(Face const& eSet)
{
  size_t nbVert=eSet.count();
  std::vector<int> eList(nbVert);
  boost::dynamic_bitset<>::size_type aRow=eSet.find_first();
  for (size_t i=0; i<nbVert; i++) {
    eList[i]=aRow;
    aRow=eSet.find_next(aRow);
  }
  return eList;
}

Face ReadFace(std::istream & is)
{
  if (!is.good()) {
    std::cerr << "ReadFace operation failed because stream is not valid\n";
    throw PermutalibException{1};
  }
  int len, eVal;
  is >> len;
  Face eFace(len);
  for (int i=0; i<len; i++) {
    is >> eVal;
    eFace[i]=eVal;
  }
  return eFace;
}



// We require x and y to be of the same size
bool operator<(Face const& x, Face const& y)
{
  size_t len=x.size();
  for (size_t i=0; i<len; i++) {
    if (x[i] == 0 && y[i] == 1)
      return true;
    if (x[i] == 1 && y[i] == 0)
      return false;
  }
  return false;
}


}
#endif
