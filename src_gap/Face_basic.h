#ifndef INCLUDE_FACE_BASIC_H
#define INCLUDE_FACE_BASIC_H

#include <iostream>
#include "exception.h"

#include <bitset>
#include <boost/dynamic_bitset.hpp>

typedef boost::dynamic_bitset<> Face;

#include <bitset>
#include <boost/dynamic_bitset.hpp>

namespace permutalib {


typedef boost::dynamic_bitset<> Face;

std::vector<int> FaceToVector(Face const& eSet)
{
  int nbVert=eSet.count();
  std::vector<int> eList(nbVert);
  int aRow=eSet.find_first();
  for (int i=0; i<nbVert; i++) {
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
  int len=x.size();
  for (int i=0; i<len; i++) {
    if (x[i] == 0 && y[i] == 1)
      return true;
    if (x[i] == 1 && y[i] == 0)
      return false;
  }
  return false;
}


}
#endif
