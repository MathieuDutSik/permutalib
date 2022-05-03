// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GAP_GRAPHICFUNCTIONALITY_H_
#define SRC_GAP_GRAPHICFUNCTIONALITY_H_

#include "COMB_Vectors.h"
#include <limits>
#include <vector>

namespace permutalib {

struct GraphSparseImmutable {
public:
  GraphSparseImmutable() = delete;
  GraphSparseImmutable(size_t const &_nbVert,
                       std::vector<size_t> const &_ListStart,
                       std::vector<size_t> const &_ListListAdj)
      : nbVert(_nbVert), ListStart(_ListStart), ListListAdj(_ListListAdj) {
    HasVertexColor = false;
  }
  ~GraphSparseImmutable() {}
  GraphSparseImmutable(std::vector<size_t> const &ListEdge,
                       size_t const &_nbVert)
      : nbVert(_nbVert) {
    HasVertexColor = false;
    size_t nbEdge = ListEdge.size() / 2;
    std::vector<size_t> ListDeg(nbVert, 0);
    for (size_t iEdge = 0; iEdge < nbEdge; iEdge++)
      for (size_t i = 0; i < 2; i++) {
        size_t eVert = ListEdge[2 * iEdge + i];
        ListDeg[eVert]++;
      }
    ListStart.resize(nbVert + 1, 0);
    size_t nbAdj = 0;
    for (size_t iVert = 0; iVert < nbVert; iVert++) {
      size_t eDeg = ListDeg[iVert];
      ListStart[iVert + 1] = ListStart[iVert] + eDeg;
      nbAdj += eDeg;
    }
    std::vector<size_t> ListPos(nbVert, 0);
    ListListAdj.resize(nbAdj, 0);
    for (size_t iEdge = 0; iEdge < nbEdge; iEdge++)
      for (size_t i = 0; i < 2; i++) {
        size_t j = 1 - i;
        size_t eVert = ListEdge[2 * iEdge + i];
        size_t fVert = ListEdge[2 * iEdge + j];
        size_t eStart = ListStart[eVert];
        size_t pos = ListPos[eVert];
        ListListAdj[eStart + pos] = fVert;
        ListPos[eVert]++;
      }
  }
  GraphSparseImmutable(GraphSparseImmutable const &eG) {
    nbVert = eG.GetNbVert();
    ListStart = eG.GetListStart();
    ListListAdj = eG.GetListListAdj();
    HasVertexColor = eG.GetHasVertexColor();
    ListVertexColor = eG.GetListVertexColor();
  }
  GraphSparseImmutable operator=(GraphSparseImmutable const &eG) {
    nbVert = eG.GetNbVert();
    ListStart = eG.GetListStart();
    ListListAdj = eG.GetListListAdj();
    HasVertexColor = eG.GetHasVertexColor();
    ListVertexColor = eG.GetListVertexColor();
    return *this;
  }
  // lighter stuff
  size_t GetNbVert() const { return nbVert; }
  std::vector<size_t> GetListListAdj() const { return ListListAdj; }
  std::vector<size_t> GetListStart() const { return ListStart; }
  bool GetHasVertexColor() const { return HasVertexColor; }
  std::vector<size_t> GetListVertexColor() const { return ListVertexColor; }
  //
  void SetHasColor(bool const &TheVal) {
    if (TheVal == HasVertexColor)
      return;
    HasVertexColor = TheVal;
    if (TheVal)
      ListVertexColor = std::vector<size_t>(nbVert);
    if (!TheVal)
      ListVertexColor.clear();
  }
  void SetColor(size_t const &iVert, size_t const &eColor) {
    ListVertexColor[iVert] = eColor;
  }
  std::vector<size_t> Adjacency(size_t const &iVert) const {
    size_t eStart = ListStart[iVert];
    size_t eEnd = ListStart[iVert + 1];
    std::vector<size_t> TheRet;
    for (size_t i = eStart; i < eEnd; i++)
      TheRet.push_back(ListListAdj[i]);
    return TheRet;
  }
  bool IsAdjacent(size_t const &iVert, size_t const &jVert) const {
    size_t eStart = ListStart[iVert];
    size_t eEnd = ListStart[iVert + 1];
    for (size_t i = eStart; i < eEnd; i++)
      if (ListListAdj[i] == jVert)
        return true;
    return false;
  }
  size_t GetColor(size_t const &iVert) const {
    if (!HasVertexColor) {
      std::cerr << "Call to GetColor while HasVertexColor=false\n";
      throw PermutalibException{1};
    }
    return ListVertexColor[iVert];
  }

private:
  size_t nbVert;
  std::vector<size_t> ListStart;
  std::vector<size_t> ListListAdj;
  bool HasVertexColor;
  std::vector<size_t> ListVertexColor;
};

template <typename Tgr>
std::vector<size_t> ConnectedComponents_vector(Tgr const &GR) {
  size_t nbVert = GR.GetNbVert();
  size_t miss_val = std::numeric_limits<size_t>::max();
  std::vector<size_t> ListStatus(nbVert, miss_val);
  size_t iStatus = 0;
  auto Assignment = [&](size_t const &iVert) -> void {
    std::vector<size_t> TheSet{iVert};
    ListStatus[iVert] = iStatus;
    while (true) {
      std::vector<size_t> NewSet;
      for (size_t &eVal : TheSet) {
        for (size_t &fVal : GR.Adjacency(eVal)) {
          if (ListStatus[fVal] == miss_val) {
            ListStatus[fVal] = iStatus;
            NewSet.push_back(fVal);
          }
        }
      }
      if (NewSet.size() == 0)
        break;
      TheSet = NewSet;
    }
    iStatus++;
  };
  for (size_t iVert = 0; iVert < nbVert; iVert++)
    if (ListStatus[iVert] == miss_val)
      Assignment(iVert);
  return ListStatus;
}

template <typename Tgr>
std::vector<std::vector<size_t>> ConnectedComponents_set(Tgr const &GR) {
  size_t nbVert = GR.GetNbVert();
  std::vector<size_t> ListStatus = ConnectedComponents_vector(GR);
  size_t nbConn = VectorMax(ListStatus) + 1;
  std::vector<std::vector<size_t>> ListConn(nbConn);
  for (size_t iVert = 0; iVert < nbVert; iVert++) {
    size_t iStatus = ListStatus[iVert];
    ListConn[iStatus].push_back(iVert);
  }
  return ListConn;
}

// clang-format off
}  // namespace permutalib
#endif  // SRC_GAP_GRAPHICFUNCTIONALITY_H_
// clang-format on
