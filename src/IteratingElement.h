// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GAP_ITERATING_ELEMENT_H_
#define SRC_GAP_ITERATING_ELEMENT_H_

/*
  Functionalities for iterating over all elements
  and generating all the elements if needed.
 */
#include "StabChainMain.h"
#include "PermGroup.h"

namespace permutalib {

template<typename Telt, typename Tidx_label>
struct IteratorType {
private:
  std::vector<StabChain<Telt, Tidx_label>> ListS;
  std::vector<size_t> ListPos;
  std::vector<size_t> ListSiz;
  std::vector<Telt> ListRes;
public:
IteratorType(std::vector<StabChain<Telt, Tidx_label>> ListS,
             std::vector<size_t> ListPos, std::vector<size_t> ListSiz,
             std::vector<Telt> ListRes)
  : ListS(ListS), ListPos(ListPos), ListSiz(ListSiz), ListRes(ListRes) {}
  const Telt &operator*() const { return ListRes[0]; }
  void IterIncrease() {
    using Tidx = typename Telt::Tidx;
    Tidx n = ListS[0]->comm->n;
    size_t len = ListPos.size();
    for (size_t i = 0; i < len; i++) {
      if (ListPos[i] < ListSiz[i] - 1) {
        ListPos[i]++;
        for (size_t j = 0; j < i; j++)
          ListPos[j] = 0;
        Telt elt(n);
        Tidx bpt = ListS[i]->orbit[0];
        Tidx img = ListS[i]->orbit[ListPos[i]];
        Tidx img_work = img;
        while (true) {
          if (img_work == bpt)
            break;
          Tidx_label idx = ListS[i]->transversal[img_work];
          elt *= ListS[i]->comm->labels[idx];
          img_work = PowAct(img, elt);
        }
        if (i != len - 1)
          elt *= ListRes[i + 1];
        for (size_t j = 0; j <= i; j++)
          ListRes[j] = elt;
        return;
      }
    }
    // This is the case of a END iterator
    for (size_t i = 0; i < len; i++)
      ListPos[i] = ListSiz[i];
  }
  IteratorType &operator++() {
    IterIncrease();
    return *this;
  }
  IteratorType operator++(int) {
    IteratorType tmp = *this;
    IterIncrease();
    return tmp;
  }
  bool operator!=(const IteratorType &x) const {
    for (size_t i = 0; i < ListPos.size(); i++)
      if (ListPos[i] != x.ListPos[i])
        return true;
    return false;
  }
  bool operator==(const IteratorType &x) const {
    for (size_t i = 0; i < ListPos.size(); i++)
      if (ListPos[i] != x.ListPos[i])
        return false;
    return true;
  }
};


template<typename Telt, typename Tidx_label>
IteratorType<Telt,Tidx_label> get_begin_iterator(StabChain<Telt, Tidx_label> const& S) {
  using Tidx = typename Telt::Tidx;
  Tidx n = S->comm->n;
  std::vector<StabChain<Telt, Tidx_label>> ListS;
  std::vector<size_t> ListPos;
  std::vector<size_t> ListSiz;
  std::vector<Telt> ListRes;
  StabChain<Telt, Tidx_label> Swork = S;
  while (Swork != nullptr) {
    size_t len = Swork->orbit.size();
    if (len == 0)
      break;
    ListS.push_back(Swork);
    ListPos.push_back(0);
    ListSiz.push_back(len);
    ListRes.push_back(Telt(n));
    Swork = Swork->stabilizer;
  }
  return IteratorType<Telt,Tidx_label>(ListS, ListPos, ListSiz, ListRes);
}

template<typename Telt, typename Tidx_label>
IteratorType<Telt,Tidx_label> get_end_iterator(StabChain<Telt, Tidx_label> const& S) {
  std::vector<size_t> ListPos;
  StabChain<Telt, Tidx_label> Swork = S;
  while (Swork != nullptr) {
    ListPos.push_back(Swork->orbit.size());
    Swork = Swork->stabilizer;
  }
  return IteratorType<Telt,Tidx_label>({}, ListPos, {}, {});
}

template<typename Telt, typename Tidx_label>
std::vector<Telt> get_all_elements(StabChain<Telt, Tidx_label> const& S) {
  std::vector<Telt> l_elt;
  if (S->orbit.size() == 0) {
    l_elt.push_back(S->comm->identity);
  } else {
    IteratorType<Telt,Tidx_label> iter = get_begin_iterator(S);
    IteratorType<Telt,Tidx_label> iter_end = get_end_iterator(S);
    while (iter != iter_end) {
      Telt const& eElt = *iter;
      l_elt.push_back(eElt);
      iter++;
    }
  }
  return l_elt;
}

template<typename Telt, typename Tidx_label>
Face exhaustive_minimum_face_orbit(StabChain<Telt, Tidx_label> const& S, Face const& f) {
  IteratorType<Telt,Tidx_label> iter = get_begin_iterator(S);
  IteratorType<Telt,Tidx_label> iter_end = get_end_iterator(S);
  Face f_minimum = f;
  while (iter != iter_end) {
    Telt const& eElt = *iter;
    Face f_img = OnSets(f, eElt);
    if (f_img < f_minimum)
      f_minimum = f_img;
    iter++;
  }
  return f_minimum;
}

// clang-format off
}  // namespace permutalib
#endif  // SRC_GAP_GROUP_H_
// clang-format on
