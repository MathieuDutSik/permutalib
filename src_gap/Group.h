#ifndef DEFINE_PERMUTALIB_GROUP_H
#define DEFINE_PERMUTALIB_GROUP_H


#include "StabChainMain.h"
#include "stbcbckt.h"
#include "nsi.h"
#include "Properties.h"
#include "NormalStructure.h"
#include <map>


namespace permutalib {



template<typename Telt>
Telt RandomElement(const std::vector<Telt>& LGen, const typename Telt::Tidx& n)
{
  size_t len = rand() % 100;
  size_t n_gen = LGen.size();
  Telt eElt(n);
  for (size_t iIter=0; iIter<len; iIter++) {
    size_t pos = size_t(rand()) % n_gen;
    eElt *= LGen[pos];
  }
  return eElt;
}





template<typename Telt_inp, typename Tint_inp>
struct Group {
public:
  // dependent types
  using Telt = Telt_inp;
  using Tidx = typename Telt::Tidx;
  using Tint = Tint_inp;
  using Tidx_label = uint16_t;
  // constructors
  Group(const StabChain<Telt,Tidx_label>& _S) : S(_S), size_tint(Order<Telt,Tidx_label,Tint>(_S)) {
  }
  Group(const std::vector<Telt>& LGen, const Tidx& n) {
#ifdef DEBUG_STABCHAINMAIN
    std::cerr << "CPP Beginning of MinimalStabChain\n";
#endif
    StabChainOptions<Tint,Tidx> options = GetStandardOptions<Tint,Tidx>(n);
#ifdef DEBUG_STABCHAINMAIN
    std::cerr << "CPP Before StabChainOp_listgen\n";
#endif
    S = StabChainOp_listgen<Telt,Tidx_label,Tint>(LGen, options);
    UnbindCycles(S);
    size_tint = Order<Telt,Tidx_label,Tint>(S);
  }
  Group(const Tidx& n) : Group({}, n) {
  }
  Group() : Group(0) {
  }
  Group(Group<Telt,Tint> && G) : S(std::move(G.S)), size_tint(G.size_tint) {
  }
  Group(const Group<Telt,Tint>& G) : S(G.S), size_tint(G.size_tint) {
  } // The S is a shared_ptr so copy is fine.
  Group<Telt,Tint>& operator=(const Group<Telt,Tint>& G) {
    // The S is a shared_ptr so copy is fine.
    S = G.S;
    size_tint = G.size_tint;
    return *this;
  }
  // Basic getters
  std::vector<Telt> GeneratorsOfGroup() const {
    return Kernel_GeneratorsOfGroup(S);
  }
  std::string GapString() const {
    std::vector<Telt> LGen = Kernel_GeneratorsOfGroup(S);
    if (LGen.size() == 0) {
      return "Group(())";
    } else {
      return "Group(" + GapStringTVector(LGen) + ")";
    }
  }
  Tint size() const {
    return size_tint;
  }
  std::map<Tidx, int> factor_size() const {
    return FactorsSizeStabChain(S);
  }
  Tidx n_act() const {
    return S->comm->n;
  }
  // operation
  Group<Telt,Tint> GroupConjugate(const Telt& x) const {
    std::vector<Telt> LGen;
    Telt xInv =~x;
    for (auto & eGen : Kernel_GeneratorsOfGroup(S)) {
      Telt eGenCj = xInv * eGen * x;
      LGen.emplace_back(eGenCj);
    }
    return Group<Telt,Tint>(LGen, S->comm->n);
  }
  // Action on points or sets
  Group<Telt,Tint> Stabilizer_OnPoints(const Tidx& x) const {
    return Group(Kernel_Stabilizer_OnPoints<Telt,Tidx_label,Tint>(S, x));
  }
  std::pair<bool,Telt> RepresentativeAction_OnPoints(const Tidx& x1, const Tidx& x2) const {
    return Kernel_RepresentativeAction_OnPoints<Telt,Tidx_label,Tint>(S, x1, x2);
  }
  Group<Telt,Tint> Stabilizer_OnSets(const Face& f) const {
    return Group(Kernel_Stabilizer_OnSets<Telt,Tidx_label,Tint>(S, f));
  }
  std::pair<bool,Telt> RepresentativeAction_OnSets(const Face& f1, const Face& f2) const {
    return Kernel_RepresentativeAction_OnSets<Telt,Tidx_label,Tint>(S, f1, f2);
  }
  bool operator==(const Group& g) const {
    return EqualityTest(S, g.S);
  }
  Face CanonicalImage(const Face& f) const {
    return Kernel_CanonicalImage<Telt,Tidx_label,Tint>(S, f);
  }
  Telt rand() const {
    return RandomElement(Kernel_GeneratorsOfGroup(S), S->comm->n);
  }
  bool IsCommutative() const {
    return Kernel_IsCommutative(S);
  }
  bool IsTransitive() const {
    return Kernel_IsTransitive(S);
  }
  bool IsPrimitive() const {
    return Kernel_IsPrimitive(S);
  }
  bool IsCyclic() const {
    return Kernel_IsCyclic<Telt,Tidx_label,Tint>(S);
  }
  std::vector<Telt> SmallGeneratingSet() const {
    return Kernel_SmallGeneratingSet<Telt,Tidx_label,Tint>(S);
  }
  // Normal structure
  bool IsNormalSubgroup(const Group<Telt,Tint>& U) const {
    return Kernel_IsNormalSubgroup(S, U.S);
  }
  Group<Telt,Tint> NormalClosure(const Group<Telt,Tint>& H) const {
    return Group<Telt,Tint>(Kernel_NormalClosure<Telt,Tidx_label,Tint>(S, H.S));
  }
  Group<Telt,Tint> DerivedSubgroup() const {
    return Group<Telt,Tint>(Kernel_DerivedSubgroup<Telt,Tidx_label,Tint>(S));
  }
  Group<Telt,Tint> CentreSubgroup() const {
    return Group<Telt,Tint>(Kernel_CentreSubgroup<Telt,Tidx_label,Tint>(S));
  }
  Group<Telt,Tint> Centralizer_elt(const Telt& x) const {
    return Group<Telt,Tint>(Kernel_Centralizer_elt<Telt,Tidx_label,Tint>(S, x));
  }
  Group<Telt,Tint> Centralizer_grp(const Group<Telt,Tint>& H) const {
    return Group<Telt,Tint>(Kernel_Centralizer_grp<Telt,Tidx_label,Tint>(S, H.S));
  }
  bool isin(const Telt& x) const {
    return IsElementInStabChain(S, x);
  }
private:
  struct IteratorType {
  private:
    std::vector<StabChain<Telt,Tidx_label>> ListS;
    std::vector<size_t> ListPos;
    std::vector<size_t> ListSiz;
    std::vector<Telt> ListRes;
  public:
    IteratorType(std::vector<StabChain<Telt,Tidx_label>> ListS,
                 std::vector<size_t> ListPos, std::vector<size_t> ListSiz, std::vector<Telt> ListRes) : ListS(ListS), ListPos(ListPos), ListSiz(ListSiz), ListRes(ListRes) {
    }
    const Telt& operator*() const {
      return ListRes[0];
    }
    void IterIncrease() {
      Tidx n = ListS[0]->comm->n;
      size_t len = ListPos.size();
      for (size_t i=0; i<len; i++) {
        if (ListPos[i] < ListSiz[i] - 1) {
          ListPos[i]++;
          for (size_t j=0; j<i; j++)
            ListPos[j] = 0;
          Telt elt(n);
          Tidx bpt = ListS[i]->orbit[0];
          Tidx img = ListS[i]->orbit[ListPos[i]];
          Tidx img_work = img;
          while(true) {
            if (img_work == bpt)
              break;
            Tidx_label idx = ListS[i]->transversal[img_work];
            elt *= ListS[i]->comm->labels[idx];
            img_work = PowAct(img, elt);
          }
          if (i != len - 1)
            elt *= ListRes[i+1];
          for (size_t j=0; j<=i; j++)
            ListRes[j] = elt;
          return;
        }
      }
      // This is the case of a END iterator
      for (size_t i=0; i<len; i++)
        ListPos[i] = ListSiz[i];
    }
    IteratorType& operator++() {
      IterIncrease();
      return *this;
    }
    IteratorType operator++(int) {
      IteratorType tmp = *this;
      IterIncrease();
      return tmp;
    }
    bool operator!=(const IteratorType& x) const {
      for (size_t i=0; i<ListPos.size(); i++)
        if (ListPos[i] != x.ListPos[i])
          return true;
      return false;
    }
    bool operator==(const IteratorType& x) const {
      for (size_t i=0; i<ListPos.size(); i++)
        if (ListPos[i] != x.ListPos[i])
          return false;
      return true;
    }
  };
public:
  using iterator=IteratorType;
  using const_iterator=IteratorType;
  const_iterator begin() const {
    Tidx n = S->comm->n;
    std::vector<StabChain<Telt,Tidx_label>> ListS;
    std::vector<size_t> ListPos;
    std::vector<size_t> ListSiz;
    std::vector<Telt> ListRes;
    StabChain<Telt,Tidx_label> Swork = S;
    while(Swork != nullptr) {
      size_t len = Swork->orbit.size();
      if (len == 0)
        break;
      ListS.push_back(Swork);
      ListPos.push_back(0);
      ListSiz.push_back(len);
      ListRes.push_back(Telt(n));
      Swork = Swork->stabilizer;
    }
    return IteratorType(ListS, ListPos, ListSiz, ListRes);
  }
  const_iterator end() const {
    std::vector<size_t> ListPos;
    StabChain<Telt,Tidx_label> Swork = S;
    while(Swork != nullptr) {
      ListPos.push_back(Swork->orbit.size());
      Swork = Swork->stabilizer;
    }
    return IteratorType({}, ListPos, {}, {});
  }
private:
  StabChain<Telt,Tidx_label> S;
  Tint size_tint;
};


template<typename Tgroup>
std::ostream& operator<<(std::ostream& os, const Tgroup& grp)
{
  boost::archive::text_oarchive oa(os);
  os << grp;
  return os;
}


template<typename Tgroup>
std::istream& operator>>(std::istream& is, Tgroup& grp)
{
  boost::archive::text_iarchive ia(is);
  is >> grp;
  return is;
}


}
#endif
