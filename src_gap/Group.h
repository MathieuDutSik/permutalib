// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GAP_GROUP_H_
#define SRC_GAP_GROUP_H_

// clang-format off
#include "StabChainMain.h"
#include "stbcbckt.h"
#include "nsi.h"
#include "Properties.h"
#include "NormalStructure.h"
#include "PermutationElt.h"
#include "AscendingChains_and_Cosets.h"
// clang-format on
#include <map>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <boost/archive/tmpdir.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/utility.hpp>

/*
  The Group class is far too rigid for us.
  ---We cannot handle different sizes occurring in some algorithm like
  Centering. This would require having the SingleSidedPerm in the code directly.
  ---Should we have a notion of trivial group? NO.
     ---PRO: Simpler constructor. But giving the number of elements is never qn
  issue really.
     ---CON: A complications of the code.
     ---CON: When returning an equivalence, what to do? All kinds of type
  problem will show up.
  ---We need to have constructor using identity element. This will be needed for
  the construction of the Kernel and the PreImages. So, that part of the changes
  is not under discussions.
  ---The API of the GAP with ListMatrGens, ListPermGens should be used as well
  in this case.
  ---When implementing (for example MyMatrix), we need to have a custom type
  that allows for
     ---isIdentity
     ---Default constructor that gives the right entry (e.g. the
  IdentityMat<T>(dim))
     ---The inverse
     ---The product, the *= and other operators.
     So, for example the MyMatrix will have to be contained in some Singleton
  class. But is that ok? The dimension n is only known dynamically, so cannot be
  part of the template parameter of the class.

  ---For computing PreImage, we have to use the Sifted permutation.
  ---For the Kernel, when we work only with permutations, things are clear. But
  for general case, we have serious thinking to do.

  Kernel computation:
  ---The command for getting the stab chain is StabChainStrong.
  ---The ExtendStabChain depends on ChangeStabChain and others.
  ---Maybe better is to encode our own function for finding that Kernel.
  ---In all objectivity, Schreier lemma provide a solution to our problem.
  However, there are several issues:
      ---The generating set is large.
      ---It requires us to work with the right cosets, which we do not have
  right now as functionality.
      ---On the contrary, the algorithm for finding the StrongStabChain, works
  with orbits and do not require this. So, further thinking is needed.

  Thinking:
  ---It seems sure that we cannot use the StabChainStrong algorithm. This is
  because that algorithm eventually has to find a base point. And precisely, we
  will not find a base point.
  ---When building the Kernel, if a stabchain algorithm can be applied, GAP
  would probably use it. So, it makes sense to look at the code.
     ---The algorithm is pretty complicated
     ---It uses the KernelOfMultiplicativeGeneralMapping and then
  CoKernelOfMultiplicativeGeneralMapping
     ---Two algorithms are used: NormalClosure   and   CoKernelGensPermHom.
     ---NormalClosure (code is in grp.gi) depends on the testing that an element
  belongs to the group.
  ---NormalClosure algorithm dependence on testing membership. This is actually
  an expensive algorithm. So, maybe the algorithm of GAP is not adequate for us.

 */

namespace permutalib {

template <typename Telt>
Telt RandomElement(const std::vector<Telt> &LGen, const Telt &id) {
  size_t len = random() % 100;
  size_t n_gen = LGen.size();
  Telt eElt = id;
  for (size_t iIter = 0; iIter < len; iIter++) {
    size_t pos = size_t(random()) % n_gen;
    eElt *= LGen[pos];
  }
  return eElt;
}

template <typename Telt_inp, typename Tint_inp> struct Group {
public:
  // dependent types
  using Telt = Telt_inp;
  using Tidx = typename Telt::Tidx;
  using Tint = Tint_inp;
  using Tidx_label = uint16_t;
  // constructors
  Group(const StabChain<Telt, Tidx_label> &_S)
      : S(_S), size_tint(Order<Telt, Tidx_label, Tint>(_S)) {}
  Group(const std::vector<Telt> &LGen, const Telt &id) {
#ifdef DEBUG_STABCHAINMAIN
    std::cerr << "CPP Beginning of MinimalStabChain\n";
#endif
    StabChainOptions<Tint, Telt> options = GetStandardOptions<Tint, Telt>(id);
#ifdef DEBUG_STABCHAINMAIN
    std::cerr << "CPP Before StabChainOp_listgen\n";
#endif
    S = StabChainOp_listgen<Telt, Tidx_label, Tint>(LGen, options);
    UnbindCycles(S);
    size_tint = Order<Telt, Tidx_label, Tint>(S);
  }
  Group(const Tidx &n) : Group({}, n) {}
  Group() : Group(0) {}
  Group(Group<Telt, Tint> &&G) : S(std::move(G.S)), size_tint(G.size_tint) {}
  Group(const Group<Telt, Tint> &G) : S(G.S), size_tint(G.size_tint) {}
  Group<Telt, Tint> &operator=(const Group<Telt, Tint> &G) {
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
  Tint size() const { return size_tint; }
  std::map<Tidx, int> factor_size() const { return FactorsSizeStabChain(S); }
  Tidx n_act() const { return S->comm->n; }
  // operation
  Group<Telt, Tint> GroupConjugate(const Telt &x) const {
    std::vector<Telt> LGen;
    Telt xInv = ~x;
    for (auto &eGen : Kernel_GeneratorsOfGroup(S)) {
      Telt eGenCj = xInv * eGen * x;
      LGen.emplace_back(eGenCj);
    }
    return Group<Telt, Tint>(LGen, S->comm->n);
  }
  // Action on points or sets
  Group<Telt, Tint> Stabilizer_OnPoints(const Tidx &x) const {
    return Group(Kernel_Stabilizer_OnPoints<Telt, Tidx_label, Tint>(S, x));
  }
  std::optional<Telt> RepresentativeAction_OnPoints(const Tidx &x1,
                                                    const Tidx &x2) const {
    return Kernel_RepresentativeAction_OnPoints<Telt, Tidx_label, Tint>(S, x1,
                                                                        x2);
  }
  Group<Telt, Tint> Stabilizer_OnSets(const Face &f) const {
    return Group(Kernel_Stabilizer_OnSets<Telt, Tidx_label, Tint>(S, f));
  }
  std::optional<Telt> RepresentativeAction_OnSets(const Face &f1,
                                                  const Face &f2) const {
    return Kernel_RepresentativeAction_OnSets<Telt, Tidx_label, Tint>(S, f1,
                                                                      f2);
  }
  bool operator==(const Group &g) const { return EqualityTest(S, g.S); }
  bool operator!=(const Group &g) const { return !EqualityTest(S, g.S); }
  Face CanonicalImage(const Face &f) const {
    return Kernel_CanonicalImage<Telt, Tidx_label, Tint>(S, f);
  }
  Telt rand() const {
    return RandomElement(Kernel_GeneratorsOfGroup(S), S->comm->identity);
  }
  Telt random() const { return rand(); }
  bool IsCommutative() const { return Kernel_IsCommutative(S); }
  bool IsTransitive() const { return Kernel_IsTransitive(S); }
  bool IsPrimitive() const { return Kernel_IsPrimitive(S); }
  bool IsCyclic() const { return Kernel_IsCyclic<Telt, Tidx_label, Tint>(S); }
  std::vector<Telt> SmallGeneratingSet() const {
    return Kernel_SmallGeneratingSet<Telt, Tidx_label, Tint>(S);
  }
  // Normal structure
  bool IsNormalSubgroup(const Group<Telt, Tint> &U) const {
    return Kernel_IsNormalSubgroup(S, U.S);
  }
  Group<Telt, Tint> NormalClosure(const Group<Telt, Tint> &H) const {
    return Group<Telt, Tint>(
        Kernel_NormalClosure<Telt, Tidx_label, Tint>(S, H.S));
  }
  Group<Telt, Tint> DerivedSubgroup() const {
    return Group<Telt, Tint>(Kernel_DerivedSubgroup<Telt, Tidx_label, Tint>(S));
  }
  Group<Telt, Tint> CentreSubgroup() const {
    return Group<Telt, Tint>(Kernel_CentreSubgroup<Telt, Tidx_label, Tint>(S));
  }
  Group<Telt, Tint> Centralizer_elt(const Telt &x) const {
    return Group<Telt, Tint>(
        Kernel_Centralizer_elt<Telt, Tidx_label, Tint>(S, x));
  }
  Group<Telt, Tint> Centralizer_grp(const Group<Telt, Tint> &H) const {
    return Group<Telt, Tint>(
        Kernel_Centralizer_grp<Telt, Tidx_label, Tint>(S, H.S));
  }
  std::vector<BlockDecomposition<Tidx>> GetSequenceBlockDecomposition() const {
    return ComputeSequenceBlockDecomposition(Kernel_GeneratorsOfGroup(S),
                                             S->comm->identity.size());
  }
  std::vector<Group<Telt, Tint>> GetAscendingChain() const {
    Telt id = S->comm->identity;
    std::vector<StabChain<Telt, Tidx_label>> l_stab =
        Kernel_AscendingChain<Telt, Tidx_label, Tint>(S);
    std::vector<Group<Telt, Tint>> l_grp;
    for (auto &e_s : l_stab) {
      Group<Telt, Tint> eGRP(e_s);
      l_grp.push_back(eGRP);
    }
    return l_grp;
  }
  Group<Telt, Tint> Intersection(Group<Telt, Tint> const &H) const {
    return Group<Telt, Tint>(
        Kernel_Intersection<Telt, Tidx_label, Tint>(S, H.S));
  }
  const Telt &get_identity() const { return S->comm->identity; }
  bool isin(const Telt &x) const { return IsElementInStabChain(S, x); }
  Telt Sift(Telt const &x) const { return SiftedPermutation(S, x); }

private:
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

public:
  using iterator = IteratorType;
  using const_iterator = IteratorType;
  using value_type = Telt;
  const_iterator begin() const {
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
    return IteratorType(ListS, ListPos, ListSiz, ListRes);
  }
  const_iterator end() const {
    std::vector<size_t> ListPos;
    StabChain<Telt, Tidx_label> Swork = S;
    while (Swork != nullptr) {
      ListPos.push_back(Swork->orbit.size());
      Swork = Swork->stabilizer;
    }
    return IteratorType({}, ListPos, {}, {});
  }

public:
  StabChain<Telt, Tidx_label> S;
  Tint size_tint;
};

template <typename Tgroup, typename TeltMatr, typename Tobj, typename Fop>
std::vector<TeltMatr>
PreImageSubgroupAction(std::vector<TeltMatr> const &ListMatrGens,
                       std::vector<typename Tgroup::Telt> const &ListPermGens,
                       TeltMatr const &id_matr, Tgroup const &stab,
                       Tobj const &x, Fop const &f_op) {
  using TeltPerm = typename Tgroup::Telt;
  using Tidx = typename TeltPerm::Tidx;
  using Telt = std::pair<TeltMatr, TeltPerm>;
  auto f_prod = [](Telt const &x, Telt const &y) -> Telt {
    return {x.first * y.first, x.second * y.second};
  };
  //
  TeltPerm id_perm = stab.get_identity();
  Tidx len = id_perm.size();
  Tgroup GRP(ListPermGens, len);
  auto f_act = [&](Tobj const &x, Telt const &u) -> Tobj {
    return f_op(x, u.second);
  };
  Telt id{id_matr, id_perm};
  std::vector<Telt> ListGens;
  for (size_t iGen = 0; iGen < ListMatrGens.size(); iGen++) {
    ListGens.push_back({ListMatrGens[iGen], ListPermGens[iGen]});
  }
  std::vector<std::pair<Tobj, Telt>> ListPair =
      OrbitPairEltRepr(ListGens, id, x, f_prod, f_act);
  std::unordered_map<Tobj, Telt> map;
  //  std::cerr << "f=" << f << "\n";
  for (auto &kv : ListPair) {
    map[kv.first] = kv.second;
    //    std::cerr << "kv.first=" << kv.first << "kv.second.second=" <<
    //    kv.second.second << "\n";
  }
  size_t nCoset = ListPair.size();
  //
  // We are using the Schreier lemma
  // See https://en.wikipedia.org/wiki/Schreier%27s_lemma
  //
  // The formula is the following for right coset. A right coset is something
  // like H r. If H is a subgroup of finite index and r1, ...., rN such that G =
  // H r1  \cup  .....  \cup  H rN If we are doing the action on the right and H
  // is the stabilizer then we get O = x G = { x r1 , ..... , x rN }
  //
  //
  //
  //
  std::unordered_set<TeltMatr> SetMatrGens;
  size_t nGen = ListMatrGens.size();
  std::cerr << "nCoset=" << nCoset << " |ListMatrGens|=" << nGen << "\n";
  for (size_t iCoset = 0; iCoset < nCoset; iCoset++) {
    Tobj const &x_cos = ListPair[iCoset].first;
    TeltMatr const &eGenMatr = ListPair[iCoset].second.first;
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
    TeltPerm const &eGenPerm = ListPair[iCoset].second.second;
#endif
    for (size_t iGen = 0; iGen < nGen; iGen++) {
      TeltMatr const &eGenMatr_B = ListMatrGens[iGen];
      TeltPerm const &eGenPerm_B = ListPermGens[iGen];
      Tobj x_img = f_op(x_cos, eGenPerm_B);
      Telt const &eElt = map[x_img];
      TeltMatr eGenMatr_new = eGenMatr * eGenMatr_B * Inverse(eElt.first);
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
      TeltPerm eGenPerm_new = eGenPerm * eGenPerm_B * Inverse(eElt.second);
      Tobj x_test = f_op(x, eGenPerm_new);
      if (x_test != x) {
        std::cerr << "iGen=" << iGen << " / " << nGen << "  iCoset=" << iCoset
                  << " / " << nCoset << "\n";
        std::cerr << "x_test=" << x_test << " x=" << x << "\n";
        std::cerr << "eGenPerm=" << eGenPerm << "\n";
        std::cerr << "eElt.second=" << eElt.second << "\n";
        std::cerr << "eGenPerm_new=" << eGenPerm_new << "\n";
        throw PermutalibException{1};
      }
#endif
      if (!IsIdentity(eGenMatr_new))
        SetMatrGens.insert(eGenMatr_new);
    }
  }
  std::cerr << "|SetMatrGens|=" << SetMatrGens.size() << "\n";
  std::vector<TeltMatr> ListMatrGens_ret;
  for (auto &eGen : SetMatrGens)
    ListMatrGens_ret.push_back(eGen);
  return ListMatrGens_ret;
}

template <typename TeltPerm, typename TeltMatr, typename Tint>
std::vector<TeltMatr>
StabilizerMatrixPermSubset(std::vector<TeltMatr> const &ListMatrGens,
                           std::vector<TeltPerm> const &ListPermGens,
                           TeltMatr const &id_matr, Face const &f) {
  using Tgroup = Group<TeltPerm, Tint>;
  using Tidx = typename TeltPerm::Tidx;
  using Tobj = Face;

  Tidx len = f.size();
  Tgroup GRP(ListPermGens, len);
  Tgroup stab = GRP.Stabilizer_OnSets(f);
  auto f_op = [&](Tobj const &x, TeltPerm const &u) -> Tobj {
    return OnSets(x, u);
  };
  return PreImageSubgroupAction<Tgroup, TeltMatr, Tobj, decltype(f_op)>(
      ListMatrGens, ListPermGens, id_matr, stab, f, f_op);
}

template <typename Tgroup, typename TeltMatr>
std::vector<TeltMatr>
PreImageSubgroup(std::vector<TeltMatr> const &ListMatrGens,
                 std::vector<typename Tgroup::Telt> const &ListPermGens,
                 TeltMatr const &id_matr, Tgroup const &eGRP) {
  using Telt = typename Tgroup::Telt;
  using Tobj = Telt;
  auto f_op = [&](Tobj const &x, Telt const &u) -> Tobj {
    return eGRP.Sift(Inverse(u) * x);
  };
  Telt id = eGRP.get_identity();
  return PreImageSubgroupAction<Tgroup, TeltMatr, Tobj, decltype(f_op)>(
      ListMatrGens, ListPermGens, id_matr, eGRP, id, f_op);
}

template <typename TeltPerm, typename TeltMatr, typename Tint>
std::optional<TeltMatr>
RepresentativeActionMatrixPermSubset(std::vector<TeltMatr> const &ListMatrGens,
                                     std::vector<TeltPerm> const &ListPermGens,
                                     TeltMatr const &id_matr, Face const &f1,
                                     Face const &f2) {
#ifdef DEBUG_REPRESENTATIVE_ACTION_MATRIX_PERM_SUBSET
  std::cerr << "Beginning of RepresentativeActionMatrixPermSubset\n";
#endif
  using Tidx = typename TeltPerm::Tidx;
  using Tgroup = Group<TeltPerm, Tint>;
  //
  Tidx len = f1.size();
  Tgroup GRP(ListPermGens, len);
#ifdef DEBUG_REPRESENTATIVE_ACTION_MATRIX_PERM_SUBSET
  std::cerr << "We have GRP\n";
  std::cerr << "GRP=" << GRP.GapString() << "\n";
#endif
  std::optional<TeltPerm> opt = GRP.RepresentativeAction_OnSets(f1, f2);
#ifdef DEBUG_REPRESENTATIVE_ACTION_MATRIX_PERM_SUBSET
  std::cerr << "We have opt\n";
#endif
  if (!opt)
    return {};
  // If we allow ourselves to compute RepresentativeAction
  Tgroup TheStab = GRP.Stabilizer_OnSets(f1);
  Tint OrbitSize = GRP.order() / TheStab.order();
  Tint CritSize = 10000;
  TeltPerm const &elt = *opt;
  size_t nGen = ListPermGens.size();
  auto f_mapping_elt=[&]() -> TeltMatr {
    //
    TeltPerm id_perm(len);
    using Tseq = SequenceType<true>;
    using Telt = PermutationElt<Tidx, Tseq>;
    using TgroupB = Group<Telt, Tint>;
    Tseq id_seq;
    Telt ePair(elt.getListVal(), id_seq);
    Telt idB(id_perm.getListVal(), id_seq);
    std::vector<Telt> ListGensB;
    std::vector<TeltMatr> ListMatrGens_inv(nGen);
    for (size_t iGen = 0; iGen < nGen; iGen++) {
      std::vector<int64_t> ListIdx{int64_t(iGen) + 1};
      Tseq e_seq(ListIdx);
      Telt fPair(ListPermGens[iGen].getListVal(), e_seq);
      ListGensB.push_back(fPair);
      ListMatrGens_inv[iGen] = Inverse(ListMatrGens[iGen]);
    }
#ifdef DEBUG_REPRESENTATIVE_ACTION_MATRIX_PERM_SUBSET
    std::cerr << "We have ListGensB\n";
#endif
    TgroupB GRP_B(ListGensB, idB);
#ifdef DEBUG_REPRESENTATIVE_ACTION_MATRIX_PERM_SUBSET
    std::cerr << "We have GRP_B\n";
#endif
    Telt res = GRP_B.Sift(ePair);
#ifdef DEBUG_REPRESENTATIVE_ACTION_MATRIX_PERM_SUBSET
    std::cerr << "We have res\n";
#endif
    //  NicePrint("res", res);
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
    //  std::cerr << "Doing the check\n";
    std::vector<Tidx> const &V = res.getListVal();
    for (Tidx u = 0; u < len; u++) {
      if (V[u] != u) {
        std::cerr << "The permutation residue is not the idenity at u=" << u
                  << "\n";
        throw PermutalibException{1};
      }
    }
#endif
    Tseq ret_seq = Inverse(res.getElt());
#ifdef DEBUG_REPRESENTATIVE_ACTION_MATRIX_PERM_SUBSET
    std::cerr << "We have ret_seq\n";
#endif
    TeltMatr ret_matr = id_matr;
    const std::vector<int64_t>& ListIdx = ret_seq.getVect();
    for (auto & eIdx : ListIdx) {
      if (eIdx > 0) {
        size_t iGen = eIdx - 1;
        ret_matr *= ListMatrGens[iGen];
      } else {
        size_t iGen = (-eIdx) - 1;
        ret_matr *= ListMatrGens_inv[iGen];
      }
    }
#ifdef DEBUG_REPRESENTATIVE_ACTION_MATRIX_PERM_SUBSET
    std::cerr << "We have ret_matr\n";
#endif
    //  std::cerr << "Returning from RepresentativeActionMatrixPermSubset\n";
    return ret_matr;
  };
  auto f_build_orbit=[&]() -> TeltMatr {
    size_t miss_val = std::numeric_limits<size_t>::max();
    std::unordered_set<Face> set;
    std::vector<std::tuple<Face,size_t,size_t>> l_x_iorig_igen;
    auto f_insert=[&](Face const& x, size_t iOrig, size_t iGen) -> void {
      if (set.count(x) == 1)
        return;
      set.insert(x);
      std::tuple<Face,size_t,size_t> tuple{x,iOrig,iGen};
      l_x_iorig_igen.push_back(tuple);
    };
    auto f_get_elt=[&](size_t const& pos) -> TeltMatr {
      size_t curr_pos = pos;
      std::vector<size_t> ListIGen;
      while (true) {
        if (curr_pos == 0)
          break;
        size_t iOrig = std::get<1>(l_x_iorig_igen[curr_pos]);
        size_t iGen = std::get<2>(l_x_iorig_igen[curr_pos]);
        ListIGen.push_back(iGen);
        curr_pos = iOrig;
      }
      TeltMatr ret_matr = id_matr;
      size_t len = ListIGen.size();
      for (size_t i=0; i<len; i++) {
        size_t iGen = ListIGen[len - 1 - i];
        ret_matr *= ListMatrGens[iGen];
      }
      return ret_matr;
    };
    f_insert(f1, miss_val, miss_val);
    size_t n_done = 0;
    while(true) {
      size_t len = l_x_iorig_igen.size();
      if (n_done == len)
        break;
      for (size_t iOrig=n_done; iOrig<len; iOrig++) {
        Face const& f = std::get<0>(l_x_iorig_igen[iOrig]);
        for (size_t iGen=0; iGen<nGen; iGen++) {
          Face f_img = FaceAct(f, ListPermGens[iGen]);
          f_insert(f_img, iOrig, iGen);
          if (TestEqual(f_img, f2)) {
            size_t pos = l_x_iorig_igen.size() - 1;
            return f_get_elt(pos);
          }
        }
      }
      n_done = len;
    }
    std::cerr << "We should never reach that stage\n";
    throw PermutalibException{1};
  };

  if (OrbitSize < CritSize) {
    return f_build_orbit();
  } else {
    return f_mapping_elt();
  }
}

// clang-format off
}  // namespace permutalib
// clang-format on

namespace boost::serialization {

template <class Archive, typename Telt, typename Tint>
inline void load(Archive &ar, permutalib::Group<Telt, Tint> &val,
                 [[maybe_unused]] const unsigned int version) {
  using Tidx = typename Telt::Tidx;
  Tidx n_act;
  size_t n_gen;
  ar &make_nvp("n_act", n_act);
  ar &make_nvp("n_gen", n_gen);
  std::vector<Telt> LGen;
  LGen.reserve(n_gen);
  for (size_t i_gen = 0; i_gen < n_gen; i_gen++) {
    std::vector<Tidx> eList(n_act);
    for (Tidx i = 0; i < n_act; i++) {
      ar &make_nvp("val", eList[i]);
    }
    Telt eGen(eList);
    LGen.emplace_back(std::move(eGen));
  }
  val = permutalib::Group<Telt, Tint>(LGen, n_act);
}

template <class Archive, typename Telt, typename Tint>
inline void save(Archive &ar, permutalib::Group<Telt, Tint> const &val,
                 [[maybe_unused]] const unsigned int version) {
  using Tidx = typename Telt::Tidx;
  Tidx n_act = val.n_act();
  ar &make_nvp("n_act", n_act);
  std::vector<Telt> LGen = val.GeneratorsOfGroup();
  size_t n_gen = LGen.size();
  ar &make_nvp("n_gen", n_gen);
  for (size_t i_gen = 0; i_gen < n_gen; i_gen++) {
    const Telt &eGen = LGen[i_gen];
    for (Tidx i = 0; i < n_act; i++) {
      Tidx pnt = eGen.at(i);
      ar &make_nvp("val", pnt);
    }
  }
}

template <class Archive, typename Telt, typename Tint>
inline void serialize(Archive &ar, permutalib::Group<Telt, Tint> &val,
                      const unsigned int version) {
  split_free(ar, val, version);
}

template <typename Telt, typename Tint>
std::ostream &operator<<(std::ostream &os,
                         const permutalib::Group<Telt, Tint> &grp) {
  boost::archive::text_oarchive oa(os);
  oa << grp;
  return os;
}

template <typename Telt, typename Tint>
std::istream &operator>>(std::istream &is, permutalib::Group<Telt, Tint> &grp) {
  boost::archive::text_iarchive ia(is);
  ia >> grp;
  return is;
}

// clang-format off
}  // namespace permutalib
#endif  // SRC_GAP_GROUP_H_
// clang-format on
