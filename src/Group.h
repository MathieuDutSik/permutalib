// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GAP_GROUP_H_
#define SRC_GAP_GROUP_H_

// clang-format off
#include "StabChainMain.h"
#include "IteratingElement.h"
#include "stbcbckt.h"
#include "nsi.h"
#include "Properties.h"
#include "NormalStructure.h"
#include "PermutationElt.h"
#include "AscendingChains_and_Cosets.h"
// clang-format on
#include <map>
#include <fstream>
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

#ifdef DEBUG
#define DEBUG_STABCHAINMAIN
#define DEBUG_REPRESENTATIVE_ACTION_MATRIX_PERM_SUBSET
#endif

#ifdef TIMINGS
#define TIMINGS_PRE_IMAGE
#endif

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

template <typename Telt, typename Tint> struct Group;

template <typename Telt, typename Tint>
struct KernelDoubleCosetComputer {
private:
  using Tidx_label = uint16_t;
  InnerDoubleCosetComputer<Telt, Tidx_label, Tint> inner;
  bool option;
public:
  KernelDoubleCosetComputer(StabChain<Telt,Tidx_label> const& G, StabChain<Telt,Tidx_label> const& U, bool const& _option) : inner(G, U), option(_option) {
  }
  std::vector<Telt> double_cosets(Group<Telt,Tint> const& U_or_V) const {
    if (option) {
      return inner.double_cosets(U_or_V.stab_chain());
    } else {
      std::vector<Telt> l_cos;
      for (auto & eCos: inner.double_cosets(U_or_V.stab_chain())) {
        l_cos.push_back(Inverse(eCos));
      }
      return l_cos;
    }
  }
  std::vector<KernelDccEntry<Telt>> double_cosets_and_stabilizers(Group<Telt,Tint> const& V) const {
    if (!option) {
      std::cerr << "GRP: The function can only be used on the V side\n";
      throw PermutalibException{1};
    }
    return inner.double_cosets_and_stabilizers(V.stab_chain());
  }
};

template <typename Telt_inp, typename Tint_inp> struct Group {
public:
  // dependent types
  using Telt = Telt_inp;
  using Tidx = typename Telt::Tidx;
  using Tint = Tint_inp;
  using Tidx_label = uint16_t;
  using RightCosets = KernelRightCosets<Telt,Tidx_label,Tint>;
  using LeftCosets = KernelLeftCosets<Telt,Tidx_label,Tint>;
  using DoubleCosetComputer = KernelDoubleCosetComputer<Telt,Tint>;
  using DccEntry = KernelDccEntry<Telt>;
private:
  StabChain<Telt, Tidx_label> S;
  Tint size_tint;
  bool use_store_canonic;
  void set_exhaustive_canonic() {
    if (size_tint < 2000)
      use_store_canonic = true;
    else
      use_store_canonic = false;
  }
  mutable std::vector<Telt> l_group_elt;
  void compute_all_element() const {
    l_group_elt = get_all_elements(S);
  }
  mutable std::optional<std::vector<Telt>> SmallGenSet;
public:
  // constructors
  Group(const StabChain<Telt, Tidx_label> &_S)
      : S(_S), size_tint(Order<Telt, Tidx_label, Tint>(_S)) {
    set_exhaustive_canonic();
  }
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
    set_exhaustive_canonic();
  }
  Group(const Tidx &n) : Group({}, n) {
    set_exhaustive_canonic();
  }
  Group() : Group(0) {
    set_exhaustive_canonic();
  }
  Group(Group<Telt, Tint> &&G) : S(std::move(G.S)), size_tint(G.size_tint) {
    set_exhaustive_canonic();
  }
  Group(const Group<Telt, Tint> &G) : S(G.S), size_tint(G.size_tint) {
    set_exhaustive_canonic();
  }
  // Standard operators
  bool operator==(const Group &g) const { return EqualityTest(S, g.S); }
  bool operator!=(const Group &g) const { return !EqualityTest(S, g.S); }
  Group<Telt, Tint> &operator=(const Group<Telt, Tint> &G) {
    // The S is a shared_ptr so copy is fine.
    S = G.S;
    size_tint = G.size_tint;
    set_exhaustive_canonic();
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
  std::string PythonString() const {
    std::vector<Telt> LGen = Kernel_GeneratorsOfGroup(S);
    std::string str_ret = "[";
    Tidx n_act = S->comm->n;
    bool IsFirst = true;
    for (auto &eGen : LGen) {
      if (!IsFirst) {
        str_ret += ",";
      }
      str_ret += "[";
      for (Tidx i=0; i<n_act; i++) {
        Tidx pnt = eGen.at(i);
        size_t pnt_s = static_cast<size_t>(pnt);
        str_ret += std::to_string(pnt_s);
      }
      str_ret += "]";
    }
    str_ret += "]";
    return str_ret;
  }
  Tint size() const { return size_tint; }
  StabChain<Telt, Tidx_label> stab_chain() const { return S; }
  std::map<Tidx, int> factor_size() const { return FactorsSizeStabChain(S); }
  Tidx n_act() const { return S->comm->n; }
  std::vector<Telt> get_all_element() const {
    if (l_group_elt.size() == 0) {
      compute_all_element();
    }
    return l_group_elt;
  }
  // Operation and creating new groups
  Group<Telt, Tint> GroupConjugate(const Telt &x) const {
    std::vector<Telt> LGen;
    Telt xInv = ~x;
    for (auto &eGen : Kernel_GeneratorsOfGroup(S)) {
      Telt eGenCj = xInv * eGen * x;
      LGen.emplace_back(eGenCj);
    }
    return Group<Telt, Tint>(LGen, S->comm->n);
  }
  // Canonical images
  Face CanonicalImage(const Face &f) const {
    return Kernel_CanonicalImage<Telt, Tidx_label, Tint>(S, f);
  }
  Face ExhaustiveCanonicalImage(const Face &f) const {
    return exhaustive_minimum_face_orbit<Telt,Tidx_label>(S, f);
  }
  Face StoreCanonicalImage(const Face &f) const {
    if (l_group_elt.size() == 0) {
      compute_all_element();
    }
    Face f_minimum = f;
    for (auto const& eElt : l_group_elt) {
      Face f_img = OnSets(f, eElt);
      if (f_img < f_minimum)
        f_minimum = f_img;
    }
    return f_minimum;
  }
  std::pair<Face,Tint> StoreCanonicalImageOrbitSize(const Face &f) const {
    if (l_group_elt.size() == 0) {
      compute_all_element();
    }
    Tidx n = n_act();
    auto comp=[&](Face const& f1, Face const& f2) -> int8_t {
      for (Tidx i=0; i<n; i++) {
        if (f1[i] < f2[i])
          return 1;
        if (f1[i] > f2[i])
          return -1;
      }
      return 0;
    };
    Face f_minimum(n);
    for (Tidx i=0; i<n; i++)
      f_minimum[i] = 1;
    int size = 0;
    for (auto const& eElt : l_group_elt) {
      Face f_img = OnSets(f, eElt);
      int8_t test = comp(f_img, f_minimum);
      if (test == 1) {
        f_minimum = f_img;
        size = 1;
      } else {
        if (test == 0) {
          size++;
        }
      }
    }
    Tint stabSize(size);
    Tint orbitSize = size_tint / stabSize;
    return {f_minimum, orbitSize};
  }
  Face OptCanonicalImage(const Face &f) const {
    if (use_store_canonic)
      return StoreCanonicalImage(f);
    else
      return Kernel_CanonicalImage<Telt, Tidx_label, Tint>(S, f);
  }
  std::pair<Face,Group<Telt,Tint>> PairCanonicalImageSubgroupStabilizer(const Face &f) const {
    std::pair<Face,StabChain<Telt,Tidx_label>> pairCan = CanonicalImage_SubgroupStabilizer<Telt,Tidx_label,Tint>(S, f);
    return {std::move(pairCan.first), Group(std::move(pairCan.second))};
  }
  std::pair<Face,Tint> CanonicalImageOrbitSize(const Face &f) const {
    std::pair<Face,StabChain<Telt,Tidx_label>> pairCan = CanonicalImage_ConjugateStabilizer<Telt,Tidx_label,Tint>(S, f);
    Tint StabSize = Order<Telt, Tidx_label, Tint>(pairCan.second);
    Tint OrbitSize = size_tint / StabSize;
    return {std::move(pairCan.first), OrbitSize};
  }
  std::pair<Face,Tint> OptCanonicalImageOrbitSize(const Face &f) const {
    if (use_store_canonic)
      return StoreCanonicalImageOrbitSize(f);
    else
      return CanonicalImageOrbitSize(f);
  }
  Face CanonicalImageInitialTriv(const Face &f) const {
    size_t max_size = std::numeric_limits<size_t>::max();
    return Kernel_GeneralCanonicalInitialTriv<Telt, Tidx_label, Tint>(S, f, max_size).first;
  }
  Face CanonicalImageInitialTrivLimited(const Face &f, size_t const& max_size) const {
    return Kernel_GeneralCanonicalInitialTriv<Telt, Tidx_label, Tint>(S, f, max_size).first;
  }
  size_t CanonicalImageInitialTrivTreeDepth(const Face &f) const {
    size_t max_size = std::numeric_limits<size_t>::max();
    return Kernel_GeneralCanonicalInitialTriv<Telt, Tidx_label, Tint>(S, f, max_size).second;
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
  Tint OrbitSize_OnSets(const Face &f) const {
    StabChain<Telt,Tidx_label> eStab = Kernel_Stabilizer_OnSets<Telt, Tidx_label, Tint>(S, f);
    Tint OrbitSize = size_tint / Order<Telt, Tidx_label, Tint>(eStab);
    return OrbitSize;
  }
  // Random elements and subgroups
  Telt rand() const {
    // Uses the generators and move them at random.
    return RandomElement(Kernel_GeneratorsOfGroup(S), S->comm->identity);
  }
  Telt uniform_rand() const {
    // Uses the stabilizer chain in order to get random elements.
    return UniformRandomElement(S);
  }
  Group<Telt, Tint> RandomSubgroup() const {
    std::vector<Telt> LGen = UsefulRandomSubgroupGenerators(S);
    return Group<Telt, Tint>(LGen, S->comm->n);
  }
  Telt random() const { return rand(); }
  // Properties of the group
  bool IsCommutative() const { return Kernel_IsCommutative(S); }
  bool IsTransitive() const { return Kernel_IsTransitive(S); }
  bool IsPrimitive() const { return Kernel_IsPrimitive(S); }
  bool IsCyclic() const { return Kernel_IsCyclic<Telt, Tidx_label, Tint>(S); }
  std::vector<Telt> SmallGeneratingSet() const {
    if (!SmallGenSet) {
      SmallGenSet = Kernel_SmallGeneratingSet<Telt, Tidx_label, Tint>(S);
    }
    return *SmallGenSet;
  }
  // Compute cosets
  RightCosets right_cosets(const Group<Telt,Tint>& H) const {
    return KernelRightCosets<Telt,Tidx_label,Tint>(H.S, S);
  }
  LeftCosets left_cosets(const Group<Telt,Tint>& H) const {
    return KernelLeftCosets<Telt,Tidx_label,Tint>(H.S, S);
  }
  std::vector<Telt> get_all_right_cosets(const Group<Telt,Tint>& H) const {
    return enumerate_right_cosets<Telt,Tidx_label,Tint>(H.S, S);
  }
  std::vector<Telt> get_all_left_cosets(const Group<Telt,Tint>& H) const {
    return enumerate_left_cosets<Telt,Tidx_label,Tint>(H.S, S);
  }
  DoubleCosetComputer double_coset_computer_v(const Group<Telt,Tint>& U) const {
    return KernelDoubleCosetComputer<Telt,Tint>(S, U.S, true);
  }
  DoubleCosetComputer double_coset_computer_u(const Group<Telt,Tint>& V) const {
    return KernelDoubleCosetComputer<Telt,Tint>(S, V.S, false);
  }
  std::vector<Telt> double_cosets(const Group<Telt,Tint>& U, const Group<Telt,Tint>& V) const {
    if (U.size_tint > V.size_tint) {
      KernelDoubleCosetComputer<Telt,Tint> dcc_v(S, U.S, true);
      return dcc_v.double_cosets(V);
    } else {
      KernelDoubleCosetComputer<Telt,Tint> dcc_u(S, V.S, false);
      return dcc_u.double_cosets(U);
    }
  }
  std::vector<DccEntry> double_cosets_and_stabilizers(const Group<Telt,Tint>& U, const Group<Telt,Tint>& V) const {
    KernelDoubleCosetComputer<Telt,Tint> dcc_v(S, U.S, true);
    return dcc_v.double_cosets_and_stabilizers(V);
  }
  // Normal structure
  bool IsNormalSubgroup(const Group<Telt, Tint> &U) const {
    return Kernel_IsNormalSubgroup(S, U.S);
  }
  bool IsSubgroup(const Group<Telt,Tint>& U) const {
    return Kernel_IsSubgroup(S, U.S);
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
        Kernel_AscendingChainSingle<Telt, Tidx_label, Tint>(S);
    std::vector<Group<Telt, Tint>> l_grp;
    for (auto &e_s : l_stab) {
      Group<Telt, Tint> eGRP(e_s);
      l_grp.push_back(eGRP);
    }
    return l_grp;
  }
  std::vector<Group<Telt, Tint>> GetAscendingChainSubgroup(Group<Telt,Tint> const& H) const {
    Telt id = S->comm->identity;
    std::vector<StabChain<Telt, Tidx_label>> l_stab =
      Kernel_AscendingChainPair<Telt, Tidx_label, Tint>(H.S, S);
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
  //
  // The iterator business
  //
  using iterator = IteratorType<Telt,Tidx_label>;
  using const_iterator = IteratorType<Telt,Tidx_label>;
  using value_type = Telt;
  const_iterator begin() const {
    return get_begin_iterator(S);
  }
  const_iterator end() const {
    return get_end_iterator(S);
  }
};


// The function f_op needs to satisfy
// f_op(f_op(x, u), v)  =  f_op(x, u * v)
template <typename TeltPerm, typename TeltMatr, typename Tobj, typename Fop, typename FinsertGen>
std::vector<std::pair<Tobj, std::pair<TeltMatr, TeltPerm>>>
PreImageSubgroupActionGenA(std::vector<TeltMatr> const &ListMatrGens,
                           std::vector<TeltPerm> const &ListPermGens,
                           TeltMatr const &id_matr, TeltPerm const &id_perm,
                           Tobj const &x_start, Fop const &f_op,
                           FinsertGen const& f_insert_gen) {
  using Telt = std::pair<TeltMatr, TeltPerm>;
  auto f_prod = [](Telt const &x, Telt const &y) -> Telt {
    return {x.first * y.first, x.second * y.second};
  };
  //
  auto f_act = [&](Tobj const &x, Telt const &u) -> Tobj {
    return f_op(x, u.second);
  };
#ifdef TIMINGS_PRE_IMAGE
  MicrosecondTime_perm time;
#endif
  Telt id{id_matr, id_perm};
  std::vector<Telt> ListGens;
  for (size_t iGen = 0; iGen < ListMatrGens.size(); iGen++) {
    ListGens.push_back({ListMatrGens[iGen], ListPermGens[iGen]});
  }
  std::vector<std::pair<Tobj, Telt>> ListPair =
      OrbitPairEltRepr(ListGens, id, x_start, f_prod, f_act);
#ifdef TIMINGS_PRE_IMAGE
  std::cerr << "|GRP: PreImageSubgroupActionGenA, OrbitPairEltRepr|=" << time << "\n";
#endif
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
  for (auto & ePair : ListPair) {
    Tobj img = f_act(x_start, ePair.second);
    if (img != ePair.first) {
      std::cerr << "GRP: |ListPair|=" << ListPair.size() << "\n";
      std::cerr << "GRP: img does not match the element\n";
      throw PermutalibException{1};
    }
  }
# ifdef TIMINGS_PRE_IMAGE
  std::cerr << "|GRP: PreImageSubgroupActionGenA, check ListPair|=" << time << "\n";
# endif
#endif
  std::unordered_map<Tobj, Telt> map;
  for (auto &kv : ListPair) {
    map[kv.first] = kv.second;
  }
#ifdef TIMINGS_PRE_IMAGE
  std::cerr << "|GRP: PreImageSubgroupActionGenA, map|=" << time << "\n";
#endif
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
  size_t nGen = ListMatrGens.size();
#ifdef TIMINGS_PRE_IMAGE
  std::cerr << "|GRP: PreImageSubgroupActionGenA, nGen|=" << time << "\n";
#endif
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
  std::cerr << "GRP: nCoset=" << nCoset << " |ListMatrGens|=" << nGen << "\n";
#endif
  for (size_t iCoset = 0; iCoset < nCoset; iCoset++) {
    Tobj const &x_cos = ListPair[iCoset].first;
    TeltMatr const &eCosMatr = ListPair[iCoset].second.first;
    TeltPerm const &eCosPerm = ListPair[iCoset].second.second;
    for (size_t iGen = 0; iGen < nGen; iGen++) {
      TeltMatr const &eGenMatr = ListMatrGens[iGen];
      TeltPerm const &eGenPerm = ListPermGens[iGen];
      Tobj x_img = f_op(x_cos, eGenPerm);
      Telt const &eElt = map[x_img];
      TeltMatr eGenStabMatr = eCosMatr * eGenMatr * Inverse(eElt.first);
      TeltPerm eGenStabPerm = eCosPerm * eGenPerm * Inverse(eElt.second);
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
      Tobj x_test = f_op(x_start, eGenStabPerm);
      if (x_test != x_start) {
        std::cerr << "GRP: iGen=" << iGen << " / " << nGen << "  iCoset=" << iCoset
                  << " / " << nCoset << "\n";
        std::cerr << "GRP: x_test=" << x_test << " x_start=" << x_start << "\n";
        std::cerr << "GRP: eCosPerm=" << eCosPerm << "\n";
        std::cerr << "GRP: eElt.second=" << eElt.second << "\n";
        std::cerr << "GRP: eGenStabPerm=" << eGenStabPerm << "\n";
        throw PermutalibException{1};
      }
#endif
      std::pair<TeltMatr, TeltPerm> pair{eGenStabMatr, eGenStabPerm};
      f_insert_gen(pair);
    }
  }
#ifdef TIMINGS_PRE_IMAGE
  std::cerr << "|GRP: PreImageSubgroupActionGenA, f_insert oper|=" << time << "\n";
#endif
  return ListPair;
}

template <typename TeltPerm, typename TeltMatr, typename Tobj, typename Fop>
std::pair<std::vector<TeltMatr>, std::vector<std::pair<Tobj, std::pair<TeltMatr, TeltPerm>>>>
PreImageSubgroupActionGen(std::vector<TeltMatr> const &ListMatrGens,
                          std::vector<TeltPerm> const &ListPermGens,
                          TeltMatr const &id_matr, TeltPerm const &id_perm,
                          Tobj const &x_start, Fop const &f_op) {
  using Telt = std::pair<TeltMatr, TeltPerm>;
  std::unordered_set<TeltMatr> SetMatrGens;
  auto f_insert_gen=[&](std::pair<TeltMatr, TeltPerm> const& pair) -> void {
    if (!IsIdentity(pair.first)) {
      SetMatrGens.insert(pair.first);
    }
  };
  std::vector<std::pair<Tobj, Telt>> ListPair =
    PreImageSubgroupActionGenA<TeltPerm,TeltMatr,Tobj,Fop,decltype(f_insert_gen)>(ListMatrGens,
                                                                                  ListPermGens,
                                                                                  id_matr, id_perm,
                                                                                  x_start, f_op, f_insert_gen);
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
  std::cerr << "GRP: |SetMatrGens|=" << SetMatrGens.size() << "\n";
#endif
  std::vector<TeltMatr> VectMatrGens;
  for (auto &eGen : SetMatrGens) {
    VectMatrGens.push_back(eGen);
  }
  return {std::move(VectMatrGens), std::move(ListPair)};
}

template<typename TeltMatr>
struct SeqTracker {
  TeltMatr x;
  SequenceType<true> seq;
  SeqTracker() {
    std::vector<int64_t> ListIdx;
    seq = SequenceType<true>(ListIdx);
  }
  SeqTracker(TeltMatr const& _x, SequenceType<true> const& _seq) : x(_x), seq(_seq) {
  }
};

template<typename TeltMatr>
bool IsIdentity(SeqTracker<TeltMatr> const& v) {
  return IsIdentity(v.x);
}

template<typename TeltMatr>
SeqTracker<TeltMatr> operator*(SeqTracker<TeltMatr> const& v1, SeqTracker<TeltMatr> const& v2) {
  TeltMatr x = v1.x * v2.x;
  SequenceType<true> seq = v1.seq * v2.seq;
  return SeqTracker(x, seq);
}

template<typename TeltMatr>
bool operator==(SeqTracker<TeltMatr> const& v1, SeqTracker<TeltMatr> const& v2) {
  return v1.x == v2.x;
}

template<typename TeltMatr>
SeqTracker<TeltMatr> Inverse(SeqTracker<TeltMatr> const& v) {
  TeltMatr x = Inverse(v.x);
  SequenceType<true> seq = Inverse(v.seq);
  return SeqTracker(x, seq);
}

template <typename TeltPerm, typename TeltMatr, typename Tobj, typename Fop>
std::vector<TeltMatr>
PreImageSubgroupAction(std::vector<TeltMatr> const &ListMatrGens,
                       std::vector<TeltPerm> const &ListPermGens,
                       TeltMatr const &id_matr, TeltPerm const &id_perm,
                       Tobj const &x, Fop const &f_op) {
#ifdef DEBUG_PRE_IMAGE_COMPLEXITY
  using Tseq = SequenceType<true>;
  using TeltMatrComb = SeqTracker<TeltMatr>;
  std::vector<TeltMatrComb> ListMatrGensComb;
  for (size_t i_elt=0; i_elt<ListMatrGens.size(); i_elt++) {
    std::vector<int64_t> ListIdx{int64_t(i_elt) + 1};
    Tseq seq(ListIdx);
    TeltMatrComb eComb(ListMatrGens[i_elt], seq);
    ListMatrGensComb.push_back(eComb);
  }
  std::vector<int64_t> ListIdx{};
  Tseq seq(ListIdx);
  TeltMatrComb id_matrComb(id_matr, seq);
  std::pair<std::vector<TeltMatrComb>,std::vector<std::pair<Tobj, std::pair<TeltMatrComb, TeltPerm>>>> pair =
    PreImageSubgroupActionGen<TeltPerm,TeltMatrComb,Tobj,Fop>(ListMatrGensComb,
                              ListPermGens,
                              id_matrComb, id_perm,
                              x, f_op);
  size_t l1_complexity = 0;
  size_t linf_complexity = 0;
  std::vector<TeltMatr> list_ret;
  for (auto & epair: pair.first) {
    list_ret.push_back(epair.x);
    size_t comp = epair.seq.complexity();
    l1_complexity += comp;
    linf_complexity = std::max(linf_complexity, comp);
  }
  std::cerr << "permutalib::PreImageSubgroupAction, l1_complexity=" << l1_complexity << " linf_complexity=" << linf_complexity << "\n";
  return list_ret;
#else
  std::pair<std::vector<TeltMatr>,std::vector<std::pair<Tobj, std::pair<TeltMatr, TeltPerm>>>> pair =
    PreImageSubgroupActionGen<TeltPerm,TeltMatr,Tobj,Fop>(ListMatrGens,
                              ListPermGens,
                              id_matr, id_perm,
                              x, f_op);
  return pair.first;
#endif
}

template <typename Tgroup, typename TeltMatr, typename Tobj, typename Fop>
std::pair<std::vector<TeltMatr>,std::vector<TeltMatr>>
PreImageSubgroupRightCosetAction(std::vector<TeltMatr> const &ListMatrGens,
                                 std::vector<typename Tgroup::Telt> const &ListPermGens,
                                 TeltMatr const &id_matr, Tgroup const &stab,
                                 Tobj const &x, Fop const &f_op) {
  using TeltPerm = typename Tgroup::Telt;
  TeltPerm id_perm = stab.get_identity();
  std::pair<std::vector<TeltMatr>,std::vector<std::pair<Tobj, std::pair<TeltMatr, TeltPerm>>>> pair =
    PreImageSubgroupActionGen(ListMatrGens,
                              ListPermGens,
                              id_matr, id_perm,
                              x, f_op);
  // Is it Right or Left cosets? Unclear at present.
  // But we for sure really want the right cosets.
  std::vector<TeltMatr> RightCosets;
  for (auto & epair : pair.second) {
    RightCosets.emplace_back(std::move(epair.second.first));
  }
  return {std::move(pair.first), std::move(RightCosets)};
}

// We compute the right cosets which are of the form
// G = \cup_{c\in C} H c
// with C the set of right cosets.
//
// If the number of cosets is too large, we can maybe
// use an ascending chain.
template <typename Tgroup, typename TeltMatr, typename Finsert>
void PreImageSubgroupKernel(std::vector<TeltMatr> const &ListMatrGens,
                            std::vector<typename Tgroup::Telt> const &ListPermGens,
                            TeltMatr const &id_matr, Tgroup const &eGRP,
                            Finsert const& f_insert) {
  using TeltPerm = typename Tgroup::Telt;
  using Tidx = typename TeltPerm::Tidx;
  using Tobj = size_t;
#ifdef TIMINGS_PRE_IMAGE
  MicrosecondTime_perm time;
#endif
  Tidx n_act = eGRP.n_act();
  TeltPerm id_perm = eGRP.get_identity();
  Tgroup GRP_big(ListPermGens, n_act);
  std::vector<TeltPerm> l_cos = GRP_big.get_all_right_cosets(eGRP);
#ifdef TIMINGS_PRE_IMAGE
  std::cerr << "|GRP: PreImageSubgroupKernel, l_cos|=" << time << "\n";
#endif
  size_t n_cos = l_cos.size();
  std::unordered_map<TeltPerm, size_t> map;
  auto f_can=[&](TeltPerm const& x) -> TeltPerm {
    return MinimalElementCosetStabChain(eGRP.stab_chain(), x);
  };
  for (size_t i_cos=0; i_cos<n_cos; i_cos++) {
    TeltPerm const& e_cos = l_cos[i_cos];
    TeltPerm f_cos = f_can(e_cos);
    map[f_cos] = i_cos;
  }
#ifdef TIMINGS_PRE_IMAGE
  std::cerr << "|GRP: PreImageSubgroupKernel, map|=" << time << "\n";
#endif
  auto f_op = [&](size_t const &x, TeltPerm const &u) -> Tobj {
    TeltPerm prod = l_cos[x] * u;
    TeltPerm prod_can = f_can(prod);
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
    if (map.count(prod_can) == 0) {
      std::cerr << "GRP: PreImageSubgroupKernel, missing entry for f_op\n";
      throw PermutalibException{1};
    }
#endif
    size_t pos = map.at(prod_can);
    return pos;
  };
  TeltPerm id = eGRP.get_identity();
  TeltPerm id_can = f_can(id);
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
  if (map.count(id_can) == 0) {
    std::cerr << "GRP: PreImageSubgroupKernel, missing entry in creation of id_can\n";
    throw PermutalibException{1};
  }
#endif
  size_t pos_id = map.at(id_can);
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
  std::cerr << "GRP: PreImageSubgroupKernel, pos_id=" << pos_id << "\n";
  std::cerr << "GRP: PreImageSubgroupKernel, |GRP_big|=" << GRP_big.size() << " |eGRP|=" << eGRP.size() << "\n";
  //  KernelCheckRightCosets<TeltPerm, uint16_t, typename Tgroup::Tint>(GRP_big.stab_chain(), eGRP.stab_chain(), l_cos);
  auto f_map_elt=[&](TeltPerm const& u) -> TeltPerm {
    std::vector<Tidx> eList;
    for (auto & eCos : l_cos) {
      TeltPerm prod = eCos * u;
      TeltPerm prod_can = f_can(prod);
      if (map.count(prod_can) == 0) {
        std::cerr << "GRP: PreImageSubgroupKernel, missing entry in creation of ListPermGens_cos\n";
        throw PermutalibException{1};
      }
      size_t pos = map.at(prod_can);
      eList.push_back(pos);
    }
    TeltPerm ePerm(eList);
    return ePerm;
  };
  size_t n_act_s = n_act;
  size_t n_act_tot = n_act_s + n_cos;
  std::cerr << "GRP: PreImageSubgroupKernel, n_act_tot=" << n_act_tot << "\n";
  auto f_big_map_elt=[&](TeltPerm const& u) -> TeltPerm {
    std::vector<Tidx> eListBig(n_act_s + n_cos);
    TeltPerm u_img = f_map_elt(u);
    for (Tidx i_act=0; i_act<n_act; i_act++) {
      eListBig[i_act] = u.at(i_act);
    }
    for (size_t i_cos=0; i_cos<n_cos; i_cos++) {
      eListBig[n_act_s + i_cos] = n_act_s + u_img.at(i_cos);
    }
    return TeltPerm(eListBig);
  };
  std::vector<TeltPerm> ListPermGens_cos_dir;
  for (auto & ePermGen : ListPermGens) {
    ListPermGens_cos_dir.push_back(f_big_map_elt(ePermGen));
  }
  Tgroup FullGRPcos_dir(ListPermGens_cos_dir, n_act + n_cos);
  if (FullGRPcos_dir.size() != GRP_big.size()) {
    std::cerr << "GRP: PreImageSubgroupKernel, FullGRPcos_dir should be equal to GRP_big in size\n";
    throw PermutalibException{1};
  }
  for (auto & u: eGRP.GeneratorsOfGroup()) {
    TeltPerm eImg = f_map_elt(u);
    if (eImg.at(pos_id) != pos_id) {
      std::cerr << "eImg should map pos_id to pos_id\n";
      throw PermutalibException{1};
    }
  }
#endif
  auto f_insert_gen=[&](std::pair<TeltMatr, TeltPerm> const& pair) -> void {
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
    if (!eGRP.isin(pair.second)) {
      std::cerr << "pair.second should belong to the subgroup eGRP\n";
      throw PermutalibException{1};
    }
#endif
    if (!IsIdentity(pair.first)) {
      f_insert(pair);
    }
  };
#ifdef TIMINGS_PRE_IMAGE
  std::cerr << "|GRP: PreImageSubgroupKernel, pre_work|=" << time << "\n";
#endif
  PreImageSubgroupActionGenA<TeltPerm,TeltMatr,Tobj,decltype(f_op),decltype(f_insert_gen)>(ListMatrGens,
                                                                                           ListPermGens,
                                                                                           id_matr, id_perm,
                                                                                           pos_id, f_op, f_insert_gen);
}

template <typename Tgroup, typename TeltMatr>
std::vector<TeltMatr>
PreImageSubgroup(std::vector<TeltMatr> const &ListMatrGens,
                 std::vector<typename Tgroup::Telt> const &ListPermGens,
                 TeltMatr const &id_matr, Tgroup const &eGRP) {
  using TeltPerm = typename Tgroup::Telt;
  std::unordered_set<TeltMatr> SetMatrGens;
  auto f_insert=[&](std::pair<TeltMatr, TeltPerm> const& pair) -> void {
    SetMatrGens.insert(pair.first);
  };
  PreImageSubgroupKernel<Tgroup,TeltMatr,decltype(f_insert)>(ListMatrGens,
                                                             ListPermGens,
                                                             id_matr, eGRP,
                                                             f_insert);
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
  std::cerr << "GRP: |SetMatrGens|=" << SetMatrGens.size() << "\n";
#endif
  std::vector<TeltMatr> VectMatrGens;
  for (auto &eGen : SetMatrGens) {
    VectMatrGens.push_back(eGen);
  }
  return VectMatrGens;
}



template <typename Tgroup, typename TeltMatr>
std::vector<std::pair<TeltMatr, typename Tgroup::Telt>>
PreImageSubgroupTotal(std::vector<TeltMatr> const &ListMatrGens,
                      std::vector<typename Tgroup::Telt> const &ListPermGens,
                      TeltMatr const &id_matr, Tgroup const &eGRP) {
  using TeltPerm = typename Tgroup::Telt;
  using Telt = std::pair<TeltMatr, TeltPerm>;
  std::unordered_set<Telt> SetGens;
  auto f_insert=[&](Telt const& pair) -> void {
    SetGens.insert(pair);
  };
  PreImageSubgroupKernel<Tgroup,TeltMatr,decltype(f_insert)>(ListMatrGens,
                                                             ListPermGens,
                                                             id_matr, eGRP,
                                                             f_insert);
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
  std::cerr << "GRP: |SetGens|=" << SetGens.size() << "\n";
#endif
  std::vector<Telt> VectGens;
  for (auto &eGen : SetGens) {
    VectGens.push_back(eGen);
  }
  return VectGens;
}

template <typename Tgroup>
Tgroup ReadGroupFromStream(std::istream& is) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  size_t nbGen;
  int n_i;
  is >> n_i;
  is >> nbGen;
  Tidx n = Tidx(n_i);
  std::vector<Telt> LGen(nbGen);
  for (size_t iGen = 0; iGen < nbGen; iGen++) {
    std::vector<Tidx> ePermV(n);
    for (Tidx i = 0; i < n; i++) {
      int eVal_i;
      is >> eVal_i;
      Tidx eVal = Tidx(eVal_i);
      if (eVal >= n) {
        std::cerr << "GRP: Values is above range\n";
        std::cerr << "GRP: i=" << int(i) << " n=" << int(n)
                  << " eVal=" << int(eVal) << "\n";
        throw permutalib::PermutalibException{1};
      }
      ePermV[i] = eVal;
    }
    Telt ePerm(ePermV);
    LGen[iGen] = ePerm;
  }
  Telt id(n);
  return Tgroup(LGen, id);
}

template <typename Tgroup>
Tgroup ReadGroupFromFile(std::string const& eFile) {
  std::ifstream is(eFile);
  if (!is.good()) {
    std::cerr << "GRP: is stream is invalid, not possible to read eFile=" << eFile << "\n";
    throw permutalib::PermutalibException{1};
  }
  return permutalib::ReadGroupFromStream<Tgroup>(is);
}

Face ConvertStringToFace(std::string const& s) {
  size_t n = s.size();
  Face f(n);
  for (size_t i=0; i<n; i++) {
    std::string eChar = s.substr(i,1);
    if (eChar != "1" && eChar != "0") {
      std::cerr << "GRP: We have eChar=" << eChar << "\n";
      std::cerr << "GRP: Allowed values are 0 and 1\n";
      throw permutalib::PermutalibException{1};
    }
    if (eChar == "1") {
      f[i] = 1;
    }
  }
  return f;
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
  TeltPerm id_perm = stab.get_identity();
  return PreImageSubgroupAction<TeltPerm, TeltMatr, Tobj, decltype(f_op)>(
      ListMatrGens, ListPermGens, id_matr, id_perm, f, f_op);
}

template <typename TeltPerm, typename TeltMatr, typename Tint>
std::pair<std::vector<TeltMatr>,std::vector<TeltMatr>>
StabilizerRightCosetMatrixPermSubset(std::vector<TeltMatr> const &ListMatrGens,
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
  return PreImageSubgroupRightCosetAction<Tgroup, TeltMatr, Tobj, decltype(f_op)>(
      ListMatrGens, ListPermGens, id_matr, stab, f, f_op);
}

template <typename TeltPerm, typename TeltMatr, typename Tint>
std::optional<TeltMatr>
RepresentativeActionMatrixPermSubset(std::vector<TeltMatr> const &ListMatrGens,
                                     std::vector<TeltPerm> const &ListPermGens,
                                     TeltMatr const &id_matr, Face const &f1,
                                     Face const &f2) {
#ifdef DEBUG_REPRESENTATIVE_ACTION_MATRIX_PERM_SUBSET
  std::cerr << "GRP: Beginning of RepresentativeActionMatrixPermSubset\n";
#endif
  using Tidx = typename TeltPerm::Tidx;
  using Tgroup = Group<TeltPerm, Tint>;
  //
  Tidx len = f1.size();
  Tgroup GRP(ListPermGens, len);
#ifdef DEBUG_REPRESENTATIVE_ACTION_MATRIX_PERM_SUBSET
  std::cerr << "GRP: We have GRP\n";
  std::cerr << "GRP: GRP=" << GRP.GapString() << "\n";
#endif
  std::optional<TeltPerm> opt = GRP.RepresentativeAction_OnSets(f1, f2);
#ifdef DEBUG_REPRESENTATIVE_ACTION_MATRIX_PERM_SUBSET
  std::cerr << "GRP: We have opt\n";
#endif
  if (!opt)
    return {};
  // If we allow ourselves to compute RepresentativeAction
  Tgroup TheStab = GRP.Stabilizer_OnSets(f1);
  Tint OrbitSize = GRP.size() / TheStab.size();
  Tint CritSize = 10000;
#ifdef DEBUG_REPRESENTATIVE_ACTION_MATRIX_PERM_SUBSET
  std::cerr << "GRP: We have OrbitSize=" << OrbitSize << " CritSize=" << CritSize << "\n";
#endif
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
    std::cerr << "GRP: We have ListGensB\n";
#endif
    TgroupB GRP_B(ListGensB, idB);
#ifdef DEBUG_REPRESENTATIVE_ACTION_MATRIX_PERM_SUBSET
    std::cerr << "GRP: We have GRP_B\n";
#endif
    Telt res = GRP_B.Sift(ePair);
#ifdef DEBUG_REPRESENTATIVE_ACTION_MATRIX_PERM_SUBSET
    std::cerr << "GRP: We have res\n";
#endif
    //  NicePrint("res", res);
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
    std::vector<Tidx> const &V = res.getListVal();
    for (Tidx u = 0; u < len; u++) {
      if (V[u] != u) {
        std::cerr << "GRP: The permutation residue is not the identity at u=" << u
                  << "\n";
        throw PermutalibException{1};
      }
    }
#endif
    Tseq ret_seq = Inverse(res.getElt());
#ifdef DEBUG_REPRESENTATIVE_ACTION_MATRIX_PERM_SUBSET
    std::cerr << "GRP: We have ret_seq\n";
#endif
    TeltMatr ret_matr = id_matr;
#ifdef DEBUG_PRE_IMAGE_COMPLEXITY
    std::cerr << "GRP: RepresentativeActionMatrixPermSubset, comp(ret_seq)=" << ret_seq.complexity() << "\n";
#endif
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
    std::cerr << "GRP: We have ret_matr\n";
#endif
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
#ifdef DEBUG_REPRESENTATIVE_ACTION_MATRIX_PERM_SUBSET
      std::cerr << "GRP: Beginning of f_get_elt\n";
#endif
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
#ifdef DEBUG_REPRESENTATIVE_ACTION_MATRIX_PERM_SUBSET
      std::cerr << "GRP: End of f_get_elt\n";
#endif
      return ret_matr;
    };
    f_insert(f1, miss_val, miss_val);
    size_t n_done = 0;
    while(true) {
      size_t len = l_x_iorig_igen.size();
#ifdef DEBUG_REPRESENTATIVE_ACTION_MATRIX_PERM_SUBSET
      std::cerr << "GRP: n_done=" << n_done << " len=" << len << "\n";
#endif
      if (n_done == len)
        break;
      for (size_t iOrig=n_done; iOrig<len; iOrig++) {
        // Copy is needed since the std::vector is growing
        Face f = std::get<0>(l_x_iorig_igen[iOrig]);
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
    std::cerr << "GRP: We should never reach that stage\n";
    throw PermutalibException{1};
  };

  if (OrbitSize < CritSize) {
    return f_build_orbit();
  } else {
    return f_mapping_elt();
  }
}

template <typename TeltPerm, typename TeltMatr, typename Tint>
struct PreImagerElement {
private:
  using Tidx = typename TeltPerm::Tidx;
  using Tgroup = Group<TeltPerm, Tint>;
  using Tseq = SequenceType<true>;
  using Telt = PermutationElt<Tidx, Tseq>;
  using TgroupB = Group<Telt, Tint>;
  TeltMatr id_matr;
  std::vector<TeltMatr> ListMatrGens;
  std::vector<TeltMatr> ListMatrGens_inv;
  std::optional<TgroupB> opt;
public:
  PreImagerElement(std::vector<TeltMatr> const &_ListMatrGens,
                   std::vector<TeltPerm> const &ListPermGens,
                   TeltMatr const &_id_matr) : id_matr(_id_matr), ListMatrGens(_ListMatrGens) {
    size_t nGen = ListMatrGens.size();
    if (nGen > 0) {
      Tseq id_seq;
      Tidx len = ListPermGens[0].size();
      TeltPerm id_perm(len);
      Telt idB(id_perm.getListVal(), id_seq);
      std::vector<Telt> ListGensB;
      for (size_t iGen = 0; iGen < nGen; iGen++) {
        std::vector<int64_t> ListIdx{int64_t(iGen) + 1};
        Tseq e_seq(ListIdx);
        Telt fPair(ListPermGens[iGen].getListVal(), e_seq);
        ListGensB.push_back(fPair);
        ListMatrGens_inv.push_back(Inverse(ListMatrGens[iGen]));
      }
      TgroupB GRP_B(ListGensB, idB);
      opt = GRP_B;
    }
  }
  std::optional<TeltMatr> get_preimage(TeltPerm const& elt) const {
    if (opt) {
      TgroupB const& GRP_B = *opt;
      Tseq id_seq;
      Telt ePair(elt.getListVal(), id_seq);
      Telt res = GRP_B.Sift(ePair);
      if (!res.isIdentity()) {
        return {};
      }
      Tseq ret_seq = Inverse(res.getElt());
#ifdef DEBUG_PRE_IMAGE_COMPLEXITY
      std::cerr << "GRP: get_preimage, comp(ret_seq)=" << ret_seq.complexity() << "\n";
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
      return ret_matr;
    } else {
      if (elt.isIdentity()) {
        return id_matr;
      } else {
        return {};
      }
    }
  }


};





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

// clang-format off
}  // boost::serialization
// clang-format on

namespace permutalib {

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




namespace std {
  template <typename TeltMatr> struct hash<permutalib::SeqTracker<TeltMatr>> {
    std::size_t operator()(const permutalib::SeqTracker<TeltMatr> &v) const {
      return std::hash<TeltMatr>()(v.x);
    }
  };
// clang-format off
}  // namespace std
// clang-format on



#endif  // SRC_GAP_GROUP_H_
// clang-format on
