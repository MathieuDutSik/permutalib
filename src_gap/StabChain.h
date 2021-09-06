#ifndef DEFINE_PERMUTALIB_STAB_CHAIN_H
#define DEFINE_PERMUTALIB_STAB_CHAIN_H

/*
  Version of the code used. It must be from Summer 2016.
 */


/*
  Questions on the GAP code:
  ---Why QuickInverseRepresentative is not used more if it is the same but maybe quicker
     than InverseRepresentative

  Not needed for the time being:
  ---MembershipTestKnownBase needed for PCGS
  ---MinimalElementCosetStabChain for coset representatives.
  ---ListStabChain since we are using a different design in our code.
  ---OrbitStabChain which seems used by no code at all.
  ---knownBase seems strange. We should be able to work completely without it,
     since I never set up the knownBase in the first place.
  ---

  For what it seems:
  ---The transversal entry are permutation (when assigned) and all included in the
     S.labels
  ---The S.labels are all identical from the top to bottom.
  ---Significantly, this seems also the case of S.identity.
  ---The attribute "relativeOrders" is related to pcgs and not needed here.
*/

#include <vector>
#include <memory>

#include "GapPrint.h"
#include "factorize.h"
#include "PermGroup.h"
#include "list.h"
#include "COMB_Vectors.h"

// Needed for the comparison with
#ifdef SYNCHRONIZED_DEBUG_GAP478
# define DEBUG_STABCHAIN
# define DEBUG_CHANGE_STAB_CHAIN
#endif

// Other debugging.(but not currently in the gap stuff)
#undef DEBUG_ADD_GEN_SCH
#undef DEBUG_INV_REP



namespace permutalib {

std::string GetIntTypeNature(int const& val)
{
  if (val == int_reducedm1)
    return "reduced-1";
  if (val == int_false)
    return "false";
  if (val == int_true)
    return "true";
  if (val == int_fail)
    return "fail";
  if (val == int_int)
    return "int";
  if (val == int_perm)
    return "perm";
  if (val == int_group)
    return "group";
  if (val == int_stablev)
    return "stab+lev";
  return "Not allowed value";
}


template<typename Telt>
typename Telt::Tidx OnPoints(typename Telt::Tidx const& val, Telt const& g)
{
  return PowAct(val, g);
}


template<typename T>
void AssignationVectorGapStyle(std::vector<T> & eVect, size_t const& pos, T const& val)
{
  size_t siz=eVect.size();
  if (pos < siz) {
    eVect[pos] = val;
    return;
  }
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
  if (pos != siz) {
    std::cerr << "Assignation leaves gap in the vector. Not allowed\n";
    throw PermutalibException{1};
  }
#endif
  eVect.push_back(val);
}


template<typename T, typename Telt>
std::vector<T> PermutedAct(std::vector<T> const& V, Telt const& g)
{
  using Tidx=typename Telt::Tidx;
  Tidx len=Tidx(V.size());
  std::vector<T> Vret(len);
  for (Tidx i=0; i<len; i++) {
    Tidx iImg=g.at(i);
    Vret[iImg] = V[i];
  }
  return Vret;
}


template<typename Telt>
size_t GetLabelIndex(std::vector<Telt> & labels, Telt const& u)
{
  size_t nbLabel=labels.size();
  for (size_t iLabel=0; iLabel<nbLabel; iLabel++)
    if (labels[iLabel] == u)
      return iLabel;
  labels.push_back(u);
  return nbLabel;
}

template<typename Telt>
size_t GetLabelIndex_const(std::vector<Telt> const& labels, Telt const& u)
{
  size_t nbLabel=labels.size();
  for (size_t iLabel=0; iLabel<nbLabel; iLabel++)
    if (labels[iLabel] == u)
      return iLabel;
  return std::numeric_limits<size_t>::max();
}

// The labels are put on top since they are all identical.
// The stabilizers can be added in any way:
// ---At the top of the chain
// ---Removed in the middle
// Therefore we need a data set that allows us to do that and
// native std::vector are not adequate.
// We choose to use the same structure of shared_ptr as the GAP
// code.


template<typename Telt>
struct CommonStabInfo {
  typename Telt::Tidx n;
  Telt identity;
  std::vector<Telt> labels;
};

template<typename Telt, typename Tidx_label>
struct StabLevel {
  std::vector<Tidx_label> transversal;
  std::vector<typename Telt::Tidx> orbit;
  std::vector<Tidx_label> genlabels; // Not used in Random algorithm
  std::vector<int8_t> cycles; // We need a vector because if we use a "Face" dynamic bitset then
                              // we cannot extend. On the other hand if we use a std::vector<bool>
                              // we are exposed to the problems of this
  bool IsBoundCycle;
  // entry below are specific to the random algorithms:
  std::vector<Telt> treegen;
  std::vector<Telt> treegeninv;
  std::vector<Telt> aux;
  int treedepth;
  int diam;
  //
  std::shared_ptr<CommonStabInfo<Telt>> comm;
  std::shared_ptr<StabLevel<Telt,Tidx_label>> stabilizer;
};

template<typename Telt, typename Tidx_label>
using StabChain = std::shared_ptr<StabLevel<Telt,Tidx_label>>;
// other possible entries:
// transimages, genimages, labelimages, idimage


template<typename Telt, typename Tidx_label>
bool IsIdenticalObj(StabChain<Telt,Tidx_label> const& S1, StabChain<Telt,Tidx_label> const& S2)
{
  return S1 == S2;
}





template<typename Telt, typename Tidx_label>
StabChain<Telt,Tidx_label> StabChainGenerators(std::vector<Telt> const& generators, typename Telt::Tidx const& n, Telt const& id)
{
  std::shared_ptr<CommonStabInfo<Telt>> comm = std::make_shared<CommonStabInfo<Telt>>(CommonStabInfo<Telt>({n, id, generators}));
  //
  std::vector<Tidx_label> transversal(n, std::numeric_limits<Tidx_label>::max());
  std::vector<typename Telt::Tidx> orbit;
  std::vector<Tidx_label> genlabels;
  for (Tidx_label igen=0; igen<Tidx_label(generators.size()); igen++)
    genlabels.push_back(igen);
  std::vector<int8_t> cycles;
  bool IsBoundCycle = false;
  std::vector<Telt> treegen;
  std::vector<Telt> treegeninv;
  std::vector<Telt> aux;
  int treedepth = 0;
  int diam = 0;
  return std::make_shared<StabLevel<Telt,Tidx_label>>(StabLevel<Telt,Tidx_label>({std::move(transversal), std::move(orbit), std::move(genlabels), std::move(cycles), IsBoundCycle, std::move(treegen), std::move(treegeninv), std::move(aux), treedepth, diam, comm, nullptr}));
}


//
// Printing functionality
//


template<typename T>
void PrintVectDebug(std::string const& str, std::vector<T> const& V)
{
  std::cerr << str << " =";
  for (auto & eVal : V)
    std::cerr << " " << eVal;
  std::cerr << "\n";
}


template<typename Telt, typename Tidx_label>
void PrintStabChainTransversals(StabChain<Telt,Tidx_label> const& S)
{
  using Tidx=typename Telt::Tidx;
  StabChain<Telt,Tidx_label> Swork = S;
  Tidx n=Swork->comm->n;
  size_t iLevel=0;
  while (Swork != nullptr) {
    std::vector<std::optional<Telt>> V(n);
    for (Tidx i=0; i<n; i++) {
      Tidx_label eVal = Swork->transversal[i];
      if (eVal == std::numeric_limits<Tidx_label>::max())
        V[i] = {};
      else
        V[i] = Swork->comm->labels[eVal];
    }
    //
    std::cerr << "CPP i=" << iLevel << " transversal=" << GapStringMissingTVector(V) << "\n";
    Swork = Swork->stabilizer;
    iLevel++;
  }
}


template<typename Telt, typename Tidx_label>
void PrintStabChainOrbits(StabChain<Telt,Tidx_label> const& S)
{
  StabChain<Telt,Tidx_label> Swork = S;
  size_t iLevel=0;
  while (Swork != nullptr) {
    std::cerr << "CPP i=" << iLevel << " orbit=" << GapStringIntVector(Swork->orbit) << "\n";
    Swork = Swork->stabilizer;
    iLevel++;
  }
}


template<typename Telt, typename Tidx_label>
std::string GetListStabCommPartition(std::vector<StabChain<Telt,Tidx_label>> const& ListS)
{
  size_t len = ListS.size();
  Face Status(len);
  std::vector<std::string> ListStr;
  for (size_t i=0; i<len; i++) {
    if (Status[i] == 0) {
      std::vector<int> LVal;
      auto ptr = ListS[i]->comm;
      for (size_t j=0; j<len; j++) {
        if (ListS[j]->comm == ptr) {
          LVal.push_back(int(j));
          Status[j] = 1;
        }
      }
      std::string estr = GapStringIntVector(LVal);
      ListStr.emplace_back(std::move(estr));
    }
  }
  return GapStringTVector(ListStr);
}


template<typename Telt, typename Tidx_label>
void PrintStabChain(StabChain<Telt,Tidx_label> const& S)
{
  using Tidx=typename Telt::Tidx;
  StabChain<Telt,Tidx_label> Sptr = S;
  Tidx n = Sptr->comm->n;
  std::cerr << "CPP Reference Partition=" << GetListStabCommPartition(ListStabChain(S)) << "\n";
  size_t iLevel = 0;
  while (Sptr != nullptr) {
    std::cerr << "CPP iLev=" << iLevel << "\n";
    std::string strTransversal = "[ ]";
    if (Sptr->transversal.size() > 0) {
      if (Tidx(Sptr->transversal.size()) != n) {
        std::cerr << "Sptr->transversal should be of length 0 or n=" << n << "\n";
        throw PermutalibException{1};
      }
      std::vector<std::optional<Telt>> V(n);
      for (Tidx i=0; i<n; i++) {
        Tidx_label eVal = Sptr->transversal[i];
        if (eVal == std::numeric_limits<Tidx_label>::max())
          V[i] = {};
        else
          V[i] = Sptr->comm->labels[eVal];
      }
      strTransversal = GapStringMissingTVector(V);
    }
    //
    std::cerr << "CPP   orbit=" << GapStringIntVector(Sptr->orbit) << "\n";
    std::cerr << "CPP   transversal=" << strTransversal << "\n";
    if (Sptr->IsBoundCycle) {
      std::cerr << "CPP   cycles=" << GapStringBoolVectorB(Sptr->cycles) << "\n";
    } else {
      std::cerr << "CPP   No cycles\n";
    }
    Sptr = Sptr->stabilizer;
    iLevel++;
  }
}


template<typename Telt, typename Tidx_label>
void PrintListStabCommPartition(std::string const& mesg, std::vector<StabChain<Telt,Tidx_label>> const& ListS)
{
  std::cerr << mesg << " ListStabCommPartition=" << GetListStabCommPartition(ListS) << "\n";
}


template<typename Telt>
std::string perm_to_string(Telt const& eVal)
{
  std::ostringstream os;
  os << eVal;
  std::string str=os.str();
  return str;
}


template<typename Telt, typename Tidx_label>
std::string GetStringExpressionOfStabChain(StabChain<Telt,Tidx_label> const& eRec)
{
  std::string strRet="record_";
  StabChain<Telt,Tidx_label> eStab = eRec;
  while (eStab != nullptr) {
    strRet += "orbit_[";
    for (auto & eVal : eStab->orbit)
      strRet += " " + std::to_string(eVal+1);
    strRet += "]";
    strRet += "_transversal_";
    for (auto & eVal : eStab->transversal) {
      if (eVal == std::numeric_limits<Tidx_label>::max())
        strRet += " " + std::to_string(eVal);
      else
        strRet += " " + perm_to_string(eStab->comm->labels[eVal]);
    }
    eStab = eStab->stabilizer;
  }
  return strRet;
}


template<typename Telt, typename Tidx_label>
std::ostream& operator<<(std::ostream& os, StabChain<Telt,Tidx_label> const& Stot)
{
  os << "CPP n=" << Stot->comm->n << " identity=" << GapStyleString(Stot->comm->identity) << "\n";
  os << "CPP labels=[ ";
  bool IsFirst=true;
  for (auto & eLabel : Stot->comm->labels) {
    if (!IsFirst)
      os << ", ";
    IsFirst=false;
    os << GapStyleString(eLabel);
  }
  os << " ]\n";
  size_t iLev=0;
  StabChain<Telt,Tidx_label> Sptr = Stot;
  while (Sptr!= nullptr) {
    os << "CPP iLev=" << iLev << "\n";
    os << "CPP  orbit =";
    for (auto & eVal : Sptr->orbit) {
      os << " " << eVal+1;
    }
    os << "\n";
    os << "CPP  transversal =";
    for (auto & eVal : Sptr->transversal) {
      if (eVal == std::numeric_limits<Tidx_label>::max())
	os << " " << eVal;
      else
	os << " " << Sptr->comm->labels[eVal];
    }
    os << "\n";
    Sptr = Sptr->stabilizer;
    iLev++;
  }
  return os;
}


//
// Non printing stuff
//


template<typename Telt, typename Tidx_label>
void UnbindCycles(StabChain<Telt,Tidx_label> const& S)
{
  StabChain<Telt,Tidx_label> Sptr = S;
  while (Sptr != nullptr) {
    Sptr->cycles.clear();
    Sptr->IsBoundCycle = false;
    Sptr = Sptr->stabilizer;
  }
}


template<typename Telt, typename Tidx_label>
int GetStabilizerDepth(StabChain<Telt,Tidx_label> const& S1)
{
  StabChain<Telt,Tidx_label> S2 = S1;
  int dep = 0;
  while (S2 != nullptr) {
    dep++;
    S2 = S2->stabilizer;
  }
  return dep;
}


template<typename Telt, typename Tidx_label>
bool HasStabStab(StabChain<Telt,Tidx_label> const& S)
{
  if (S == nullptr)
    return false;
  if (S->stabilizer == nullptr)
    return false;
  return true;
}


template<typename Telt, typename Tidx_label>
StabChain<Telt,Tidx_label> ShallowCopy(StabChain<Telt,Tidx_label> const& S)
{
  StabLevel<Telt,Tidx_label> Sret = *S;
  Sret.comm = S->comm;
  Sret.stabilizer = S->stabilizer;
  return std::make_shared<StabLevel<Telt,Tidx_label>>(Sret);
}



template<typename Telt, typename Tidx_label>
StabChain<Telt,Tidx_label> StructuralCopy(StabChain<Telt,Tidx_label> const& S)
{
  if (S == nullptr)
    return nullptr;
  using Tcomm=std::shared_ptr<CommonStabInfo<Telt>>;
  std::vector<Tcomm> ListLabels;
  std::vector<Tcomm> ListLabelsImg;
  auto get_comm=[&](Tcomm const& e_comm) -> Tcomm {
     for (size_t i=0; i<ListLabels.size(); i++)
       if (ListLabels[i] == e_comm)
         return ListLabelsImg[i];
     Tcomm comm_new = std::make_shared<CommonStabInfo<Telt>>(*e_comm);
     ListLabels.push_back(e_comm);
     ListLabelsImg.push_back(comm_new);
     return comm_new;
  };
  StabChain<Telt,Tidx_label> Sptr = S;
  std::vector<std::shared_ptr<StabLevel<Telt,Tidx_label>>> ListPtr;
  while (Sptr != nullptr) {
    ListPtr.push_back(Sptr);
    Sptr = Sptr->stabilizer;
  }
  size_t len = ListPtr.size();
  StabChain<Telt,Tidx_label> Sret = nullptr;
  for (size_t i=0; i<len; i++) {
    size_t j = len - 1 - i;
    std::shared_ptr<StabLevel<Telt,Tidx_label>> S2 = std::make_shared<StabLevel<Telt,Tidx_label>>(*ListPtr[j]);
    S2->stabilizer = Sret;
    S2->comm = get_comm(ListPtr[j]->comm);
    Sret = S2;
  }
#ifdef DEBUG
  std::string str1 = GetStringExpressionOfStabChain(S);
  std::string str2 = GetStringExpressionOfStabChain(Sret);
  if (str1 != str2) {
    std::cerr << "We fail to have equality of entries\n";
    std::cerr << "str1=" << str1 << "\n";
    std::cerr << "str2=" << str2 << "\n";
    throw PermutalibException{1};
  }
#endif
  return Sret;
}


template<typename Telt, typename Tidx_label>
StabChain<Telt,Tidx_label> CopyStabChain(StabChain<Telt,Tidx_label> const& S)
{
  return StructuralCopy(S);
}


template<typename Telt, typename Tidx_label>
StabLevel<Telt,Tidx_label> EmptyStabLevel(std::shared_ptr<CommonStabInfo<Telt>> const& comm)
{
  std::vector<Tidx_label> transversal(comm->n, std::numeric_limits<Tidx_label>::max());
  std::vector<typename Telt::Tidx> orbit;
  std::vector<Tidx_label> genlabels;
  std::vector<int8_t> cycles;
  bool IsBoundCycle = false;
  std::vector<Telt> treegen;
  std::vector<Telt> treegeninv;
  std::vector<Telt> aux;
  int treedepth = 0;
  int diam = 0;
  return {std::move(transversal), std::move(orbit), std::move(genlabels), std::move(cycles), IsBoundCycle, std::move(treegen), std::move(treegeninv), std::move(aux), treedepth, diam, comm, nullptr};
}



template<typename Telt, typename Tidx_label>
StabChain<Telt,Tidx_label> EmptyStabChain(typename Telt::Tidx const& n)
{
  Telt id(n);
  std::shared_ptr<CommonStabInfo<Telt>> comm = std::make_shared<CommonStabInfo<Telt>>(CommonStabInfo<Telt>({n, id, {id}}));
  return std::make_shared<StabLevel<Telt,Tidx_label>>(EmptyStabLevel<Telt,Tidx_label>(comm));
}


template<typename Telt, typename Tidx_label>
StabChain<Telt,Tidx_label> EmptyStabChainPlusNode(typename Telt::Tidx const& n, typename Telt::Tidx const& bas)
{
  StabChain<Telt,Tidx_label> S = EmptyStabChain<Telt,Tidx_label>(n);
  InitializeSchreierTree(S, bas);
  return S;
}


template<typename Telt, typename Tidx_label>
StabChain<Telt,Tidx_label> EmptyStabChainPlusCommon(std::shared_ptr<CommonStabInfo<Telt>> const& comm)
{
  StabChain<Telt,Tidx_label> S = std::make_shared<StabLevel<Telt,Tidx_label>>(EmptyStabLevel<Telt,Tidx_label>(comm));
  return S;
}


template<typename Telt, typename Tidx_label>
StabChain<Telt,Tidx_label> EmptyStabChainPlusCommonPlusNode(std::shared_ptr<CommonStabInfo<Telt>> const& comm, typename Telt::Tidx const& bas)
{
  StabChain<Telt,Tidx_label> S = std::make_shared<StabLevel<Telt,Tidx_label>>(EmptyStabLevel<Telt,Tidx_label>(comm));
  InitializeSchreierTree(S, bas);
  return S;
}



template<typename Telt, typename Tidx_label>
typename Telt::Tidx BasePoint(StabChain<Telt,Tidx_label> const& S)
{
  using Tidx=typename Telt::Tidx;
  if (S == nullptr)
    return std::numeric_limits<Tidx>::max();
  if (S->orbit.size() == 0)
    return std::numeric_limits<Tidx>::max();
  return S->orbit[0];
}





template<typename Telt, typename Tidx_label>
void RemoveStabChain(StabChain<Telt,Tidx_label> & Stot)
{
  using Tidx=typename Telt::Tidx;
  Stot->stabilizer = nullptr;
  Stot->genlabels.clear();
  Stot->orbit.clear();
  Stot->transversal.clear();
  Stot->cycles.clear();
  Stot->IsBoundCycle = false;
  Stot->treegen.clear();
  Stot->treegeninv.clear();
  Stot->aux.clear();
  //
  Telt id = Stot->comm->identity;
  Tidx n = Stot->comm->n;
  Stot->comm = std::make_shared<CommonStabInfo<Telt>>(CommonStabInfo<Telt>({n, id, {id}}));
}


template<typename Telt, typename Tidx_label>
bool IsInBasicOrbit(StabChain<Telt,Tidx_label> const& S, typename Telt::Tidx const& pnt)
{
  Tidx_label eVal=S->transversal[pnt];
  if (eVal == std::numeric_limits<Tidx_label>::max())
    return false;
  return true;
}


template<typename Telt, typename Tidx_label>
Telt InverseRepresentative(StabChain<Telt,Tidx_label> const& S, typename Telt::Tidx const& pnt)
{
  using Tidx=typename Telt::Tidx;
  Tidx bpt=S->orbit[0];
  Telt rep=S->comm->identity;
  Tidx pntw=pnt;
#ifdef DEBUG_INV_REP
  std::cerr << "CPP INVREP bpt=" << int(bpt+1) << " pnt=" << int(pntw+1) << "\n";
#endif
  while(pntw != bpt) {
    Tidx_label idx=S->transversal[pntw];
    const Telt& te=S->comm->labels[idx];
#ifdef DEBUG_INV_REP
    std::cerr << "CPP INVREP te=" << te << "\n";
#endif
    pntw = PowAct(pntw, te);
#ifdef DEBUG_INV_REP
    std::cerr << "CPP INVREP   pnt=" << int(pntw+1) << "\n";
#endif
    rep *= te;
#ifdef DEBUG_INV_REP
    std::cerr << "CPP INVREP   rep=" << rep << "\n";
#endif
  }
#ifdef DEBUG_INV_REP
  std::cerr << "CPP INVREP return rep=" << rep << "\n";
#endif
  return rep;
}


template<typename Telt, typename Tidx_label>
std::vector<Telt> InverseRepresentativeWord(StabChain<Telt,Tidx_label> const& S, typename Telt::Tidx const& pnt)
{
  using Tidx=typename Telt::Tidx;
  Tidx bpt=S->orbit[0];
  std::vector<Telt> word;
  Tidx pntw=pnt;
  while(pntw != bpt) {
    Tidx_label idx=S->transversal[pntw];
    const Telt& te=S->comm->labels[idx];
    pntw = PowAct(pntw, te);
    word.push_back(te);
  }
  return word;
}


template<typename Telt, typename Tidx_label>
Telt SiftedPermutation(StabChain<Telt,Tidx_label> const& S, Telt const& g)
{
  using Tidx=typename Telt::Tidx;
  Telt gW = g;
  StabChain<Telt,Tidx_label> Sptr = S;
  while(true) {
    if (Sptr->stabilizer == nullptr || gW.isIdentity())
      return gW;
    Tidx bpt = Sptr->orbit[0];
    Tidx img = PowAct(bpt, gW);
    if (Sptr->transversal[img] == std::numeric_limits<Tidx_label>::max())
      return gW;
    while(img != bpt) {
      Tidx_label idx = Sptr->transversal[img];
      gW *= Sptr->comm->labels[idx];
      img = PowAct(bpt, gW);
    }
    Sptr = Sptr->stabilizer;
  }
  return gW;
}


template<typename Telt, typename Tidx_label>
bool IsElementInStabChain(StabChain<Telt,Tidx_label> const& S, Telt const& g)
{
  Telt res = SiftedPermutation(S, g);
  return res.isIdentity();
}



// Testing that S1 is a subset of S2
template<typename Telt, typename Tidx_label>
bool InclusionTest(const StabChain<Telt,Tidx_label>& S1, const StabChain<Telt,Tidx_label>& S2)
{
  for (auto & eGen : Kernel_GeneratorsOfGroup(S1)) {
    Telt res = SiftedPermutation(S2, eGen);
    if (!res.isIdentity())
      return false;
  }
  return true;
}


template<typename Telt, typename Tidx_label>
bool EqualityTest(const StabChain<Telt,Tidx_label>& S1, const StabChain<Telt,Tidx_label>& S2)
{
  if (!InclusionTest(S1,S2))
    return false;
  if (!InclusionTest(S2,S1))
    return false;
  return true;
}




template<typename Telt, typename Tidx_label>
std::vector<typename Telt::Tidx> BaseStabChain(StabChain<Telt,Tidx_label> const& S)
{
  using Tidx=typename Telt::Tidx;
  StabChain<Telt,Tidx_label> Sptr = S;
  std::vector<Tidx> base;
  while (true) {
    if (Sptr == nullptr || Sptr->orbit.size() == 0)
      break;
    base.push_back(Sptr->orbit[0]);
    Sptr = Sptr->stabilizer;
  }
  return base;
}



template<typename Telt, typename Tidx_label, typename Tint>
Tint SizeStabChain(StabChain<Telt,Tidx_label> const& S)
{
  Tint size=1;
  StabChain<Telt,Tidx_label> Sptr = S;
  while (Sptr != nullptr) {
    Tint siz_i = Tint(Sptr->orbit.size());
    if (siz_i == 0)
      break;
    size *= siz_i;
    Sptr = Sptr->stabilizer;
  }
  return size;
}




template<typename Telt, typename Tidx_label>
std::map<typename Telt::Tidx, int> FactorsSizeStabChain(StabChain<Telt,Tidx_label> const& S)
{
  using Tidx = typename Telt::Tidx;
  std::map<Tidx, int> MapPrimeMult;
  StabChain<Telt,Tidx_label> Sptr = S;
  while (Sptr != nullptr) {
    Tidx siz = Tidx(Sptr->orbit.size());
    if (siz == 0)
      break;
    std::vector<Tidx> V = factorize(siz);
    for (auto & eVal : V)
      MapPrimeMult[eVal]++;
    Sptr = Sptr->stabilizer;
  }
  return MapPrimeMult;
}




template<typename Telt, typename Tidx_label>
std::vector<Telt> StrongGeneratorsStabChain(StabChain<Telt,Tidx_label> const& S)
{
  StabChain<Telt,Tidx_label> Sptr = S;
  std::unordered_set<Telt> sgs_set;
  while (Sptr != nullptr) {
    const std::vector<Telt>& labels = Sptr->comm->labels;
    for (auto & pos : Sptr->genlabels)
      sgs_set.insert(labels[pos]);
    Sptr = Sptr->stabilizer;
  }
  std::vector<Telt> sgs;
  sgs.insert(sgs.end(), sgs_set.begin(), sgs_set.end());
  return sgs;
}



template<typename Telt, typename Tidx_label>
std::vector<Telt> GeneratorsStabChain(StabChain<Telt,Tidx_label> const& S)
{
  std::vector<Telt> sgs;
  sgs.reserve(S->genlabels.size());
  for (auto & ePos : S->genlabels) {
#ifdef DEBUG_STABCHAIN
    std::cerr << "DEBUG ePos=" << ePos << "\n";
#endif
    sgs.push_back(S->comm->labels[ePos]);
  }
  return sgs;
}



template<typename Telt, typename Tidx_label>
std::vector<Telt> Kernel_GeneratorsOfGroup(StabChain<Telt,Tidx_label> const& S)
{
  return StrongGeneratorsStabChain(S);
}



template<typename Telt, typename Tidx_label>
Telt LargestElementStabChain(StabChain<Telt,Tidx_label> const& S)
{
  using Tidx=typename Telt::Tidx;
  StabChain<Telt,Tidx_label> Sptr = S;
  Telt rep=Sptr->identity;
  while (Sptr != nullptr) {
    if (Sptr->genlabels.size() == 0)
      break;
    Tidx pnt=Sptr->orbit[0];
    Tidx min=0;
    Tidx val=0;
    for (auto & i : Sptr->orbit) {
      Tidx img=PowAct(i, rep);
      if (img > val) {
	min=i;
	val=img;
      }
    }
    while (pnt != min) {
      Tidx_label idx=Sptr->transversal[min];
      const Telt& gen=Sptr->comm->labels[idx];
      rep = LeftQuotient(gen, rep);
      min = PowAct(min, gen);
    }
    Sptr = Sptr->stabilizer;
  }
  return rep;
}



// is base is empty then this just replaces the IsBound(options.base)
template<typename Tint, typename Tidx>
struct StabChainOptions {
  Tidx n;
  std::vector<Tidx> base;
  std::vector<Tidx> knownBase;
  int random;
  bool reduced;
  Tint size;
  Tint limit;
};


template<typename Tint, typename Tidx>
StabChainOptions<Tint,Tidx> GetStandardOptions(Tidx const& n)
{
  std::vector<Tidx> base;
  std::vector<Tidx> knownBase;
  int random = 1000;
  bool reduced=true;
  Tint size=0;
  Tint limit=0;
  return {n, std::move(base), std::move(knownBase), random, reduced, size, limit};
}


template<typename Telt>
bool IsTrivial_ListGen(std::vector<Telt> const& LGen)
{
  for (auto & eElt : LGen)
    if (!eElt.isIdentity())
      return false;
  return true;
}

template<typename Telt, typename Tidx_label>
std::vector<typename Telt::Tidx> MovedPoints(StabChain<Telt,Tidx_label> const& S)
{
  using Tidx = typename Telt::Tidx;
  std::unordered_set<Telt> LGen;
  StabChain<Telt,Tidx_label> Sptr = S;
  while (Sptr != nullptr) {
    for (auto & eIdx : Sptr->genlabels) {
      Telt eGen = S->comm->labels[eIdx];
      LGen.insert(eGen);
    }
    Sptr = Sptr->stabilizer;
  }
  auto IsMoved=[&](Tidx const& ePt) -> bool {
    for (auto & eGen : LGen)
      if (eGen.at(ePt) != ePt)
	return true;
    return false;
  };
  Tidx n=S->comm->n;
  std::vector<Tidx> LMoved;
  LMoved.reserve(n);
  for (Tidx i=0; i<n; i++)
    if (IsMoved(i))
      LMoved.push_back(i);
  return LMoved;
}


template<typename Telt, typename Tidx_label>
bool IsTrivial(StabChain<Telt,Tidx_label> const& G)
{
  StabChain<Telt,Tidx_label> Sptr = G;
  while (Sptr != nullptr) {
    for (auto & eIdx : Sptr->genlabels)
      if (!Sptr->comm->labels[eIdx].isIdentity())
        return false;
    Sptr = Sptr->stabilizer;
  }
  return true;
}






template<typename Telt>
typename Telt::Tidx LargestMovedPoint(std::vector<Telt> const& LGen)
{
  using Tidx = typename Telt::Tidx;
  if (LGen.size() == 0)
    return 0;
  Tidx n=LGen[0].size();
  Tidx eMov=0;
  for (auto & eGen : LGen)
    for (Tidx u=1; u<n; u++) {
      Tidx v=eGen.at(u);
      if (u != v) {
        if (u > eMov)
          eMov = u;
      }
    }
  eMov++;
  return eMov;
}


template<typename Telt, typename Tidx_label>
void InitializeSchreierTree(StabChain<Telt,Tidx_label> & S, typename Telt::Tidx const& pnt)
{
  using Tidx=typename Telt::Tidx;
  Tidx n=S->comm->n;
  //
  S->orbit = {pnt};
  //
  std::vector<Tidx_label> transversal(n, std::numeric_limits<Tidx_label>::max());
  transversal[pnt] = 0;
  S->transversal = transversal;
}


template<typename Telt, typename Tidx_label>
void InsertTrivialStabilizer(StabChain<Telt,Tidx_label> & Stot, typename Telt::Tidx const& pnt)
{
#ifdef DEBUG_STABCHAIN
  std::cerr << "CPP begin InsertTrivialStabilizer\n";
#endif
  Stot->stabilizer = ShallowCopy(Stot);
  Stot->transversal = {};
  Stot->orbit = {};
  Stot->genlabels = Stot->stabilizer->genlabels;
  Stot->cycles = {};
  Stot->IsBoundCycle = false;
  Stot->treegen = {};
  Stot->treegen = {};
  Stot->aux = {};
  Stot->treedepth = -1;
  Stot->diam = -1;
  InitializeSchreierTree(Stot, pnt);
}


template<typename Telt, typename Tidx_label, typename Tint>
StabChain<Telt,Tidx_label> StabChainOp_trivial_group(StabChain<Telt,Tidx_label> const& Stot, StabChainOptions<Tint, typename Telt::Tidx> const& options)
{
  using Tidx = typename Telt::Tidx;
#ifdef DEBUG
  std::cerr << "CPP Call to StabChainOp (trivial group)\n";
#endif
  StabChain<Telt,Tidx_label> S = EmptyStabChain<Telt,Tidx_label>(Stot->comm->n);
  if (options.base.size() > 0 && !options.reduced) {
    StabChain<Telt,Tidx_label> T = S;
    for (Tidx const& pnt : options.base) {
      InsertTrivialStabilizer( T, pnt );
      T = T->stabilizer;
    }
  }
  return S;
}



template<typename Telt, typename Tidx_label>
std::vector<StabChain<Telt,Tidx_label>> ListStabChain(StabChain<Telt,Tidx_label> const& S)
{
  std::vector<StabChain<Telt,Tidx_label>> ListStab;
  StabChain<Telt,Tidx_label> Sptr = S;
  while (Sptr != nullptr) {
    ListStab.push_back(Sptr);
    Sptr = Sptr->stabilizer;
  }
  return ListStab;
}



template<typename Telt, typename Tidx_label>
StabChain<Telt,Tidx_label> StabChainBaseStrongGenerators(std::vector<typename Telt::Tidx> const& base, std::vector<Telt> const& sgs)
{
  using Tidx=typename Telt::Tidx;
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
  if (sgs.size() == 0) {
    std::cerr << "sgs is empty. StabChainBaseStrongGenerators is broken in that case\n";
    throw PermutalibException{1};
  }
#endif
  Tidx n=sgs[0].size();
  StabChain<Telt,Tidx_label> S = EmptyStabChain<Telt,Tidx_label>(n);
  size_t nbGen = sgs.size();
  Face status(nbGen);
  size_t basSiz=base.size();
  for (size_t iBas=0; iBas<basSiz; iBas++) {
    Tidx pnt=base[iBas];
    std::vector<Telt> sgsFilt;
    for (size_t i=0; i<nbGen; i++)
      if (status[i] == 0)
	sgsFilt.push_back(sgs[i]);
    InsertTrivialStabilizer(S, pnt);
#ifdef DEBUG
    std::cerr << "CPP Before call to AddGeneratorsExtendSchreierTree from StabChainStrongGenerators\n";
#endif
    AddGeneratorsExtendSchreierTree(S, sgsFilt);
    for (size_t iGen=0; iGen<nbGen; iGen++)
      if (status[iGen] == 0 && PowAct(pnt, sgs[iGen]) != pnt)
	status[iGen] = 1;
  }
  return S;
}





template<typename Telt, typename Tidx_label>
void AddGeneratorsExtendSchreierTree(StabChain<Telt,Tidx_label> & S, std::vector<Telt> const& newgens)
{
  using Tidx=typename Telt::Tidx;
#ifdef DEBUG_ADD_GEN_SCH
  std::cerr << "CPP AGEST : Beginning of AddGeneratorsExtendSchreierTree\n";
  //  std::cerr << "CPP AGEST 1: genlabels=" << GapStringIntVector(S->genlabels) << "\n";
  StabChain<Telt,Tidx_label> Swrite = S;
  int idxwrt=0;
  while (Swrite != nullptr) {
    //    std::cerr << "DEBUG idxwrt=" << idxwrt << " |labels|=" << S->comm->labels.size() << "\n";
    idxwrt++;
    Swrite = Swrite->stabilizer;
  }
#endif
  size_t nbLabel = S->comm->labels.size();
  Face old=BlistList_direct(nbLabel, S->genlabels);
  old[0]=true;
  Face ald=old;
#ifdef DEBUG_ADD_GEN_SCH
  std::cerr << "CPP AGEST newgens=" << GapStringTVector(newgens) << "\n";
  std::cerr << "DEBUG AGEST 1: old=" << GapStringBoolVector(old) << "\n";
  std::cerr << "DEBUG AGEST 1: ald=" << GapStringBoolVector(ald) << "\n";
  std::cerr << "DEBUG AGEST labels=" << GapStringTVector(S->comm->labels) << "\n";
  std::cerr << "DEBUG AGEST 2: genlabels=" << GapStringIntVector(S->genlabels) << "\n";
#endif
  for (auto & gen : newgens) {
    Tidx_label pos = PositionVect_ui<Telt,Tidx_label>(S->comm->labels, gen);
    if (pos == std::numeric_limits<Tidx_label>::max()) {
      S->comm->labels.push_back(gen);
      old.push_back(false);
      ald.push_back(true);
      Tidx_label posG=Tidx_label(S->comm->labels.size()) - 1;
#ifdef DEBUG_ADD_GEN_SCH
      std::cerr << "CPP AGEST  genlabels insert 1:\n";
#endif
      S->genlabels.push_back(posG);
    } else {
      if (!ald[pos]) {
#ifdef DEBUG_ADD_GEN_SCH
        std::cerr << "CPP AGEST  genlabels insert 2:\n";
#endif
	S->genlabels.push_back(pos);
      }
    }
  }
#ifdef DEBUG_ADD_GEN_SCH
  std::cerr << "DEBUG AGEST 2: old=" << GapStringBoolVector(old) << "\n";
  std::cerr << "DEBUG AGEST 2: ald=" << GapStringBoolVector(ald) << "\n";
#endif

  size_t len = S->orbit.size();
#ifdef DEBUG_ADD_GEN_SCH
  std::cerr << "CPP AGEST len=" << len << "\n";
  std::cerr << "CPP AGEST S->orbit=" << GapStringIntVector(S->orbit) << "\n";
#endif
  size_t i=0;
#ifdef DEBUG_ADD_GEN_SCH
  //  std::cerr << "XXX ELIMINATE begin\n";
#endif
  if (S->IsBoundCycle) {
#ifdef DEBUG_ADD_GEN_SCH
    std::cerr << "CPP AGEST Cycles S.cycles=" << GapStringBoolVectorB(S->cycles) << "\n";
#endif
    while (i < S->orbit.size()) {
#ifdef DEBUG_ADD_GEN_SCH
      std::cerr << "CPP   AGEST i=" << int(i+1) << " |cycles|=" << S->cycles.size() << "\n";
#endif
      for (Tidx_label& j : S->genlabels) {
	if (i >= len || old[j] == 0) {
	  Tidx img=SlashAct(S->orbit[i], S->comm->labels[j]);
#ifdef DEBUG_ADD_GEN_SCH
	  std::cerr << "CPP     AGEST img=" << int(img+1) << " g=" << S->comm->labels[j] << "\n";
          //	  std::cerr << "DEBUG |S->transversal|=" << S->transversal.size() << " img=" << img << "\n";
#endif
	  if (S->transversal[img] != std::numeric_limits<Tidx_label>::max()) {
#ifdef DEBUG_ADD_GEN_SCH
	    std::cerr << "CPP       |S->cycles|=" << S->cycles.size() << " i=" << int(i+1) << "\n";
#endif
            AssignationVectorGapStyle(S->cycles, i, int8_t(true));
            //            S->cycles[i] = true;
#ifdef DEBUG_ADD_GEN_SCH
	    std::cerr << "CPP       AGEST assign true\n";
#endif
	  } else {
	    S->transversal[img] = j;
#ifdef DEBUG_ADD_GEN_SCH
	    std::cerr << "CPP       AGEST S.transversal[img]=" << S->comm->labels[j] << "\n";
#endif
	    S->orbit.push_back(img);
#ifdef DEBUG_ADD_GEN_SCH
	    std::cerr << "CPP       AGEST insert img\n";
#endif
	    S->cycles.push_back(false);
	  }
	}
      }
      i++;
    }
  } else {
#ifdef DEBUG_ADD_GEN_SCH
    std::cerr << "CPP AGEST No Cycles\n";
#endif
    while (i < S->orbit.size()) {
      for (Tidx_label& j : S->genlabels) {
	if (i >= len || old[j] == 0) {
	  Tidx img=SlashAct(S->orbit[i], S->comm->labels[j]);
	  if (S->transversal[img] == std::numeric_limits<Tidx_label>::max()) {
	    S->transversal[img] = j;
	    S->orbit.push_back(img);
	  }
	}
      }
      i++;
    }
  }
#ifdef DEBUG_ADD_GEN_SCH
  //  std::cerr << "XXX ELIMINATE end\n";
#endif
}



template<typename Telt, typename Tidx_label>
void ChooseNextBasePoint(StabChain<Telt,Tidx_label> & S, std::vector<typename Telt::Tidx> const& base, std::vector<Telt> const& newgens)
{
  using Tidx=typename Telt::Tidx;
#ifdef DEBUG_STABCHAIN
  std::cerr << "CPP base = " << GapStringIntVector(base) << "\n";
#endif
  auto IsFullyStable=[&](Tidx const& eBas) -> bool {
    for (auto & eGen : newgens) {
      if (PowAct(eBas, eGen) != eBas)
	return false;
    }
    return true;
  };
  size_t i = 0;
  size_t len=base.size();
  while(true) {
    if (i == len)
      break;
    Tidx eBas=base[i];
    if (!IsFullyStable(eBas))
      break;
    i++;
  }
  Tidx pnt;
  if (i < len) {
    pnt = base[i];
  } else {
    pnt = SmallestMovedPoint(newgens);
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
    if (pnt == std::numeric_limits<Tidx>::max()) {
      std::cerr << "The SmallestMovePoint return the maximum value\n";
      std::cerr << "It likely means that there is no smallest moved points.\n";
      std::cerr << "Please debug\n";
      throw PermutalibException{1};
    }
#endif
  }
  Tidx bpt;
  size_t pos;
  size_t miss_val = std::numeric_limits<size_t>::max();
  if (S->orbit.size() > 0) {
    bpt = S->orbit[0];
    pos = PositionVect_ui<Tidx,size_t>(base, bpt);
  } else {
    bpt = std::numeric_limits<Tidx>::max();
    pos = miss_val;
  }
#ifdef DEBUG_STABCHAIN
  std::cerr << "CPP pnt=" << PosFail_to_string(pnt) << " bpt=" << ConstrainedIntInfinity_to_string(bpt, S->comm->n) << " pos=" << PosFail_to_string(pos) << "\n";
#endif
  if ((pos != miss_val && i < pos) || (pos == miss_val && i < base.size()) || (pos == miss_val && pnt < bpt)) {
#ifdef DEBUG_STABCHAIN
    std::cerr << "CPP matching test\n";
#endif
    InsertTrivialStabilizer(S, pnt);
    if (S->stabilizer->IsBoundCycle) {
      std::vector<int8_t> eFace = {0};
      S->cycles = eFace;
      S->IsBoundCycle = true;
#ifdef DEBUG_STABCHAIN
      std::cerr << "CPP   Initializing cycles\n";
#endif
    }
  }
#ifdef DEBUG_STABCHAIN
  std::cerr << "CPP Exiting ChooseNextBasePoint\n";
#endif
}



template<typename Telt, typename Tidx_label, typename Tint>
void StabChainStrong(StabChain<Telt,Tidx_label> & S, std::vector<Telt> const& newgens, StabChainOptions<Tint, typename Telt::Tidx> const& options)
{
  using Tidx = typename Telt::Tidx;
#ifdef DEBUG_STABCHAIN
  std::cerr << "CPP Begin StabChainStrong : newgens=" << GapStringTVector(newgens) << "\n";
  std::cerr << "DEBUG |S->genlabels|=" << S->genlabels.size() << "\n";
  std::cerr << "CPP StabChainStrong 1: genlabels=" << GapStringIntVector(S->genlabels) << "\n";
#endif
  ChooseNextBasePoint(S, options.base, newgens);
#ifdef DEBUG_STABCHAIN
  std::cerr << "CPP StabChainStrong 2: genlabels=" << GapStringIntVector(S->genlabels) << "\n";
#endif

  Tidx pnt = S->orbit[0];
  Tidx len = Tidx(S->orbit.size());
  Tidx old = Tidx(S->genlabels.size());
#ifdef DEBUG_STABCHAIN
  std::cerr << "CPP Before AddGeneratorsExtendSchreierTree\n";
#endif
  AddGeneratorsExtendSchreierTree(S, newgens);

  //# If a new generator fixes the base point, put it into the stabilizer.
#ifdef DEBUG_STABCHAIN
  std::cerr << "CPP newgens=" << GapStringTVector(newgens) << "\n";
#endif
  for (auto & eGen : newgens) {
#ifdef DEBUG_STABCHAIN
    std::cerr << "CPP eGen=" << eGen << " eGen=" << GapStyleString(eGen) << "\n";
#endif
    if (!eGen.isIdentity() && PowAct(pnt, eGen) == pnt) {
#ifdef DEBUG_STABCHAIN
      std::cerr << "CPP   1: Calling StabChainStrong with eGen=" << GapStyleString(eGen) << "\n";
#endif
      StabChainStrong(S->stabilizer, {eGen}, options);
    }
  }

  // # Compute the Schreier generators (seems to work better backwards).
  std::vector<Tidx> pnts = ClosedInterval<Tidx>(0, Tidx(S->orbit.size()));
  if (S->IsBoundCycle)
    pnts = ListBlist(pnts, S->cycles);
#ifdef DEBUG_STABCHAIN
  std::cerr << "CPP pnts = " << GapStringIntVector(pnts) << "\n";
  std::cerr << "CPP Usecycle=" << S->IsBoundCycle << "\n";
  if (S->IsBoundCycle) {
    std::cerr << "CPP cycles=" << GapStringBoolVectorB(S->cycles) << "\n";
  }
#endif
  Tidx gen1=0;
#ifdef DEBUG_STABCHAIN
  std::cerr << "CPP StabChainStrong O=" << GapStringIntVector(S->orbit) << "\n";
#endif
  for (const Tidx & i : Reversed(pnts)) {
    Tidx p=S->orbit[i];
#ifdef DEBUG_STABCHAIN
    std::cerr << "CPP StabChainStrong i=" << int(i+1) << " p=" << int(p+1) << "\n";
#endif
    Telt rep=InverseRepresentative(S, p);
    if (i < len)
      gen1 = old;
#ifdef DEBUG_STABCHAIN
    std::cerr << "CPP StabChainStrong gen1=" << int(gen1+1) << " rep=" << rep << "\n";
#endif
    Tidx end_seq = Tidx(S->genlabels.size());
    for (Tidx j=gen1; j<end_seq; j++) {
      const Telt& g = S->comm->labels[ S->genlabels[j] ];
#ifdef DEBUG_STABCHAIN
      std::cerr << "CPP StabChainStrong   j=" << int(j+1) << " g=" << g << "\n";
#endif
      if (S->transversal[ SlashAct(p, g) ] != S->genlabels[j]) {
        Telt sch = SiftedPermutation(S, Inverse(g*rep));
#ifdef DEBUG_STABCHAIN
	std::cerr << "CPP sch=" << sch << " g=" << g << " rep=" << rep << "\n";
#endif
	if (!sch.isIdentity())
	  StabChainStrong(S->stabilizer, {sch}, options);
      }
    }
  }
#ifdef DEBUG_STABCHAIN
  std::cerr << "CPP exiting StabChainStrong\n";
#endif
}


template<typename Telt, typename Tidx_label, typename Tint>
void ClosureGroup_options(StabChain<Telt,Tidx_label> & S, Telt const& g, StabChainOptions<Tint, typename Telt::Tidx> const& options)
{
  Telt sch = SiftedPermutation(S, g);
  if (!sch.isIdentity())
    StabChainStrong(S, {sch}, options);
}


template<typename Telt, typename Tidx_label, typename Tint>
void ClosureGroup(StabChain<Telt,Tidx_label> & S, Telt const& g)
{
  using Tidx = typename Telt::Tidx;
  Tidx n = S->comm->n;
  StabChainOptions<Tint,Tidx> options = GetStandardOptions<Tint,Tidx>(n);
  ClosureGroup_options(S, g, options);
}



template<typename Telt, typename Tidx_label>
bool StabChainForcePoint(StabChain<Telt,Tidx_label> & Stot, typename Telt::Tidx const& pnt)
{
#ifdef DEBUG_STABCHAIN
  std::cerr << "CPP Beginning of StabChainForcePoint pnt=" << int(pnt+1) << "\n";
  PrintStabChain(Stot);
  std::cerr << "DEBUG |Stot->transversal|=" << Stot->transversal.size() << "\n";
#endif
  if (Stot->transversal.size() == 0 || Stot->transversal[pnt] == std::numeric_limits<Tidx_label>::max()) {
#ifdef DEBUG_STABCHAIN
    std::cerr << "CPP Matching the first test\n";
#endif
    if (IsFixedStabilizer(Stot, pnt )) {
#ifdef DEBUG_STABCHAIN
      std::cerr << "CPP Matching the second test\n";
#endif
      InsertTrivialStabilizer(Stot, pnt);
    } else {
      if (!StabChainForcePoint(Stot->stabilizer, pnt) || !StabChainSwap(Stot)) {
#ifdef DEBUG_STABCHAIN
        std::cerr << "CPP StabChainForcePoint, return false\n";
#endif
	return false;
      }
    }
  }
#ifdef DEBUG_STABCHAIN
  std::cerr << "CPP StabChainForcePoint, return true\n";
  PrintStabChain(Stot);
#endif
  return true;
}


template<typename Telt, typename Tidx_label>
std::vector<Telt> GetListGenerators(StabChain<Telt,Tidx_label> const& Stot)
{
  size_t len = Stot->genlabels.size();
  std::vector<Telt> LGens(len);
  for (size_t i=0; i<len; i++)
    LGens[i] = Stot->comm->labels[Stot->genlabels[i]];
  return LGens;
}




template<typename Telt, typename Tidx_label>
bool StabChainSwap(StabChain<Telt,Tidx_label> & Stot)
{
  using Tidx=typename Telt::Tidx;
#ifdef DEBUG_STABCHAIN
  std::cerr << "CPP Beginning of StabChainSwap\n";
  PrintStabChain(Stot);
#endif
  Tidx n = Stot->comm->n;
  Tidx a = Stot->orbit[0];
  Tidx b = Stot->stabilizer->orbit[0];
#ifdef DEBUG_STABCHAIN
  std::cerr << "CPP StabChainSwap a=" << int(a+1) << " b=" << int(b+1) << "\n";
#endif
  //
  std::vector<Telt> LGens = GetListGenerators(Stot);
  // We have some missing entries in the S.generators.
  // Apparently, the S.generators contains generators that are not used.
#ifdef DEBUG_STABCHAIN
  std::cerr << "CPP LGens=" << GapStringTVector(LGens) << "\n";
#endif
  //
  //  StabChain<Telt> Ttot = EmptyStabChainPlusNode<Telt>(n, b);
  StabChain<Telt,Tidx_label> Ttot = EmptyStabChainPlusCommonPlusNode<Telt,Tidx_label>(Stot->comm, b);
  AddGeneratorsExtendSchreierTree(Ttot, LGens);
#ifdef DEBUG_STABCHAIN
  std::cerr << "CPP StabChainSwap : after first AGEST\n";
#endif
  //
  StabChain<Telt,Tidx_label> Tstab = EmptyStabChainPlusNode<Telt,Tidx_label>(n, a);
  if (Stot->stabilizer != nullptr) {
    if (Stot->stabilizer->stabilizer != nullptr) {
      std::vector<Telt> LGensB = GetListGenerators(Stot->stabilizer->stabilizer);
#ifdef DEBUG_STABCHAIN
      std::cerr << "CPP StabChainSwap : before second AGEST gens=" << GapStringTVector(LGensB) << "\n";
#endif
      AddGeneratorsExtendSchreierTree(Tstab, LGensB);
#ifdef DEBUG_STABCHAIN
      std::cerr << "CPP StabChainSwap : after second AGEST\n";
#endif
    }
  }
#ifdef DEBUG_STABCHAIN
  std::cerr << "CPP After AddGeneratorsExtendSchreierTree 1 : Tstab=\n";
  PrintStabChain(Tstab);
#endif
  //
  size_t ind = 0;
  size_t len = Stot->orbit.size() * Stot->stabilizer->orbit.size() / Ttot->orbit.size();
#ifdef DEBUG_STABCHAIN
  std::cerr << "CPP StabChainSwap |Tstab->orbit|=" << int(Tstab->orbit.size()) << " len=" << int(len) << "\n";
#endif
  while (Tstab->orbit.size() < len) {
#ifdef DEBUG_STABCHAIN
    std::cerr << "CPP beginning of loop\n";
#endif
    Tidx pnt;
    while(true) {
      ind++;
      if (ind >= Stot->orbit.size())
	return false;
#ifdef DEBUG_STABCHAIN
      std::cerr << "CPP |orbit|=" << Stot->orbit.size() << " ind=" << int(ind+1) << "\n";
#endif
      pnt = Stot->orbit[ind];
#ifdef DEBUG_STABCHAIN
      std::cerr << "CPP ind=" << int(ind+1) << " pnt=" << int(pnt+1) << "\n";
      std::cerr << "DEBUG |Tstab->transversal|=" << Tstab->transversal.size() << " pnt=" << pnt << "\n";
#endif
      if (Tstab->transversal[pnt] == std::numeric_limits<Tidx_label>::max())
	break;
    }
#ifdef DEBUG_STABCHAIN
    std::cerr << "CPP ind=" << int(ind+1) << " pnt=" << int(pnt+1) << "\n";
#endif
    Tidx img = b;
    Tidx i = pnt;
    while (i != a) {
      Tidx_label posGen=Stot->transversal[i];
      const Telt& x = Stot->comm->labels[posGen];
      img = PowAct(img, x);
      i = PowAct(i, x);
    }
#ifdef DEBUG_STABCHAIN
    std::cerr << "CPP i=" << int(i+1) << " img=" << int(img+1) << "\n";
#endif
    if (Stot->stabilizer->transversal[img] != std::numeric_limits<Tidx_label>::max()) {
#ifdef DEBUG_STABCHAIN
      std::cerr << "CPP Determining gen, step 1\n";
#endif
      Telt gen = Stot->comm->identity;
      while (true) {
        Tidx pnt_gen = PowAct(pnt,gen);
        if (pnt_gen == a)
          break;
#ifdef DEBUG_STABCHAIN
        std::cerr << "CPP pnt^gen=" << int(pnt_gen + 1) << "\n";
#endif
	Tidx_label posGen=Stot->transversal[pnt_gen];
#ifdef DEBUG_STABCHAIN
        std::cerr << "DEBUG 1 posGen=" << posGen << "\n";
#endif
	gen *= Stot->comm->labels[posGen];
      }
#ifdef DEBUG_STABCHAIN
      std::cerr << "CPP Determining gen, step 2 gen=" << gen << "\n";
#endif
      while (true) {
        Tidx b_gen = PowAct(b, gen);
        if (b_gen == b)
          break;
#ifdef DEBUG_STABCHAIN
        std::cerr << "CPP b^gen=" << int(b_gen+1) << "\n";
#endif
	Tidx_label posGen=Stot->stabilizer->transversal[b_gen];
#ifdef DEBUG_STABCHAIN
        std::cerr << "DEBUG 2 posGen=" << posGen << "\n";
#endif
	gen *= Stot->stabilizer->comm->labels[posGen];
      }
#ifdef DEBUG_STABCHAIN
      std::cerr << "CPP Determining gen, step 3 gen=" << gen << "\n";
#endif
      AddGeneratorsExtendSchreierTree(Tstab, {gen});
#ifdef DEBUG_STABCHAIN
      std::cerr << "CPP After AddGeneratorsExtendSchreierTree 2 : Tstab=\n";
      PrintStabChain(Tstab);
#endif
    }
  }
#ifdef DEBUG_STABCHAIN
  std::cerr << "CPP After while loop\n";
  std::cerr << "CPP T=\n";
  PrintStabChain(Ttot);
  std::cerr << "CPP Tstab=\n";
  PrintStabChain(Tstab);
#endif
  auto MapAtLevel=[&](StabChain<Telt,Tidx_label> & Swork, StabChain<Telt,Tidx_label> const& insStab) -> void {
    Swork->genlabels = insStab->genlabels;
    Swork->comm = insStab->comm;
    Swork->orbit = insStab->orbit;
    Swork->transversal = insStab->transversal;
  };
  MapAtLevel(Stot, Ttot);
#ifdef DEBUG_STABCHAIN
  PrintStabChain(Stot);
  std::cerr << "CPP StabChainSwap 1: |orbit|=" << Tstab->orbit.size() << "\n";
#endif
  if (Tstab->orbit.size() == 1) {
    Stot->stabilizer = Stot->stabilizer->stabilizer;
#ifdef DEBUG_STABCHAIN
    std::cerr << "CPP StabChainSwap 2:\n";
#endif
  } else {
    MapAtLevel(Stot->stabilizer, Tstab);
#ifdef DEBUG_STABCHAIN
    std::cerr << "CPP StabChainSwap 3:\n";
#endif
  }
#ifdef DEBUG_STABCHAIN
  PrintStabChain(Stot);
#endif
  return true;
}





// maybe use std::map<T, T> instead
template<typename T, typename Thom>
T LabsLims(T const& lab, Thom hom, std::vector<T> & labs, std::vector<T> & lims)
{
  size_t pos=PositionVect_ui<T,size_t>(labs, lab);
  if (pos == std::numeric_limits<size_t>::max()) {
    labs.push_back(lab);
    T img=hom(lab);
    lims.push_back(img);
  }
  return lims[pos];
}



// The original ConjugateStabChain code
// It is simplified from the original one with us being in the case
// IsPerm(hom) and IsPerm(map).
template<typename Telt, typename Tidx_label, typename Fhom, typename Fmap, typename Fcond>
StabChain<Telt,Tidx_label> ConjugateStabChain(StabChain<Telt,Tidx_label> & Stot, StabChain<Telt,Tidx_label> & Ttot, const Fhom& hom, const Fmap& map, const Fcond& cond)
{
  using Tidx=typename Telt::Tidx;
  Tidx n  = Stot->comm->n;
  Telt id = Stot->comm->identity;
  StabChain<Telt,Tidx_label> Sptr = Stot;
  StabChain<Telt,Tidx_label> Tptr = Ttot;

  using Tcomm=std::shared_ptr<CommonStabInfo<Telt>>;
  std::vector<Tcomm> ListLabels;
  std::vector<Tcomm> ListLabelsImg;

  auto get_comm=[&](Tcomm const& e_comm) -> Tcomm {
    for (size_t iLabel=0; iLabel<ListLabels.size(); iLabel++) {
      if (e_comm == ListLabels[iLabel])
        return ListLabelsImg[iLabel];
    }
    std::vector<Telt> labels = ListT(e_comm->labels, hom);
    Tcomm comm_new = std::make_shared<CommonStabInfo<Telt>>(CommonStabInfo<Telt>({n, id, labels}));
    ListLabels.push_back(e_comm);
    ListLabelsImg.push_back(comm_new);
    return comm_new;
  };

  std::vector<int> genlabels;
  while (cond(Sptr)) {
    if (Sptr->transversal.size() > 0) {
      std::vector<Tidx_label> NewTransversal(n);
      for (Tidx i=0; i<n; i++) {
	Tidx iImg=map(i);
	NewTransversal[iImg] = Sptr->transversal[i];
      }
      Sptr->transversal = NewTransversal;
    }
    Tptr->comm = get_comm(Sptr->comm);
    Tptr->genlabels = Sptr->genlabels;
    Tptr->orbit = ListT(Sptr->orbit, map);

    // Going to the next level.
    Sptr = Sptr->stabilizer;
    if (Tptr->stabilizer == nullptr)
      Tptr->stabilizer = EmptyStabChainPlusCommon<Telt,Tidx_label>(Tptr->comm);
    Tptr = Tptr->stabilizer;
  }
  //
  // Mapping the labels that showed up.
  //
  return Tptr;
}



template<typename Telt, typename Tidx_label>
void ConjugateStabChain_Element(StabChain<Telt,Tidx_label> & Stot, Telt const& cnj)
{
  using Tidx=typename Telt::Tidx;
  Telt cnjInv=~cnj;
  auto hom=[&](Telt const& x) -> Telt {
    return cnjInv*x*cnj;
  };
  auto map=[&](Tidx const& x) -> Tidx {
    return cnj.at(x);
  };
  auto cond=[](StabChain<Telt,Tidx_label> const& S) -> bool {
    return S->stabilizer != nullptr;
  };
  (void)ConjugateStabChain(Stot, Stot, hom, map, cond);
}


template<typename Telt, typename Tidx_label>
std::string PrintTopOrbit(StabChain<Telt,Tidx_label> const& S)
{
  if (S == nullptr) {
    return "unset";
  }
  size_t len=S->orbit.size();
  std::string str = "[ ";
  for (size_t u=0; u<len; u++) {
    if (u>0)
      str += ", ";
    str += std::to_string(S->orbit[u]+1);
  }
  str += " ]";
  return str;
}



// value of reduced
//  reduced = -1 corresponds to reduced = -1 in GAP code
//  reduced = 0  corresponds to reduced = false in GAP code
//  reduced = 1  corresponds to reduced = true in GAP code
template<typename Telt, typename Tidx_label>
bool ChangeStabChain(StabChain<Telt,Tidx_label> & Gptr, std::vector<typename Telt::Tidx> const& base, int const& reduced)
{
  using Tidx=typename Telt::Tidx;
#ifdef DEBUG_CHANGE_STAB_CHAIN
  std::string strG_orig=GetStringExpressionOfStabChain(Gptr);
  std::string strG_current=strG_orig;
  std::cerr << "CPP Beginning ChangeStabChain, GetStabilizerDepth = " << GetStabilizerDepth(Gptr) << "\n";
#endif
  Telt cnj = Gptr->comm->identity;
  StabChain<Telt,Tidx_label> Sptr = Gptr;
#ifdef DEBUG_CHANGE_STAB_CHAIN
  std::string strS_current=GetStringExpressionOfStabChain(Sptr);
  std::cerr << "CPP we have strS_current\n";

  auto KeyUpdating=[&](std::string const& str) {
    std::cerr << "CPP Gptr at str=" << str << "\n";
    PrintStabChain(Gptr);
    std::cerr << "CPP Sptr at str=" << str << "\n";
    PrintStabChain(Sptr);
  };
  KeyUpdating("Begin ChangeStabChain");
#endif
  std::vector<Tidx> newBase;
  size_t i=0;
  size_t basSiz=base.size();
#ifdef DEBUG_CHANGE_STAB_CHAIN
  std::cerr << "CPP ChangeStabChain base = " << GapStringIntVector(base) << "\n";
  std::cerr << "CPP ChangeStabChain 1 orbit=" << PrintTopOrbit(Gptr) << "\n";
#endif
  while (HasStabStab(Sptr) || i < basSiz) {
#ifdef DEBUG_CHANGE_STAB_CHAIN
    KeyUpdating("Before BasePoint");
#endif
    Tidx old=BasePoint(Sptr);
#ifdef DEBUG_CHANGE_STAB_CHAIN
    std::cerr << "CPP ChangeStabChain old=" << PosFalse_to_string(old) << " i=" << int(i+1) << " |base|=" << basSiz << "\n";
    KeyUpdating("After BasePoint");
#endif
    if (Sptr->genlabels.size() == 0 && (reduced == int_true || i >= basSiz)) {
#ifdef DEBUG_CHANGE_STAB_CHAIN
      std::cerr << "CPP Before RemoveStabChain\n";
      KeyUpdating("Before RemoveStabChain");
#endif
      RemoveStabChain(Sptr);
#ifdef DEBUG_CHANGE_STAB_CHAIN
      KeyUpdating("After RemoveStabChain");
#endif
      i = basSiz;
    } else if (i < basSiz) {
      Tidx newpnt = SlashAct(base[i], cnj);
#ifdef DEBUG_CHANGE_STAB_CHAIN
      std::cerr << "CPP newpnt=" << int(newpnt+1) << "\n";
#endif
      i++;
      if (reduced == int_reducedm1) {
	newBase.push_back(newpnt);
	if (newpnt != old) {
	  if (IsFixedStabilizer(Sptr, newpnt)) {
	    InsertTrivialStabilizer(Sptr, newpnt);
#ifdef DEBUG_CHANGE_STAB_CHAIN
            KeyUpdating("After InsertTrivialStabilizer");
#endif
	  }
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
	  else {
	    std::cerr << "CPP <base> must be an extension of base of <G>\n";
	    throw PermutalibException{1};
	  }
#endif
	}
	Sptr = Sptr->stabilizer;
#ifdef DEBUG_CHANGE_STAB_CHAIN
        KeyUpdating("After S:=S.stabilizer 1");
#endif
      } else if (reduced == int_false || !IsFixedStabilizer(Sptr, newpnt )) {
	if (Sptr->stabilizer != nullptr) {
#ifdef DEBUG_CHANGE_STAB_CHAIN
          KeyUpdating("Before StabChainForcePoint");
#endif
	  if (!StabChainForcePoint(Sptr, newpnt)) {
#ifdef DEBUG_CHANGE_STAB_CHAIN
            KeyUpdating("After StabChainForcePoint, return false");
#endif
	    return false;
          }
#ifdef DEBUG_CHANGE_STAB_CHAIN
          KeyUpdating("After StabChainForcePoint, return true");
          std::cerr << "CPP 1: cnj=" << cnj << "\n";
#endif
	  cnj = LeftQuotient(InverseRepresentative(Sptr, newpnt), cnj);
#ifdef DEBUG_CHANGE_STAB_CHAIN
          std::cerr << "CPP 2: cnj=" << cnj << "\n";
#endif
	} else {
	  InsertTrivialStabilizer(Sptr, newpnt);
#ifdef DEBUG_CHANGE_STAB_CHAIN
          KeyUpdating("After InsertTrivialStabilizer");
#endif
	}
	newBase.push_back(Sptr->orbit[0]);
	Sptr = Sptr->stabilizer;
#ifdef DEBUG_CHANGE_STAB_CHAIN
        KeyUpdating("After S:=S.stabilizer 2");
#endif
      }
    } else if (PositionVect_ui<Tidx,Tidx_label>(newBase, old) != std::numeric_limits<Tidx_label>::max() || (reduced == int_true && Sptr->orbit.size() == 1)) {
#ifdef DEBUG_CHANGE_STAB_CHAIN
      std::cerr << "CPP Stabilizer shift in ChangeStabChain\n";
#endif
      Sptr->comm = Sptr->stabilizer->comm;
      Sptr->genlabels = Sptr->stabilizer->genlabels;
#ifdef DEBUG_CHANGE_STAB_CHAIN
      KeyUpdating("MissCase 1");
#endif
      if (Sptr->stabilizer->orbit.size() > 0) {
        Sptr->orbit = Sptr->stabilizer->orbit;
        Sptr->transversal = Sptr->stabilizer->transversal;
      } else {
        Sptr->orbit.clear();
        Sptr->transversal.clear();
      }
#ifdef DEBUG_CHANGE_STAB_CHAIN
      KeyUpdating("MissCase 2");
#endif
      if (Sptr->stabilizer->stabilizer != nullptr) {
        Sptr->stabilizer = Sptr->stabilizer->stabilizer;
      } else {
        Sptr->stabilizer = nullptr;
      }
#ifdef DEBUG_CHANGE_STAB_CHAIN
      KeyUpdating("MissCase 3");
#endif
    } else {
      newBase.push_back(old);
      Sptr = Sptr->stabilizer;
#ifdef DEBUG_CHANGE_STAB_CHAIN
      KeyUpdating("After S:=S.stabilizer 3");
#endif
    }
  }
#ifdef DEBUG_CHANGE_STAB_CHAIN
  std::cerr << "CPP LEAVE i=" << int(i+1) << " |base|=" << basSiz << "\n";
  KeyUpdating("After the loop");
  std::cerr << "CPP Before ConjugateStabChain cnj=" << cnj << "\n";
#endif
  if (!cnj.isIdentity())
    ConjugateStabChain_Element(Gptr, cnj);
#ifdef DEBUG_CHANGE_STAB_CHAIN
  KeyUpdating("After ConjugateStabChain");
  std::cerr << "CPP Leaving ChangeStabChain\n";
#endif
  return true;
}

template<typename Telt, typename Tidx_label>
bool ExtendStabChain(StabChain<Telt,Tidx_label> & Stot, std::vector<typename Telt::Tidx> const& base)
{
#ifdef DEBUG_STABCHAIN
  std::cerr << "CPP Beginning of ExtendStabChain\n";
#endif
  return ChangeStabChain(Stot, base, int_reducedm1);
}


template<typename Telt, typename Tidx_label>
bool ReduceStabChain(StabChain<Telt,Tidx_label> & Stot)
{
  return ChangeStabChain(Stot, {}, int_true);
}




template<typename Telt, typename Tidx_label>
bool TestEqualityStabChain(StabChain<Telt,Tidx_label> const& L, StabChain<Telt,Tidx_label> const& R)
{
  StabChain<Telt,Tidx_label> Lptr = L;
  StabChain<Telt,Tidx_label> Rptr = R;
  //#define DEBUG_EQUALITY
  while(true) {
    if (Lptr == nullptr && Rptr != nullptr) {
#ifdef DEBUG_EQUALITY
      std::cerr << "TestEqualityStabChain 1\n";
#endif
      return false;
    }
    if (Lptr != nullptr && Rptr == nullptr) {
#ifdef DEBUG_EQUALITY
      std::cerr << "TestEqualityStabChain 2\n";
#endif
      return false;
    }
    if (Lptr == nullptr)
      break;
    if (Lptr->orbit != Rptr->orbit) {
#ifdef DEBUG_EQUALITY
      std::cerr << "TestEqualityStabChain 3\n";
#endif
      return false;
    }
    if (Lptr->transversal.size() != Rptr->transversal.size()) {
#ifdef DEBUG_EQUALITY
      std::cerr << "TestEqualityStabChain 4\n";
#endif
      return false;
    }
    size_t lenL=Lptr->transversal.size();
    Tidx_label miss_val = std::numeric_limits<Tidx_label>::max();
    for (size_t u=0; u<lenL; u++) {
      if (Lptr->transversal[u] == miss_val && Rptr->transversal[u] != miss_val) {
#ifdef DEBUG_EQUALITY
        std::cerr << "TestEqualityStabChain 5\n";
#endif
	return false;
      }
      if (Lptr->transversal[u] != miss_val && Rptr->transversal[u] == miss_val) {
#ifdef DEBUG_EQUALITY
        std::cerr << "TestEqualityStabChain 6\n";
#endif
	return false;
      }
      if (Lptr->transversal[u] != miss_val) {
	Tidx_label idxL=Lptr->transversal[u];
	Tidx_label idxR=Rptr->transversal[u];
	const Telt& permL=Lptr->comm->labels[idxL];
	const Telt& permR=Rptr->comm->labels[idxR];
	if (permL != permR) {
#ifdef DEBUG_EQUALITY
          std::cerr << "TestEqualityStabChain 7\n";
#endif
	  return false;
        }
      }
    }
    Rptr = Rptr->stabilizer;
    Lptr = Lptr->stabilizer;
  }
#ifdef DEBUG_EQUALITY
  std::cerr << "TestEqualityStabChain 8\n";
#endif
  return true;
}



template<typename Telt, typename Tidx_label>
void SetStabChainFromLevel(std::vector<StabChain<Telt,Tidx_label>> & R_list,
                           std::vector<StabChain<Telt,Tidx_label>> const& L_list,
                           size_t const& levbegin, size_t const& levend)
{
  for (size_t iLev=levbegin; iLev<levend; iLev++)
    R_list[iLev] = StructuralCopy(L_list[iLev]);
}




template<typename Telt, typename Tidx_label>
bool IsFixedStabilizer(StabChain<Telt,Tidx_label> const& S, typename Telt::Tidx const& pnt)
{
  for (const Tidx_label & posGen : S->genlabels)
    if (pnt != PowAct(pnt, S->comm->labels[posGen]))
      return false;
  return true;
}


template<typename Telt, typename Tidx_label>
Telt MinimalElementCosetStabChain(StabChain<Telt,Tidx_label> const& Stot, Telt const& g)
{
  using Tidx=typename Telt::Tidx;
  Telt gRet=g;
  StabChain<Telt,Tidx_label> Sptr = Stot;
  while (Sptr != nullptr) {
    if (Sptr->genlabels.size() == 0)
      return gRet;
    Tidx pMin=std::numeric_limits<Tidx>::max();
    for (auto & i : Sptr->orbit) {
      Tidx a=PowAct(i,gRet);
      if (a < pMin)
	pMin = a;
    }
    Tidx bp=Sptr->orbit[0];
    Tidx pp=SlashAct(pMin, gRet);
    while (bp != pp) {
      Tidx_label pos=Sptr->transversal[pp];
      gRet = LeftQuotient(Sptr->comm->labels[pos], gRet);
      pp = SlashAct(pMin, gRet);
    }
    Sptr = Sptr->stabilizer;
  }
  return gRet;
}






template<typename Tret, typename Telt, typename Tidx_label, typename F_pt, typename F_elt>
StabChain<Tret,Tidx_label> HomomorphismMapping(StabChain<Telt,Tidx_label> const& Stot, F_pt f_pt, F_elt f_elt)
{
  //  using Tidx = typename Telt::Tidx;
  using Tret_idx = typename Tret::Tidx;
  Tret idMap = f_elt(Stot->comm->identity);
  Tret_idx nMap=idMap.size();
  auto fVector =[&](std::vector<Telt> const& V) -> std::vector<Tret> {
    size_t len = V.size();
    std::vector<Tret> Vret(len);
    for (size_t i=0; i<len; i++)
      Vret[i] = f_elt(V[i]);
    return Vret;
  };



  std::vector<std::pair< std::shared_ptr<CommonStabInfo<Telt>>, std::shared_ptr<CommonStabInfo<Tret>> >> ListPairComm;
  auto get_mapped_comm=[&](const std::shared_ptr<CommonStabInfo<Telt>>& comm_in) -> std::shared_ptr<CommonStabInfo<Tret>> {
    for (auto & e_pair : ListPairComm)
      if (e_pair.first == comm_in)
        return e_pair.second;
    std::vector<Tret> labelsMap = fVector(comm_in->labels);
    std::shared_ptr<CommonStabInfo<Tret>> comm_out = std::make_shared<CommonStabInfo<Tret>>(CommonStabInfo<Tret>({nMap, idMap, std::move(labelsMap)}));
    ListPairComm.push_back({comm_in, comm_out});
    return comm_out;
  };




  std::vector<StabChain<Tret,Tidx_label>> ListLevel;

  StabChain<Telt,Tidx_label> Sptr = Stot;
  while (Sptr != nullptr) {
    std::vector<Tret_idx> orbit;
    for (auto & eVal : Sptr->orbit)
      orbit.push_back(f_pt(eVal));
    StabChain<Tret,Tidx_label> Swork = std::make_shared<StabLevel<Telt,Tidx_label>>(StabLevel<Telt,Tidx_label>({Sptr->transversal, orbit, Sptr->genlabels, Sptr->cycles, Sptr->IsBoundCycle, fVector(Sptr->treegen), fVector(Sptr->treegeninv), fVector(Sptr->aux), Sptr->treedepth, Sptr->diam, get_mapped_comm(Sptr->comm), nullptr}));
    ListLevel.push_back(Swork);
    Sptr = Sptr->stabilizer;
  }
  size_t n_lev = ListLevel.size();
  for (size_t i=0; i<n_lev-1; i++)
    ListLevel[i]->stabilizer = ListLevel[i+1];
  return ListLevel[0];
}










}


#endif
