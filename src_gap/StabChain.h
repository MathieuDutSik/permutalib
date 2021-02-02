#ifndef DEFINE_STAB_CHAIN
#define DEFINE_STAB_CHAIN

#define DEBUG_PERMUTALIB

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
#include "PermGroup.h"



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



template<typename T>
void AssignationVectorGapStyle(std::vector<T> & eVect, int const& pos, T const& val)
{
  int siz=eVect.size();
  if (pos < siz) {
    eVect[pos] = val;
    return;
  }
#ifdef DEBUG
  if (pos != siz) {
    std::cerr << "Assignation leaves gap in the vector. Not allowed\n";
    throw TerminalException{1};
  }
#endif
  eVect.push_back(val);
}


template<typename T>
void PrintVectDebug(std::string const& str, std::vector<T> const& V)
{
  std::cerr << str << " =";
  for (auto & eVal : V)
    std::cerr << " " << eVal;
  std::cerr << "\n";
}



template<typename T, typename Telt>
std::vector<T> PermutedAct(std::vector<T> const& V, Telt const& g)
{
  int len=V.size();
  std::vector<T> Vret(len);
  for (int i=0; i<len; i++) {
    int iImg=g.at(i);
    Vret[iImg] = V[i];
  }
  return Vret;
}


template<typename Telt>
int GetLabelIndex(std::vector<Telt> & labels, Telt const& u)
{
  int nbLabel=labels.size();
  for (int iLabel=0; iLabel<nbLabel; iLabel++)
    if (labels[iLabel] == u)
      return iLabel;
  labels.push_back(u);
  return nbLabel;
}

template<typename Telt>
int GetLabelIndex_const(std::vector<Telt> const& labels, Telt const& u)
{
  int nbLabel=labels.size();
  for (int iLabel=0; iLabel<nbLabel; iLabel++)
    if (labels[iLabel] == u)
      return iLabel;
  return -1;
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
  int n;
  Telt identity;
  std::vector<Telt> labels;
};

template<typename Telt>
struct StabLevel {
  std::vector<int> transversal;
  std::vector<int> orbit;
  std::vector<int> genlabels; // Not used in Random algorithm
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
  std::shared_ptr<StabLevel<Telt>> stabilizer;
};

template<typename Telt>
using StabChain = std::shared_ptr<StabLevel<Telt>>;
// other possible entries:
// transimages, genimages, labelimages, idimage


template<typename Telt>
bool IsIdenticalObj(StabChain<Telt> const& S1, StabChain<Telt> const& S2)
{
  return S1 == S2;
}





template<typename Telt>
StabChain<Telt> StabChainGenerators(std::vector<Telt> const& generators, int const& n, Telt const& id)
{
  std::shared_ptr<CommonStabInfo<Telt>> comm = std::make_shared<CommonStabInfo<Telt>>(CommonStabInfo<Telt>({n, id, generators}));
  //
  std::vector<int> transversal = std::vector<int>(n, -1);
  std::vector<int> orbit;
  std::vector<int> genlabels;
  for (int igen=0; igen<int(generators.size()); igen++)
    genlabels.push_back(igen);
  std::vector<int8_t> cycles;
  bool IsBoundCycle = false;
  std::vector<Telt> treegen;
  std::vector<Telt> treegeninv;
  std::vector<Telt> aux;
  int treedepth = 0;
  int diam = 0;
  return std::make_shared<StabLevel<Telt>>(StabLevel<Telt>({transversal, orbit, genlabels, cycles, IsBoundCycle, treegen, treegeninv, aux, treedepth, diam, comm, nullptr}));
}






template<typename Telt>
void PrintStabChainTransversals(StabChain<Telt> const& S)
{
  StabChain<Telt> Swork = S;
  int n=Swork->comm->n;
  int iLevel=0;
  while (Swork != nullptr) {
    std::vector<std::optional<Telt>> V(n);
    for (int i=0; i<n; i++) {
      int eVal = Swork->transversal[i];
      if (eVal == -1)
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

template<typename Telt>
void UnbindCycles(StabChain<Telt> const& S)
{
  StabChain<Telt> Swork = S;
  while (Swork != nullptr) {
    Swork->cycles.clear();
    Swork->IsBoundCycle = false;
    Swork = Swork->stabilizer;
  }
}



template<typename Telt>
void PrintStabChainOrbits(StabChain<Telt> const& S)
{
  StabChain<Telt> Swork = S;
  int iLevel=0;
  while (Swork != nullptr) {
    std::cerr << "CPP i=" << iLevel << " orbit=" << GapStringIntVector(Swork->orbit) << "\n";
    Swork = Swork->stabilizer;
    iLevel++;
  }
}


template<typename Telt>
std::string GetListStabCommPartition(std::vector<StabChain<Telt>> const& ListS)
{
  int len = ListS.size();
  std::vector<int> Status(len,0);
  std::vector<std::string> ListStr;
  for (int i=0; i<len; i++) {
    if (Status[i] == 0) {
      std::vector<int> LVal;
      auto ptr = ListS[i]->comm;
      for (int j=0; j<len; j++) {
        if (ListS[j]->comm == ptr) {
          LVal.push_back(j);
          Status[j] = 1;
        }
      }
      std::string estr = GapStringIntVector(LVal);
      ListStr.push_back(estr);
    }
  }
  return GapStringTVector(ListStr);
}


template<typename Telt>
void PrintStabChain(StabChain<Telt> const& S)
{
  StabChain<Telt> Swork = S;
  int n = Swork->comm->n;
  std::cerr << "CPP Partition=" << GetListStabCommPartition(ListStabChain(S)) << "\n";
  int iLevel = 0;
  while (Swork != nullptr) {
    std::cerr << "CPP iLev=" << iLevel << "\n";
    std::string strTransversal = "[ ]";
    if (Swork->transversal.size() > 0) {
      if (int(Swork->transversal.size()) != n) {
        std::cerr << "Swork->transversal should be of length 0 or n=" << n << "\n";
        throw TerminalException{1};
      }
      std::vector<std::optional<Telt>> V(n);
      for (int i=0; i<n; i++) {
        int eVal = Swork->transversal[i];
        if (eVal == -1)
          V[i] = {};
        else
          V[i] = Swork->comm->labels[eVal];
      }
      strTransversal = GapStringMissingTVector(V);
    }
    //
    std::cerr << "CPP   orbit=" << GapStringIntVector(Swork->orbit) << "\n";
    std::cerr << "CPP   transversal=" << strTransversal << "\n";
    //    std::cerr << "XXX ELIMINATE begin\n";
    if (Swork->IsBoundCycle) {
      std::cerr << "CPP   cycles=" << GapStringBoolVectorB(Swork->cycles) << "\n";
    } else {
      std::cerr << "CPP   No cycles\n";
    }
    //    std::cerr << "XXX ELIMINATE end\n";
    Swork = Swork->stabilizer;
    iLevel++;
  }
}


template<typename Telt>
void PrintListStabCommPartition(std::string const& mesg, std::vector<StabChain<Telt>> const& ListS)
{
  std::cerr << mesg << " ListStabCommPartition=" << GetListStabCommPartition(ListS) << "\n";
}



template<typename Telt>
int GetStabilizerDepth(StabChain<Telt> const& S1)
{
  StabChain<Telt> S2 = S1;
  int dep = 0;
  while(true) {
    if (S2 == nullptr)
      break;
    dep++;
    S2 = S2->stabilizer;
  }
  return dep;
}


template<typename Telt>
std::string perm_to_string(Telt const& eVal)
{
  std::ostringstream os;
  os << eVal;
  std::string str=os.str();
  return str;
}


template<typename Telt>
std::string GetStringExpressionOfStabChain(StabChain<Telt> const& eRec)
{
  std::string strRet="record_";
  StabChain<Telt> eStab = eRec;
  //srd::cerr << "GetStringExpressionOfStabChain, step 1\n";
  while(true) {
    //srd::cerr << "GetStringExpressionOfStabChain, step 2\n";
    if (eStab == nullptr)
      break;
    //srd::cerr << "GetStringExpressionOfStabChain, step 3\n";
    strRet += "orbit_[";
    for (auto & eVal : eStab->orbit)
      strRet += " " + std::to_string(eVal+1);
    strRet += "]";
    //srd::cerr << "GetStringExpressionOfStabChain, step 4\n";
    strRet += "_transversal_";
    for (auto & eVal : eStab->transversal) {
      if (eVal == -1)
        strRet += " " + std::to_string(eVal);
      else
        strRet += " " + perm_to_string(eStab->comm->labels[eVal]);
    }
    //srd::cerr << "GetStringExpressionOfStabChain, step 5\n";
    eStab = eStab->stabilizer;
    //srd::cerr << "GetStringExpressionOfStabChain, step 6\n";
  }
  //srd::cerr << "GetStringExpressionOfStabChain, step 7\n";
  return strRet;
}




template<typename Telt>
std::ostream& operator<<(std::ostream& os, StabChain<Telt> const& Stot)
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
  int iLev=0;
  StabChain<Telt> Sptr = Stot;
  while(true) {
    if (Sptr == nullptr)
      break;
    os << "CPP iLev=" << iLev << "\n";
    os << "CPP  orbit =";
    for (auto & eVal : Sptr->orbit) {
      os << " " << eVal+1;
    }
    os << "\n";
    os << "CPP  transversal =";
    for (auto & eVal : Sptr->transversal) {
      if (eVal == -1)
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


template<typename Telt>
StabChain<Telt> ShallowCopy(StabChain<Telt> const& S)
{
  StabLevel<Telt> Sret = *S;
  Sret.comm = S->comm;
  Sret.stabilizer = S->stabilizer;
  return std::make_shared<StabLevel<Telt>>(Sret);
}



template<typename Telt>
StabChain<Telt> StructuralCopy(StabChain<Telt> const& S)
{
  if (S == nullptr)
    return nullptr;
  using Tcomm=std::shared_ptr<CommonStabInfo<Telt>>;
  std::vector<Tcomm> ListLabels;
  std::vector<Tcomm> ListLabelsImg;
  auto get_comm=[&](Tcomm const& e_comm) -> Tcomm {
     for (size_t i=0; i<ListLabels.size(); i++) {
       if (ListLabels[i] == e_comm)
         return ListLabelsImg[i];
     }
     Tcomm comm_new = std::make_shared<CommonStabInfo<Telt>>(*e_comm);
     ListLabels.push_back(e_comm);
     ListLabelsImg.push_back(comm_new);
     return comm_new;
  };
  StabChain<Telt> Sptr = S;
  std::vector<std::shared_ptr<StabLevel<Telt>>> ListPtr;
  while(true) {
    if (Sptr == nullptr)
      break;
    ListPtr.push_back(Sptr);
    Sptr = Sptr->stabilizer;
  }
  size_t len = ListPtr.size();
  StabChain<Telt> Sret = nullptr;
  for (size_t i=0; i<len; i++) {
    size_t j = len - 1 - i;
    std::shared_ptr<StabLevel<Telt>> S2 = std::make_shared<StabLevel<Telt>>(*ListPtr[j]);
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
    throw TerminalException{1};
  }
#endif
  return Sret;
}



template<typename Telt>
StabChain<Telt> RestrictedStabChain(StabChain<Telt> const& Stot, int const& eLev)
{
  int nbLevel=Stot.stabilizer.size();
  std::vector<StabLevel<Telt>> stabilizerRed;
  for (int uLev=eLev; uLev<nbLevel; uLev++)
    stabilizerRed.push_back(Stot.stabilizer[uLev]);
  return {Stot.n, Stot.identity, Stot.labels, stabilizerRed};
}


template<typename Telt>
StabLevel<Telt> EmptyStabLevel(std::shared_ptr<CommonStabInfo<Telt>> const& comm)
{
  std::vector<int> transversal = std::vector<int>(comm->n, -1);
  std::vector<int> orbit;
  std::vector<int> genlabels;
  std::vector<int8_t> cycles;
  bool IsBoundCycle = false;
  std::vector<Telt> treegen;
  std::vector<Telt> treegeninv;
  std::vector<Telt> aux;
  int treedepth = 0;
  int diam = 0;
  return {transversal, orbit, genlabels, cycles, IsBoundCycle, treegen, treegeninv, aux, treedepth, diam, comm, nullptr};
}



template<typename Telt>
StabChain<Telt> EmptyStabChain(int const& n)
{
  Telt id(n);
  std::shared_ptr<CommonStabInfo<Telt>> comm = std::make_shared<CommonStabInfo<Telt>>(CommonStabInfo<Telt>({n, id, {id}}));
  return std::make_shared<StabLevel<Telt>>(EmptyStabLevel<Telt>(comm));
}


template<typename Telt>
StabChain<Telt> EmptyStabChainPlusNode(int const& n, int const& bas)
{
  StabChain<Telt> S = EmptyStabChain<Telt>(n);
  InitializeSchreierTree(S, bas);
  return S;
}


template<typename Telt>
StabChain<Telt> EmptyStabChainPlusCommon(std::shared_ptr<CommonStabInfo<Telt>> const& comm)
{
  StabChain<Telt> S = std::make_shared<StabLevel<Telt>>(EmptyStabLevel<Telt>(comm));
  return S;
}


template<typename Telt>
StabChain<Telt> EmptyStabChainPlusCommonPlusNode(std::shared_ptr<CommonStabInfo<Telt>> const& comm, int const& bas)
{
  StabChain<Telt> S = std::make_shared<StabLevel<Telt>>(EmptyStabLevel<Telt>(comm));
  InitializeSchreierTree(S, bas);
  return S;
}







template<typename Telt>
int BasePoint(StabChain<Telt> const& S)
{
  if (S == nullptr)
    return -1;
  if (S->orbit.size() == 0)
    return -1;
  return S->orbit[0];
}





// almost certainly wrong code below.
template<typename Telt>
void RemoveStabChain(StabChain<Telt> & Stot)
{
  Stot->stabilizer = nullptr;
  Stot->genlabels.clear();
  Stot->orbit.clear();
  Stot->transversal.clear();
  Stot->cycles.clear();
  Stot->IsBoundCycle = false;
  Stot->treegen.clear();
  Stot->treegeninv.clear();
  Stot->aux.clear();
}


template<typename Telt>
bool IsInBasicOrbit(StabChain<Telt> const& S, int const& pnt)
{
  int eVal=S->transversal[pnt];
  if (eVal == -1)
    return false;
  return true;
}

// This correspond to the code of
// pnt^g in GAP code.
/*
template<typename Telt>
int PowAct(int const& pnt, Telt const& g)
{
  return g.at(pnt);
}
*/

template<typename Telt>
Telt InverseRepresentative(StabChain<Telt> const& S, int const& pnt)
{
  int bpt=S->orbit[0];
  Telt rep=S->comm->identity;
  int pntw=pnt;
#undef DEBUG_INV_REP
#ifdef DEBUG_INV_REP
  std::cerr << "CPP INVREP bpt=" << (bpt+1) << " pnt=" << (pntw+1) << "\n";
#endif
  while(pntw != bpt) {
    int idx=S->transversal[pntw];
    Telt te=S->comm->labels[idx];
#ifdef DEBUG_INV_REP
    std::cerr << "CPP INVREP te=" << te << "\n";
#endif
    pntw=PowAct(pntw, te);
#ifdef DEBUG_INV_REP
    std::cerr << "CPP INVREP   pnt=" << (pntw+1) << "\n";
#endif
    rep = rep * te;
#ifdef DEBUG_INV_REP
    std::cerr << "CPP INVREP   rep=" << rep << "\n";
#endif
  }
#ifdef DEBUG_INV_REP
  std::cerr << "CPP INVREP return rep=" << rep << "\n";
#endif
  return rep;
}


template<typename Telt>
std::vector<Telt> InverseRepresentativeWord(StabChain<Telt> const& S, int const& pnt)
{
  int bpt=S->orbit[0];
  std::vector<Telt> word;
  int pntw=pnt;
  while(pntw != bpt) {
    int idx=S->transversal[pntw];
    Telt te=S->comm->labels[idx];
    pntw = PowAct(pntw, te);
    word.push_back(te);
  }
  return word;
}


template<typename Telt>
Telt SiftedPermutation(StabChain<Telt> const& S, Telt const& g)
{
  Telt gW=g;
  //  std::cerr << "CPP Before shared pointer copy\n";
  StabChain<Telt> Sptr = S;
  //  std::cerr << "CPP After shared pointer copy\n";
  while(true) {
    if (Sptr->stabilizer == nullptr || gW.isIdentity())
      return gW;
    int bpt = Sptr->orbit[0];
    int img = PowAct(bpt, gW);
    if (Sptr->transversal[img] == -1)
      return gW;
    while(true) {
      if (img == bpt)
	break;
      int idx = Sptr->transversal[img];
      gW = gW * Sptr->comm->labels[idx];
      img = PowAct(bpt, gW);
    }
    Sptr = Sptr->stabilizer;
  }
  return gW;
}


template<typename Telt>
std::vector<int> BaseStabChain(StabChain<Telt> const& S)
{
  StabChain<Telt> Sptr = S;
  std::vector<int> base;
  while(true) {
    if (Sptr == nullptr)
      break;
    base.push_back(Sptr->orbit[0]);
    Sptr = Sptr->stabilizer;
  }
  return base;
}



template<typename Telt, typename Tint>
Tint SizeStabChain(StabChain<Telt> const& S)
{
  Tint size=1;
  StabChain<Telt> Sptr = S;
  while(true) {
    if (Sptr == nullptr)
      break;
    int siz=Sptr->orbit.size();
    if (siz == 0)
      break;
    Tint siz_i = siz;
    size *= siz_i;
    Sptr = Sptr->stabilizer;
  }
  return size;
}

template<typename Telt>
std::vector<Telt> StrongGeneratorsStabChain(StabChain<Telt> const& S)
{
  StabChain<Telt> Sptr = S;
  std::set<Telt> sgs_set;
  while(true) {
    if (Sptr == nullptr)
      break;
    std::vector<Telt> const& labels = Sptr->comm->labels;
    for (auto & pos : Sptr->genlabels)
      sgs_set.insert(labels[pos]);
    Sptr = Sptr->stabilizer;
  }
  std::vector<Telt> sgs;
  for (auto & ePos : sgs_set)
    sgs.push_back(ePos);
  return sgs;
}



template<typename Telt>
std::vector<Telt> GeneratorsStabChain(StabChain<Telt> const& S)
{
  std::vector<Telt> sgs;
  for (auto & ePos : S->genlabels) {
    std::cerr << "DEBUG ePos=" << ePos << "\n";
    sgs.push_back(S->comm->labels[ePos]);
  }
  return sgs;
}





template<typename Telt>
std::vector<Telt> GeneratorsOfGroup(StabChain<Telt> const& S)
{
  return StrongGeneratorsStabChain(S);
}





template<typename Telt>
Telt LargestElementStabChain(StabChain<Telt> const& S)
{
  StabChain<Telt> Sptr = S;

  Telt rep=Sptr->identity;
  while(true) {
    if (Sptr == nullptr)
      break;
    if (Sptr->genlabels.size() == 0)
      break;
    int pnt=Sptr->orbit[0];
    int min=0;
    int val=0;
    for (auto & i : Sptr->orbit) {
      int img=PowAct(i, rep);
      if (img > val) {
	min=i;
	val=img;
      }
    }
    while(true) {
      if (pnt == min)
	break;
      int idx=Sptr->transversal[min];
      Telt gen=Sptr->comm->labels[idx];
      rep=LeftQuotient(gen,rep);
      min=PowAct(min, gen);
    }
    Sptr = Sptr->stabilizer;
  }
  return rep;
}


/*
template<typename Telt>
std::vector<Telt> ElementsStabChain(StabChain<Telt> const& Stot)
{
  std::vector<Telt> elms;
  auto LevelIncrease=[&](StabLevel<Telt> const& eLev) -> void {
    std::vector<Telt> NewElms;
    for (auto & pnt : eLev.orbit) {
      Telt rep=eLev.identity;
      while (PowAct(eLev.orbit[0], rep) != pnt) {
	int jpt=SlashAct(pnt, rep);
	int idx=Stot.stabilizer[eLev].transversal[jpt];
	rep=LeftQuotient(Stot->comm->labels[idx], rep);
      }
      for (auto & eStb : elms)
	NewElms.push_back(eStb * rep);
    }
    elms=NewElms;
  };
  elms={Stot.stabilizer[0].identity};
  int len=Stot.stabilizer.size();
  for (int iLev=0; iLev<len; iLev++) {
    int jLev=len-1-iLev;
    LevelIncrease(Stot.stabilizer[jLev]);
  }
  return elms;
}
*/




// is base is empty then this just replaces the IsBound(options.base)
template<typename Tint>
struct StabChainOptions {
  int n;
  std::vector<int> base;
  std::vector<int> knownBase;
  int random;
  bool reduced;
  Tint size;
  Tint limit;
};


template<typename Tint>
StabChainOptions<Tint> GetStandardOptions(int const& n)
{
  std::vector<int> base;
  std::vector<int> knownBase;
  int random = 1000;
  bool reduced=true;
  Tint size=0;
  Tint limit=0;
  return {n, base, knownBase, random, reduced, size, limit};
}


template<typename Telt>
bool IsTrivial_ListGen(std::vector<Telt> const& LGen)
{
  for (auto & eElt : LGen)
    if (!eElt.isIdentity())
      return false;
  return true;
}

template<typename Telt>
std::vector<int> MovedPoints(StabChain<Telt> const& S)
{
  std::set<int> LIdx;
  StabChain<Telt> Sptr = S;
  while(true) {
    if (Sptr == nullptr)
      break;
    for (auto & eIdx : Sptr->genlabels)
      LIdx.insert(eIdx);
    Sptr = Sptr->stabilizer;
  }
  auto IsMoved=[&](int const& ePt) -> bool {
    for (auto & eIdx : LIdx) {
      if (S->comm->labels[eIdx].at(ePt) != ePt)
	return true;
    }
    return false;
  };
  std::vector<int> LMoved;
  int n=S->comm->n;
  for (int i=0; i<n; i++)
    if (IsMoved(i))
      LMoved.push_back(i);
  return LMoved;
}


template<typename Telt>
bool IsTrivial(StabChain<Telt> const& G)
{
  StabChain<Telt> Sptr = G;
  while(true) {
    if (Sptr == nullptr)
      break;
    for (auto & eIdx : Sptr->genlabels)
      if (!Sptr->comm->labels[eIdx].isIdentity())
        return false;
    Sptr = Sptr->stabilizer;
  }
  return true;
}






template<typename Telt>
int LargestMovedPoint(std::vector<Telt> const& LGen)
{
  if (LGen.size() == 0)
    return -1;
  int n=LGen[0].size();
  std::vector<int> Status(n, 1);
  for (auto & eGen : LGen) {
    for (int u=0; u<n; u++) {
      int v=eGen.at(u);
      if (u != v)
	Status[u]=0;
    }
  }
  int eMov=0;
  for (int u=0; u<n; u++)
    if (Status[u] == 0)
      eMov=u;
  eMov++;
  return eMov;
}

template<typename Telt>
void InitializeSchreierTree(StabChain<Telt> & S, int const& pnt)
{
  //  std::cerr << "Beginning of InitializeSchreierTree\n";
  int n=S->comm->n;
  //
  S->orbit = {pnt};
  //
  std::vector<int> transversal(n, -1);
  transversal[pnt] = 0;
  S->transversal = transversal;
  //  std::cerr << "   Ending of InitializeSchreierTree\n";
}

template<typename Telt>
void InsertTrivialStabilizer(StabChain<Telt> & Stot, int const& pnt)
{
  std::cerr << "CPP begin InsertTrivialStabilizer\n";
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


template<typename Telt, typename Tint>
StabChain<Telt> StabChainOp_trivial_group(StabChain<Telt> const& Stot, StabChainOptions<Tint> const& options)
{
  std::cerr << "CPP Call to StabChainOp (trivial group)\n";
  StabChain<Telt> S = EmptyStabChain<Telt>(Stot->comm->n);
  if (options.base.size() > 0 && !options.reduced) {
    StabChain<Telt> T = S;
    for (int const& pnt : options.base) {
      InsertTrivialStabilizer( T, pnt );
      T = T->stabilizer;
    }
  }
  return S;
}




template<typename Telt>
std::vector<StabChain<Telt>> ListStabChain(StabChain<Telt> const& S)
{
  std::vector<StabChain<Telt>> ListStab;
  StabChain<Telt> Sptr = S;
  while(true) {
    if (Sptr == nullptr)
      break;
    ListStab.push_back(Sptr);
    Sptr = Sptr->stabilizer;
  }
  return ListStab;
}




template<typename Telt>
StabChain<Telt> StabChainBaseStrongGenerators(std::vector<int> const& base, std::vector<Telt> const& sgs)
{
#ifdef DEBUG_PERMUTALIB
  if (sgs.size() == 0) {
    std::cerr << "sgs is empty. StabChainBaseStrongGenerators is broken in that case\n";
    throw TerminalException{1};
  }
#endif
  int n=sgs[0].size();
  StabChain<Telt> S = EmptyStabChain<Telt>(n);
  int nbGen=sgs.size();
  Face status(nbGen);
  for (int i=0; i<nbGen; i++)
    status[i] = 1;
  int basSiz=base.size();
  for (int iBas=0; iBas<basSiz; iBas++) {
    int pnt=base[iBas];
    std::vector<Telt> sgsFilt;
    for (int i=0; i<nbGen; i++)
      if (status[i] == 1)
	sgsFilt.push_back(sgs[i]);
    InsertTrivialStabilizer(S, pnt);
    std::cerr << "CPP Before call to AddGeneratorsExtendSchreierTree from StabChainStrongGenerators\n";
    AddGeneratorsExtendSchreierTree(S, iBas, sgsFilt);
    for (int i=0; i<nbGen; i++)
      if (status[i] == 1 && PowAct(pnt, sgs[i]) != pnt)
	status[i]=0;
  }
  return S;
}

template<typename T>
std::vector<T> SortVector(std::vector<T> const& f)
{
  std::vector<T> RetF = f;
  sort(RetF.begin(), RetF.end(),
       [](T const& x, T const& y) -> bool {
         if (x<y)
           return true;
         return false;
       });
  return RetF;
}





template<typename Telt>
void AddGeneratorsExtendSchreierTree(StabChain<Telt> & S, std::vector<Telt> const& newgens)
{
#define DEBUG_ADD_GEN_SCH
#ifdef DEBUG_ADD_GEN_SCH
  std::cerr << "CPP AGEST : Beginning of AddGeneratorsExtendSchreierTree\n";
  //  std::cerr << "CPP AGEST 1: genlabels=" << GapStringIntVector(S->genlabels) << "\n";
  StabChain<Telt> Swrite = S;
  int idxwrt=0;
  while(Swrite != nullptr) {
    //    std::cerr << "DEBUG idxwrt=" << idxwrt << " |labels|=" << S->comm->labels.size() << "\n";
    idxwrt++;
    Swrite = Swrite->stabilizer;
  }
#endif
  int nbLabel=S->comm->labels.size();
  std::vector<int> ListAtt(nbLabel);
  for (int i=0; i<nbLabel; i++)
    ListAtt[i]=i;
  Face old=BlistList(ListAtt, S->genlabels);
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
    int pos = PositionVect(S->comm->labels, gen);
    if (pos == -1) {
      S->comm->labels.push_back(gen);
      old.push_back(false);
      ald.push_back(true);
      int posG=S->comm->labels.size() - 1;
#ifdef DEBUG_ADD_GEN_SCH
      std::cerr << "CPP AGEST  genlabels insert 1:\n";
      //      std::cerr << "CPP AGEST  genlabels insert 1: pos=" << (posG+1) << "\n";
      //std::cerr << "CPP AGEST  genlabels insert X: pos=" << (posG+1) << "\n";
#endif
      S->genlabels.push_back(posG);
    } else {
      if (!ald[pos]) {
#ifdef DEBUG_ADD_GEN_SCH
        std::cerr << "CPP AGEST  genlabels insert 2:\n";
        //        std::cerr << "CPP AGEST  genlabels insert 2: pos=" << (pos+1) << "\n";
	//std::cerr << "CPP AGEST  genlabels insert X: pos=" << (pos+1) << "\n";
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
      std::cerr << "CPP   AGEST i=" << (i+1) << " |cycles|=" << S->cycles.size() << "\n";
#endif
      for (int& j : S->genlabels) {
	if (i >= len || old[j] == 0) {
	  int img=SlashAct(S->orbit[i], S->comm->labels[j]);
#ifdef DEBUG_ADD_GEN_SCH
	  std::cerr << "CPP     AGEST img=" << (img+1) << " g=" << S->comm->labels[j] << "\n";
          //	  std::cerr << "DEBUG |S->transversal|=" << S->transversal.size() << " img=" << img << "\n";
#endif
	  if (S->transversal[img] != -1) {
#ifdef DEBUG_ADD_GEN_SCH
	    std::cerr << "CPP       |S->cycles|=" << S->cycles.size() << " i=" << (i+1) << "\n";
#endif
            AssignationVectorGapStyle(S->cycles, i, int8_t(true));
            //            S->cycles[i] = true;
#ifdef DEBUG_ADD_GEN_SCH
	    std::cerr << "CPP       AGEST assign true\n";
#endif
	  } else {
	    S->transversal[img]=j;
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
      for (int& j : S->genlabels) {
	if (i >= len || old[j] == 0) {
	  int img=SlashAct(S->orbit[i], S->comm->labels[j]);
	  if (S->transversal[img] == -1) {
	    S->transversal[img]=j;
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



template<typename Telt>
void ChooseNextBasePoint(StabChain<Telt> & S, std::vector<int> const& base, std::vector<Telt> const& newgens)
{
  std::cerr << "CPP base = " << GapStringIntVector(base) << "\n";
  auto IsFullyStable=[&](int const& eBas) -> bool {
    for (auto & eGen : newgens) {
      if (PowAct(eBas, eGen) != eBas)
	return false;
    }
    return true;
  };
  int i = 0;
  int len=base.size();
  while(true) {
    if (i == len)
      break;
    int eBas=base[i];
    if (!IsFullyStable(eBas))
      break;
    i++;
  }
  int pnt;
  if (i < len)
    pnt=base[i];
  else
    pnt=SmallestMovedPoint(newgens);
  int bpt, pos;
  if (S->orbit.size() > 0) {
    bpt = S->orbit[0];
    pos = PositionVect(base, bpt);
  }
  else {
    bpt = S->comm->n + 444; // value in GAP is infinity
    pos = -1;
  }
  std::cerr << "CPP pnt=" << PosFail_to_string(pnt) << " bpt=" << ConstrainedIntInfinity_to_string(bpt, S->comm->n) << " pos=" << PosFail_to_string(pos) << "\n";
  if ((pos != -1 && i < pos) || (pos == -1 && i<int(base.size())) || (pos == -1 && pnt < bpt)) {
    std::cerr << "CPP matching test\n";
    //    std::cerr << "CPP Before InsertTrivialStabilizer S=" << S << "\n";
    InsertTrivialStabilizer(S, pnt);
    //    std::cerr << "CPP After InsertTrivialStabilizer S=" << S << "\n";
    if (S->stabilizer->IsBoundCycle) {
      std::vector<int8_t> eFace = {0};
      S->cycles = eFace;
      S->IsBoundCycle=true;
      std::cerr << "CPP   Initializing cycles\n";
    }
  }
  std::cerr << "CPP Exiting ChooseNextBasePoint\n";
}



template<typename Telt, typename Tint>
void StabChainStrong(StabChain<Telt> & S, std::vector<Telt> const& newgens, StabChainOptions<Tint> const& options)
{
  std::cerr << "CPP Begin StabChainStrong : newgens=" << GapStringTVector(newgens) << "\n";
  std::cerr << "CPP StabChainStrong 1: genlabels=" << GapStringIntVector(S->genlabels) << "\n";
  ChooseNextBasePoint(S, options.base, newgens);
  std::cerr << "CPP StabChainStrong 2: genlabels=" << GapStringIntVector(S->genlabels) << "\n";

  int pnt = S->orbit[0];
  int len = S->orbit.size();
  int old = S->genlabels.size();
  std::cerr << "CPP Before AddGeneratorsExtendSchreierTree\n";
  AddGeneratorsExtendSchreierTree(S, newgens);

  //# If a new generator fixes the base point, put it into the stabilizer.
  std::cerr << "CPP newgens=" << GapStringTVector(newgens) << "\n";
  for (auto & eGen : newgens) {
    std::cerr << "CPP eGen=" << eGen << " eGen=" << GapStyleString(eGen) << "\n";
    if (!eGen.isIdentity() && PowAct(pnt, eGen) == pnt) {
      std::cerr << "CPP   1: Calling StabChainStrong with eGen=" << GapStyleString(eGen) << "\n";
      StabChainStrong(S->stabilizer, {eGen}, options);
    }
  }

  // # Compute the Schreier generators (seems to work better backwards).
  std::vector<int> pnts = ClosedInterval(0, S->orbit.size());
  if (S->IsBoundCycle)
    pnts = ListBlist(pnts, S->cycles);
  std::cerr << "CPP pnts = " << GapStringIntVector(pnts) << "\n";
  std::cerr << "CPP Usecycle=" << S->IsBoundCycle << "\n";
  if (S->IsBoundCycle) {
    std::cerr << "CPP cycles=" << GapStringBoolVectorB(S->cycles) << "\n";
  }
  int gen1=0;
  std::cerr << "CPP StabChainStrong O=" << GapStringIntVector(S->orbit) << "\n";
  for (int& i : Reversed(pnts)) {
    int p=S->orbit[i];
    std::cerr << "CPP StabChainStrong i=" << (i+1) << " p=" << (p+1) << "\n";
    Telt rep=InverseRepresentative(S, p);
    if (i < len)
      gen1=old;
    std::cerr << "CPP StabChainStrong gen1=" << (gen1+1) << " rep=" << rep << "\n";
    for (int & j : ClosedInterval(gen1, S->genlabels.size())) {
      Telt g = S->comm->labels[ S->genlabels[j] ];
      std::cerr << "CPP StabChainStrong   j=" << (j+1) << " g=" << g << "\n";
      if (S->transversal[ SlashAct(p, g) ] != S->genlabels[j]) {
        Telt sch = SiftedPermutation(S, Inverse(g*rep));
	std::cerr << "CPP sch=" << sch << " g=" << g << " rep=" << rep << "\n";
	if (!sch.isIdentity()) {
	  StabChainStrong(S->stabilizer, {sch}, options );
	}
      }
    }
  }
  std::cerr << "CPP exiting StabChainStrong\n";
}


template<typename Telt>
bool StabChainForcePoint(StabChain<Telt> & Stot, int const& pnt)
{
  std::cerr << "CPP Beginning of StabChainForcePoint pnt=" << (pnt+1) << "\n";
  PrintStabChain(Stot);
  std::cerr << "DEBUG |Stot->transversal|=" << Stot->transversal.size() << "\n";
  if (Stot->transversal.size() == 0 || Stot->transversal[pnt] == -1) {
    std::cerr << "CPP Matching the first test\n";
    if (IsFixedStabilizer(Stot, pnt )) {
      std::cerr << "CPP Matching the second test\n";
      InsertTrivialStabilizer(Stot, pnt);
    }
    else {
      if (!StabChainForcePoint(Stot->stabilizer, pnt) || !StabChainSwap(Stot)) {
        std::cerr << "CPP StabChainForcePoint, return false\n";
	return false;
      }
    }
  }
  std::cerr << "CPP StabChainForcePoint, return true\n";
  PrintStabChain(Stot);
  return true;
}


template<typename Telt>
std::vector<Telt> GetListGenerators(StabChain<Telt> const& Stot)
{
  std::vector<Telt> LGens;
  for (auto & posGen : Stot->genlabels)
    LGens.push_back(Stot->comm->labels[posGen]);
  return LGens;
}




template<typename Telt>
bool StabChainSwap(StabChain<Telt> & Stot)
{
  std::cerr << "CPP Beginning of StabChainSwap\n";
  PrintStabChain(Stot);
  int n = Stot->comm->n;
  int a = Stot->orbit[0];
  int b = Stot->stabilizer->orbit[0];
  std::cerr << "CPP StabChainSwap a=" << (a+1) << " b=" << (b+1) << "\n";
  //
  std::vector<Telt> LGens = GetListGenerators(Stot);
  // We have some missing entries in the S.generators.
  // Apparently, the S.generators contains generators that are not used.
  std::cerr << "CPP LGens=" << GapStringTVector(LGens) << "\n";
  //
  //  StabChain<Telt> Ttot = EmptyStabChainPlusNode<Telt>(n, b);
  StabChain<Telt> Ttot = EmptyStabChainPlusCommonPlusNode(Stot->comm, b);
  AddGeneratorsExtendSchreierTree(Ttot, LGens);
  std::cerr << "CPP StabChainSwap : after first AGEST\n";
  //
  StabChain<Telt> Tstab = EmptyStabChainPlusNode<Telt>(n, a);
  if (Stot->stabilizer != nullptr) {
    if (Stot->stabilizer->stabilizer != nullptr) {
      std::vector<Telt> LGensB = GetListGenerators(Stot->stabilizer->stabilizer);
      std::cerr << "CPP StabChainSwap : before second AGEST gens=" << GapStringTVector(LGensB) << "\n";
      AddGeneratorsExtendSchreierTree(Tstab, LGensB);
      std::cerr << "CPP StabChainSwap : after second AGEST\n";
    }
  }
  std::cerr << "CPP After AddGeneratorsExtendSchreierTree 1 : Tstab=\n";
  PrintStabChain(Tstab);
  //
  size_t ind = 0;
  size_t len = Stot->orbit.size() * Stot->stabilizer->orbit.size() / Ttot->orbit.size();
  std::cerr << "CPP StabChainSwap |Tstab->orbit|=" << int(Tstab->orbit.size()) << " len=" << int(len) << "\n";
  while (Tstab->orbit.size() < len) {
    std::cerr << "CPP beginning of loop\n";
    int pnt;
    while(true) {
      ind++;
      if (ind >= Stot->orbit.size())
	return false;
      std::cerr << "CPP |orbit|=" << Stot->orbit.size() << " ind=" << (ind+1) << "\n";
      pnt = Stot->orbit[ind];
      std::cerr << "CPP ind=" << (ind+1) << " pnt=" << (pnt+1) << "\n";
      std::cerr << "DEBUG |Tstab->transversal|=" << Tstab->transversal.size() << " pnt=" << pnt << "\n";
      if (Tstab->transversal[pnt] == -1)
	break;
    }
    std::cerr << "CPP ind=" << (ind+1) << " pnt=" << (pnt+1) << "\n";
    int img = b;
    int i = pnt;
    while (i != a) {
      int posGen=Stot->transversal[i];
      img = PowAct(img, Stot->comm->labels[posGen]);
      i = PowAct(i, Stot->comm->labels[posGen]);
    }
    std::cerr << "CPP i=" << (i+1) << " img=" << (img+1) << "\n";
    if (Stot->stabilizer->transversal[img] != -1) {
      std::cerr << "CPP Determining gen, step 1\n";
      Telt gen = Stot->comm->identity;
      while (PowAct(pnt, gen) != a) {
        std::cerr << "CPP pnt^gen=" << (PowAct(pnt, gen)+1) << "\n";
	int posGen=Stot->transversal[PowAct(pnt, gen)];
        std::cerr << "DEBUG posGen=" << posGen << "\n";
	gen = gen * Stot->comm->labels[posGen];
      }
      std::cerr << "CPP Determining gen, step 2 gen=" << gen << "\n";
      while (PowAct(b, gen) != b) {
        std::cerr << "CPP b^gen=" << (PowAct(b, gen)+1) << "\n";
	int posGen=Stot->stabilizer->transversal[PowAct(b, gen)];
        std::cerr << "DEBUG posGen=" << posGen << "\n";
	gen = gen * Stot->stabilizer->comm->labels[posGen];
      }
      std::cerr << "CPP Determining gen, step 3 gen=" << gen << "\n";
      AddGeneratorsExtendSchreierTree(Tstab, {gen});
      std::cerr << "CPP After AddGeneratorsExtendSchreierTree 2 : Tstab=\n";
      PrintStabChain(Tstab);
    }
  }
  std::cerr << "CPP After while loop\n";
  std::cerr << "CPP T=\n";
  PrintStabChain(Ttot);
  std::cerr << "CPP Tstab=\n";
  PrintStabChain(Tstab);
  auto MapAtLevel=[&](StabChain<Telt> & Swork, StabChain<Telt> const& insStab) -> void {
    Swork->genlabels = insStab->genlabels;
    Swork->comm = insStab->comm;
    Swork->orbit = insStab->orbit;
    Swork->transversal = insStab->transversal;
  };
  MapAtLevel(Stot, Ttot);
  PrintStabChain(Stot);
  std::cerr << "CPP StabChainSwap 1: |orbit|=" << Tstab->orbit.size() << "\n";
  if (Tstab->orbit.size() == 1) {
    Stot->stabilizer = Stot->stabilizer->stabilizer;
    std::cerr << "CPP StabChainSwap 2:\n";
  }  else {
    MapAtLevel(Stot->stabilizer, Tstab);
    std::cerr << "CPP StabChainSwap 3:\n";
  }
  PrintStabChain(Stot);
  return true;
}





// maybe use std::map<T, T> instead
template<typename T>
T LabsLims(T const& lab, std::function<T(T const&)> const& hom, std::vector<T> & labs, std::vector<T> & lims)
{
  int pos=PositionVect(labs, lab);
  if (pos == -1) {
    int pos=labs.size();
    labs.push_back(lab);
    T img=hom(lab);
    lims.push_back(img);
  }
  return lims[pos];
}



// The original ConjugateStabChain code
// It is simplified from the original one with us being in the case
// IsPerm(hom) and IsPerm(map).
template<typename Telt, typename Fhom, typename Fmap, typename Fcond>
StabChain<Telt> ConjugateStabChain(StabChain<Telt> & Stot, StabChain<Telt> & Ttot, Fhom const& hom, Fmap const& map, Fcond const& cond)
{
  int n=Stot->comm->n;
  Telt id=Stot->comm->identity;
  StabChain<Telt> Sptr = Stot;
  StabChain<Telt> Tptr = Ttot;

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
      std::vector<int> NewTransversal(n);
      for (int i=0; i<n; i++) {
	int iImg=map(i);
	int eVal=Sptr->transversal[i];
	NewTransversal[iImg] = eVal;
      }
      Sptr->transversal=NewTransversal;
    }
    Tptr->comm = get_comm(Sptr->comm);
    Tptr->genlabels = Sptr->genlabels;
    Tptr->orbit = ListT(Sptr->orbit, map);

    // Going to the next level.
    Sptr = Sptr->stabilizer;
    if (Tptr->stabilizer == nullptr)
      Tptr->stabilizer = EmptyStabChainPlusCommon<Telt>(Tptr->comm);
    Tptr = Tptr->stabilizer;
  }
  //
  // Mapping the labels that showed up.
  //
  return Tptr;
}



template<typename Telt>
void ConjugateStabChain_Element(StabChain<Telt> & Stot, Telt const& cnj)
{
  Telt cnjInv=~cnj;
  auto hom=[&](Telt const& x) -> Telt {
    Telt retV = cnjInv*x*cnj;
    return retV;
  };
  auto map=[&](int const& x) -> int {
    return cnj.at(x);
  };
  auto cond=[](StabChain<Telt> const& S) -> bool {
    if (S->stabilizer == nullptr)
      return false;
    return true;
  };
  (void)ConjugateStabChain(Stot, Stot, hom, map, cond);
}


template<typename Telt>
std::string PrintTopOrbit(StabChain<Telt> const& S)
{
  if (S == nullptr) {
    return "unset";
  }
  int len=S->orbit.size();
  std::string str = "[ ";
  for (int u=0; u<len; u++) {
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
template<typename Telt>
bool ChangeStabChain(StabChain<Telt> & Gptr, std::vector<int> const& base, int const& reduced)
{
#define DEBUG_CHANGE_STAB_CHAIN
#ifdef DEBUG_CHANGE_STAB_CHAIN
  std::string strG_orig=GetStringExpressionOfStabChain(Gptr);
  std::string strG_current=strG_orig;
  std::cerr << "CPP Beginning ChangeStabChain, GetStabilizerDepth = " << GetStabilizerDepth(Gptr) << "\n";
#endif
  Telt cnj = Gptr->comm->identity;
  StabChain<Telt> Sptr = Gptr;
#ifdef DEBUG_CHANGE_STAB_CHAIN
  std::string strS_current=GetStringExpressionOfStabChain(Sptr);
  std::cerr << "CPP we have strS_current\n";

  int idx=0;
  auto KeyUpdating=[&](std::string const& str) {
    bool DoPrint;
    idx++;
    std::string strGloc=GetStringExpressionOfStabChain(Gptr);
    std::string strSloc=GetStringExpressionOfStabChain(Sptr);
    DoPrint=false;
    if (DoPrint) {
      std::cerr << "CPP KU At " << idx << " of " << str << " dep(G)/dep(S)=" << GetStabilizerDepth(Gptr) << "/" << GetStabilizerDepth(Sptr) << "\n";
      std::cerr << "CPP KU At step " << idx << " of " << str << "\n";
      std::cerr << "CPP sgs(G) = " << GapStringTVector(SortVector(StrongGeneratorsStabChain(Gptr))) << "\n";
      std::cerr << "CPP sgs(S) = " << GapStringTVector(SortVector(StrongGeneratorsStabChain(Sptr))) << "\n";
      if (strG_current == strGloc) {
        std::cerr << "CPP   KU At step " << idx << " of " << str << " no change of G\n";
      }
      else {
        std::cerr << "CPP   KU At step " << idx << " of " << str << " CHANGE of G\n";
        strG_current=strGloc;
      }
      if (strS_current == strSloc) {
        std::cerr << "CPP   KU At step " << idx << " of " << str << " no change of S\n";
      }
      else {
        std::cerr << "CPP   KU At step " << idx << " of " << str << " CHANGE of S\n";
        strS_current=strSloc;
      }
    }
    std::cerr << "CPP Gptr at str=" << str << "\n";
    PrintStabChain(Gptr);
    std::cerr << "CPP Sptr at str=" << str << "\n";
    PrintStabChain(Sptr);
  };
#endif
  std::vector<int> newBase;
  size_t i=0;
  size_t basSiz=base.size();
#ifdef DEBUG_CHANGE_STAB_CHAIN
  std::cerr << "CPP ChangeStabChain base = " << GapStringIntVector(base) << "\n";
  std::cerr << "CPP ChangeStabChain 1 orbit=" << PrintTopOrbit(Gptr) << "\n";
#endif
  while (GetStabilizerDepth(Sptr) > 1 || i < basSiz) {
#ifdef DEBUG_CHANGE_STAB_CHAIN
    std::cerr << "CPP GetStabilizerDepth(S)=" << GetStabilizerDepth(Sptr) << " GetStabilizerDepth(G)=" << GetStabilizerDepth(Gptr) << "\n";
#endif
    int old=BasePoint(Sptr);
#ifdef DEBUG_CHANGE_STAB_CHAIN
    std::cerr << "CPP ChangeStabChain old=" << PosFalse_to_string(old) << " i=" << (i+1) << " |base|=" << basSiz << "\n";
    KeyUpdating("After BasePoint");
#endif
    //    std::cerr << "eLev=" << eLev << "\n";
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
      int newpnt = SlashAct(base[i], cnj);
#ifdef DEBUG_CHANGE_STAB_CHAIN
      std::cerr << "CPP newpnt=" << (newpnt+1) << "\n";
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
#ifdef DEBUG_PERMUTALIB
	  else {
	    std::cerr << "CPP <base> must be an extension of base of <G>\n";
	    throw TerminalException{1};
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
    } else if (PositionVect(newBase, old) != -1 || (reduced == int_true && Sptr->orbit.size() == 1)) {
#ifdef DEBUG_CHANGE_STAB_CHAIN
      std::cerr << "CPP Stabilizer shift in ChangeStabChain\n";
#endif

      Sptr->genlabels = Sptr->stabilizer->genlabels;
      if (Sptr->stabilizer->orbit.size() > 0) {
        Sptr->orbit = Sptr->stabilizer->orbit;
        Sptr->transversal = Sptr->stabilizer->transversal;
      } else {
        Sptr->orbit.clear();
        Sptr->transversal.clear();
      }
      if (Sptr->stabilizer->stabilizer != nullptr) {
        Sptr->stabilizer = Sptr->stabilizer->stabilizer;
      } else {
        Sptr->stabilizer = nullptr;
      }
    } else {
      newBase.push_back(old);
      Sptr = Sptr->stabilizer;
#ifdef DEBUG_CHANGE_STAB_CHAIN
      KeyUpdating("After S:=S.stabilizer 3");
#endif
    }
  }
#ifdef DEBUG_CHANGE_STAB_CHAIN
  std::cerr << "CPP LEAVE i=" << (i+1) << " |base|=" << basSiz << "\n";
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

template<typename Telt>
bool ExtendStabChain(StabChain<Telt> & Stot, std::vector<int> const& base)
{
  std::cerr << "CPP Beginning of ExtendStabChain\n";
  return ChangeStabChain(Stot, base, int_reducedm1);
}


template<typename Telt>
bool ReduceStabChain(StabChain<Telt> & Stot)
{
  return ChangeStabChain(Stot, {}, int_true);
}




template<typename Telt>
bool TestEqualityStabChain(StabChain<Telt> const& L, StabChain<Telt> const& R)
{
  StabChain<Telt> Lptr = L;
  StabChain<Telt> Rptr = R;
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
    int lenL=Lptr->transversal.size();
    for (int u=0; u<lenL; u++) {
      if (Lptr->transversal[u] == -1 && Rptr->transversal[u] != -1) {
#ifdef DEBUG_EQUALITY
        std::cerr << "TestEqualityStabChain 5\n";
#endif
	return false;
      }
      if (Lptr->transversal[u] != -1 && Rptr->transversal[u] == -1) {
#ifdef DEBUG_EQUALITY
        std::cerr << "TestEqualityStabChain 6\n";
#endif
	return false;
      }
      if (Lptr->transversal[u] != -1) {
	int idxL=Lptr->transversal[u];
	int idxR=Rptr->transversal[u];
	Telt permL=Lptr->comm->labels[idxL];
	Telt permR=Rptr->comm->labels[idxR];
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


template<typename Telt>
void SetStabChainFromLevel(std::vector<StabChain<Telt>> & R_list, std::vector<StabChain<Telt>> const& L_list,
                           int const& levbegin, int const& levend)
{
  for (int iLev=levbegin; iLev<levend; iLev++) {
    R_list[iLev] = StructuralCopy(L_list[iLev]);
  }
}




template<typename Telt>
bool IsFixedStabilizer(StabChain<Telt> const& S, int const& pnt)
{
  for (auto & posGen : S->genlabels) {
    if (pnt != PowAct(pnt, S->comm->labels[posGen]))
      return false;
  }
  return true;
}


template<typename Telt>
Telt MinimalElementCosetStabChain(StabChain<Telt> const& Stot, Telt const& g)
{
  Telt gRet=g;
  StabChain<Telt> Sptr = Stot;
  while(true) {
    if (Sptr == nullptr)
      break;
    if (Sptr->genlabels.size() == 0)
      return gRet;
    int pMin=Sptr->comm->n + 1;
    for (auto & i : Sptr->orbit) {
      int a=PowAct(i,gRet);
      if (a < pMin)
	pMin=a;
    }
    int bp=Sptr->orbit[0];
    int pp=SlashAct(pMin, gRet);
    while (bp != pp) {
      int pos=Sptr->transversal[pp];
      gRet=LeftQuotient(Sptr->comm->labels[pos], gRet);
      pp = SlashAct(pMin, gRet);
    }
    Sptr = Sptr->stabilizer;
  }
  return gRet;
}






template<typename Telt, typename Tret>
StabChain<Tret> HomomorphismMapping(StabChain<Telt> const& Stot, std::function<Tret(Telt const&)> const& f)
{
  Tret idMap = f(Stot->comm->identity);
  int nMap=idMap.size();
  auto fVector =[&](std::vector<Telt> const& V) -> std::vector<Tret> {
    std::vector<Tret> Vret;
    for (auto & eElt : V)
      Vret.push_back(f(eElt));
    return Vret;
  };
  std::vector<Tret> labelsMap = fVector(Stot->comm->labels);
  std::shared_ptr<CommonStabInfo<Tret>> comm = std::make_shared<CommonStabInfo<Tret>>(CommonStabInfo<Tret>({nMap, idMap, labelsMap}));

  StabChain<Telt> Sptr = Stot;
  StabChain<Tret> Swork = nullptr;
  StabChain<Tret> Sreturn = nullptr;

  while(true) {
    if (Sptr == nullptr)
      break;
    Swork = std::make_shared<StabLevel<Telt>>(StabLevel<Telt>({Sptr->transversal, Sptr->orbit, Sptr->genlabels, Sptr->cycles, Sptr->IsBoundCycle, fVector(Sptr->treegen), fVector(Sptr->treegeninv), fVector(Sptr->aux), Sptr->treedepth, Sptr->diam, comm, nullptr}));
    if (Sreturn == nullptr)
      Sreturn = Swork;
    Sptr = Sptr->stabilizer;
    Swork = Swork->stabilizer;
  }
  return Sreturn;
}




}


#endif
