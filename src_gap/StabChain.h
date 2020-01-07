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
  ---StabChainForcePoint seems not to be used
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

static const int int_reducedm1 = -1;
static const int int_false = 0;
static const int int_true = 1;
static const int int_fail = 2;
static const int int_int  = 3;
static const int int_perm = 4;
static const int int_group = 5;
static const int int_stablev = 6;


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
  bool UseCycle;
  std::vector<Telt> labels;
};

template<typename Telt>
struct StabLevel {
  std::vector<int> transversal;
  std::vector<int> orbit;
  std::vector<int> genlabels; // Not used in Random algorithm
  Face cycles;
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
void PrintStabChainTransversals(StabChain<Telt> const& S)
{
  StabChain<Telt> Swork = S;
  int n=Swork->comm->n;
  int iLevel=0;
  while(Swork != nullptr) {
    std::vector<std::optional<Telt>> V(n);
    for (int i=0; i<n; i++) {
      int eVal = Swork->transversal[i];
      if (eVal == -1)
        V[i] = {};
      else
        V[i] = Swork->comm->labels[eVal];
    }
    //
    std::cerr << "CPP i=" << iLevel << " " << GapStringMissingTVector(V) << "\n";
    Swork = Swork->stabilizer;
    iLevel++;
  }
}


/*
template<typename Telt>
struct StabChain {
  int n;
  Telt identity;
  bool UseCycle;
  std::vector<Telt> labels;
  std::vector<StabLevel<Telt>> stabilizer;
};
*/

template<typename Telt>
int GetStabilizerDepth(StabChain<Telt> const& Sptr)
{
  StabChain<Telt> Ssec = Sptr;
  int dep = 0;
  while(true) {
    if (Ssec == nullptr)
      break;
    dep++;
    Ssec = Ssec->stabilizer;
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
    /*
    os << "  genlabels =";
    for (auto & eVal : Sptr->genlabels)
      os << " " << eVal;
    os << "\n";
    */
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
  std::shared_ptr<CommonStabInfo<Telt>> comm_new = std::make_shared<CommonStabInfo<Telt>>(*(S->comm));
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
    S2->comm = comm_new;
    Sret = S2;
  }
  std::string str1 = GetStringExpressionOfStabChain(S);
  std::string str2 = GetStringExpressionOfStabChain(Sret);
  if (str1 != str2) {
    std::cerr << "We fail to have equality of entries\n";
    std::cerr << "str1=" << str1 << "\n";
    std::cerr << "str2=" << str2 << "\n";
    throw TerminalException{1};
  }
  return Sret;
}



template<typename Telt>
StabChain<Telt> RestrictedStabChain(StabChain<Telt> const& Stot, int const& eLev)
{
  int nbLevel=Stot.stabilizer.size();
  std::vector<StabLevel<Telt>> stabilizerRed;
  for (int uLev=eLev; uLev<nbLevel; uLev++)
    stabilizerRed.push_back(Stot.stabilizer[uLev]);
  return {Stot.n, Stot.identity, Stot.UseCycle, Stot.labels, stabilizerRed};
}


template<typename Telt>
StabLevel<Telt> EmptyStabLevel(std::shared_ptr<CommonStabInfo<Telt>> const& comm)
{
  std::vector<int> transversal = std::vector<int>(comm->n, -1);
  std::vector<int> orbit;
  std::vector<int> genlabels;
  Face cycles;
  std::vector<Telt> treegen;
  std::vector<Telt> treegeninv;
  std::vector<Telt> aux;
  int treedepth = 0;
  int diam = 0;
  return {transversal, orbit, genlabels, cycles, treegen, treegeninv, aux, treedepth, diam, comm, nullptr};
}



template<typename Telt>
StabChain<Telt> EmptyStabChain(int const& n)
{
  Telt id(n);
  std::shared_ptr<CommonStabInfo<Telt>> comm = std::make_shared<CommonStabInfo<Telt>>(CommonStabInfo<Telt>({n, id, false, {id}}));
  return std::make_shared<StabLevel<Telt>>(EmptyStabLevel<Telt>(comm));
}


template<typename Telt>
StabChain<Telt> EmptyStabChainPlusNode(int const& n, int const& bas)
{
  StabChain<Telt> S = EmptyStabChain<Telt>(n);
  InitializeSchreierTree(S, bas);
  //  S->orbit.push_back(bas);
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
  Stot->stabilizer == nullptr;
  Stot->genlabels.clear();
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
  std::set<int> sgs_set;
  StabChain<Telt> Sptr = S;
  while(true) {
    if (Sptr == nullptr)
      break;
    int siz = Sptr->genlabels.size();
    if (siz == 0)
      break;
    for (auto & pos : Sptr->genlabels)
      sgs_set.insert(pos);
    Sptr = Sptr->stabilizer;
  }
  std::vector<Telt> sgs;
  for (auto & ePos : sgs_set)
    sgs.push_back(S->comm->labels[ePos]);
  return sgs;
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
  std::set<int> LIdx;
  StabChain<Telt> Sptr = G;
  while(true) {
    if (Sptr == nullptr)
      break;
    for (auto & eIdx : Sptr->genlabels)
      LIdx.insert(eIdx);
    Sptr = Sptr->stabilizer;
  }
  for (auto & eIdx : LIdx)
    if (!G->comm->labels[eIdx].isIdentity())
      return false;
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
  Stot->stabilizer = ShallowCopy(Stot);
  Stot->transversal = {};
  Stot->orbit = {};
  Stot->genlabels = Stot->stabilizer->genlabels;
  Stot->cycles = Face(Stot->comm->n);
  Stot->treegen = {};
  Stot->treegen = {};
  Stot->aux = {};
  Stot->treedepth = -1;
  Stot->diam = -1;
  InitializeSchreierTree(Stot, pnt);
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

std::string GapStringIntVector(std::vector<int> const& f)
{
  std::string str;
  str += "[ ";
  int len=f.size();
  for (int i=0; i<len; i++) {
    if (i>0)
      str += ", ";
    str += std::to_string(f[i]+1);
  }
  str += " ]";
  return str;
}

std::string GapStringBoolVector(Face const& f)
{
  std::string str;
  str += "[ ";
  int len=f.size();
  for (int i=0; i<len; i++) {
    if (i>0)
      str += ", ";
    if (f[i])
      str += "true";
    else
      str += "false";
  }
  str += " ]";
  return str;
}

template<typename T>
std::string GapStringTVector(std::vector<T> const& f)
{
  std::ostringstream os;
  os << "[ ";
  int len=f.size();
  for (int i=0; i<len; i++) {
    if (i>0)
      os << ", ";
    os << f[i];
  }
  os << " ]";
  return os.str();
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
  std::cerr << "CPP AGEST 1: genlabels=" << GapStringIntVector(S->genlabels) << "\n";
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
  std::cerr << "CPP AGEST 1: old=" << GapStringBoolVector(old) << "\n";
  std::cerr << "CPP AGEST 1: ald=" << GapStringBoolVector(ald) << "\n";
  std::cerr << "CPP AGEST labels=" << GapStringTVector(S->comm->labels) << "\n";
  std::cerr << "CPP AGEST 2: genlabels=" << GapStringIntVector(S->genlabels) << "\n";
#endif
  for (auto & gen : newgens) {
    int pos = PositionVect(S->comm->labels, gen);
    if (pos == -1) {
      S->comm->labels.push_back(gen);
      old.push_back(false);
      ald.push_back(true);
      int posG=S->comm->labels.size() - 1;
#ifdef DEBUG_ADD_GEN_SCH
      std::cerr << "CPP AGEST  genlabels insert 1: pos=" << (posG+1) << "\n";
#endif
      S->genlabels.push_back(posG);
    }
    else {
      if (!ald[pos]) {
#ifdef DEBUG_ADD_GEN_SCH
	std::cerr << "CPP AGEST  genlabels insert 2: pos=" << (pos+1) << "\n";
#endif
	S->genlabels.push_back(pos);
      }
    }
  }
#ifdef DEBUG_ADD_GEN_SCH
  std::cerr << "CPP AGEST 2: old=" << GapStringBoolVector(old) << "\n";
  std::cerr << "CPP AGEST 2: ald=" << GapStringBoolVector(ald) << "\n";
#endif

  int len = S->orbit.size();
  int i=0;
  if (S->comm->UseCycle) {
#ifdef DEBUG_ADD_GEN_SCH
    std::cerr << "CPP AGEST Cycles len=" << len << "\n";
#endif
    while (i < int(S->orbit.size())) {
#ifdef DEBUG_ADD_GEN_SCH
      std::cerr << "CPP   AGEST i=" << (i +1) << "\n";
#endif
      for (int& j : S->genlabels) {
	if (i > len-1 || old[j] == 0) {
	  int img=SlashAct(S->orbit[i], S->comm->labels[j]);
#ifdef DEBUG_ADD_GEN_SCH
	  std::cerr << "CPP     AGEST img=" << (img +1) << " g=" << S->comm->labels[j] << "\n";
#endif
	  if (S->transversal[img] != -1) {
	    S->cycles[i]=true;
#ifdef DEBUG_ADD_GEN_SCH
	    std::cerr << "CPP       AGEST assign true\n";
#endif
	  }
	  else {
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
  }
  else {
#ifdef DEBUG_ADD_GEN_SCH
    std::cerr << "CPP AGEST No Cycles len=" << len << "\n";
#endif
    while (i < int(S->orbit.size())) {
      for (int& j : S->genlabels) {
	if (i > len || old[j] == 0) {
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
  std::cerr << "CPP BPT/POS bpt=" << (bpt+1) << " pos=" << pos << "\n";
  if ((pos != -1 && i < pos) || (pos == -1 && i<int(base.size())) || (pos == -1 && pnt < bpt)) {
    std::cerr << "CPP pnt=" << (pnt+1) << " bpt=" << (bpt+1) << " pos=" << pos << "\n";
    //    std::cerr << "CPP Before InsertTrivialStabilizer S=" << S << "\n";
    InsertTrivialStabilizer(S, pnt);
    //    std::cerr << "CPP After InsertTrivialStabilizer S=" << S << "\n";
    if (S->comm->UseCycle) {
      Face eFace(1);
      eFace[0] = 0;
      S->cycles = eFace;
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
  for (auto & eGen : newgens)
    if (eGen.isIdentity() == false && PowAct(pnt, eGen) == pnt) {
      std::cerr << "CPP   1: Calling StabChainStrong with eGen=" << GapStyleString(eGen) << "\n";
      StabChainStrong(S->stabilizer, {eGen}, options);
    }

  // # Compute the Schreier generators (seems to work better backwards).
  std::vector<int> pnts = ClosedInterval(0, S->orbit.size());
  if (S->comm->UseCycle)
    pnts = ListBlist(pnts, S->cycles);
  std::cerr << "CPP pnts = " << GapStringIntVector(pnts) << "\n";
  std::cerr << "CPP Usecycle=" << S->comm->UseCycle << "\n";
  if (S->comm->UseCycle) {
    std::cerr << "CPP cycles=" << GapStringBoolVector(S->cycles) << "\n";
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
}


template<typename Telt>
bool StabChainForcePoint(StabChain<Telt> & Stot, int const& pnt)
{
  std::cerr << "CPP Beginning of StabChainForcePoint\n";
  if (Stot->transversal[pnt] == -1) {
    std::cerr << "CPP Matching the first test\n";
    if (IsFixedStabilizer(Stot, pnt )) {
      std::cerr << "CPP Matching the second test\n";
      InsertTrivialStabilizer(Stot, pnt);
    }
    else {
      if (!StabChainForcePoint(Stot->stabilizer, pnt) || !StabChainSwap(Stot))
	return false;
    }
  }
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
  int n=Stot->comm->n;
  int a = Stot->orbit[0];
  int b = Stot->stabilizer->orbit[0];
  //
  std::vector<Telt> LGens = GetListGenerators(Stot);
  //
  StabChain<Telt> Ttot = EmptyStabChainPlusNode<Telt>(n, b);
  AddGeneratorsExtendSchreierTree(Ttot, LGens);
  //
  StabChain<Telt> Tstab = EmptyStabChainPlusNode<Telt>(n, a);
  if (Tstab->stabilizer != nullptr) {
    if (Tstab->stabilizer->stabilizer != nullptr) {
      std::vector<Telt> LGensB = GetListGenerators(Stot->stabilizer->stabilizer);
      AddGeneratorsExtendSchreierTree(Tstab, LGensB);
    }
  }
  //
  int ind = 0;
  int len = Stot->orbit.size() * Stot->stabilizer->orbit.size() / Ttot->orbit.size();
  while (int(Tstab->orbit.size()) < len) {
    int pnt;
    while(true) {
      ind++;
      if (ind > int(Stot->orbit.size()))
	return false;
      pnt = Stot->orbit[ind];
      if (Tstab->transversal[pnt] == -1)
	break;
    }
    int img = b;
    int i = pnt;
    while (i != a) {
      int posGen=Stot->transversal[i];
      img = PowAct(img, Stot->comm->labels[posGen]);
      i = PowAct(i, Stot->comm->labels[posGen]);
    }
    if (Stot->stabilizer->transversal[img] != -1) {
      Telt gen = Stot->comm->identity;
      while (PowAct(pnt, gen) != a) {
	int posGen=Stot->transversal[PowAct(pnt, gen)];
	gen = gen * Stot->comm->labels[posGen];
      }
      while (PowAct(b, gen) != b) {
	int posGen=Stot->stabilizer->transversal[PowAct(pnt, gen)];
	gen = gen * Stot->comm->labels[posGen];
      }
      AddGeneratorsExtendSchreierTree(Tstab, {gen});
    }
  }
  auto MappingIndex=[&](StabChain<Telt> const& Wtot, int const& idx) -> int {
    if (idx == -1)
      return -1;
    Telt eElt=Wtot->comm->labels[idx];
    return GetLabelIndex(Stot->comm->labels, eElt);
  };
  auto MapAtLevel=[&](StabChain<Telt> const& Wtot) -> void {
    Stot->genlabels.clear();
    for (int const& posGen : Wtot->genlabels) {
      int posGenMap=MappingIndex(Wtot, posGen);
      Stot->genlabels.push_back(posGenMap);
    }
    Stot->orbit = Wtot->orbit;
    for (int u=0; u<n; u++) {
      int idx=Wtot->transversal[u];
      int idxMap=MappingIndex(Wtot, idx);
      Stot->transversal[u] = idxMap;
    }
  };
  MapAtLevel(Ttot);
  if (Tstab->orbit.size() == 1)
    Stot->stabilizer = Stot->stabilizer->stabilizer;
  else
    MapAtLevel(Tstab->stabilizer);
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

// We significantly change the functionality
// The function ConjugateStabChain is used only for the conjugation in the part of the code
// that interest us.
// The action that correspond to the code as written are
// OnTuples(orbit, map) with map=cnj. This is the normal action on points.
// OnTuples(labels, hom) with hom=cnj. This is the action by conjugacy:
// The action are therefore of the
template<typename Telt>
void ConjugateStabChain(StabChain<Telt> & Stot, Telt const& cnj)
{
#undef DEBUG_CONJ_STAB_CHAIN
  int n=Stot->comm->n;
  Telt cnjInv=~cnj;
  auto hom=[&](Telt const& x) -> Telt {
#ifdef DEBUG_CONJ_STAB_CHAIN
    std::cerr << "hom: cnj=" << cnj << " x=" << x << " cnjInv=" << cnjInv << "\n";
#endif
    Telt retV = cnjInv*x*cnj;
#ifdef DEBUG_CONJ_STAB_CHAIN
    std::cerr << "retV=" << retV << "\n";
#endif
    return retV;
  };
  auto map=[&](int const& x) -> int {
    return cnj.at(x);
  };
  StabChain<Telt> Sptr = Stot;
#ifdef DEBUG_CONJ_STAB_CHAIN
  std::cerr << "ConjugateStabChain, step 1\n";
#endif
  std::unordered_map<int,int> MappedTrans;
  while(true) {
#ifdef DEBUG_CONJ_STAB_CHAIN
    std::cerr << "ConjugateStabChain, step 2\n";
#endif
    if (Sptr == nullptr)
      break;
#ifdef DEBUG_CONJ_STAB_CHAIN
    std::cerr << "ConjugateStabChain, step 3\n";
#endif
    if (Sptr->transversal.size() > 0) {
      std::vector<int> NewTransversal(n);
      //      std::cerr << "Before loop n=" << n << "\n";
      for (int i=0; i<n; i++) {
	int iImg=cnj.at(i);
	//	std::cerr << "i=" << i << " iImg=" << iImg << "\n";
	int eVal=Sptr->transversal[i];
	//	std::cerr << "eVal=" << eVal << "\n";
	NewTransversal[iImg] = eVal;
        if (eVal != -1)
          MappedTrans[eVal] = -1;
	//	std::cerr << "After assignation\n";
      }
      //      std::cerr << " After loop\n";
      Sptr->transversal=NewTransversal;
    }
    
#ifdef DEBUG_CONJ_STAB_CHAIN
    std::cerr << "ConjugateStabChain, step 4\n";
#endif
    Sptr->treegen = ListT(Sptr->treegen, hom);
    Sptr->treegeninv = ListT(Sptr->treegeninv, hom);
    Sptr->aux = ListT(Sptr->aux, hom);
    Sptr->orbit = ListT(Sptr->orbit, map);
    Sptr = Sptr->stabilizer;
#ifdef DEBUG_CONJ_STAB_CHAIN
    std::cerr << "ConjugateStabChain, step 9\n";
#endif
  }
#ifdef DEBUG_CONJ_STAB_CHAIN
  std::cerr << "ConjugateStabChain, step 10\n";
#endif
  //
  // Mapping the labels that showed up.
  //
  for (auto & ePair : MappedTrans) {
    Telt img = hom(Stot->comm->labels[ePair.first]);
    int pos = PositionVect(Stot->comm->labels, img);
    if (pos == -1) {
      pos = Stot->comm->labels.size();
      Stot->comm->labels.push_back(img);
    }
    ePair.second = pos;
  }
  //
  // Now remapping the labels that occurred.
  //
  Sptr = Stot;
  while(true) {
    if (Sptr == nullptr)
      break;
    for (size_t i=0; i<Sptr->transversal.size(); i++) {
      int eVal = Sptr->transversal[i];
      if (eVal != -1)
        Sptr->transversal[i] = MappedTrans[eVal];
    }
    for (size_t i=0; i<Sptr->genlabels.size(); i++) {
      int eVal = Sptr->genlabels[i];
      Sptr->genlabels[i] = MappedTrans[eVal];
    }
    Sptr = Sptr->stabilizer;
  }
#ifdef DEBUG_CONJ_STAB_CHAIN
  std::cerr << "ConjugateStabChain, step 11\n";
#endif
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


std::string ConvertIntFalse(int const& pos)
{
  if (pos == -1)
    return "false";
  return std::to_string(pos+1);
}



// value of reduced
//  reduced = -1 corresponds to reduced = -1 in GAP code
//  reduced = 0  corresponds to reduced = false in GAP code
//  reduced = 1  corresponds to reduced = true in GAP code
template<typename Telt>
bool ChangeStabChain(StabChain<Telt> & Gptr, std::vector<int> const& base, int const& reduced)
{
  std::string strG_orig=GetStringExpressionOfStabChain(Gptr);
  std::string strG_current=strG_orig;
#define DEBUG_CHANGE_STAB_CHAIN
#ifdef DEBUG_CHANGE_STAB_CHAIN
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
  };
#endif
  std::vector<int> newBase;
  int i=0;
  int basSiz=base.size();
#ifdef DEBUG_CHANGE_STAB_CHAIN
  std::cerr << "CPP ChangeStabChain base = " << GapStringIntVector(base) << "\n";
  std::cerr << "CPP ChangeStabChain 1 orbit=" << PrintTopOrbit(Gptr) << "\n";
#endif
  while (GetStabilizerDepth(Sptr) > 1 || i<basSiz) {
#ifdef DEBUG_CHANGE_STAB_CHAIN
    std::cerr << "CPP GetStabilizerDepth(S)=" << GetStabilizerDepth(Sptr) << " GetStabilizerDepth(G)=" << GetStabilizerDepth(Gptr) << "\n";
#endif
    int old=BasePoint(Sptr);
#ifdef DEBUG_CHANGE_STAB_CHAIN
    std::cerr << "CPP ChangeStabChain old=" << ConvertIntFalse(old) << " i=" << (i+1) << " |base|=" << basSiz << "\n";
    KeyUpdating("After BasePoint");
#endif
    //    std::cerr << "eLev=" << eLev << "\n";
    if (Sptr->genlabels.size() == 0 && (reduced == int_true || i >= basSiz)) {
#ifdef DEBUG_CHANGE_STAB_CHAIN
      std::cerr << "CPP Before RemoveStabChain\n";
#endif
      int dep1=GetStabilizerDepth(Sptr);
      RemoveStabChain(Sptr);
      int dep2=GetStabilizerDepth(Sptr);
#ifdef DEBUG_CHANGE_STAB_CHAIN
      std::cerr << "CPP RemoveStabChain dep1=" << dep1 << " dep2=" << dep2 << "\n";
      KeyUpdating("After RemoveStabChain");
#endif
      i = basSiz;
    }
    else if (i < basSiz) {
      int newpnt = SlashAct(base[i], cnj);
      i++;
      if (reduced == int_reducedm1) {
	newBase.push_back(newpnt);
	if (newpnt != old) {
	  if (IsFixedStabilizer(Sptr, newpnt)) {
            int dep1=GetStabilizerDepth(Sptr);
	    InsertTrivialStabilizer(Sptr, newpnt);
            int dep2=GetStabilizerDepth(Sptr);
#ifdef DEBUG_CHANGE_STAB_CHAIN
            std::cerr << "CPP InsertTrivialStabilizer1 dep1=" << dep1 << " dep2=" << dep2 << "\n";
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
      }
      else if (reduced == int_false || !IsFixedStabilizer(Sptr, newpnt )) {
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
	}
	else {
          int dep1=GetStabilizerDepth(Sptr);
	  InsertTrivialStabilizer(Sptr, newpnt);
          int dep2=GetStabilizerDepth(Sptr);
#ifdef DEBUG_CHANGE_STAB_CHAIN
          std::cerr << "CPP InsertTrivialStabilizer2 dep1=" << dep1 << " dep2=" << dep2 << "\n";
          KeyUpdating("After InsertTrivialStabilizer");
#endif
	}
	newBase.push_back(Sptr->orbit[0]);
	Sptr = Sptr->stabilizer;
#ifdef DEBUG_CHANGE_STAB_CHAIN
        KeyUpdating("After S:=S.stabilizer 2");
#endif
      }
    }
    else if (PositionVect(newBase, old) != -1 || (reduced == int_true && Sptr->orbit.size() == 1)) {
#ifdef DEBUG_CHANGE_STAB_CHAIN
      std::cerr << "CPP Stabilizer shift in ChangeStabChain\n";
#endif
      int dep1=GetStabilizerDepth(Sptr);
      Sptr->stabilizer = Sptr->stabilizer->stabilizer;
      int dep2=GetStabilizerDepth(Sptr);
#ifdef DEBUG_CHANGE_STAB_CHAIN
      std::cerr << "CPP Manual Removal dep1=" << dep1 << " dep2=" << dep2 << "\n";
#endif
    }
    else {
      newBase.push_back(old);
      Sptr = Sptr->stabilizer;
#ifdef DEBUG_CHANGE_STAB_CHAIN
      KeyUpdating("After S:=S.stabilizer 3");
#endif
    }
  }
#ifdef DEBUG_CHANGE_STAB_CHAIN
  std::cerr << "CPP LEAVE GetStabilizerDepth(S)=" << GetStabilizerDepth(Sptr) << " i=" << (i+1) << " |base|=" << basSiz << "\n";
  std::cerr << "CPP Ending ChangeStabChain, GetStabilizerDepth(G) = " << GetStabilizerDepth(Gptr) << "\n";
  std::cerr << "CPP Ending ChangeStabChain, GetStabilizerDepth(S) = " << GetStabilizerDepth(Sptr) << "\n";
  std::cerr << "CPP sgs(G) = " << GapStringTVector(SortVector(StrongGeneratorsStabChain(Gptr))) << "\n";
  std::cerr << "CPP sgs(S) = " << GapStringTVector(SortVector(StrongGeneratorsStabChain(Sptr))) << "\n";
  std::cerr << "CPP ChangeStabChain 2 orbit=" << PrintTopOrbit(Gptr) << "\n";
  std::cerr << "CPP Before ConjugateStabChain cnj=" << cnj << "\n";
#endif
  if (!cnj.isIdentity())
    ConjugateStabChain(Gptr, cnj);
#ifdef DEBUG_CHANGE_STAB_CHAIN
  std::cerr << "CPP ChangeStabChain 3 orbit=" << PrintTopOrbit(Gptr) << "\n";
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
  while(true) {
    if (Lptr == nullptr && Rptr != nullptr)
      return false;
    if (Lptr != nullptr && Rptr == nullptr)
      return false;
    if (Lptr == nullptr)
      break;
    if (Lptr->orbit != R->orbit)
      return false;
    if (Lptr->transversal.size() != Rptr->transversal.size())
      return false;
    int lenL=Lptr->transversal.size();
    for (int u=0; u<lenL; u++) {
      if (Lptr->transversal[u] == -1 && Rptr->transversal[u] != -1)
	return false;
      if (Lptr->transversal[u] != -1 && Rptr->transversal[u] == -1)
	return false;
      if (Lptr->transversal[u] != -1) {
	int idxL=Lptr->transversal[u];
	int idxR=Rptr->transversal[u];
	Telt permL=L->comm->labels[idxL];
	Telt permR=L->comm->labels[idxR];
	if (permL != permR)
	  return false;
      }
    }
    Rptr = Rptr->stabilizer;
    Lptr = Lptr->stabilizer;
  }
  return true;
}


template<typename Telt>
void SetStabChainFromLevel(StabChain<Telt> & R, StabChain<Telt> const& L, int const& lev)
{
  StabChain<Telt> Rptr = R;
  StabChain<Telt> Lptr = L;
  int n=Lptr->comm->n;
  while(true) {
#ifdef DEBUG_PERMUTALIB
    if ((Rptr == nullptr && Lptr != nullptr) || (Rptr != nullptr && Lptr == nullptr)) {
      std::cerr << "Inconsistency\n";
      throw TerminalException{1};
    }
#endif
    if (Rptr == nullptr)
      break;
    Rptr->orbit = Lptr->orbit;
    for (int i=0; i<n; i++) {
      int idx=Lptr->transversal[i];
      int posPerm;
      if (idx == -1) {
	posPerm=-1;
      }
      else {
	Telt ePerm=Lptr->comm->labels[idx];
	posPerm = GetLabelIndex(Rptr->comm->labels, ePerm);
      }
      Rptr->transversal[i] = posPerm;
    }
    Rptr->genlabels.clear();
    for (int eVal : Lptr->genlabels) {
      Telt ePerm=Lptr->comm->labels[eVal];
      int pos = GetLabelIndex(Rptr->comm->labels, ePerm);
      Rptr->genlabels.push_back(pos);
    }
    Rptr = Rptr->stabilizer;
    Lptr = Lptr->stabilizer;
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
      gRet=LeftQuotient(Stot->comm->labels[pos], gRet);
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
  std::shared_ptr<CommonStabInfo<Tret>> comm = std::make_shared<CommonStabInfo<Tret>>(CommonStabInfo<Tret>({nMap, idMap, Stot->comm->UseCycle, labelsMap}));

  StabChain<Telt> Sptr = Stot;
  StabChain<Tret> Swork = nullptr;
  StabChain<Tret> Sreturn = nullptr;

  while(true) {
    if (Sptr == nullptr)
      break;
    Swork = std::make_shared<StabLevel<Telt>>(StabLevel<Telt>({Sptr->transversal, Sptr->orbit, Sptr->genlabels, Sptr->cycles, fVector(Sptr->treegen), fVector(Sptr->treegeninv), fVector(Sptr->aux), Sptr->treedepth, Sptr->diam, comm, nullptr}));
    if (Sreturn == nullptr)
      Sreturn = Swork;
    Sptr = Sptr->stabilizer;
    Swork = Swork->stabilizer;
  }
  return Sreturn;
}




}


#endif
