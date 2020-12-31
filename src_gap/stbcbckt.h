#ifndef DEFINE_STBCBCKT_H
#define DEFINE_STBCBCKT_H

#include "StabChainMain.h"
#include "partition.h"
//#include "COMB_Combinatorics.h"
#include "Combinatorics.h"
#include "COMB_Vectors.h"
#include "plus_infinity.h"

/*
#############################################################################
##
#W  stbcbckt.gi                 GAP library                    Heiko Theißen
##
##
#Y  Copyright (C)  1997,  Lehrstuhl D für Mathematik,  RWTH Aachen, Germany
#Y  (C) 1998 School Math and Comp. Sci., University of St Andrews, Scotland
#Y  Copyright (C) 2002 The GAP Group
##
*/

namespace permutalib {

template<typename Telt>
struct permPlusBool {
  int status; // values in {int_false, int_true, int_perm}
  Telt val;
};

template<typename Telt>
struct StabChainPlusLev {
  int status; // possible values in {int_false, int_true, int_int, int_stablev}
  int value_int;
  StabChain<Telt> Stot;
};


template<typename Telt>
StabChainPlusLev<Telt> StructuralCopy(StabChainPlusLev<Telt> const& S)
{
  return {S.status, S.value_int, StructuralCopy(S.Stot)};
}


std::string STRING_VectInt(std::vector<int> const& V)
{
  if (V.size() == 0)
    return "[  ]";
  std::string str = "[";
  for (size_t i=0; i<V.size(); i++) {
    if (i>0)
      str += ",";
    str += " " + std::to_string(V[i]);
  }
  str += " ]";
  return str;
}


  // The ExtendedT in gap seems to be passing by value for img entries that gets
  // modified
template<typename Telt>
permPlusBool<Telt> ExtendedT(Telt const& t, int const& pnt, int img, int const& simg, StabChainPlusLev<Telt> const& S)
{
  std::cerr << "CPP ExtendedT sgs(S.Stot)=" << GapStringTVectorB(SortVector(StrongGeneratorsStabChain(S.Stot))) << "\n";
  if (simg == -1)
    img = SlashAct(img, t);
  else
    img = simg;
  // If <G> fixes <pnt>, nothing more can  be changed, so test whether <pnt>
  // = <img>.
  int bpt = BasePoint(S.Stot);
  std::cerr << "CPP img=" << (img+1) << " bpt=" << PosFalse_to_string(bpt) << " pnt=" << (pnt+1) << "\n";
  if (bpt != pnt) {
    std::cerr << "CPP Case bpt != pnt\n";
    if (pnt != img) {
      std::cerr << "CPP ExtendedT, return false 1\n";
      return {int_false, {}};
    }
    else {
      return {int_perm, t};
    }
  }
  if (S.Stot->transversal[img] == -1) {
    std::cerr << "CPP ExtendedT, return false 2\n";
    return {int_false, {}};
  }
  std::cerr << "CPP Final case t=" << t << "\n";
  std::cerr << "CPP sgs(S.Stot)=" << GapStringTVector(SortVector(StrongGeneratorsStabChain(S.Stot))) << "\n";
  return {int_perm, LeftQuotient(InverseRepresentative(S.Stot, img), t)};
}


template<typename Telt>
bool IsBool(StabChainPlusLev<Telt> const& S)
{
  if (S.status == int_true || S.status == int_false)
    return true;
  return false;
}

template<typename Telt>
int BasePoint(StabChainPlusLev<Telt> const& S)
{
  return BasePoint(S.Stot);
}

struct singStrat {
  int p;
  int s;
  int i;
};


struct Refinement {
public:
  Refinement(int const& val1, int const& val2) {
    nature = 0;
    inputProcessfix = {val1, val2};
  }
  Refinement(Partition const& ePart, std::vector<singStrat> const& strat) {
    nature = 1;
    inputIntersection = {ePart, strat};
  }
  int nature; // 0 for PROCESSFIX, 1 for INTERSECTION
  std::pair<int,int> inputProcessfix;
  std::pair<Partition,std::vector<singStrat>> inputIntersection;
};

bool IsInsertableRefinement(Refinement const& eRfm)
{
  if (eRfm.nature == 0)
    return true;
  if (eRfm.nature == 1) {
    if (eRfm.inputIntersection.second.size() == 0)
      return false;
    return true;
  }
#ifdef DEBUG
  std::cerr << "We should never reach that stage\n";
  throw TerminalException{1};
#endif
  return true;
}


// The underscore nature of a function can be seen in stbcbckt top.
// Since we did not implement all the algorithms of stbcbckt, the value is always 0.
int UnderscoreNature(int const& nature)
{
  return 0;
}




template<typename Telt>
struct dataType {
  Partition P;
};


// The rbase if the main data type in the work
// critical thing is about the "lev" array in rbase.lev[d]
// and also rbase.lev2[d]
// lev is initially empty in EmptyRBase
// rbase.lev is expanded in RegisterRBasePoint
// rbase.lev[d] and rbase.lev2[d] values are NOT changed over the operations.
// Therefore if they are assigned as shared pointer from rbase.level
// a priori the value change with rbase.level and rbase.level2
// rbase.orbit is used in the computation (Wrong sizes to be resolved)
//
// What is apparent in runs is that there is a link between rbase.lev[d] and rbase.level
// We have Add(rbase.lev, rbase.level)
// Later rbase.level is changed and the value of rbase.lev[d] accordingly as well. This proves the
// existence of a link which is otherwise hard to find explicitly.
//
// The critical thing is therefore to find out when the .level is changed.
// We have to find out the formalism underlying the use of lev[d].
//
// FACTS:
// ---The rbase changes over time.
// ---rbase.lev[1].stabilizer is not the same as rbase.lev[2]
// ---rbase.lev is changed at
//    if not ProcessFixpoint( rbase, pnt )  then
//    ProcessFixpoint( rbase, pnt );
//    strat := StratMeetPartition( rbase, P, O ); # Which also uses ProcessFixpoint
//    Therefore ProcessFixpoint seems to be the critical entry
//      and the function used appears to be ChangeStabChain
// ---After the KeyUpdating the rbase.level is NOT the same as rbase.lev[last]
//    This is because we have the increment    rbase.level := rbase.level.stabilizer;
// ---ChangeStabChain is very complicated. We have S:=G, G not modified but
//      G has been altered by the operation. But S is not the same as G in the end.
//      So, the operations are just very complicated to follow.
//         Three operations (may more) on S modify also G:
//           StabChainForcePoint, InsertTrivialStabilizer and RemoveStabChain.
//      Only one operation seems to make S different from G:
//          The operations S:=S.stabilizer are operations that makes S different from G.
//    Therefore the operations of ChangeStabChain appears to be understood.
//    It can clearly be programmed down though that is reasonably non-trivial.
//    But S does not show up later on.
// ---So we need to find an adequate mechanism for dealing with the problem
//    of pointers.
template<typename Telt>
struct rbaseType {
  std::vector<int> domain;
  std::vector<int> base;
  std::vector<int> where;
  //
  dataType<Telt> data;
  std::vector<std::vector<int>> fix;
  //
  std::vector<std::vector<Refinement>> rfm;
  Partition partition;
  std::vector<StabChainPlusLev<Telt>> lev;
  StabChainPlusLev<Telt> level;
  //
  std::vector<StabChainPlusLev<Telt>> lev2;
  StabChainPlusLev<Telt> level2;
  //
  std::vector<std::string> levkey;
};


template<typename Telt>
void KeyUpdatingRbase(std::string const& str, rbaseType<Telt> & rbase)
{
  bool DoPrint=false;
  std::vector<std::string> ListKey;
  for (auto & x : rbase.lev)
    ListKey.push_back(GetStringExpressionOfStabChain(x.Stot));
  //
  if (DoPrint) {
    size_t len = rbase.lev.size();
    std::string strO = "[ ";
    for (size_t u=0; u<len; u++) {
      if (u>0)
        strO += ", ";
      strO += PrintTopOrbit(rbase.lev[u].Stot);
    }
    strO += " ]";
    std::cerr << "CPP KUR: at " << str << "\n";
    std::cerr << "CPP   Lorbit=" << strO << "\n";
    bool test = ListKey[len-1] == GetStringExpressionOfStabChain(rbase.level.Stot);
    std::cerr << "CPP KUR: at " << str << " test_equality=" << test << "\n";
    for (size_t i=0; i<len; i++)
      if (i<rbase.levkey.size())
        if (rbase.levkey[i] != ListKey[i])
          std::cerr << "CPP  KUR: Change of key at i=" << (i+1) << "\n";
    rbase.levkey = ListKey;
  }
}

template<typename Telt>
std::string ListOrbitOfRbaseLEV(rbaseType<Telt> const& rbase)
{
  std::string str = "[ ";
  int sizLev=rbase.lev.size();
  for (int iLev=0; iLev<sizLev; iLev++) {
    if (iLev > 0)
      str += ", ";
    std::vector<int> eOrb=rbase.lev[iLev].Stot->orbit;
    int len=eOrb.size();
    str += "[ ";
    for (int u=0; u<len; u++) {
      if (u>0)
	str += ", ";
      str += std::to_string(eOrb[u]);
    }
    str += " ]";
  }
  str += " ]";
  return str;
}


template<typename Telt>
void PrintRBaseLevel(rbaseType<Telt> const& rbase, std::string const& str)
{
  if (rbase.level.status == int_int) {
    std::cerr << str << " PRBL rbase.level, integer : " << rbase.level.value_int << "\n";
  }
  else {
    if (rbase.level.status == int_stablev) {
      int len=rbase.lev.size();
      std::cerr << str << " |rbase.lev|=" << len << "\n";
      for (int eD=0; eD<len; eD++) {
        std::cerr << "CPP rbase.lev[" << (eD+1) << "]=" << GapStringTVector(SortVector(StrongGeneratorsStabChain(rbase.lev[eD].Stot))) << "\n";
        PrintStabChainTransversals(rbase.lev[eD].Stot);
        PrintStabChainOrbits(rbase.lev[eD].Stot);
      }
      std::cerr << "CPP rbase.level=" << GapStringTVector(SortVector(StrongGeneratorsStabChain(rbase.level.Stot))) << "\n";
      std::cerr << str << " PRBL rbase.level, record, |genlabels|=" << rbase.level.Stot->genlabels.size() << "\n";
      std::cerr << str << " PRBL orbit=" << PrintTopOrbit(rbase.level.Stot) << "\n";
    }
    else {
      std::cerr << str << " PRBL rbase.level=" << GetIntTypeNature(rbase.level.status) << "\n";
    }
  }
}



template<typename Telt>
bool ProcessFixpoint_rbase(rbaseType<Telt> & rbase, int const& pnt)
{
  std::cerr << "CPP ProcessFixpoint_rbase beginning\n";
  if (rbase.level2.status != int_true && rbase.level2.status != int_false) {
    std::cerr << "CPP Before ChangeStabChain level2\n";
    ChangeStabChain(rbase.level2.Stot, {pnt}, int_true);
    PrintRBaseLevel(rbase, "CPP After CSC level2");
    std::cerr << "CPP After ChangeStabChain level2\n";
    if (BasePoint(rbase.level2) == pnt) {
      std::cerr << "CPP Going to stabilizer of level2\n";
      rbase.level2.Stot = rbase.level2.Stot->stabilizer;
    }
  }
  if (rbase.level.status == int_int) {
    rbase.level.value_int--;
  }
  else {
    std::cerr << "CPP Before ChangeStabChain level\n";
    ChangeStabChain(rbase.level.Stot, {pnt}, int_true);
    PrintRBaseLevel(rbase, "CPP After CSC level");
    std::cerr << "CPP After ChangeStabChain level\n";
    if (BasePoint(rbase.level) == pnt) {
      std::cerr << "CPP Going to stabilizer of level\n";
      rbase.level.Stot = rbase.level.Stot->stabilizer;
    }
    else {
      std::cerr << "CPP returning false\n";
      return false;
    }
  }
  return true;
}



template<typename Telt>
struct imageType {
  int depth;
  dataType<Telt> data;
  Partition partition;
  permPlusBool<Telt> perm;
  StabChainPlusLev<Telt> level;
  std::vector<int> bimg;
  //
  permPlusBool<Telt> perm2;
  StabChainPlusLev<Telt> level2;
};


template<typename Telt>
bool ProcessFixpoint_image(imageType<Telt> & image, int const& pnt, int & img, int const& simg)
{
  if (image.perm.status != int_true) {
    std::cerr << "CPP PFI  sgs(level)=" << GapStringTVector(SortVector(StrongGeneratorsStabChain(image.level.Stot))) << "\n";
    std::cerr << "CPP Case image.perm.status = true\n";
    std::cerr << "CPP Before ExtendedT img=" << (img+1) << "\n";
    permPlusBool<Telt> t = ExtendedT(image.perm.val, pnt, img, simg, image.level);
    std::cerr << "CPP After ExtendedT img=" << (img+1) << "\n";
    if (t.status == int_false) {
      std::cerr << "CPP Returning false 1\n";
      return false;
    }
    else {
      if (BasePoint(image.level ) == pnt)
        image.level.Stot = image.level.Stot->stabilizer;
    }
    image.perm = t;
  }
  if (image.level2.status != int_false) {
    std::cerr << "CPP PFI sgs(level2)=" << GapStringTVector(SortVector(StrongGeneratorsStabChain(image.level2.Stot))) << "\n";
    std::cerr << "CPP Case image.perm.status = false\n";
    permPlusBool<Telt> t = ExtendedT(image.perm2.val, pnt, img, -1, image.level2);
    if (t.status == int_false) {
      std::cerr << "CPP Returning false 2\n";
      return false;
    }
    else {
      if (BasePoint(image.level2 ) == pnt)
        image.level2.Stot = image.level2.Stot->stabilizer;
    }
    image.perm2 = t;
  }
  return true;
}


template<typename Telt>
bool IsTrivialRBase(rbaseType<Telt> const& rbase)
{
  std::cerr << "CPP IsTrivialRBase : IsInt()=";
  if (rbase.level.status == int_int)
    std::cerr << "true  value_int=" << rbase.level.value_int;
  else
    std::cerr << "false";
  if (rbase.level.status == int_int) {
  }
  std::cerr << "\n";
  //
  std::cerr << "CPP IsTrivialRBase : stab=";
  if (rbase.level.status == int_stablev) {
    std::cerr << "true  |genlabels|=" << rbase.level.Stot->genlabels.size();
  }
  else {
    std::cerr << "false";
  }
  std::cerr << "\n";
  //
  if (rbase.level.status == int_int) {
    if (rbase.level.value_int <= 1) {
      std::cerr << "CPP IsTrivialRBase, leaving at case 1 with True\n";
      return true;
    }
  }
  if (rbase.level.status == int_stablev) {
    if (rbase.level.Stot->genlabels.size() == 0) {
      std::cerr << "CPP IsTrivialRBase, leaving at case 2 with True\n";
      return true;
    }
  }
  std::cerr << "CPP IsTrivialRBase, leaving at case 3 with False\n";
  return false;
}



template<typename Telt>
rbaseType<Telt> EmptyRBase(std::vector<StabChain<Telt>> const& G, bool const& IsId, std::vector<int> const& Omega, Partition const& P)
{
  rbaseType<Telt> rbase;
  rbase.domain = Omega;
  rbase.base = {};
  rbase.where = {};
  rbase.rfm = {};
  rbase.partition = P;
  rbase.lev = {};
  if (G.size() == 2) {
    int n = GetNumberPoint(P);
    if (IsId) {
      rbase.level2.status = int_true;
      rbase.level2.Stot = EmptyStabChain<Telt>(n);
    }
    else {
      rbase.level2 = {int_stablev, -555, G[1]};
      std::cerr << "CPP rbase Before bool print\n";
      std::cerr << "CPP bool=" << (rbase.level2.Stot->cycles.size() > 0) << "\n";
      std::cerr << "CPP rbase After bool print\n";
      rbase.lev2 = {};
    }
  }
  else {
    rbase.level2.status = int_false;
  }
  rbase.level = {int_stablev, -666, G[0]};
  for (auto & pnt : Fixcells(P))
    ProcessFixpoint_rbase(rbase, pnt);
  return rbase;
}




template<typename Telt>
bool MeetPartitionStrat(rbaseType<Telt> const& rbase, imageType<Telt> & image, Partition const& S, Telt const& g, std::vector<singStrat> const& strat)
{
  if (strat.size() == 0)
    return false;
  for (auto & pRec : strat) {
    int eFix=FixpointCellNo(image.partition, pRec.i);
    std::cerr << "CPP ProcessFixpoint_image, Case MeetPartitionStrat\n";
    if ((pRec.p == -1 && !ProcessFixpoint_image(image, pRec.s, eFix, -1)) ||
	(pRec.p != -1 && SplitCell_Partition(image.partition, pRec.p, S, pRec.s, g, pRec.i ) != pRec.i))
      return false;
  }
  return true;
}

/*
#############################################################################
##
#F  StratMeetPartition( <rbase>, <P>, <S>, <g> )  . construct a meet strategy
##
##  Entries in <strat> have the following meaning:
##    [p,s,i] (p<>0) means that `0 < |P[p]\cap S[s]/g| = i < |P[p]|',
##            i.e., a new cell with <i> points was appended to <P>
##                  (and these <i> have been taken out of `P[p]'),
##    [0,a,p] means that fixpoint <a> was mapped to fixpoint in `P[p]',
##            i.e., `P[p]' has become a one-point cell.
##
*/
template<typename Telt>
std::vector<singStrat> StratMeetPartition(rbaseType<Telt> & rbase, Partition & P, Partition const& S, Telt const& g)
{
  //  std::cerr << "StratMeetPartition begin P\n";
  //  RawPrintPartition(P);
  //  std::cerr << "StratMeetPartition begin S\n";
  //  RawPrintPartition(S);
  //  std::cerr << "Now working\n";
  std::vector<singStrat> strat;
  std::vector<int> cellsP = P.cellno;
  if (!g.isIdentity()) {
    for (int i=0; i<NumberCells(P); i++) {
      std::vector<int> cell = Cell( P, i );
      for (auto & eVal : cell) {
        int img=PowAct(eVal, g);
	cellsP[img] = i;
      }
    }
  }
  //  PrintVectDebug("P.cellno=", P.cellno);
  //  PrintVectDebug("cellsP=", cellsP);
  // If <S> is just a set, it is interpreted as partition ( <S>|<S>^compl ).
  int nrcells = NumberCells(S) - 1;

  for (int s=0; s<nrcells; s++) {
    // now split with cell number s of S.
    //    std::cerr << "s=" << s << "\n";
    std::vector<int> p=Cell(S, s);
    //    PrintVectDebug("p", p);

    std::vector<int> p2;
    for (auto & eVal : p)
      p2.push_back(cellsP[eVal]);
    //    PrintVectDebug("p2", p2);
    CollectedResult<int> p3=Collected(p2);
    std::vector<int> splits;
    for (int h=0; h<int(p3.LVal.size()); h++) {
      // a cell will split iff it contains more points than are in the s-cell
      //      std::cerr << "h=" << h << " mult=" << p3.LMult[h] << " val=" << p3.LVal[h] << "\n";
      //      std::cerr << "Before if test\n";
      if (P.lengths[p3.LVal[h]] > p3.LMult[h])
        splits.push_back(p3.LVal[h]);
    }
    for (auto & pVal : splits) {
      // Last argument true means that the cell will split.
      int i = SplitCell_Partition(P, pVal, S, s, g, true);
      if (!g.isIdentity()) {
	std::vector<int> cell = Cell(P, NumberCells(P));
	for (auto & eVal : cell) {
	  int img=PowAct(eVal, g);
	  cellsP[img] = NumberCells(P);
	}
      }
      strat.push_back({pVal, s, i});
      // If  we have one  or two  new fixpoints, put  them  into the base.
      if (i == 0) {
        int pnt = FixpointCellNo(P, NumberCells(P));
	ProcessFixpoint_rbase(rbase, pnt);
	strat.push_back({-1, pnt, NumberCells(P)});
	if (IsTrivialRBase(rbase))
	  return strat;
      }
      if (P.lengths[pVal] == 1) {
        int pnt = FixpointCellNo(P, pVal);
	ProcessFixpoint_rbase(rbase, pnt);
	strat.push_back({-1, pnt, pVal});
	if (IsTrivialRBase(rbase))
	  return strat;
      }
    }
  }
  return strat;
}

template<typename Telt>
void AddRefinement(rbaseType<Telt> & rbase, int const& pos, Refinement const& eRfm)
{
  std::cerr << "CPP beginning of AddRefinement\n";
  if (IsInsertableRefinement(eRfm)) {
    std::cerr << "CPP Doing RFM insertion\n";
    rbase.rfm[pos].push_back(eRfm);
  }
  for (size_t i=0; i<rbase.rfm.size(); i++) {
    std::cerr << "CPP i=" << (i+1) << " |rbase.rfm[i]|=" << rbase.rfm[i].size() << "\n";
  }
}


template<typename Telt>
void RegisterRBasePoint(Partition & P, rbaseType<Telt> & rbase, int const& pnt, Telt const& TheId)
{
  if (rbase.level2.status != int_true && rbase.level2.status != int_false) {
    std::cerr << "CPP Inserting rbase.level2 into rbase.lev2\n";
    std::cerr << "CPP rbase.level2.status=" << GetIntTypeNature(rbase.level2.status) << "\n";
    rbase.lev2.push_back(rbase.level2);
  }
  PrintRBaseLevel(rbase, "CPP RegisterRBasePoint 1");
  rbase.lev.push_back(rbase.level);
  rbase.base.push_back(pnt);
  KeyUpdatingRbase("RegisterRBasePoint 1", rbase);
  int k = IsolatePoint(P, pnt);
  NicePrintPartition("CPP After IsolatePoint P", P);
  KeyUpdatingRbase("RegisterRBasePoint 1.1", rbase);
  if (!ProcessFixpoint_rbase(rbase, pnt)) {
    std::cerr << "CPP INFO: Warning R-base point is already fixed\n";
  }
  //  rbase.lev.push_back(rbase.level);
  KeyUpdatingRbase("RegisterRBasePoint 1.2", rbase);
  PrintRBaseLevel(rbase, "CPP RegisterRBasePoint 2");
  rbase.where.push_back(k);
  int len=rbase.rfm.size();
  rbase.rfm.push_back({});
  std::cerr << "CPP Before P.lengths test k=" << (k+1) << " len=" << rbase.rfm.size() << "\n";
  KeyUpdatingRbase("RegisterRBasePoint 1.3", rbase);
  if (P.lengths[k] == 1) {
    std::cerr << "CPP Matching P.lengths test\n";
    int pnt = FixpointCellNo(P, k);
    std::cerr << "CPP Section P.lengths after FixpointCellNo pnt=" << (pnt+1) << "\n";
    PrintRBaseLevel(rbase, "CPP RegisterRBasePoint 2.1");
    ProcessFixpoint_rbase(rbase, pnt);
    PrintRBaseLevel(rbase, "CPP RegisterRBasePoint 2.2");
    KeyUpdatingRbase("RegisterRBasePoint 1.4", rbase);
    std::cerr << "CPP Section P.lengths after ProcessFixpoint_rbase\n";
    AddRefinement(rbase, len, Refinement({pnt,k}));
    std::cerr << "CPP After AddRefinement 1\n";
    KeyUpdatingRbase("RegisterRBasePoint 1.5", rbase);
  }
  PrintRBaseLevel(rbase, "CPP RegisterRBasePoint 3");
  KeyUpdatingRbase("RegisterRBasePoint 2", rbase);
  if (rbase.level2.status != int_false) {
    std::cerr << "CPP Matching the ! false test\n";
    auto MainInsert=[&](StabChainPlusLev<Telt> const& lev) -> void {
      if (lev.status != int_int) {
	std::vector<Telt> LGen = StrongGeneratorsStabChain(lev.Stot);
	std::cerr << "CPP StrongGeneratorsStabChain(lev) = " << GapStringTVector(SortVector(LGen)) << "\n";
	Partition O = OrbitsPartition(LGen, lev.Stot->comm->n, rbase.domain);
	NicePrintPartition("CPP Before StratMeetPartition O", O);
        KeyUpdatingRbase("RegisterRBasePoint 2.1", rbase);
	std::vector<singStrat> strat = StratMeetPartition(rbase, P, O, TheId);
        KeyUpdatingRbase("RegisterRBasePoint 2.2", rbase);
        AddRefinement(rbase, len, Refinement({O,strat}));
        std::cerr << "CPP After AddRefinement 2\n";
      }
    };
    if (rbase.level2.status == int_true) {
      std::cerr << "CPP Before call to MainInsert(level)\n";
      MainInsert(rbase.level);
    }
    else {
      std::cerr << "CPP Before call to MainInsert(level2)\n";
      MainInsert(rbase.level2);
    }
  }
  //  rbase.lev.push_back(rbase.level);
  KeyUpdatingRbase("RegisterRBasePoint 3", rbase);
}




template<typename Telt>
void NextRBasePoint(Partition & P, rbaseType<Telt> & rbase, Telt const& TheId)
{
  //  std::cerr << "Working with NextRBasePoint rbase.level2.status=" << GetIntTypeNature(rbase.level2.status) << "\n";
  //  RawPrintPartition(P);
  std::vector<int> lens = P.lengths;
  std::vector<int> order = ClosedInterval(0, NumberCells(P));
  std::cerr << "CPP lens=[ ";
  for (size_t i=0; i<lens.size(); i++) {
    if (i>0)
      std::cerr << ", ";
    std::cerr << lens[i];
  }
  std::cerr << " ]\n";
  //  std::cerr << "Before SortParallel\n";
  //  PrintVectDebug(" lens", lens);
  //  PrintVectDebug("order", order);
  //
  SortParallel(lens, order);
  //  std::cerr << "After SortParallel\n";
  //  PrintVectDebug(" lens", lens);
  //  PrintVectDebug("order", order);
  //


  int k = PositionProperty(lens, [](int const& x) -> int {return x != 1;});
  //  std::cerr << "Starting at k=" << k << "\n";
  int l = -1;
  if (rbase.level.status == int_int) {
    l = 0;
  }
  else {
    while (true) {
      //      std::cerr << "Before PositionProperty operation len[k]=" << lens[k] << "\n";
      l = PositionProperty(ClosedInterval(0, lens[k]), [&](int const& i) -> bool {
	  return !IsFixedStabilizer(rbase.level.Stot, P.points[i+P.firsts[order[k]]]);});
      //      std::cerr << "At k=" << k << " found l=" << l << "\n";
      if (l != -1)
	break;
      k++;
    }
  }
  //  std::cerr << "k=" << k << " l=" << l << "\n";
  //  std::cerr << "order[k]=" << order[k] << "\n";
  //  std::cerr << "P.firsts[order[k]]=" << P.firsts[order[k]] << "\n";
  int p = P.points[ P.firsts[ order[k] ] + l ];
  std::cerr << "CPP p=" << (p+1) << "\n";
  NicePrintPartition("CPP Before RegisterRBasePoint P", P);
  PrintRBaseLevel(rbase, "CPP Before RegisterRBasePoint");
  RegisterRBasePoint(P, rbase, p, TheId);
}

template<typename Telt>
bool Refinements_ProcessFixpoint(rbaseType<Telt> & rbase, imageType<Telt> & image, int const& pnt, int const& cellnum)
{
 int img = FixpointCellNo(image.partition, cellnum);
 std::cerr << "CPP ProcessFixpoint_image, Case Refinements_ProcessFixpoint\n";
 return ProcessFixpoint_image(image, pnt, img, -1);
}


template<typename Telt>
bool Refinements_Intersection(rbaseType<Telt> & rbase, imageType<Telt> & image, Partition const& Q, std::vector<singStrat> const& strat)
{
  Telt t;
  if (image.level2.status == int_false) {
    t = image.perm.val;
  }
  else {
    t = image.perm2.val;
  }
  Telt tinv =Inverse(t);
  return MeetPartitionStrat(rbase, image, Q, t, strat);
}

// The function RRefine is doing the computation using CallFuncList
// It processes a number of refinement strategies.
// The functions Refinements used return only booleans
//
template<typename Telt>
int RRefine(rbaseType<Telt> & rbase, imageType<Telt> & image, bool const& uscore)
{
  std::cerr << "CPP uscore=" << uscore << "\n";
  auto BoolToInt=[&](bool const& val) -> int {
    if (val)
      return int_true;
    return int_false;
  };
  auto Evaluation=[&](Refinement const& eRef) -> bool {
    if (eRef.nature == int_false)
      return Refinements_ProcessFixpoint(rbase, image, eRef.inputProcessfix.first, eRef.inputProcessfix.second);
    if (eRef.nature == int_true)
      return Refinements_Intersection(rbase, image, eRef.inputIntersection.first, eRef.inputIntersection.second);
    return true;
  };
  if (!uscore) {
    for (auto & Rf : rbase.rfm[image.depth]) {
      std::cerr << "CPP Doing one CallFuncList 1\n";
      bool t = Evaluation(Rf);
      if (!t) {
	return int_fail;
      }
      else {
	if (!t) {
	  return BoolToInt(t);
	}
      }
    }
    return int_true;
  }
  else {
    for (auto & Rf : rbase.rfm[image.depth]) {
      std::cerr << "CPP Doing one CallFuncList 2\n";
      if (UnderscoreNature(Rf.nature)) {
	bool t = Evaluation(Rf);
	if (!t) {
	  return int_fail;
	}
	else {
	  if (!t) {
	    return BoolToInt(t);
	  }
	}
      }
    }
    return int_true;
  }
}


template<typename Telt>
bool PBIsMinimal(std::vector<int> const& range, int const& a, int const& b, StabChain<Telt> const& S)
{
  if (IsInBasicOrbit(S, b)) {
    for (auto & pVal : S->orbit) {
      if (a > pVal)
        return false;
    }
    return true;
  }
  if (b < a)
    return false;
  if (IsFixedStabilizer(S, b))
    return true;

  std::vector<int> orb{b};
  int pos=0;
  Face old = BlistList(range, orb);
  while(true) {
    int siz=orb.size();
    if (pos == siz)
      break;
    for (int i=pos; i<siz; i++) {
      int pnt=orb[i];
      for (auto & lVal : S->genlabels) {
        int img = PowAct(pnt, S->comm->labels[lVal]);
        if (!old[img]) {
          if (img < a)
            return false;
          old[img]=true;
          orb.push_back(img);
        }
      }
    }
    pos=siz;
  }
  return true;
}

template<typename Telt>
void SubtractBlistOrbitStabChain(Face & blist, std::vector<Telt> const& LGen, int const& pnt_in)
{
  std::vector<int> orb{pnt_in};
  blist[pnt_in]=false;
  int pos=0;
  //  int PrevPos=0;
  while(true) {
    int siz = orb.size();
    if (pos == siz) {
      break;
    }
    for (int ePos=pos; ePos<siz; ePos++) {
      int pnt=orb[ePos];
      for (auto& eGen : LGen) {
        int img = PowAct(pnt, eGen);
        if (blist[img]) {
          blist[img]=false;
          orb.push_back(img);
        }
      }
    }
    pos=siz;
  }
}


template<typename Telt>
struct ResultPBT {
  int nature; // Allowed values in {int_group, int_fail, int_perm}.
              // int_group for group
              // int_fail for fail
              // int_perm for equivalence element
  StabChain<Telt> stab;
  Telt res;
};


template<typename Telt>
Telt MappingPermListList(int const& n, std::vector<int> const& src, std::vector<int> const& dst)
{
  std::vector<int> ListImage(n);
  Face StatusSrc(n);
  Face StatusDst(n);
  for (int i=0; i<n; i++) {
    StatusSrc[i]=1;
    StatusDst[i]=1;
  }
  int len = src.size();
  for (int i=0; i<len; i++) {
    ListImage[src[i]] = dst[i];
    StatusSrc[src[i]] = 0;
    StatusDst[dst[i]] = 0;
  }
  int sizRemain = n - len;
  boost::dynamic_bitset<>::size_type posSrc=StatusSrc.find_first();
  boost::dynamic_bitset<>::size_type posDst=StatusDst.find_first();
  for (int u=0; u<sizRemain; u++) {
    ListImage[posSrc] = posDst;
    posSrc=StatusSrc.find_next(posSrc);
    posDst=StatusDst.find_next(posDst);
  }
  return Telt(ListImage);
}


template<typename T>
void AssignationVectorGapStyle(std::vector<T> & eVect, int const& pos, T const& val)
{
  int siz=eVect.size();
  if (pos < siz)
    eVect[pos] = val;
#ifdef DEBUG
  if (pos != siz) {
    std::cerr << "Assignation leaves gap in the vector. Not allowed\n";
    throw TerminalException{1};
  }
#endif
  eVect.push_back(val);
}

std::string GetStringGAP(Face const& f)
{
  std::string str = "[ ";
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

void PrintVectorORB(std::string const& str, std::vector<Face> const& eV)
{
  int len=0;
  int siz=eV.size();
  if (siz > 0)
    len = eV[0].size();
  std::cerr << "Printing std:vector<Face> : " << str << " |eV|=" << siz << " len=" << len << "\n";
  for (int i=0; i<siz; i++) {
    std::cerr << "i=" << i << " v=" << GetStringGAP(eV[i]) << "\n";
  }
}



template<typename Telt, typename Tint>
ResultPBT<Telt> PartitionBacktrack(StabChain<Telt> const& G, std::function<bool(Telt const&)> const& Pr, bool const& repr, rbaseType<Telt> & rbase, dataType<Telt> const& data, StabChain<Telt> & L, StabChain<Telt> & R)
{
  int n=G->comm->n;
  std::cerr << "CPP PartitionBacktrack step 1\n";
  std::cerr << "CPP L=\n";
  PrintStabChain(L);
  std::cerr << "CPP R=\n";
  PrintStabChain(R);
  std::cerr << "CPP INIT sgs(G)=" << GapStringTVector(SortVector(StrongGeneratorsStabChain(G))) << "\n";
  std::cerr << "CPP INIT sgs(L)=" << GapStringTVector(SortVector(StrongGeneratorsStabChain(L))) << "\n";
  std::cerr << "CPP INIT sgs(R)=" << GapStringTVector(SortVector(StrongGeneratorsStabChain(R))) << "\n";
  imageType<Telt> image;
  std::vector<Face> orB; // backup of <orb>. We take a single entry. Not sure it is correct
  int nrback;
  std::vector<Face> orb;
  std::vector<std::vector<int>> org; // intersected (mapped) basic orbits of <G>
  Tplusinfinity<int> blen(true, 0);
  int dd, branch; // branch is level where $Lstab\ne Rstab$ starts
  std::vector<int> range;    // range for construction of <orb>
  Partition oldcel;       // old value of <image.partition.cellno>
  std::vector<int> oldcel_cellno;
  std::vector<StabChain<Telt>> L_list, R_list;
  std::function<permPlusBool<Telt>(int const&,bool const&)> PBEnumerate = [&](int const& d, bool const & wasTriv) -> permPlusBool<Telt> {
    std::cerr << "CPP PBEnumerate, step 1, d=" << (d+1) << " wasTriv=" << wasTriv << "\n";
    //    PrintVectorORB("orb", orb);
    //    PrintVectorORB("orB", orB);
    permPlusBool<Telt> oldprm, oldprm2;
    int a;                // current R-base point
    permPlusBool<Telt> t; // group element constructed, to be handed upwards
    int m;                // initial number of candidates in <orb>
    int max;              // maximal number of candidates still needed
    boost::dynamic_bitset<>::size_type b;        // image of base point currently being considered

    if (image.perm.status == int_false) {
      std::cerr << "CPP PBEnumerate, EXIT 1, image.perm.status=" << GetIntTypeNature(image.perm.status) << "\n";
      return {int_fail, {}};
    }
    image.depth = d;
    std::cerr << "CPP PBEnumerate, step 2\n";
    PrintRBaseLevel(rbase, "CPP Step 2");
    //    PrintVectorORB("orb", orb);
    //    PrintVectorORB("orB", orB);

    // Store the original values of <image.*>.
    int undoto = NumberCells(image.partition);
    if (image.perm.status == int_true) {
      oldcel = image.partition;
    }
    else {
      oldcel_cellno = image.partition.cellno;
      oldprm = image.perm;
    }
    std::cerr << "CPP PBEnumerate, step 3\n";
    if (image.level2.status != int_false)
      oldprm2 = image.perm2;
    else
      oldprm2.status = int_false;
    std::cerr << "CPP PBEnumerate, step 4 d=" << (d+1) << " |rbase.base|=" << rbase.base.size() << "\n";
    PrintRBaseLevel(rbase, "CPP Step 4");
    // Recursion comes to an end  if all base  points have been prescribed
    // images.
    if (d >= int(rbase.base.size())) {
      std::cerr << "CPP Matching d > Length(rbase.base) test\n";
      if (IsTrivialRBase(rbase)) {
	blen = rbase.base.size();
	std::cerr << "CPP IsTrivialRBase matching test blen=" << blen << " wasTriv=" << wasTriv << "\n";
	// Do     not  add the   identity    element  in the  subgroup
	// construction.
	if (wasTriv) {
	  std::cerr << "CPP wasTriv Critical, step 1\n";
	  // In the subgroup case, assign to  <L> and <R> stabilizer
	  // chains when the R-base is complete.
	  StabChainOptions<Tint> options = GetStandardOptions<Tint>(n);
	  options.base = rbase.base;
	  options.reduced = false;
	  std::cerr << "CPP Before computation of ListStabChain Order(L)=" << Order<Telt,mpz_class>(L) << "\n";
          std::cerr << "CPP sgs(L)=" << GapStringTVector(SortVector(StrongGeneratorsStabChain(L))) << " base=" << GapStringIntVector(rbase.base) << "\n";
	  L_list = ListStabChain(StabChainOp_stabchain_nofalse<Telt,Tint>(L, options));
	  std::cerr << "CPP ListStabChain |L|=" << L_list.size() << "\n";
	  std::cerr << "CPP wasTriv Critical, step 2\n";
	  R_list = L_list;
	  std::cerr << "CPP wasTriv Critical, step 3\n";
	  std::cerr << "CPP PBEnumerate, EXIT 2 |L|=" << L_list.size() << "\n";
	  return {int_fail,{}};
	}
	else {
	  permPlusBool<Telt> prm;
	  if (image.perm.status == int_true)
	    prm = {int_perm, MappingPermListList<Telt>(n, rbase.fix[rbase.base.size()-1], Fixcells(image.partition))};
	  else
	    prm = image.perm;
	  if (image.level2.status != int_false) {
	    if (SiftedPermutation(image.level2.Stot, prm.val * Inverse(image.perm2.val)).isIdentity()) {
	      std::cerr << "CPP PBEnumerate, EXIT 3 |L|=" << L_list.size() << "\n";
	      return prm;
	    }
	  }
	  else {
	    if (Pr(prm.val)) {
	      std::cerr << "CPP PBEnumerate, EXIT 4 |L|=" << L_list.size() << "\n";
	      return {int_perm, prm.val};
	    }
	  }
	  std::cerr << "CPP PBEnumerate, EXIT 5 |L|=" << L_list.size() << "\n";
	  return {int_fail, {}};
	}
	// Construct the   next refinement  level. This  also  initializes
	// <image.partition> for the case ``image = base point''.
      }
      else {
	std::cerr << "CPP Not matching IsTrivialRBase test\n";
        PrintRBaseLevel(rbase, "CPP Before NextRBasePoint");
	NextRBasePoint(rbase.partition, rbase, G->comm->identity);
	PrintRBaseLevel(rbase, "CPP After NextRBasePoint");
	if (image.perm.status == int_true)
	  rbase.fix.push_back(Fixcells(rbase.partition));
	std::cerr << "CPP After Fixcells insert\n";
	std::vector<int> eNewF(range.size(), 0);
	org.push_back(eNewF);
	if (repr) {
	  // In  the representative  case,  change  the   stabilizer
	  // chains of <L> and <R>.
          std::cerr << "CPP Before ChangeStabChain L_list[d]\n";
	  ChangeStabChain(L_list[d], {rbase.base[d]}, int_false);
          PrintStabChainTransversals(L_list[d]);
          std::cerr << "CPP After ChangeStabChain L_list[d]\n";
	  //	  L[ d + 1 ] := L[ d ].stabilizer;
          std::cerr << "CPP Before ChangeStabChain R_list[d]\n";
	  ChangeStabChain(R_list[d], {rbase.base[d]}, int_false);
          PrintStabChainTransversals(R_list[d]);
          std::cerr << "CPP After ChangeStabChain R_list[d]\n";
	  //	  R[ d + 1 ] := R[ d ].stabilizer;
	}
      }
    }
    a = rbase.base[d];
    std::cerr << "CPP PBEnumerate, step 5\n";
    PrintRBaseLevel(rbase, "CPP Step 5");

    // Intersect  the current cell of <P>  with  the mapped basic orbit of
    // <G> (and also with the one of <H> in the intersection case).
    if (image.perm.status == int_true) {
      AssignationVectorGapStyle(orb, d, BlistList(range, Cell(oldcel, rbase.where[d]) ));
      if (image.level2.status != int_false) {
	b = orb[d].find_first();
	while (b != boost::dynamic_bitset<>::npos) {
	  if (!IsInBasicOrbit(rbase.lev2[d].Stot, SlashAct(b, image.perm2.val)))
	    orb[d][b] = false;
	  b = orb[d].find_next(b);
	}
      }
      std::cerr << "CPP ORB: Case image.perm=true d=" << (d+1) << " orb[d]=" << GetStringGAP(orb[d]) << "\n";
    }
    else {
      std::cerr << "CPP image.perm<>true orb=" << GapStringListBoolVector(orb) << "\n";
      AssignationVectorGapStyle(orb, d, BlistList(range, {}));
      // line below needs to be checked.
      std::cerr << "CPP ORB: Before pVal loop d=" << (d+1) << " orb[d]=" << GetStringGAP(orb[d]) << "\n";
      std::cerr << "CPP RBASE: List(...) = [ ";
      bool IsF=true;
      for (auto & eRec : rbase.lev) {
        if (!IsF)
          std::cerr << ", ";
        IsF=false;
	std::cerr << PrintTopOrbit(eRec.Stot);
      }
      std::cerr << " ]\n";
      //      std::cerr << "CPP d=" << d << " |rbase.lev|=" << rbase.lev.size() << "\n";
      for (auto & pVal : rbase.lev[d].Stot->orbit) {
	b = PowAct(pVal, image.perm.val);
	std::cerr << "CPP pVal=" << (pVal+1) << " b=" << (b+1) << "\n";
	if (oldcel_cellno[b] == rbase.where[d]) {
	  bool DoOper=false;
	  if (image.level2.status == int_false)
	    DoOper=true;
	  if (!DoOper) {
            //	    std::cerr << "CPP d=" << d << " |rbase.lev2|=" << rbase.lev2.size() << "\n";
	    if (IsInBasicOrbit(rbase.lev2[d].Stot, SlashAct(b,image.perm2.val)))
	      DoOper=true;
	  }
	  if (DoOper) {
	    orb[d][b] = true;
	    org[d][b] = pVal;
	  }
	}
      }
      std::cerr << "CPP ORB: After pVal loop d=" << (d+1) << " orb[d]=" << GetStringGAP(orb[d]) << "\n";
      //      std::cerr << "After pVal loop\n";
    }
    std::cerr << "CPP PBEnumerate, step 6\n";
    PrintRBaseLevel(rbase, "CPP Step 6");
    if (d == 0 && ForAll(G->comm->labels, [&](Telt const& x){return PowAct(a, x) == a;})) {
      orb[d][a]=true; // ensure a is a possible image (can happen if acting on permutations with more points)
      std::cerr << "CPP ORB: After assignation d=" << (d+1) << " orb[d]=" << GetStringGAP(orb[d]) << "\n";
    }
    AssignationVectorGapStyle(orB, d, orb[d]);
    std::cerr << "CPP PBEnumerate, step 7, wasTriv=" << wasTriv << "\n";
    PrintRBaseLevel(rbase, "CPP Step 7");
    
    // Loop  over the candidate images  for the  current base point. First
    // the special case image = base up to current level.
    if (wasTriv) {
      std::cerr << "CPP wasTriv Critical, step 4\n";
      AssignationVectorGapStyle(image.bimg, d, a);
      std::cerr << "CPP wasTriv Critical, step 5\n";
      // Refinements that start with '_' must be executed even when base
      // = image since they modify image.data, etc.
      RRefine(rbase, image, true);
      PrintRBaseLevel(rbase, "CPP After RRefine");
      // Recursion.
      PBEnumerate(d + 1, true);
      std::cerr << "CPP wasTriv Critical, step 6 d=" << (d+1) << "\n";
      image.depth = d;
      // Now we  can  remove  the  entire   <R>-orbit of <a>  from   the
      // candidate list.
      std::cerr << "CPP ORB 1: Before subtract d=" << (d+1) << " orb[d]=" << GetStringGAP(orb[d]) << "\n";
      std::cerr << "CPP |L|=" << L_list.size() << "\n";
      SubtractBlist(orb[d], BlistList(range, L_list[d]->orbit));
      std::cerr << "CPP ORB 1: After subtract d=" << (d+1) << " orb[d]=" << GetStringGAP(orb[d]) << "\n";
    }
    std::cerr << "CPP PBEnumerate, step 8\n";

    // Only the early points of the orbit have to be considered.
    //    PrintVectorORB("orb", orb);
    //    PrintVectorORB("orB", orB);
    m = SizeBlist( orB[d] );
    if (m < int(L_list[d]->orbit.size()) ) {
      std::cerr << "CPP PBEnumerate, EXIT 6 |L|=" << L_list.size() << "\n";
      return {int_fail,{}};
    }
    max = PositionNthTrueBlist(orB[d], m - L_list[d]->orbit.size());
    std::cerr << "CPP PBEnumerate, step 9\n";
    std::cerr << "CPP wasTriv=" << wasTriv << " a=" << (a+1) << " max=" << (max+1) << "\n";
    //    PrintVectorORB("orb", orb);
    //    PrintVectorORB("orB", orB);

    if (wasTriv && a > max) {
      m--;
      std::cerr << "CPP Before test m=" << m << " Length(L[d].orbit)=" << L_list[d]->orbit.size() << "\n";
      if (m < int(L_list[d]->orbit.size()) ) {
	std::cerr << "CPP PBEnumerate, EXIT 7 |L|=" << L_list.size() << "\n";
	return {int_fail,{}};
      }
      max = PositionNthTrueBlist( orB[d], m - L_list[d]->orbit.size());
    }
    std::cerr << "CPP PBEnumerate, step 10\n";
    // Now the other possible images.
    b = orb[d].find_first();
    while (b != boost::dynamic_bitset<>::npos) {
      int b_int = int(b);
      std::cerr << "CPP b=" << (b+1) << " b_int=" << (b_int+1) << "\n";
      // Try to prune the node with prop 8(ii) of Leon paper.
      if (!repr && !wasTriv) {
	dd = branch;
	while (dd < d) {
	  if (IsInBasicOrbit(L_list[dd], a) && !PBIsMinimal(range, R_list[dd]->orbit[0], b_int, R_list[d]))
	    dd = d + 1;
	  else
	    dd = dd + 1;
	}
      }
      else {
	dd = d;
      }
      std::cerr << "CPP dd=" << (dd+1) << " d=" << (d+1) << "\n";
      if (dd == d) {
        std::cerr << "CPP equality dd=d\n";
	// Undo the  changes made to  <image.partition>, <image.level>
	// and <image.perm>.
	for (int i=undoto+1; i<NumberCells(image.partition); i++)
	  UndoRefinement(image.partition);
	if (image.perm.status != int_true) {
          std::cerr << "CPP assignation image.level\n";
	  image.level = rbase.lev[d];
	  image.perm = oldprm;
	}
	if (image.level2.status != int_false) {
          std::cerr << "CPP assignation image.level2\n";
	  image.level2 = rbase.lev2[d];
	  image.perm2  = oldprm2;
	}
	// If <b> could not be prescribed as image for  <a>, or if the
	// refinement was impossible, give up for this image.
        std::cerr << "CPP Before AssignationVectorGapStyle b_int=" << (b_int+1) << "\n";
	AssignationVectorGapStyle(image.bimg, d, b_int);
        std::cerr << "CPP Before IsolatePoint b_int=" << (b_int+1) << "\n";
	IsolatePoint( image.partition, b_int );
        std::cerr << "CPP ProcessFixpoint_image, Case PartitionBacktrack 1\n";
        std::cerr << "CPP Before ProcessFixpoint_image b_int=" << (b_int+1) << "\n";
	bool val = ProcessFixpoint_image(image, a, b_int, org[d][b_int]);
	std::cerr << "CPP a=" << (a+1) << " b=" << (b_int+1) << " org[d][b]=" << (org[d][b_int]+1) << " val=" << val << "\n";
	if (val)
	  t.status = RRefine(rbase, image, false);
	else
	  t.status = int_fail;
        std::cerr << "CPP After assignment of t. t.status=" << GapStringBool(t.status) << "\n";

	if (t.status != int_fail) {
          std::cerr << "CPP case of not fail\n";
	  // Subgroup case, base <> image   at current level:   <R>,
	  //   which until now is identical to  <L>, must be changed
	  //   without affecting <L>, so take a copy.
          std::cerr << "CPP wasTriv=" << wasTriv << " d=" << (d+1) << "\n";
          //          std::cerr << "CPP L[d]=";
          //          PrintStabChainTransversals(L_list[d]);
          //          std::cerr << "CPP R[d]=";
          //          PrintStabChainTransversals(R_list[d]);
          std::cerr << "CPP TestEquality=" << TestEqualityStabChain(L_list[d], R_list[d]) << "\n";
	  if (wasTriv && TestEqualityStabChain(L_list[d], R_list[d])) {
            std::cerr << "CPP Assigning R from d\n";
	    SetStabChainFromLevel(R, L, d);
	    branch = d;
	  }
          std::cerr << "CPP After wasTriv test\n";
          std::cerr << "CPP d=" << (d+1) << " blen=" << blen << "\n";
	  if (2 * (d+1) <= blen) {
            std::cerr << "CPP Before ChangeStabChain R_list[d] 2\n";
	    ChangeStabChain(R_list[d], {b_int}, int_false);
            std::cerr << "CPP After ChangeStabChain R_list[d] 2\n";
	    //	    R[ d + 1 ] = R[ d ].stabilizer;
	  }
	  else {
            std::cerr << "CPP Beginning else case\n";
	    std::vector<Telt> LGen = StrongGeneratorsStabChain( R);
            std::cerr << "CPP First generating step done\n";
	    std::vector<Telt> LGenB = Filtered(LGen, [&](Telt const& gen) -> bool {return PowAct(b_int, gen) == b_int;});
            std::cerr << "CPP |LGenB|=" << LGenB.size() << "\n";
	    //	    R[ d + 1 ] := rec( generators := Filtered( R[ d + 1 ], gen -> b ^ gen = b ) );
	    int largMov=LargestMovedPoint(LGenB);
	    StabChainOptions<Tint> options = GetStandardOptions<Tint>(n);
	    options.base = ClosedInterval(0, largMov);
            std::cerr << "XXX ELIMINATE begin\n";
	    R_list[d+1] = StabChainOp_listgen(LGenB, options);
            std::cerr << "XXX ELIMINATE end\n";
	  }
	}
        std::cerr << "CPP t step 2\n";
	//	PrintVectorORB("orb", orb);
	//	PrintVectorORB("orB", orB);

	// Recursion.
	if (t.status == int_true) {
	  t = PBEnumerate(d + 1, false);
	  nrback++;
	  image.depth = d;
	}
        std::cerr << "CPP t step 3\n";

	// If   <t>   =   fail, either   the   recursive   call  was
	//   unsuccessful,  or all new  elements   have been added  to
	//   levels  below  the current one   (this happens if  base =
	//   image up to current level).
	if (t.status != int_fail) {
          std::cerr << "CPP Matching t<>fail\n";
	  // Representative case, element found: Return it.
	  // Subgroup case, base <> image  before current level:  We
	  //   need  only find  a representative  because we already
	  //   know the stabilizer of <L> at an earlier level.
	  if (repr || !wasTriv) {
	    std::cerr << "CPP PBEnumerate, EXIT 8 |L|=" << L_list.size() << "\n";
	    return t;
	  }
	  else {
	    // Subgroup case, base  <> image at current level: Enlarge
	    //   <L>    with  <t>. Decrease <max>     according to the
	    //   enlarged <L>. Reset <R> to the enlarged <L>.
	    //	    for (int dd=0; dd<d; dd++)
	    //	      AddGeneratorsExtendSchreierTree( L[ dd ], {t});
            // It is a little bit unclear why the loop was removed and a single call to
            // AGEST with L_list[dd].
            for (int dd=0; dd<=d; dd++) {
              std::cerr << "CPP Before AGEST dd=" << (dd+1) << "\n";
              AddGeneratorsExtendSchreierTree(L_list[dd], {t.val});
            }
	    if (m < int(L_list[d]->orbit.size())) {
	      std::cerr << "CPP PBEnumerate, EXIT 9\n";
	      return {int_fail,{}};
	    }
	    max = PositionNthTrueBlist( orB[d], m - L_list[d]->orbit.size());
	    SetStabChainFromLevel(R, L, d);
	  }
	}
        std::cerr << "CPP t step 4\n";

	// Now  we can remove the   entire <R>-orbit  of <b> from  the
	// candidate list.
	std::cerr << "CPP ORB 2: Before subtract d=" << (d+1) << " orb[d]=" << GetStringGAP(orb[d]) << "\n";
	if (R_list[d]->transversal[b] != -1)
	  SubtractBlist(orb[d], BlistList(range, R_list[d]->orbit));
	else
	  SubtractBlistOrbitStabChain(orb[d], StrongGeneratorsStabChain(R), b_int);
	std::cerr << "CPP ORB 2: After subtract d=" << (d+1) << " orb[d]=" << GetStringGAP(orb[d]) << "\n";
	b = orb[d].find_next(b);
	std::cerr << "CPP End of the loop. Now b=" << PosFail_to_string(b) << "\n";
      }

    }
    std::cerr << "CPP PBEnumerate, step 11, EXIT 10 |L|=" << L_list.size() << "\n";
    return {int_fail, {}};
  };

  nrback=0; // count the number of times we jumped up

  // Trivial cases first.
  if (IsTrivial(G)) {
    if (!repr)
      return {int_group, G, {}};
    if (Pr(G->comm->identity))
      return {int_perm, {}, G->comm->identity};
    else
      return {int_fail, {}, {}};
  }

  // Construct the <image>.
  image.data=data;
  image.depth=1;
  if (repr) {
    image.partition = data.P;
  }
  else {
    image.partition = rbase.partition;
  }
  if (IsBool(rbase.level2)) {
    image.level2 = {int_false, -777, {}};
  }
  else {
    image.level2 = rbase.level2;
    image.perm2  = {int_perm, G->comm->identity};
  }

  // If  <Pr> is  function,   multiply  permutations. Otherwise, keep   them
  // factorized.
  image.perm = {int_perm, G->comm->identity};
  //  image.level = rbase.chain;

  if (repr) {
    // In the representative case, map the  fixpoints of the partitions at
    // the root of the search tree.
    if (rbase.partition.lengths != image.partition.lengths) {
      image.perm.status = int_false;
    }
    else {
      std::vector<int> fix  = Fixcells(rbase.partition);
      std::vector<int> fixP = Fixcells(image.partition);
      for (int i=0; i<int(fix.size()); i++) {
        std::cerr << "CPP ProcessFixpoint_image, Case PartitionBacktrack 2\n";
	ProcessFixpoint_image(image, fix[i], fixP[i], -1);
      }
    }
    // In   the representative case,   assign  to <L>  and <R>  stabilizer
    // chains.
    //    L := ListStabChain( CopyStabChain( StabChainMutable( L ) ) );
    //    R := ListStabChain( CopyStabChain( StabChainMutable( R ) ) );
  }

  int lenD=rbase.domain[rbase.domain.size()-1];
  for (int i=0; i<=lenD; i++)
    range.push_back(i);
  permPlusBool<Telt> rep = PBEnumerate(0, !repr);
  if (!repr) {
    return {int_group, L, {}};
  }
  else {
    if (rep.status == int_perm)
      return {int_perm, {}, rep.val};
    else
      return {int_fail, {}, {}};
  }
}

bool IsSubset(Face const& f, std::vector<int> const& g)
{
  for (auto & eVal : g)
    if (f[eVal] == 0)
      return false;
  return true;
}


Face Difference_face(Face const& Phi, std::vector<int> const& Omega)
{
  Face fRet = Phi;
  for (auto & eVal : Omega)
    fRet[eVal] = 0;
  return fRet;
}


template<typename Telt>
Face OnSets(Face const& f, Telt const& g)
{
  int n=f.size();
  Face fRet(n);
  boost::dynamic_bitset<>::size_type b = f.find_first();
  while (b != boost::dynamic_bitset<>::npos) {
    int eImg=g.at(b);
    fRet[eImg]=1;
    b = f.find_next(b);
  }
  return fRet;
}


template<typename Telt, typename Tint>
ResultPBT<Telt> RepOpSetsPermGroup(StabChain<Telt> const& G, bool const& repr, Face const& Phi, Face const& Psi)
{
  std::cerr << "CPP Beginning of RepOpSetsPermGroup\n";
  PrintStabChain(G);
  std::cerr << "CPP UseCycle=" << (G->cycles.size() > 0) << "\n";
  std::cerr << "CPP After bool print\n";
  int n=G->comm->n;
  std::vector<int> Omega = MovedPoints(G);
  std::cerr << "CPP n=" << n << " Omega=" << STRING_VectInt(Omega) << "\n";
  if (repr && Phi.size() != Psi.size())
    return {int_fail, {}, {}};
  if (IsSubset(Phi, Omega) || ForAll(Omega, [&](int const &p) -> bool {return !Phi[p];})) {
    if (repr) {
      if (Difference_face(Phi, Omega) != Difference_face(Psi, Omega))
	return {int_fail, {}, {}};
      else
	return {int_perm, {}, G->comm->identity};
    }
    else {
      return {int_group, G, {}};
    }
  }
  else {
    if (repr && (IsSubset(Psi, Omega) || ForAll(Omega, [&](int const& p) -> bool {return !Psi[p];})))
      return {int_fail, {}, {}};
  }
  auto GetPartitionFromPair=[&](Face const& Ph) -> Partition {
    std::vector<int> IntVect;
    std::vector<int> DiffVect;
    for (auto & eVal : Omega) {
      if (Ph[eVal] == 1)
	IntVect.push_back(eVal);
      if (Ph[eVal] == 0)
	DiffVect.push_back(eVal);
    }
    return GetPartition({IntVect, DiffVect});
  };


  Partition P = GetPartitionFromPair(Phi);
  Partition Q = GetPartitionFromPair(Psi);


  auto GetSubgroup=[&](Face const& Ph) -> StabChain<Telt> {
    std::vector<Telt> LGen = StrongGeneratorsStabChain(G);
    std::cerr << "CPP GetSubgroup, |LGen|=" << LGen.size() << "\n";
    std::cerr << "CPP GetSubgroup, LGen=";
    for (auto & eGen : LGen)
      std::cerr << " " << eGen;
    std::cerr << "\n";
    std::vector<Telt> sgs=Filtered(StrongGeneratorsStabChain(G), [&](Telt const& g)->bool{return OnSets(Ph, g) == Ph;});
    std::cerr << "CPP sgs=";
    for (auto & eGen : sgs)
      std::cerr << " " << eGen;
    std::cerr << "\n";
    return MinimalStabChain<Telt,Tint>(sgs, n);
  };


  StabChain<Telt> L = GetSubgroup(Phi);
  StabChain<Telt> R;
  std::cerr << "CPP repr=" << repr << "\n";
  if (repr)
    R = GetSubgroup(Psi);
  else
    R = L;
  std::cerr << "CPP Orders: |R|=" << SizeStabChain<Telt,Tint>(R) << " |L|=" << SizeStabChain<Telt,Tint>(L) << "\n";
  rbaseType<Telt> rbase = EmptyRBase<Telt>({G, G}, true, Omega, P);
  std::cerr << "CPP RepOpSetsPermGroup rbase.level2.status=" << GetIntTypeNature(rbase.level2.status) << "\n";
  std::vector<int> Phi_vect = FaceToVector(Phi);
  std::function<bool(Telt const&)> Pr=[&](Telt const& gen) -> bool {
    for (auto & i : Phi_vect) {
      int iImg=gen.at(i);
      if (Psi[iImg] == 0)
	return false;
    }
    return true;
  };
  std::cerr << "CPP Before call to PartitionBacktrack\n";
  std::cerr << "CPP bool=" << (rbase.level2.Stot->cycles.size() > 0) << "\n";
  std::cerr << "CPP After bool print\n";
  return PartitionBacktrack<Telt,Tint>( G, Pr, repr, rbase, {Q}, L, R );
}

template<typename Telt,typename Tint>
StabChain<Telt> Stabilizer_OnSets(StabChain<Telt> const& G, Face const& Phi)
{
  std::cerr << "CPP Beginning of Stabilizer_OnSets\n";
  bool repr=false;
  return RepOpSetsPermGroup<Telt,Tint>(G, repr, Phi, Phi).stab;
}



template<typename Telt,typename Tint>
std::pair<bool,Telt> RepresentativeAction_OnSets(StabChain<Telt> const& G, Face const& f1, Face const& f2)
{
  std::cerr << "CPP Beginning of RepresentativeAction_OnSets\n";
  bool repr=true;
  ResultPBT<Telt> eRec = RepOpSetsPermGroup<Telt,Tint>(G, repr, f1, f2);
  if (eRec.nature == int_fail)
    return {false, {}};
  return {true, eRec.res};
}




}

#endif
