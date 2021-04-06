#ifndef DEFINE_PERMUTALIB_STBCBCKT_H
#define DEFINE_PERMUTALIB_STBCBCKT_H

#include "StabChainMain.h"
#include "partition.h"
#include "Combinatorics.h"
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

#ifdef SYNCHRONIZED_DEBUG_GAP478
# define DEBUG_STBCBCKT
#endif

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


// The ExtendedT in gap seems to be passing by value for img entries that gets
// modified
template<typename Telt>
permPlusBool<Telt> ExtendedT(Telt const& t, typename Telt::Tidx const& pnt, typename Telt::Tidx img, typename Telt::Tidx const& simg, StabChainPlusLev<Telt> const& S)
{
  using Tidx=typename Telt::Tidx;
#ifdef DEBUG_STBCBCKT
  std::cerr << "CPP ExtendedT sgs(S.Stot)=" << GapStringTVectorB(SortVector(StrongGeneratorsStabChain(S.Stot))) << "\n";
#endif
  if (simg == std::numeric_limits<Tidx>::max())
    img = SlashAct(img, t);
  else
    img = simg;
  // If <G> fixes <pnt>, nothing more can  be changed, so test whether <pnt>
  // = <img>.
  Tidx bpt = BasePoint(S.Stot);
#ifdef DEBUG_STBCBCKT
  std::cerr << "CPP img=" << int(img+1) << " bpt=" << PosFalse_to_string(bpt) << " pnt=" << int(pnt+1) << "\n";
#endif
  if (bpt != pnt) {
#ifdef DEBUG_STBCBCKT
    std::cerr << "CPP Case bpt != pnt\n";
#endif
    if (pnt != img) {
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP ExtendedT, return false 1\n";
#endif
      return {int_false, {}};
    } else {
      return {int_perm, t};
    }
  }
  if (S.Stot->transversal[img] == -1) {
#ifdef DEBUG_STBCBCKT
    std::cerr << "CPP ExtendedT, return false 2\n";
#endif
    return {int_false, {}};
  }
#ifdef DEBUG_STBCBCKT
  std::cerr << "CPP Final case t=" << t << "\n";
  std::cerr << "CPP sgs(S.Stot)=" << GapStringTVector(SortVector(StrongGeneratorsStabChain(S.Stot))) << "\n";
#endif
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


template<typename Tidx>
struct singStrat {
  Tidx p;
  Tidx s;
  Tidx i;
};


template<typename Tidx>
struct Refinement {
public:
  Refinement(Tidx const& val1, Tidx const& val2) {
    nature = 0;
    inputProcessfix = {val1, val2};
  }
  Refinement(Partition<Tidx> const& ePart, std::vector<singStrat<Tidx>> const& strat) {
    nature = 1;
    inputIntersection = {ePart, strat};
  }
  int nature; // 0 for PROCESSFIX, 1 for INTERSECTION
  std::pair<Tidx,Tidx> inputProcessfix;
  std::pair<Partition<Tidx>,std::vector<singStrat<Tidx>>> inputIntersection;
};


template<typename Tidx>
bool IsInsertableRefinement(Refinement<Tidx> const& eRfm)
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
  throw PermutalibException{1};
#endif
  return true;
}


// The underscore nature of a function can be seen in stbcbckt top.
// Since we did not implement all the algorithms of stbcbckt, the value is always 0.
int UnderscoreNature(int const& nature)
{
  return 0;
}




template<typename Tidx>
struct dataType {
  Partition<Tidx>& P;
  dataType(Partition<Tidx>& _P) : P(_P)
  {
  }
  dataType<Tidx> operator=(dataType<Tidx>& data)
  {
    return dataType(data.P);
  }
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
  std::vector<typename Telt::Tidx> domain;
  std::vector<typename Telt::Tidx> base;
  std::vector<typename Telt::Tidx> where;
  //
  StabChain<Telt> chain;
  std::vector<std::vector<typename Telt::Tidx>> fix;
  //
  std::vector<std::vector<Refinement<typename Telt::Tidx>>> rfm;
  Partition<typename Telt::Tidx> partition;
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
    ListKey.emplace_back(GetStringExpressionOfStabChain(x.Stot));
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
          std::cerr << "CPP  KUR: Change of key at i=" << int(i+1) << "\n";
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
  } else {
    if (rbase.level.status == int_stablev) {
      int len=rbase.lev.size();
      std::cerr << str << " |rbase.lev|=" << len << "\n";
      for (int eD=0; eD<len; eD++) {
        PrintStabChain(rbase.lev[eD].Stot);
        std::cerr << "CPP rbase.lev[" << int(eD+1) << "]=" << GapStringTVector(SortVector(StrongGeneratorsStabChain(rbase.lev[eD].Stot))) << "\n";
      }
      std::cerr << "CPP rbase.level=" << GapStringTVector(SortVector(StrongGeneratorsStabChain(rbase.level.Stot))) << "\n";
      std::cerr << str << " PRBL rbase.level, record, |genlabels|=" << rbase.level.Stot->genlabels.size() << "\n";
      std::cerr << str << " PRBL orbit=" << PrintTopOrbit(rbase.level.Stot) << "\n";
    } else {
      std::cerr << str << " PRBL rbase.level=" << GetIntTypeNature(rbase.level.status) << "\n";
    }
  }
}



template<typename Telt>
bool ProcessFixpoint_rbase(rbaseType<Telt> & rbase, typename Telt::Tidx const& pnt)
{
#ifdef DEBUG_STBCBCKT
  std::cerr << "CPP ProcessFixpoint_rbase beginning pnt=" << int(pnt+1) << "\n";
#endif
  if (rbase.level2.status != int_true && rbase.level2.status != int_false) {
#ifdef DEBUG_STBCBCKT
    std::cerr << "CPP Before ChangeStabChain level2\n";
#endif
    ChangeStabChain(rbase.level2.Stot, {pnt}, int_true);
#ifdef DEBUG_STBCBCKT
    PrintRBaseLevel(rbase, "CPP After CSC level2");
    std::cerr << "CPP After ChangeStabChain level2\n";
#endif
    if (BasePoint(rbase.level2) == pnt) {
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP Going to stabilizer of level2\n";
#endif
      rbase.level2.Stot = rbase.level2.Stot->stabilizer;
    }
  }
  if (rbase.level.status == int_int) {
    rbase.level.value_int--;
  } else {
#ifdef DEBUG_STBCBCKT
    std::cerr << "CPP Before ChangeStabChain level\n";
#endif
    ChangeStabChain(rbase.level.Stot, {pnt}, int_true);
#ifdef DEBUG_STBCBCKT
    PrintRBaseLevel(rbase, "CPP After CSC level");
    std::cerr << "CPP After ChangeStabChain level\n";
#endif
    if (BasePoint(rbase.level) == pnt) {
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP Going to stabilizer of level\n";
#endif
      rbase.level.Stot = rbase.level.Stot->stabilizer;
    } else {
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP returning false\n";
#endif
      return false;
    }
  }
  return true;
}



template<typename Telt>
struct imageType {
  int depth;
  Partition<typename Telt::Tidx>& partition;
  permPlusBool<Telt> perm;
  StabChainPlusLev<Telt> level;
  std::vector<typename Telt::Tidx> bimg;
  //
  permPlusBool<Telt> perm2;
  StabChainPlusLev<Telt> level2;
  imageType(Partition<typename Telt::Tidx>& _partition) : partition(_partition)
  {
  }
};


template<typename Telt>
bool ProcessFixpoint_image(imageType<Telt> & image, typename Telt::Tidx const& pnt, typename Telt::Tidx & img, typename Telt::Tidx const& simg)
{
  using Tidx=typename Telt::Tidx;
  if (image.perm.status != int_true) {
#ifdef DEBUG_STBCBCKT
    std::cerr << "CPP PFI sgs(level)=" << GapStringTVector(SortVector(StrongGeneratorsStabChain(image.level.Stot))) << "\n";
    std::cerr << "CPP Case image.perm.status = true\n";
    std::cerr << "CPP Before ExtendedT img=" << int(img+1) << "\n";
#endif
    permPlusBool<Telt> t = ExtendedT(image.perm.val, pnt, img, simg, image.level);
#ifdef DEBUG_STBCBCKT
    std::cerr << "CPP After ExtendedT img=" << int(img+1) << "\n";
#endif
    if (t.status == int_false) {
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP Returning false 1\n";
#endif
      return false;
    } else {
      if (BasePoint(image.level ) == pnt)
        image.level.Stot = image.level.Stot->stabilizer;
    }
    image.perm = t;
  }
  if (image.level2.status != int_false) {
#ifdef DEBUG_STBCBCKT
    std::cerr << "CPP PFI sgs(level2)=" << GapStringTVector(SortVector(StrongGeneratorsStabChain(image.level2.Stot))) << "\n";
    std::cerr << "CPP Case image.perm.status = false\n";
#endif
    permPlusBool<Telt> t = ExtendedT(image.perm2.val, pnt, img, std::numeric_limits<Tidx>::max(), image.level2);
    if (t.status == int_false) {
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP Returning false 2\n";
#endif
      return false;
    } else {
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
#ifdef DEBUG_STBCBCKT
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
  } else {
    std::cerr << "false";
  }
  std::cerr << "\n";
#endif
  //
  if (rbase.level.status == int_int) {
    if (rbase.level.value_int <= 1) {
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP IsTrivialRBase, leaving at case 1 with True\n";
#endif
      return true;
    }
  }
  if (rbase.level.status == int_stablev) {
    if (rbase.level.Stot->genlabels.size() == 0) {
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP IsTrivialRBase, leaving at case 2 with True\n";
#endif
      return true;
    }
  }
#ifdef DEBUG_STBCBCKT
  std::cerr << "CPP IsTrivialRBase, leaving at case 3 with False\n";
#endif
  return false;
}



template<typename Telt>
rbaseType<Telt> EmptyRBase(std::vector<StabChain<Telt>> const& G, bool const& IsId, std::vector<typename Telt::Tidx> const& Omega, Partition<typename Telt::Tidx> const& P)
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
    } else {
      rbase.level2 = {int_stablev, -555, G[1]};
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP rbase Before bool print\n";
      std::cerr << "CPP bool=" << rbase.level2.Stot->IsBoundCycle << "\n";
      std::cerr << "CPP rbase After bool print\n";
#endif
      rbase.lev2 = {};
    }
  } else {
    rbase.level2.status = int_false;
  }
  rbase.chain = CopyStabChain(G[0]);
  rbase.level = {int_stablev, -666, rbase.chain};
  for (auto & pnt : Fixcells(P)) {
#ifdef DEBUG_STBCBCKT
    std::cerr << "CPP Fixcells call ProcessFixpoint_rbase\n";
#endif
    ProcessFixpoint_rbase(rbase, pnt);
  }
  return rbase;
}




template<typename Telt>
bool MeetPartitionStrat(rbaseType<Telt> const& rbase, imageType<Telt> & image, Partition<typename Telt::Tidx> const& S, Telt const& g, std::vector<singStrat<typename Telt::Tidx>> const& strat)
{
  using Tidx=typename Telt::Tidx;
#ifdef DEBUG_STBCBCKT
  std::cerr << "CPP Running MeetPartitionStrat\n";
#endif
  if (strat.size() == 0)
    return false;
  for (auto & pRec : strat) {
#ifdef DEBUG_STBCBCKT
    std::cerr << "CPP ProcessFixpoint_image, Case MeetPartitionStrat\n";
#endif
    if (pRec.p == std::numeric_limits<Tidx>::max()) {
      Tidx eFix = FixpointCellNo(image.partition, pRec.i);
      if (!ProcessFixpoint_image(image, pRec.s, eFix, std::numeric_limits<Tidx>::max()))
        return false;
    }
    if (pRec.p != std::numeric_limits<Tidx>::max() && SplitCell_Partition(image.partition, pRec.p, S, pRec.s, g, pRec.i ) != pRec.i)
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
std::vector<singStrat<typename Telt::Tidx>> StratMeetPartition(rbaseType<Telt> & rbase, Partition<typename Telt::Tidx> & P, Partition<typename Telt::Tidx> const& S, Telt const& g)
{
  using Tidx=typename Telt::Tidx;
  std::vector<singStrat<Tidx>> strat;
  std::vector<Tidx> cellsP = P.cellno;
  if (!g.isIdentity()) {
    for (Tidx i=0; i<NumberCells(P); i++) {
      std::vector<Tidx> cell = Cell( P, i );
      for (auto & eVal : cell) {
        Tidx img=PowAct(eVal, g);
	cellsP[img] = i;
      }
    }
  }
  // If <S> is just a set, it is interpreted as partition ( <S>|<S>^compl ).
  Tidx nrcells = NumberCells(S) - 1;

  for (Tidx s=0; s<nrcells; s++) {
    // now split with cell number s of S.
    std::vector<Tidx> p=Cell(S, s);

    std::vector<int> p2;
    for (auto & eVal : p)
      p2.push_back(cellsP[eVal]);
    //    PrintVectDebug("p2", p2);
    CollectedResult<int> p3=Collected(p2);
    std::vector<Tidx> splits;
    for (int h=0; h<int(p3.LVal.size()); h++) {
      // a cell will split iff it contains more points than are in the s-cell
      if (P.lengths[p3.LVal[h]] > p3.LMult[h])
        splits.push_back(p3.LVal[h]);
    }
    for (auto & pVal : splits) {
      // Last argument true means that the cell will split.
      Tidx i = SplitCell_Partition(P, pVal, S, s, g, std::numeric_limits<Tidx>::max());
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP g=" << g << " i=" << i << "\n";
#endif
      if (!g.isIdentity()) {
	std::vector<Tidx> cell = Cell(P, NumberCells(P));
	for (auto & eVal : cell) {
	  Tidx img=PowAct(eVal, g);
	  cellsP[img] = NumberCells(P);
	}
      }
      strat.push_back({pVal, s, i});
      // If  we have one  or two  new fixpoints, put  them  into the base.
      if (i == 1) {
        Tidx iPart = NumberCells(P) - 1;
#ifdef DEBUG_STBCBCKT
        std::cerr << "CPP NumberCells=" << int(iPart+1) << "\n";
#endif
        Tidx pnt = FixpointCellNo(P, iPart);
#ifdef DEBUG_STBCBCKT
        std::cerr << "CPP FixpointCellNo - NumberCells\n";
#endif
	ProcessFixpoint_rbase(rbase, pnt);
	strat.push_back({std::numeric_limits<Tidx>::max(), pnt, iPart});
	if (IsTrivialRBase(rbase))
	  return strat;
      }
      if (P.lengths[pVal] == 1) {
        Tidx pnt = FixpointCellNo(P, pVal);
#ifdef DEBUG_STBCBCKT
        std::cerr << "CPP FixpointCellNo - pVal\n";
#endif
	ProcessFixpoint_rbase(rbase, pnt);
	strat.push_back({std::numeric_limits<Tidx>::max(), pnt, pVal});
	if (IsTrivialRBase(rbase))
	  return strat;
      }
    }
  }
  return strat;
}

template<typename Telt>
void AddRefinement(rbaseType<Telt> & rbase, int const& pos, Refinement<typename Telt::Tidx> const& eRfm)
{
#ifdef DEBUG_STBCBCKT
  std::cerr << "CPP beginning of AddRefinement\n";
#endif
  if (IsInsertableRefinement(eRfm)) {
#ifdef DEBUG_STBCBCKT
    std::cerr << "CPP Doing RFM insertion\n";
#endif
    rbase.rfm[pos].emplace_back(eRfm);
  }
#ifdef DEBUG_STBCBCKT
  for (size_t i=0; i<rbase.rfm.size(); i++) {
    std::cerr << "CPP i=" << int(i+1) << " |rbase.rfm[i]|=" << rbase.rfm[i].size() << "\n";
  }
#endif
}


template<typename Telt>
void RegisterRBasePoint(Partition<typename Telt::Tidx> & P, rbaseType<Telt> & rbase, int const& pnt, Telt const& TheId)
{
  using Tidx=typename Telt::Tidx;
  if (rbase.level2.status != int_true && rbase.level2.status != int_false) {
#ifdef DEBUG_STBCBCKT
    std::cerr << "CPP Inserting rbase.level2 into rbase.lev2\n";
    std::cerr << "CPP rbase.level2.status=" << GetIntTypeNature(rbase.level2.status) << "\n";
#endif
    rbase.lev2.push_back(rbase.level2);
  }
#ifdef DEBUG_STBCBCKT
  PrintRBaseLevel(rbase, "CPP RegisterRBasePoint 1");
#endif
  rbase.lev.push_back(rbase.level);
  rbase.base.push_back(pnt);
#ifdef DEBUG_STBCBCKT
  KeyUpdatingRbase("RegisterRBasePoint 1", rbase);
#endif
  Tidx k = IsolatePoint<Tidx>(P, pnt);
#ifdef DEBUG_STBCBCKT
  NicePrintPartition("CPP After IsolatePoint P", P);
  KeyUpdatingRbase("RegisterRBasePoint 1.1", rbase);
#endif
  if (!ProcessFixpoint_rbase(rbase, pnt)) {
#ifdef DEBUG_STBCBCKT
    std::cerr << "CPP INFO: Warning R-base point is already fixed\n";
#endif
  }
  //  rbase.lev.push_back(rbase.level);
#ifdef DEBUG_STBCBCKT
  KeyUpdatingRbase("RegisterRBasePoint 1.2", rbase);
  PrintRBaseLevel(rbase, "CPP RegisterRBasePoint 2");
#endif
  rbase.where.push_back(k);
  int len=rbase.rfm.size();
  rbase.rfm.push_back({});
#ifdef DEBUG_STBCBCKT
  std::cerr << "CPP Before P.lengths test k=" << int(k+1) << " len=" << rbase.rfm.size() << "\n";
  KeyUpdatingRbase("RegisterRBasePoint 1.3", rbase);
#endif
  if (P.lengths[k] == 1) {
#ifdef DEBUG_STBCBCKT
    std::cerr << "CPP Matching P.lengths test\n";
#endif
    Tidx pnt = FixpointCellNo(P, k);
#ifdef DEBUG_STBCBCKT
    std::cerr << "CPP Section P.lengths after FixpointCellNo pnt=" << int(pnt+1) << "\n";
    PrintRBaseLevel(rbase, "CPP RegisterRBasePoint 2.1");
#endif
    ProcessFixpoint_rbase(rbase, pnt);
#ifdef DEBUG_STBCBCKT
    PrintRBaseLevel(rbase, "CPP RegisterRBasePoint 2.2");
    KeyUpdatingRbase("RegisterRBasePoint 1.4", rbase);
    std::cerr << "CPP Section P.lengths after ProcessFixpoint_rbase\n";
#endif
    AddRefinement(rbase, len, Refinement<Tidx>({pnt,k}));
#ifdef DEBUG_STBCBCKT
    std::cerr << "CPP After AddRefinement 1\n";
    KeyUpdatingRbase("RegisterRBasePoint 1.5", rbase);
#endif
  }
#ifdef DEBUG_STBCBCKT
  PrintRBaseLevel(rbase, "CPP RegisterRBasePoint 3");
  KeyUpdatingRbase("RegisterRBasePoint 2", rbase);
#endif
  if (rbase.level2.status != int_false) {
#ifdef DEBUG_STBCBCKT
    std::cerr << "CPP Matching the ! false test\n";
#endif
    auto MainInsert=[&](StabChainPlusLev<Telt> const& lev) -> void {
      if (lev.status != int_int) {
	std::vector<Telt> LGenStrong = StrongGeneratorsStabChain(lev.Stot);
#ifdef DEBUG_STBCBCKT
	std::cerr << "CPP StrongGeneratorsStabChain(lev) = " << GapStringTVector(SortVector(LGenStrong)) << "\n";
#endif
        std::vector<Telt> LGen = GeneratorsStabChain(lev.Stot);
	Partition<Tidx> O = OrbitsPartition(LGen, lev.Stot->comm->n, rbase.domain);
#ifdef DEBUG_STBCBCKT
	NicePrintPartition("CPP Before StratMeetPartition O", O);
        KeyUpdatingRbase("RegisterRBasePoint 2.1", rbase);
#endif
	std::vector<singStrat<Tidx>> strat = StratMeetPartition(rbase, P, O, TheId);
#ifdef DEBUG_STBCBCKT
        KeyUpdatingRbase("RegisterRBasePoint 2.2", rbase);
#endif
        AddRefinement(rbase, len, Refinement<Tidx>({O,strat}));
#ifdef DEBUG_STBCBCKT
        std::cerr << "CPP After AddRefinement 2\n";
#endif
      }
    };
    if (rbase.level2.status == int_true) {
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP Before call to MainInsert(level)\n";
#endif
      MainInsert(rbase.level);
    } else {
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP Before call to MainInsert(level2)\n";
#endif
      MainInsert(rbase.level2);
    }
  }
  //  rbase.lev.push_back(rbase.level);
#ifdef DEBUG_STBCBCKT
  KeyUpdatingRbase("RegisterRBasePoint 3", rbase);
#endif
}




template<typename Telt>
void NextRBasePoint(Partition<typename Telt::Tidx> & P, rbaseType<Telt> & rbase, Telt const& TheId)
{
  using Tidx=typename Telt::Tidx;
  std::vector<Tidx> lens = P.lengths;
  std::vector<Tidx> order = ClosedInterval<Tidx>(0, NumberCells(P));
#ifdef DEBUG_STBCBCKT
  std::cerr << "CPP lens=[ ";
  for (size_t i=0; i<lens.size(); i++) {
    if (i>0)
      std::cerr << ", ";
    std::cerr << lens[i];
  }
  std::cerr << " ]\n";
#endif
  SortParallel_PairList(lens, order);


  int k = PositionProperty(lens, [](Tidx const& x) -> bool {return x != 1;});
  int l = -1;
  if (rbase.level.status == int_int) {
    l = 0;
  } else {
    while (true) {
      l = PositionProperty(ClosedInterval<Tidx>(0, lens[k]), [&](Tidx const& i) -> bool {
	  return !IsFixedStabilizer(rbase.level.Stot, P.points[i+P.firsts[order[k]]]);});
      if (l != -1)
	break;
      k++;
    }
  }
  int p = P.points[ P.firsts[ order[k] ] + l ];
#ifdef DEBUG_STBCBCKT
  std::cerr << "CPP p=" << int(p+1) << "\n";
  NicePrintPartition("CPP Before RegisterRBasePoint P", P);
  PrintRBaseLevel(rbase, "CPP Before RegisterRBasePoint");
#endif
  RegisterRBasePoint(P, rbase, p, TheId);
}

template<typename Telt>
bool Refinements_ProcessFixpoint(rbaseType<Telt> & rbase, imageType<Telt> & image, typename Telt::Tidx const& pnt, typename Telt::Tidx const& cellnum)
{
  using Tidx=typename Telt::Tidx;
  Tidx img = FixpointCellNo(image.partition, cellnum);
#ifdef DEBUG_STBCBCKT
  std::cerr << "CPP ProcessFixpoint_image, Case Refinements_ProcessFixpoint\n";
#endif
  return ProcessFixpoint_image(image, pnt, img, std::numeric_limits<Tidx>::max());
}


template<typename Telt>
bool Refinements_Intersection(rbaseType<Telt> & rbase, imageType<Telt> & image, Partition<typename Telt::Tidx> const& Q, std::vector<singStrat<typename Telt::Tidx>> const& strat)
{
  Telt t;
  if (image.level2.status == int_false) {
    t = image.perm.val;
  } else {
    t = image.perm2.val;
  }
  Telt tinv =Inverse(t);
  return MeetPartitionStrat(rbase, image, Q, tinv, strat);
}

// The function RRefine is doing the computation using CallFuncList
// It processes a number of refinement strategies.
// The functions Refinements used return only booleans
//
template<typename Telt>
int RRefine(rbaseType<Telt> & rbase, imageType<Telt> & image, bool const& uscore)
{
  using Tidx=typename Telt::Tidx;
#ifdef DEBUG_STBCBCKT
  std::cerr << "CPP uscore=" << uscore << "\n";
#endif
  auto BoolToInt=[&](bool const& val) -> int {
    if (val)
      return int_true;
    return int_false;
  };
  auto Evaluation=[&](Refinement<Tidx> const& eRef) -> bool {
    if (eRef.nature == int_false)
      return Refinements_ProcessFixpoint(rbase, image, eRef.inputProcessfix.first, eRef.inputProcessfix.second);
    if (eRef.nature == int_true)
      return Refinements_Intersection(rbase, image, eRef.inputIntersection.first, eRef.inputIntersection.second);
    return true;
  };
  if (!uscore) {
#ifdef DEBUG_STBCBCKT
    std::cerr << "CPP case of NOT uscore\n";
#endif
    for (auto & Rf : rbase.rfm[image.depth]) {
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP Doing one CallFuncList 1\n";
#endif
      bool t = Evaluation(Rf);
      if (!t) {
#ifdef DEBUG_STBCBCKT
        std::cerr << "CPP 1 return fail\n";
#endif
	return int_fail;
      } else {
	if (!t) {
#ifdef DEBUG_STBCBCKT
          std::cerr << "CPP 1 return t\n";
#endif
	  return BoolToInt(t);
	}
      }
    }
#ifdef DEBUG_STBCBCKT
    std::cerr << "CPP 1 return true\n";
#endif
    return int_true;
  } else {
#ifdef DEBUG_STBCBCKT
    std::cerr << "CPP case of uscore\n";
#endif
    for (auto & Rf : rbase.rfm[image.depth]) {
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP Doing one CallFuncList 2\n";
#endif
      if (UnderscoreNature(Rf.nature)) {
	bool t = Evaluation(Rf);
	if (!t) {
#ifdef DEBUG_STBCBCKT
          std::cerr << "CPP 2 return fail\n";
#endif
	  return int_fail;
	} else {
	  if (!t) {
#ifdef DEBUG_STBCBCKT
            std::cerr << "CPP 2 return t\n";
#endif
	    return BoolToInt(t);
	  }
	}
      }
    }
#ifdef DEBUG_STBCBCKT
    std::cerr << "CPP 2 return true\n";
#endif
    return int_true;
  }
}


template<typename Telt>
bool PBIsMinimal(std::vector<typename Telt::Tidx> const& range, typename Telt::Tidx const& a, typename Telt::Tidx const& b, StabChain<Telt> const& S)
{
  using Tidx=typename Telt::Tidx;
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

  std::vector<Tidx> orb{b};
  size_t pos=0;
  Face old = BlistList(range, orb);
  while(true) {
    size_t siz=orb.size();
    if (pos == siz)
      break;
    for (size_t i=pos; i<siz; i++) {
      Tidx pnt=orb[i];
      for (auto & lVal : S->genlabels) {
        Tidx img = PowAct(pnt, S->comm->labels[lVal]);
        if (!old[img]) {
          if (img < a)
            return false;
          old[img]=true;
          orb.emplace_back(img);
        }
      }
    }
    pos=siz;
  }
  return true;
}

template<typename Telt>
void SubtractBlistOrbitStabChain(Face & blist, std::vector<Telt> const& LGen, typename Telt::Tidx const& pnt_in)
{
  using Tidx=typename Telt::Tidx;
  std::vector<Tidx> orb{pnt_in};
  blist[pnt_in]=false;
  size_t pos=0;
  while(true) {
    size_t siz = orb.size();
    if (pos == siz) {
      break;
    }
    for (size_t ePos=pos; ePos<siz; ePos++) {
      Tidx pnt=orb[ePos];
      for (auto& eGen : LGen) {
        Tidx img = PowAct(pnt, eGen);
        if (blist[img]) {
          blist[img]=false;
          orb.emplace_back(img);
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
Telt MappingPermListList(typename Telt::Tidx const& n, std::vector<typename Telt::Tidx> const& src, std::vector<typename Telt::Tidx> const& dst)
{
  using Tidx = typename Telt::Tidx;
  std::vector<Tidx> ListImage(n);
  Face StatusSrc(n);
  Face StatusDst(n);
  for (Tidx i=0; i<n; i++) {
    StatusSrc[i]=1;
    StatusDst[i]=1;
  }
  Tidx len = src.size();
  for (Tidx i=0; i<len; i++) {
    ListImage[src[i]] = dst[i];
    StatusSrc[src[i]] = 0;
    StatusDst[dst[i]] = 0;
  }
  Tidx sizRemain = n - len;
  boost::dynamic_bitset<>::size_type posSrc=StatusSrc.find_first();
  boost::dynamic_bitset<>::size_type posDst=StatusDst.find_first();
  for (Tidx u=0; u<sizRemain; u++) {
    ListImage[posSrc] = posDst;
    posSrc=StatusSrc.find_next(posSrc);
    posDst=StatusDst.find_next(posDst);
  }
  return Telt(ListImage);
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

template<typename Telt>
imageType<Telt> BuildInitialImage(bool const&repr, rbaseType<Telt> & rbase, dataType<typename Telt::Tidx> & data)
{
  if (repr) {
    return imageType<Telt>(data.P);
  } else {
    return imageType<Telt>(rbase.partition);
  }
};


template<typename Telt, typename Tint>
ResultPBT<Telt> PartitionBacktrack(StabChain<Telt> const& G, std::function<bool(Telt const&)> const& Pr, bool const& repr, rbaseType<Telt> & rbase, dataType<typename Telt::Tidx> & data, StabChain<Telt> & L, StabChain<Telt> & R)
{
  using Tidx=typename Telt::Tidx;
  Tidx n=G->comm->n;
  Telt id = G->comm->identity;
#ifdef DEBUG_STBCBCKT
  std::cerr << "CPP PartitionBacktrack step 1\n";
  std::cerr << "CPP L=\n";
  PrintStabChain(L);
  std::cerr << "CPP R=\n";
  PrintStabChain(R);
  std::cerr << "CPP INIT sgs(G)=" << GapStringTVector(SortVector(StrongGeneratorsStabChain(G))) << "\n";
  std::cerr << "CPP INIT sgs(L)=" << GapStringTVector(SortVector(StrongGeneratorsStabChain(L))) << "\n";
  std::cerr << "CPP INIT sgs(R)=" << GapStringTVector(SortVector(StrongGeneratorsStabChain(R))) << "\n";
#endif
  imageType<Telt> image = BuildInitialImage(repr, rbase, data);
  std::vector<Face> orB; // backup of <orb>.
  int nrback;
  std::vector<Face> orb;
  std::vector<std::vector<Tidx>> org; // intersected (mapped) basic orbits of <G>
  Tplusinfinity<int> blen(true, 0);
  int dd, branch; // branch is level where $Lstab\ne Rstab$ starts
  std::vector<Tidx> range;    // range for construction of <orb>
  Tidx lenD=rbase.domain[rbase.domain.size()-1];
  for (Tidx i=0; i<=lenD; i++)
    range.push_back(i);
  Partition<Tidx> oldcel;       // old value of <image.partition.cellno>
  std::vector<Tidx> oldcel_cellno;
  std::vector<StabChain<Telt>> L_list, R_list;
  std::function<permPlusBool<Telt>(int const&,bool const&)> PBEnumerate = [&](int const& d, bool const & wasTriv) -> permPlusBool<Telt> {
#ifdef DEBUG_STBCBCKT
    std::cerr << "CPP PBEnumerate, step 1, d=" << int(d+1) << " wasTriv=" << wasTriv << "\n";
#endif
    permPlusBool<Telt> oldprm, oldprm2;
    Tidx a;                // current R-base point
    permPlusBool<Telt> t; // group element constructed, to be handed upwards
    int m;                // initial number of candidates in <orb>
    boost::dynamic_bitset<>::size_type max;              // maximal number of candidates still needed
    boost::dynamic_bitset<>::size_type b;        // image of base point currently being considered

    if (image.perm.status == int_false) {
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP PBEnumerate, EXIT 1, image.perm.status=" << GetIntTypeNature(image.perm.status) << "\n";
#endif
      return {int_fail, {}};
    }
    image.depth = d;
#ifdef DEBUG_STBCBCKT
    std::cerr << "CPP PBEnumerate, step 2\n";
    PrintRBaseLevel(rbase, "CPP Step 2");
#endif

    // Store the original values of <image.*>.
    Tidx undoto = NumberCells(image.partition);
    if (image.perm.status == int_true) {
      oldcel = image.partition;
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP Assigning from image.partition\n";
#endif
    } else {
      oldcel_cellno = image.partition.cellno;
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP Assigning from image.partition.cellno\n";
      std::cerr << "CPP oldcel=" << GapStringIntVector(oldcel_cellno) << "\n";
#endif
      oldprm = image.perm;
    }
#ifdef DEBUG_STBCBCKT
    std::cerr << "CPP PBEnumerate, step 3\n";
#endif
    if (image.level2.status != int_false)
      oldprm2 = image.perm2;
    else
      oldprm2.status = int_false;
#ifdef DEBUG_STBCBCKT
    std::cerr << "CPP PBEnumerate, step 4 d=" << int(d+1) << " |rbase.base|=" << rbase.base.size() << "\n";
    PrintRBaseLevel(rbase, "CPP Step 4");
    if (L_list.size() > 0) {
      std::cerr << "CPP |L|=" << L_list.size() << "\n";
      PrintListStabCommPartition("CPP Step 4", L_list);
    }
    if (R_list.size() > 0) {
      std::cerr << "CPP |R|=" << R_list.size() << "\n";
      PrintListStabCommPartition("CPP Step 4", R_list);
    }
#endif
    // Recursion comes to an end  if all base  points have been prescribed
    // images.
    if (d >= int(rbase.base.size())) {
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP Matching d > Length(rbase.base) test\n";
#endif
      if (IsTrivialRBase(rbase)) {
	blen = rbase.base.size();
#ifdef DEBUG_STBCBCKT
	std::cerr << "CPP IsTrivialRBase matching test blen=" << blen << " wasTriv=" << wasTriv << "\n";
#endif
	// Do     not  add the   identity    element  in the  subgroup
	// construction.
	if (wasTriv) {
#ifdef DEBUG_STBCBCKT
	  std::cerr << "CPP wasTriv Critical, step 1\n";
#endif
	  // In the subgroup case, assign to  <L> and <R> stabilizer
	  // chains when the R-base is complete.
	  StabChainOptions<Tint,Tidx> options = GetStandardOptions<Tint,Tidx>(n);
	  options.base = rbase.base;
	  options.reduced = false;
#ifdef DEBUG_STBCBCKT
	  std::cerr << "CPP Before computation of ListStabChain Order(L)=" << Order<Telt,Tint>(L) << "\n";
          std::cerr << "CPP sgs(L)=" << GapStringTVector(SortVector(StrongGeneratorsStabChain(L))) << " base=" << GapStringIntVector(rbase.base) << "\n";
          std::cerr << "CPP assigning L sequence\n";
#endif
	  L_list = ListStabChain(StabChainOp_stabchain_nofalse<Telt,Tint>(L, options));
#ifdef DEBUG_STBCBCKT
          PrintListStabCommPartition("CPP ListStabChain", L_list);
#endif
	  R_list = L_list; // Corresponds to R := ShallowCopy( L)
#ifdef DEBUG_STBCBCKT
	  std::cerr << "CPP PBEnumerate, EXIT 2\n";
#endif
	  return {int_fail,{}};
	} else {
	  permPlusBool<Telt> prm;
	  if (image.perm.status == int_true)
	    prm = {int_perm, MappingPermListList<Telt>(n, rbase.fix[rbase.base.size()-1], Fixcells(image.partition))};
	  else
	    prm = image.perm;
	  if (image.level2.status != int_false) {
	    if (SiftedPermutation(image.level2.Stot, prm.val * Inverse(image.perm2.val)).isIdentity()) {
#ifdef DEBUG_STBCBCKT
	      std::cerr << "CPP PBEnumerate, EXIT 3\n";
#endif
	      return prm;
	    }
	  } else {
	    if (Pr(prm.val)) {
#ifdef DEBUG_STBCBCKT
	      std::cerr << "CPP PBEnumerate, EXIT 4\n";
#endif
	      return {int_perm, prm.val};
	    }
	  }
#ifdef DEBUG_STBCBCKT
	  std::cerr << "CPP PBEnumerate, EXIT 5\n";
#endif
	  return {int_fail, {}};
	}
	// Construct the   next refinement  level. This  also  initializes
	// <image.partition> for the case ``image = base point''.
      } else {
#ifdef DEBUG_STBCBCKT
	std::cerr << "CPP Not matching IsTrivialRBase test\n";
        PrintRBaseLevel(rbase, "CPP Before NextRBasePoint");
        std::cerr << "CPP Before NextRBasePoint image.p.c=" << GapStringIntVector(image.partition.cellno) << "\n";
#endif
	NextRBasePoint(rbase.partition, rbase, id);
#ifdef DEBUG_STBCBCKT
        std::cerr << "CPP After NextRBasePoint image.p.c=" << GapStringIntVector(image.partition.cellno) << "\n";
	PrintRBaseLevel(rbase, "CPP After NextRBasePoint");
#endif
	if (image.perm.status == int_true)
	  rbase.fix.emplace_back(Fixcells(rbase.partition));
#ifdef DEBUG_STBCBCKT
	std::cerr << "CPP After Fixcells insert\n";
#endif
	std::vector<Tidx> eNewF(range.size(), 0);
	org.emplace_back(eNewF);
	if (repr) {
	  // In  the representative  case,  change  the   stabilizer
	  // chains of <L> and <R>.
#ifdef DEBUG_STBCBCKT
          std::cerr << "CPP Before ChangeStabChain L_list[d]\n";
#endif
	  ChangeStabChain(L_list[d], {rbase.base[d]}, int_false);
#ifdef DEBUG_STBCBCKT
          PrintStabChain(L_list[d]);
          std::cerr << "CPP After ChangeStabChain L_list[d]\n";
#endif
          AssignationVectorGapStyle(L_list, d+1, L_list[ d ]->stabilizer);
#ifdef DEBUG_STBCBCKT
          std::cerr << "CPP Before ChangeStabChain R_list[d]\n";
#endif
	  ChangeStabChain(R_list[d], {rbase.base[d]}, int_false);
#ifdef DEBUG_STBCBCKT
          PrintStabChain(R_list[d]);
          std::cerr << "CPP After ChangeStabChain R_list[d]\n";
#endif
          AssignationVectorGapStyle(R_list, d+1, R_list[ d ]->stabilizer);
#ifdef DEBUG_STBCBCKT
          std::cerr << "CPP L[d]=\n";
          PrintStabChain(L_list[d]);
          std::cerr << "CPP R[d]=\n";
          PrintStabChain(R_list[d]);
          std::cerr << "CPP L[d+1]=\n";
          PrintStabChain(L_list[d+1]);
          std::cerr << "CPP R[d+1]=\n";
          PrintStabChain(R_list[d+1]);
#endif
	}
      }
    }
    a = rbase.base[d];
#ifdef DEBUG_STBCBCKT
    std::cerr << "CPP PBEnumerate, step 5\n";
    PrintRBaseLevel(rbase, "CPP Step 5");
#endif

    // Intersect  the current cell of <P>  with  the mapped basic orbit of
    // <G> (and also with the one of <H> in the intersection case).
    if (image.perm.status == int_true) {
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP orb assign 1: d=" << int(d+1) << "\n";
#endif
      AssignationVectorGapStyle(orb, d, BlistList(range, Cell(oldcel, rbase.where[d]) ));
      if (image.level2.status != int_false) {
	b = orb[d].find_first();
	while (b != boost::dynamic_bitset<>::npos) {
	  if (!IsInBasicOrbit(rbase.lev2[d].Stot, SlashAct(Tidx(b), image.perm2.val)))
	    orb[d][b] = false;
	  b = orb[d].find_next(b);
	}
      }
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP ORB: Case image.perm=true d=" << int(d+1) << " orb[d]=" << GetStringGAP(orb[d]) << "\n";
#endif
    } else {
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP image.perm<>true orb=" << GapStringListBoolVector(orb) << "\n";
      std::cerr << "CPP orb assign 2: d=" << int(d+1) << "\n";
#endif
      AssignationVectorGapStyle(orb, d, BlistList(range, {}));
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP After assignation |orb|=" << orb.size() << "\n";
      // line below needs to be checked.
      std::cerr << "CPP ORB: Before pVal loop d=" << int(d+1) << " orb[d]=" << GetStringGAP(orb[d]) << "\n";
      std::cerr << "CPP RBASE: List(...) = [ ";
      bool IsF=true;
      for (auto & eRec : rbase.lev) {
        if (!IsF)
          std::cerr << ", ";
        IsF=false;
	std::cerr << PrintTopOrbit(eRec.Stot);
      }
      std::cerr << " ]\n";
#endif
      for (auto & pVal : rbase.lev[d].Stot->orbit) {
	b = PowAct(pVal, image.perm.val);
#ifdef DEBUG_STBCBCKT
	std::cerr << "CPP pVal=" << int(pVal+1) << " b=" << int(b+1) << "\n";
        std::cerr << "CPP oldcell=" << int(oldcel_cellno[b]+1) << " rbase.where=" << int(rbase.where[d]+1) << "\n";
#endif
	if (oldcel_cellno[b] == rbase.where[d]) {
	  bool DoOper=false;
	  if (image.level2.status == int_false)
	    DoOper=true;
	  if (!DoOper) {
	    if (IsInBasicOrbit(rbase.lev2[d].Stot, SlashAct(Tidx(b), image.perm2.val)))
	      DoOper=true;
	  }
	  if (DoOper) {
	    orb[d][b] = true;
	    org[d][b] = pVal;
	  }
	}
      }
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP ORB: After pVal loop d=" << int(d+1) << " orb[d]=" << GetStringGAP(orb[d]) << "\n";
#endif
    }
#ifdef DEBUG_STBCBCKT
    std::cerr << "CPP PBEnumerate, step 6 orb=" << GapStringListBoolVector(orb) << "\n";
    PrintRBaseLevel(rbase, "CPP Step 6");
#endif
    if (d == 0 && ForAll(G->comm->labels, [&](Telt const& x){return PowAct(a, x) == a;})) {
      orb[d][a]=true; // ensure a is a possible image (can happen if acting on permutations with more points)
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP ORB: After assignation d=" << int(d+1) << " orb[d]=" << GetStringGAP(orb[d]) << "\n";
#endif
    }
    AssignationVectorGapStyle(orB, d, orb[d]);
#ifdef DEBUG_STBCBCKT
    std::cerr << "CPP PBEnumerate, step 7, wasTriv=" << wasTriv << "\n";
    PrintRBaseLevel(rbase, "CPP Step 7");
#endif

    // Loop  over the candidate images  for the  current base point. First
    // the special case image = base up to current level.
    if (wasTriv) {
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP wasTriv Critical, step 4\n";
#endif
      AssignationVectorGapStyle(image.bimg, d, a);
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP wasTriv Critical, step 5\n";
      // Refinements that start with '_' must be executed even when base
      // = image since they modify image.data, etc.

      std::cerr << "CPP Before RRefine 1 rbase.p.c=" << GapStringIntVector(rbase.partition.cellno) << "\n";
      std::cerr << "CPP Before RRefine 1 image.p.c=" << GapStringIntVector(image.partition.cellno) << "\n";
#endif
      RRefine(rbase, image, true);
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP After RRefine 1 image.p.c=" << GapStringIntVector(image.partition.cellno) << "\n";
      PrintRBaseLevel(rbase, "CPP After RRefine");
#endif
      // Recursion.
      PBEnumerate(d + 1, true);
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP wasTriv Critical, step 6 d=" << int(d+1) << "\n";
#endif
      image.depth = d;
      // Now we  can  remove  the  entire   <R>-orbit of <a>  from   the
      // candidate list.
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP ORB 1: Before subtract d=" << int(d+1) << " orb[d]=" << GetStringGAP(orb[d]) << "\n";
#endif
      SubtractBlist(orb[d], BlistList(range, L_list[d]->orbit));
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP ORB 1: After subtract d=" << int(d+1) << " orb[d]=" << GetStringGAP(orb[d]) << "\n";
#endif
    }
#ifdef DEBUG_STBCBCKT
    std::cerr << "CPP PBEnumerate, step 8\n";
#endif

    // Only the early points of the orbit have to be considered.
    m = SizeBlist(orB[d] );
    if (m < int(L_list[d]->orbit.size()) ) {
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP PBEnumerate, EXIT 6\n";
#endif
      return {int_fail,{}};
    }
    max = PositionNthTrueBlist(orB[d], m - L_list[d]->orbit.size());
#ifdef DEBUG_STBCBCKT
    std::cerr << "CPP PBEnumerate, step 9\n";
    std::cerr << "CPP wasTriv=" << wasTriv << " a=" << int(a+1) << " max=" << int(max+1) << "\n";
#endif

    boost::dynamic_bitset<>::size_type a_size = a;
    if (wasTriv && a_size > max) {
      m--;
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP Before test m=" << m << " Length(L[d].orbit)=" << L_list[d]->orbit.size() << "\n";
#endif
      if (m < int(L_list[d]->orbit.size()) ) {
#ifdef DEBUG_STBCBCKT
	std::cerr << "CPP PBEnumerate, EXIT 7\n";
#endif
	return {int_fail,{}};
      }
      max = PositionNthTrueBlist( orB[d], m - L_list[d]->orbit.size());
    }
#ifdef DEBUG_STBCBCKT
    std::cerr << "CPP PBEnumerate, step 10\n";
#endif
    // Now the other possible images.
    b = orb[d].find_first();
    while (b != boost::dynamic_bitset<>::npos) {
      Tidx b_int = Tidx(b);
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP b=" << int(b+1) << " b_int=" << int(b_int+1) << " d=" << int(d+1) << "\n";
      std::cerr << "CPP |R[d].orbit|=" << R_list[d]->orbit.size() << "\n";
#endif
      // Try to prune the node with prop 8(ii) of Leon paper.
      if (!repr && !wasTriv && R_list[d]->orbit.size() > 0) {
#ifdef DEBUG_STBCBCKT
        std::cerr << "CPP matching if test\n";
#endif
	dd = branch;
	while (dd < d) {
#ifdef DEBUG_STBCBCKT
          std::cerr << "CPP while dd=" << int(dd+1) << " d=" << int(d+1) << "\n";
#endif
	  if (IsInBasicOrbit(L_list[dd], a) && !PBIsMinimal(range, R_list[dd]->orbit[0], b_int, R_list[d])) {
	    dd = d + 1;
#ifdef DEBUG_STBCBCKT
            std::cerr << "CPP first case\n";
#endif
          } else {
	    dd = dd + 1;
#ifdef DEBUG_STBCBCKT
            std::cerr << "CPP second case\n";
#endif
          }
	}
      } else {
	dd = d;
      }
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP dd=" << int(dd+1) << " d=" << int(d+1) << "\n";
      std::cerr << "CPP L[d]=\n";
      PrintStabChain(L_list[d]);
      std::cerr << "CPP R[d]=\n";
      PrintStabChain(R_list[d]);
#endif
      if (dd == d) {
#ifdef DEBUG_STBCBCKT
        std::cerr << "CPP equality dd=d undoto=" << undoto << " |image.partition|=" << NumberCells(image.partition) << "\n";
#endif
	// Undo the  changes made to  <image.partition>, <image.level>
	// and <image.perm>.
        Tidx nbCell=NumberCells(image.partition);
	for (Tidx i=undoto+1; i<=nbCell; i++) {
#ifdef DEBUG_STBCBCKT
          std::cerr << "CPP Before UndoRefinement cellno=" << GapStringIntVector(image.partition.cellno) << " i=" << i << "\n";
#endif
	  UndoRefinement(image.partition);
#ifdef DEBUG_STBCBCKT
          std::cerr << "CPP After UndoRefinement cellno=" << GapStringIntVector(image.partition.cellno) << " i=" << i << "\n";
#endif
        }
#ifdef DEBUG_STBCBCKT
        std::cerr << "CPP After UndoRefinement loop\n";
#endif
	if (image.perm.status != int_true) {
#ifdef DEBUG_STBCBCKT
          std::cerr << "CPP assignation image.level\n";
#endif
	  image.level = rbase.lev[d];
	  image.perm = oldprm;
	}
	if (image.level2.status != int_false) {
#ifdef DEBUG_STBCBCKT
          std::cerr << "CPP assignation image.level2\n";
#endif
	  image.level2 = rbase.lev2[d];
	  image.perm2  = oldprm2;
	}
	// If <b> could not be prescribed as image for  <a>, or if the
	// refinement was impossible, give up for this image.
#ifdef DEBUG_STBCBCKT
        std::cerr << "CPP Before AssignationVectorGapStyle b_int=" << int(b_int+1) << "\n";
#endif
	AssignationVectorGapStyle(image.bimg, d, b_int);
#ifdef DEBUG_STBCBCKT
        std::cerr << "CPP Before IsolatePoint b_int=" << int(b_int+1) << "\n";
        std::cerr << "CPP Before IsolatePoint cellno=" << GapStringIntVector(image.partition.cellno) << "\n";
#endif
	IsolatePoint<Tidx>( image.partition, b_int );
#ifdef DEBUG_STBCBCKT
        std::cerr << "CPP After IsolatePoint cellno=" << GapStringIntVector(image.partition.cellno) << "\n";
        std::cerr << "CPP ProcessFixpoint_image, Case PartitionBacktrack 1\n";
        std::cerr << "CPP Before ProcessFixpoint_image b_int=" << int(b_int+1) << "\n";
#endif
	bool val = ProcessFixpoint_image(image, a, b_int, org[d][b_int]);
#ifdef DEBUG_STBCBCKT
	std::cerr << "CPP a=" << int(a+1) << " b=" << int(b_int+1) << " org[d][b]=" << int(org[d][b_int]+1) << " val=" << val << "\n";
#endif
	if (val) {
#ifdef DEBUG_STBCBCKT
          std::cerr << "CPP Before RRefine 2 oldcel=" << GapStringIntVector(image.partition.cellno) << "\n";
#endif
	  t.status = RRefine(rbase, image, false);
#ifdef DEBUG_STBCBCKT
          std::cerr << "CPP After RRefine 2 oldcel=" << GapStringIntVector(image.partition.cellno) << "\n";
#endif
        } else {
#ifdef DEBUG_STBCBCKT
          std::cerr << "CPP assignation of t to fail\n";
#endif
	  t.status = int_fail;
        }
#ifdef DEBUG_STBCBCKT
        std::cerr << "CPP After assignment of t. t.status=" << GapStringTrueFalseFail(t.status) << "\n";
#endif

	if (t.status != int_fail) {
#ifdef DEBUG_STBCBCKT
          std::cerr << "CPP case of not fail\n";
#endif
	  // Subgroup case, base <> image   at current level:   <R>,
	  //   which until now is identical to  <L>, must be changed
	  //   without affecting <L>, so take a copy.
#ifdef DEBUG_STBCBCKT
          std::cerr << "CPP wasTriv=" << wasTriv << " d=" << int(d+1) << "\n";
          std::cerr << "CPP L[d]=\n";
          PrintStabChain(L_list[d]);
          std::cerr << "CPP R[d]=\n";
          PrintStabChain(R_list[d]);
          std::cerr << "CPP IsIdenticalObj=" << IsIdenticalObj(L_list[d], R_list[d]) << "\n";
#endif
	  if (wasTriv && IsIdenticalObj(L_list[d], R_list[d])) {
#ifdef DEBUG_STBCBCKT
            std::cerr << "CPP Assigning R from d\n";
#endif
	    SetStabChainFromLevel(R_list, L_list, d, rbase.base.size());
	    branch = d;
#ifdef DEBUG_STBCBCKT
            std::cerr << "CPP assignation branch=" << int(branch+1) << "\n";
#endif
	  }
#ifdef DEBUG_STBCBCKT
          std::cerr << "CPP After wasTriv test\n";
          std::cerr << "CPP d=" << int(d+1) << " blen=" << blen << "\n";
#endif
	  if (2 * (d+1) <= blen) {
#ifdef DEBUG_STBCBCKT
            std::cerr << "CPP Before ChangeStabChain R_list[d] 2\n";
            std::cerr << "CPP XXX Before ChangeStabChain R[d]=\n";
            PrintStabChain(R_list[d]);
#endif
	    ChangeStabChain(R_list[d], {b_int}, int_false);
#ifdef DEBUG_STBCBCKT
            std::cerr << "CPP XXX After ChangeStabChain R[d]=\n";
            PrintStabChain(R_list[d]);
            std::cerr << "CPP After ChangeStabChain R_list[d] 2\n";
#endif
	    R_list[ d + 1 ] = R_list[ d ]->stabilizer;
	  } else {
#ifdef DEBUG_STBCBCKT
            std::cerr << "CPP Beginning else case\n";
#endif
	    std::vector<Telt> LGen = StrongGeneratorsStabChain( R_list[d] );
#ifdef DEBUG_STBCBCKT
            std::cerr << "CPP LGen=" << GapStringTVector(LGen) << "\n";
            std::cerr << "CPP First generating step done b=" << int(b_int+1) << "\n";
#endif
	    std::vector<Telt> LGenB = Filtered(LGen, [&](Telt const& gen) -> bool {return PowAct(b_int, gen) == b_int;});
#ifdef DEBUG_STBCBCKT
            std::cerr << "CPP LGenB=" << GapStringTVector(LGenB) << "\n";
            std::cerr << "XXX ELIMINATE begin\n";
#endif
	    R_list[d+1] = StabChainGenerators(LGenB, n, id);
#ifdef DEBUG_STBCBCKT
            std::cerr << "XXX ELIMINATE end\n";
            std::cerr << "CPP After assignation R[d+1]=\n";
            PrintStabChainOrbits(R_list[d+1]);
#endif
	  }
#ifdef DEBUG_STBCBCKT
          std::cerr << "CPP R[d+1]=\n";
          PrintStabChainOrbits(R_list[d+1]);
#endif
	}
#ifdef DEBUG_STBCBCKT
        std::cerr << "CPP t step 2\n";
#endif

	// Recursion.
	if (t.status == int_true) {
	  t = PBEnumerate(d + 1, false);
#ifdef DEBUG_STBCBCKT
          std::cerr << "CPP After PBEnumerate Recursion case\n";
#endif
	  nrback++;
	  image.depth = d;
	}
#ifdef DEBUG_STBCBCKT
        std::cerr << "CPP t step 3\n";
#endif

	// If   <t>   =   fail, either   the   recursive   call  was
	//   unsuccessful,  or all new  elements   have been added  to
	//   levels  below  the current one   (this happens if  base =
	//   image up to current level).
	if (t.status != int_fail) {
#ifdef DEBUG_STBCBCKT
          std::cerr << "CPP Matching t<>fail\n";
#endif
	  // Representative case, element found: Return it.
	  // Subgroup case, base <> image  before current level:  We
	  //   need  only find  a representative  because we already
	  //   know the stabilizer of <L> at an earlier level.
	  if (repr || !wasTriv) {
#ifdef DEBUG_STBCBCKT
	    std::cerr << "CPP PBEnumerate, EXIT 8\n";
#endif
	    return t;
	  } else {
	    // Subgroup case, base  <> image at current level: Enlarge
	    //   <L>    with  <t>. Decrease <max>     according to the
	    //   enlarged <L>. Reset <R> to the enlarged <L>.
	    //	    for (int dd=0; dd<d; dd++)
	    //	      AddGeneratorsExtendSchreierTree( L[ dd ], {t});
            // It is a little bit unclear why the loop was removed and a single call to
            // AGEST with L_list[dd].
#ifdef DEBUG_STBCBCKT
            PrintListStabCommPartition("CPP AddGen", L_list);
#endif
            for (int dd=0; dd<=d; dd++) {
#ifdef DEBUG_STBCBCKT
              std::cerr << "CPP Before AGEST dd=" << int(dd+1) << "\n";
#endif
              AddGeneratorsExtendSchreierTree(L_list[dd], {t.val});
            }
	    if (m < int(L_list[d]->orbit.size())) {
#ifdef DEBUG_STBCBCKT
	      std::cerr << "CPP PBEnumerate, EXIT 9\n";
#endif
	      return {int_fail,{}};
	    }
	    max = PositionNthTrueBlist( orB[d], m - L_list[d]->orbit.size());
	    SetStabChainFromLevel(R_list, L_list, d, rbase.base.size());
            //            PrintListStabCommPartition(R_list);
#ifdef DEBUG_STBCBCKT
            PrintListStabCommPartition("CPP SetStab", L_list);
#endif
	  }
	}
#ifdef DEBUG_STBCBCKT
        std::cerr << "CPP t step 4\n";
#endif

	// Now  we can remove the   entire <R>-orbit  of <b> from  the
	// candidate list.
#ifdef DEBUG_STBCBCKT
	std::cerr << "CPP ORB 2: Before subtract d=" << int(d+1) << " orb[d]=" << GetStringGAP(orb[d]) << "\n";
#endif
	if (R_list[d]->transversal.size() > 0 && R_list[d]->transversal[b] != -1) {
#ifdef DEBUG_STBCBCKT
          std::cerr << "CPP subtract case 1\n";
#endif
	  SubtractBlist(orb[d], BlistList(range, R_list[d]->orbit));
        } else {
#ifdef DEBUG_STBCBCKT
          std::cerr << "CPP subtract case 2\n";
#endif
	  SubtractBlistOrbitStabChain(orb[d], StrongGeneratorsStabChain(R_list[d]), b_int);
        }
#ifdef DEBUG_STBCBCKT
	std::cerr << "CPP ORB 2: After subtract d=" << int(d+1) << " orb[d]=" << GetStringGAP(orb[d]) << "\n";
#endif
      }
      b = orb[d].find_next(b);
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP End of the loop. 1 Now b=" << PosFail_to_string(b) << "\n";
#endif
      if (b > max)
        b = boost::dynamic_bitset<>::npos;
#ifdef DEBUG_STBCBCKT
      std::cerr << "CPP End of the loop. 2 Now b=" << PosFail_to_string(b) << "\n";
#endif

    }
#ifdef DEBUG_STBCBCKT
    std::cerr << "CPP PBEnumerate, step 11, EXIT 10\n";
#endif
    return {int_fail, {}};
  };

  nrback=0; // count the number of times we jumped up

  // Trivial cases first.
  if (IsTrivial(G)) {
    if (!repr)
      return {int_group, G, {}};
    if (Pr(id))
      return {int_perm, {}, id};
    else
      return {int_fail, {}, {}};
  }

  // Construct the <image>.
  image.depth=1;
  if (repr) {
    image.partition = data.P;
  } else {
    image.partition = rbase.partition;
  }
  if (IsBool(rbase.level2)) {
    image.level2 = {int_false, -777, {}};
  } else {
    image.level2 = rbase.level2;
    image.perm2  = {int_perm, id};
  }

  // If  <Pr> is  function,   multiply  permutations. Otherwise, keep   them
  // factorized.
  image.perm = {int_perm, id};
  image.level = {int_stablev, -666, rbase.chain};

  if (repr) {
    // In the representative case, map the  fixpoints of the partitions at
    // the root of the search tree.
    if (rbase.partition.lengths != image.partition.lengths) {
      image.perm.status = int_false;
    } else {
      std::vector<Tidx> fix  = Fixcells(rbase.partition);
      std::vector<Tidx> fixP = Fixcells(image.partition);
      for (Tidx i=0; i<Tidx(fix.size()); i++) {
#ifdef DEBUG_STBCBCKT
        std::cerr << "CPP ProcessFixpoint_image, Case PartitionBacktrack 2 i=" << int(i+1) << " fix=" << int(fix[i]+1) << " fixP=" << int(fixP[i]+1) << "\n";
#endif
	ProcessFixpoint_image(image, fix[i], fixP[i], std::numeric_limits<Tidx>::max());
      }
    }
    L_list = ListStabChain(StructuralCopy(L));
    R_list = ListStabChain(StructuralCopy(R));
  }

  permPlusBool<Telt> rep = PBEnumerate(0, !repr);
#ifdef DEBUG_STBCBCKT
  std::cerr << "CPP After PBEnumerate Call 0, repr\n";
#endif
  if (!repr) {
    return {int_group, L_list[0], {}};
  } else {
    if (rep.status == int_perm)
      return {int_perm, {}, rep.val};
    else
      return {int_fail, {}, {}};
  }
}

template<typename Tidx>
bool IsSubset(Face const& f, std::vector<Tidx> const& g)
{
  for (auto & eVal : g)
    if (f[eVal] == 0)
      return false;
  return true;
}


template<typename Tidx>
Face Difference_face(Face const& Phi, std::vector<Tidx> const& Omega)
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
  using Tidx=typename Telt::Tidx;
#ifdef DEBUG_STBCBCKT
  std::cerr << "CPP Beginning of RepOpSetsPermGroup\n";
  PrintStabChain(G);
#endif
  int n=G->comm->n;
  std::vector<Tidx> Omega = MovedPoints(G);
  //  int n_phi = Phi.size();
  //  std::cerr << "n_phi=" << n_phi << " n=" << n << "\n";
#ifdef DEBUG_STBCBCKT
  std::cerr << "CPP Omega=" << GapStringIntVector(Omega) << "\n";
#endif
  if (repr && Phi.size() != Psi.size())
    return {int_fail, {}, {}};
  if (IsSubset(Phi, Omega) || ForAll(Omega, [&](int const &p) -> bool {return !Phi[p];})) {
    if (repr) {
      if (Difference_face(Phi, Omega) != Difference_face(Psi, Omega))
	return {int_fail, {}, {}};
      else
	return {int_perm, {}, G->comm->identity};
    } else {
      return {int_group, G, {}};
    }
  } else {
    if (repr && (IsSubset(Psi, Omega) || ForAll(Omega, [&](int const& p) -> bool {return !Psi[p];})))
      return {int_fail, {}, {}};
  }
  auto GetPartitionFromPair=[&](Face const& Ph) -> Partition<Tidx> {
    std::vector<Tidx> IntVect;
    std::vector<Tidx> DiffVect;
    for (auto & eVal : Omega) {
      if (Ph[eVal] == 1)
	IntVect.push_back(eVal);
      if (Ph[eVal] == 0)
	DiffVect.push_back(eVal);
    }
    //    std::cerr << "CPP IntVect=" << GapStringIntVector(IntVect) << " DiffVect=" << GapStringIntVector(DiffVect) << "\n";
    return GetPartition<Tidx>({IntVect, DiffVect});
  };


  Partition<Tidx> P = GetPartitionFromPair(Phi);
  Partition<Tidx> Q;
  if (repr) {
    Q = GetPartitionFromPair(Psi);
  } else {
    Q = P;
  }


  std::vector<Telt> LGen = StrongGeneratorsStabChain(G);
#ifdef DEBUG_STBCBCKT
  std::cerr << "CPP G : LGen=" << GapStringTVector(SortVector(LGen)) << "\n";
  std::cerr << "CPP repr=" << repr << "\n";
#endif

  auto GetSubgroup=[&](Face const& Ph) -> StabChain<Telt> {
    std::vector<Telt> sgs=Filtered(LGen, [&](Telt const& g)->bool{return OnSets(Ph, g) == Ph;});
#ifdef DEBUG_STBCBCKT
    std::cerr << "CPP SelectedGens=" << GapStringTVector(SortVector(sgs)) << "\n";
#endif
    return MinimalStabChain<Telt,Tint>(SortVector(sgs), n);
  };


  StabChain<Telt> L = GetSubgroup(Phi);
  StabChain<Telt> R;
  if (repr)
    R = GetSubgroup(Psi);
  else
    R = L;
#ifdef DEBUG_STBCBCKT
  std::cerr << "CPP Orders: |R|=" << SizeStabChain<Telt,Tint>(R) << " |L|=" << SizeStabChain<Telt,Tint>(L) << "\n";
#endif
  rbaseType<Telt> rbase = EmptyRBase<Telt>({G, G}, true, Omega, P);
  //#ifdef DEBUG_STBCBCKT
  //  std::cerr << "CPP RepOpSetsPermGroup rbase.level2.status=" << GetIntTypeNature(rbase.level2.status) << "\n";
  //#endif
  std::vector<int> Phi_vect = FaceToVector(Phi);
  std::function<bool(Telt const&)> Pr=[&](Telt const& gen) -> bool {
    for (auto & i : Phi_vect) {
      int iImg=gen.at(i);
      if (Psi[iImg] == 0)
	return false;
    }
    return true;
  };
#ifdef DEBUG_STBCBCKT
  std::cerr << "CPP Before call to PartitionBacktrack\n";
#endif
  dataType<Tidx> data(Q);
  return PartitionBacktrack<Telt,Tint>( G, Pr, repr, rbase, data, L, R );
}



template<typename Telt,typename Tint>
StabChain<Telt> Kernel_Stabilizer_OnSets(StabChain<Telt> const& G, Face const& Phi)
{
#ifdef DEBUG_STBCBCKT
  std::cerr << "CPP Beginning of Stabilizer_OnSets\n";
#endif
  size_t n = G->comm->n;
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
  if (Phi.size() != n) {
    std::cerr << "We should have Phi of size equal to n\n";
    throw PermutalibException{1};
  }
#endif
  bool repr=false;
  if (2 * Phi.count() > n) {
    Face PhiC(n);
    for (size_t i=0; i<n; i++)
      PhiC[i] = 1 - Phi[i];
    return RepOpSetsPermGroup<Telt,Tint>(G, repr, PhiC, PhiC).stab;
  } else {
    return RepOpSetsPermGroup<Telt,Tint>(G, repr, Phi, Phi).stab;
  }
}



template<typename Telt,typename Tint>
StabChain<Telt> Kernel_Stabilizer_OnPoints(StabChain<Telt> const& G, typename Telt::Tidx const& x)
{
  using Tidx=typename Telt::Tidx;
#ifdef DEBUG_STBCBCKT
  std::cerr << "CPP Beginning of Stabilizer_OnSets\n";
#endif
  Tidx n = G->comm->n;
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
  if (x >= n) {
    std::cerr << "We should have x < n\n";
    throw PermutalibException{1};
  }
#endif
  Face Phi(n);
  Phi[x]=1;
  bool repr=false;
  return RepOpSetsPermGroup<Telt,Tint>(G, repr, Phi, Phi).stab;
}



template<typename Telt,typename Tint>
std::pair<bool,Telt> Kernel_RepresentativeAction_OnSets(StabChain<Telt> const& G, Face const& f1, Face const& f2)
{
  size_t n = G->comm->n;
#ifdef DEBUG_STBCBCKT
  std::cerr << "CPP Beginning of RepresentativeAction_OnSets\n";
#endif
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
  if (f1.size() != n || f2.size() != n) {
    std::cerr << "We should have f1 and f2 of size equal to n\n";
    throw PermutalibException{1};
  }
#endif
  if (f1.count() != f2.count())
    return {false, {}};
  if (f1 == f2)
    return {true, G->comm->identity};
  auto Process_ResultPBT=[&](ResultPBT<Telt> const& eRec) -> std::pair<bool,Telt> {
    if (eRec.nature == int_fail)
      return {false, {}};
    return {true, eRec.res};
  };
  bool repr=true;
  // Put the false for debugging.
  if (2 * f1.count() > n && false) {
    Face f1C(n), f2C(n);
    for (size_t i=0; i<n; i++) {
      f1C[i] = 1 - f1[i];
      f2C[i] = 1 - f2[i];
    }
    ResultPBT<Telt> eRec = RepOpSetsPermGroup<Telt,Tint>(G, repr, f1C, f2C);
    return Process_ResultPBT(eRec);
  } else {
    ResultPBT<Telt> eRec = RepOpSetsPermGroup<Telt,Tint>(G, repr, f1, f2);
    return Process_ResultPBT(eRec);
  }
}

template<typename Telt,typename Tint>
std::pair<bool,Telt> Kernel_RepresentativeAction_OnPoints(StabChain<Telt> const& G, typename Telt::Tidx const& x1, typename Telt::Tidx const& x2)
{
  using Tidx=typename Telt::Tidx;
  Tidx n = G->comm->n;
#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
  if (x1 >= n || x2 >= n) {
    std::cerr << "We should have x1 <n && x2 < n. x1=" << x1 << " x2=" << x2 << " n=" << n << "\n";
    throw PermutalibException{1};
  }
#endif
  bool repr=true;
  Face f1(n), f2(n);
  f1[x1] = 1;
  f2[x2] = 1;
  ResultPBT<Telt> eRec = RepOpSetsPermGroup<Telt,Tint>(G, repr, f1, f2);
  if (eRec.nature == int_fail)
    return {false, {}};
  return {true, eRec.res};
}


}

#endif
