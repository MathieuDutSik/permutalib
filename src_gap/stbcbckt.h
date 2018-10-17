#ifndef DEFINE_STBCBCKT_H
#define DEFINE_STBCBCKT_H

#include "StabChainMain.h"
#include "partition.h"
#include "COMB_Combinatorics.h"
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
  int eLev;
};


 
template<typename Telt>
permPlusBool<Telt> ExtendedT(Telt const& t, int const& pnt, int& img, int const& simg, StabChainPlusLev<Telt> const& S)
{
  if (simg == -1)
    img = SlashAct(img, t);
  else
    img = simg;
    
  // If <G> fixes <pnt>, nothing more can  be changed, so test whether <pnt>
  // = <img>.
  int bpt = BasePoint(S.Stot, S.eLev);
  if (bpt != pnt) {
    if (pnt != img)
      return {int_false,{}};
    else
      return {int_perm, t};
  }
  if (S.Stot.stabilizer[S.eLev].transversal[img] == -1) {
    return {int_false,{}};
  }
  //      Telt u = InverseRepresentative(S.Stot, S.eLev, img);
  //      t = LeftQuotient(u, t);
  return {int_perm, LeftQuotient(InverseRepresentative(S.Stot, S.eLev, img), t)};
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
  return BasePoint(S.Stot, S.eLev);
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
};


template<typename Telt>
bool ProcessFixpoint_rbase(rbaseType<Telt> & rbase, int const& pnt)
{
  if (rbase.level2.status != int_true && rbase.level2.status != int_false) {
    ChangeStabChain(rbase.level2.Stot, rbase.level2.eLev, {pnt}, int_true);
    if (BasePoint(rbase.level2) == pnt) 
      rbase.level2.eLev++;
  }
  if (rbase.level.status == int_int) {
    rbase.level.value_int--;
  }
  else {
    ChangeStabChain(rbase.level.Stot, rbase.level.eLev, {pnt}, int_true);
    if (BasePoint(rbase.level) == pnt) {
      rbase.level.eLev++;
    }
    else {
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
    permPlusBool<Telt> t = ExtendedT(image.perm.val, pnt, img, simg, image.level);
    if (t.status == int_false)
      return false;
    else {
      if (BasePoint(image.level ) == pnt)
        image.level.eLev++;
    }
    image.perm = t;
  }
  if (image.level2.status != int_false) {
    permPlusBool<Telt> t = ExtendedT(image.perm2.val, pnt, img, -1, image.level2);
    if (t.status == int_false)
      return false;
    else {
      if (BasePoint(image.level2 ) == pnt)
        image.level2.eLev++;
    }
    image.perm2 = t;
  }
  return true;
}


template<typename Telt>
bool IsTrivialRBase(rbaseType<Telt> const& rbase)
{
  if (rbase.level.status == int_int) {
    if (rbase.level.value_int <= 1)
      return true;
  }
  if (rbase.level.status == int_stablev) {
    int eLev=rbase.level.eLev;
    if (rbase.level.Stot.stabilizer[eLev].genlabels.size() == 0)
      return true;
  }
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
    if (IsId) {
      rbase.level2.status=int_true;
      rbase.level2.Stot.UseCycle = false;
    }
    else {
      rbase.level2 = {int_stablev, -444, G[1], 0};
      std::cerr << "rbase Before bool print\n";
      std::cerr << "bool=" << rbase.level2.Stot.UseCycle << "\n";
      std::cerr << "rbase After bool print\n";
      rbase.lev2 = {};
    }
  }
  else {
    rbase.level2.status=int_false;
  }
  rbase.level = {int_stablev, -444, G[0], 0};
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
  std::vector<singStrat> strat;
  std::vector<int> cellsP;
  if (g.isIdentity()) {
    cellsP = P.cellno;
  }
  else {
    cellsP = std::vector<int>(P.cellno.size(), 0);
    for (int i=0; i<NumberCells(P); i++) {
      std::vector<int> cell = Cell( P, i );
      for (auto & eVal : cell) {
        int img=PowAct(eVal, g);
	cellsP[img] = i;
      }
    }
  }
  // If <S> is just a set, it is interpreted as partition ( <S>|<S>^compl ).
  int nrcells = NumberCells(S) - 1;

  for (int s=0; s<nrcells; s++) {
    // now split with cell number s of S.
    std::vector<int> p=Cell(S, s);
    std::vector<int> p2;
    for (auto & eVal : p)
      p2.push_back(cellsP[eVal]);
    CollectedResult<int> p3=Collected(p);
    std::vector<int> splits;
    for (int h=0; h<int(p3.LVal.size()); h++) {
      // a cell will split iff it contains more points than are in the s-cell
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
void RegisterRBasePoint(Partition & P, rbaseType<Telt> & rbase, int const& pnt, Telt const& TheId)
{
  if (rbase.level2.status != int_true && rbase.level2.status != int_false) {
    std::cerr << "Inserting rbase.level2 into rbase.lev2\n";
    std::cerr << "rbase.level2.status=" << GetIntTypeNature(rbase.level2.status) << "\n";
    rbase.lev2.push_back(rbase.level2);
  }
  rbase.lev.push_back(rbase.level);
  rbase.base.push_back(pnt);
  int k = IsolatePoint(P, pnt);
  rbase.where.push_back(k);
  int len=rbase.rfm.size();
  rbase.rfm.push_back({});
  if (P.lengths[k] == 1) {
    int pnt = FixpointCellNo(P, k);
    ProcessFixpoint_rbase(rbase, pnt);
    rbase.rfm[len].push_back(Refinement({pnt,k}));
  }
  if (rbase.level2.status != int_false) {
    auto MainInsert=[&](StabChainPlusLev<Telt> const& lev) -> void {
      if (lev.status != int_int) {
	Partition O = OrbitsPartition(StrongGeneratorsStabChain(lev.Stot, lev.eLev), rbase.domain);
	std::vector<singStrat> strat = StratMeetPartition(rbase, P, O, TheId);
	rbase.rfm[len].push_back(Refinement({O,strat}));
      }
    };
    if (rbase.level2.status == int_true)
      MainInsert(rbase.level);
    else
      MainInsert(rbase.level2);
  }
}




template<typename Telt>
void NextRBasePoint(Partition & P, rbaseType<Telt> & rbase, Telt const& TheId)
{
  std::cerr << "Working with NextRBasePoint rbase.level2.status=" << GetIntTypeNature(rbase.level2.status) << "\n";
  std::vector<int> lens = P.lengths;
  std::vector<int> order = ClosedInterval(0, NumberCells(P));
  SortParallel(lens, order);
  int k = PositionProperty(lens, [](int const& x) -> int {return x != 1;});
  int l = -1;
  while(l == -1) {
    if (rbase.level.status == int_int) {
      l = 0;
    }
    else {
      l = PositionProperty(ClosedInterval(0, lens[k]), [&](int const& i) -> int {
	  return !IsFixedStabilizer(rbase.level.Stot, rbase.level.eLev, P.points[i+P.firsts[order[k]]]);});
    }
    k++;
  }
  int p = P.points[ P.firsts[ order[k] ] + l ];
  RegisterRBasePoint(P, rbase, p, TheId);
}

template<typename Telt>
bool Refinements_ProcessFixpoint(rbaseType<Telt> & rbase, imageType<Telt> & image, int const& pnt, int const& cellnum)
{
 int img = FixpointCellNo(image.partition, cellnum);
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
bool PBIsMinimal(std::vector<int> const& range, int const& a, int const& b, StabChain<Telt> const& Stot, int const& eLev)
{
  if (IsInBasicOrbit(Stot, eLev, b)) {
    for (auto & pVal : Stot.stabilizer[eLev].orbit) {
      if (a > pVal)
        return false;
    }
    return true;
  }
  if (b < a)
    return false;
  if (IsFixedStabilizer(Stot, eLev, b))
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
      for (auto & lVal : Stot.stabilizer[eLev].genlabels) {
        int img = PowAct(pnt, Stot.labels[lVal]);
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

template<typename Telt, typename Tint>
ResultPBT<Telt> PartitionBacktrack(StabChain<Telt> const& G, std::function<bool(Telt const&)> const& Pr, bool const& repr, rbaseType<Telt> & rbase, dataType<Telt> const& data, StabChain<Telt> & L, StabChain<Telt> & R)
{
  int n=G.n;
  std::cerr << "PartitionBacktrack step 1\n";
  std::cerr << "rbase.level2.status=" << GetIntTypeNature(rbase.level2.status) << "\n";
  imageType<Telt> image;
  std::cerr << "PartitionBacktrack step 2\n";
  Face orB_sing; // backup of <orb>. We take a single entry. Not sure it is correct
  int nrback;
  std::vector<Face> orb;
  std::vector<Face> org; // intersected (mapped) basic orbits of <G>
  Tplusinfinity<int> blen(true, 0);
  int dd, branch; // branch is level where $Lstab\ne Rstab$ starts
  std::vector<int> range;    // range for construction of <orb>
  Partition oldcel;       // old value of <image.partition.cellno>
  std::vector<int> oldcel_cellno;
  std::cerr << "PartitionBacktrack step 3\n";
  std::function<permPlusBool<Telt>(int const&,bool const&)> PBEnumerate = [&](int const& d, bool const & wasTriv) -> permPlusBool<Telt> {
    std::cerr << "Beginning of PBEnumerate\n";
    permPlusBool<Telt> oldprm, oldprm2;
    int a;                // current R-base point
    permPlusBool<Telt> t; // group element constructed, to be handed upwards
    int m;                // initial number of candidates in <orb>
    int max;              // maximal number of candidates still needed
    boost::dynamic_bitset<>::size_type b;        // image of base point currently being considered

    if (image.perm.status == int_false) {
      std::cerr << "PBEnumerate, case 1, image.perm.status=" << GetIntTypeNature(image.perm.status) << "\n";
      return {int_fail, {}};
    }
    image.depth = d;

    // Store the original values of <image.*>.
    int undoto = NumberCells(image.partition);
    if (image.perm.status == int_true) {
      oldcel = image.partition;
    }
    else {
      oldcel_cellno = image.partition.cellno;
      oldprm = image.perm;
    }
    if (image.level2.status != int_false)
      oldprm2 = image.perm2;
    else
      oldprm2.status = int_false;
    // Recursion comes to an end  if all base  points have been prescribed
    // images.
    if (d >= int(rbase.base.size())) {
      if (IsTrivialRBase(rbase)) {
	blen = rbase.base.size();
	// Do     not  add the   identity    element  in the  subgroup
	// construction.
	if (wasTriv) {
	  // In the subgroup case, assign to  <L> and <R> stabilizer
	  // chains when the R-base is complete.
	  StabChainOptions<Tint> options = GetStandardOptions<Tint>(n);
	  options.base = rbase.base;
	  options.reduced = false;
	  L = StabChainOp<Telt,Tint>(StrongGeneratorsStabChain(L,0), options);
	  R = L;
	  std::cerr << "PBEnumerate, case 2\n";
	  return {int_fail,{}};
	}
	else {
	  permPlusBool<Telt> prm;
	  if (image.perm.status == int_true)
	    prm = {int_perm, MappingPermListList<Telt>(n, rbase.fix[rbase.base.size()-1], Fixcells(image.partition))};
	  else
	    prm = image.perm;
	  if (image.level2.status != int_false) {
	    if (SiftedPermutation(image.level2.Stot, image.level2.eLev, prm.val * Inverse(image.perm2.val)).isIdentity()) {
	      std::cerr << "PBEnumerate, case 3\n";
	      return prm;
	    }
	  }
	  else {
	    if (Pr(prm.val)) {
	      std::cerr << "PBEnumerate, case 4\n";
	      return {int_perm, prm.val};
	    }
	  }
	  std::cerr << "PBEnumerate, case 5\n";
	  return {int_fail, {}};
	}
	// Construct the   next refinement  level. This  also  initializes
	// <image.partition> for the case ``image = base point''.
      }
      else {
	//	if (!repr) {
	//	  oldcel = StructuralCopy( oldcel );
	//	}
	NextRBasePoint(rbase.partition, rbase, G.identity);
	if (image.perm.status == int_true)
	  rbase.fix.push_back(Fixcells(rbase.partition));
	Face eNewF(range.size());
	org.push_back(eNewF);
	if (repr) {
	  // In  the representative  case,  change  the   stabilizer
	  // chains of <L> and <R>.
	  ChangeStabChain(L, d, {rbase.base[d]}, int_false);
	  //	  L[ d + 1 ] := L[ d ].stabilizer;
	  ChangeStabChain(R, d, {rbase.base[d]}, int_false);
	  //	  R[ d + 1 ] := R[ d ].stabilizer;
	}
      }
    }
    a = rbase.base[d];
    
    // Intersect  the current cell of <P>  with  the mapped basic orbit of
    // <G> (and also with the one of <H> in the intersection case).
    if (image.perm.status == int_true) {
      AssignationVectorGapStyle(orb, d, BlistList(range, Cell(oldcel, rbase.where[d]) ));
      if (image.level2.status != int_false) {
	b = orb[d].find_first();
	while (b != boost::dynamic_bitset<>::npos) {
	  if (!IsInBasicOrbit(rbase.lev2[d].Stot, rbase.lev2[d].eLev, SlashAct(b, image.perm2.val))) 
	    orb[d][b] = false;
	  b = orb[d].find_next(b);
	}
      }
    }
    else {
      std::cerr << "Before call to BlistList d=" << d << " |orb|=" << orb.size() << "\n";
      AssignationVectorGapStyle(orb, d, BlistList(range, {}));
      std::cerr << "After  call to BlistList\n";
				// line below needs to be checked.
      for (auto & pVal : rbase.lev[d].Stot.stabilizer[0].orbit) {
	b = PowAct(pVal, image.perm.val);
	if (oldcel_cellno[b] == rbase.where[d]) {
	  bool DoOper=false;
	  if (image.level2.status == int_false)
	    DoOper=true;
	  if (!DoOper) {
	    std::cerr << "d=" << d << " |rbase.lev2|=" << rbase.lev2.size() << "\n";
	    if (IsInBasicOrbit(rbase.lev2[d].Stot, rbase.lev2[d].eLev, SlashAct(b,image.perm2.val)))
	      DoOper=true;
	  }
	  if (DoOper) {
	    orb[d][b] = true;
	    org[d][b] = pVal;
	  }
	}
      }
    }
    if (d == 1 && ForAll(G.labels, [&](Telt const& x){return PowAct(a, x) == a;})) {
      orb[d][a]=true; // ensure a is a possible image (can happen if acting on permutations with more points)
    }
    orB_sing = orb[d];
    
    // Loop  over the candidate images  for the  current base point. First
    // the special case image = base up to current level.
    if (wasTriv) {
      AssignationVectorGapStyle(image.bimg, d, a);
      // Refinements that start with '_' must be executed even when base
      // = image since they modify image.data, etc.
      RRefine(rbase, image, true);
      // Recursion.
      PBEnumerate(d + 1, true);
      image.depth = d;
      // Now we  can  remove  the  entire   <R>-orbit of <a>  from   the
      // candidate list.
      SubtractBlist(orb[d], BlistList(range, L.stabilizer[d].orbit));
    }
    
    // Only the early points of the orbit have to be considered.
    m = SizeBlist( orB_sing );
    if (m < int(L.stabilizer[d].orbit.size()) ) {
      std::cerr << "PBEnumerate, case 6\n";
      return {int_fail,{}};
    }
    max = PositionNthTrueBlist(orB_sing, m - L.stabilizer[d].orbit.size());
    
    if (wasTriv && a > max) {
      m--;
      if (m < int(L.stabilizer[d].orbit.size()) ) {
	std::cerr << "PBEnumerate, case 7\n";
	return {int_fail,{}};
      }
      max = PositionNthTrueBlist( orB_sing, m - L.stabilizer[d].orbit.size());
    }
    // Now the other possible images.
    b = orb[d].find_first();
    while (b != boost::dynamic_bitset<>::npos) {
      int b_int = int(b);
      // Try to prune the node with prop 8(ii) of Leon paper.
      if (!repr && !wasTriv) {
	dd = branch;
	while (dd < d) {
	  if (IsInBasicOrbit(L, dd, a) && !PBIsMinimal(range, R.stabilizer[dd].orbit[0], b_int, R, d))
	    dd = d + 1;
	  else
	    dd = dd + 1;
	}
      }
      else {
	dd = d;
      }
      if (dd == d) {
	// Undo the  changes made to  <image.partition>, <image.level>
	// and <image.perm>.
	for (int i=undoto+1; i<NumberCells(image.partition); i++) 
	  UndoRefinement(image.partition);
	if (image.perm.status != int_true) {
	  image.level = rbase.lev[ d ];
	  image.perm = oldprm;
	}
	if (image.level2.status != int_false) {
	  image.level2 = rbase.lev2[ d ];
	  image.perm2  = oldprm2;
	}
	// If <b> could not be prescribed as image for  <a>, or if the
	// refinement was impossible, give up for this image.
	AssignationVectorGapStyle(image.bimg, d, b_int);
	IsolatePoint( image.partition, b_int );
	
	if (ProcessFixpoint_image(image, a, b_int, org[d][b_int]))
	  t.status = RRefine(rbase, image, false);
	else
	  t.status = int_fail;
	
	if (t.status != int_fail) {
	  // Subgroup case, base <> image   at current level:   <R>,
	  //   which until now is identical to  <L>, must be changed
	  //   without affecting <L>, so take a copy.
	  if (wasTriv && TestEqualityAtLevel(L, R, d)) {
	    SetStabChainFromLevel(R, L, d);
	    branch = d;
	  }
	  if (2 * d <= blen) {
	    ChangeStabChain(R, d, {b_int}, int_false);
	    //	    R[ d + 1 ] = R[ d ].stabilizer;
	  }
	  else {
	    std::vector<Telt> LGen = StrongGeneratorsStabChain( R, d);
	    std::vector<Telt> LGenB = Filtered(LGen, [&](Telt const& gen) -> bool {return PowAct(b_int, gen) == b_int;});
	    //	    R[ d + 1 ] := rec( generators := Filtered( R[ d + 1 ], gen -> b ^ gen = b ) );
	    int largMov=LargestMovedPoint(LGenB);
	    StabChainOptions<Tint> options = GetStandardOptions<Tint>(n);
	    options.base = ClosedInterval(0, largMov);
	    StabChainStrong(R, d+1, LGenB, options);
	  }
	}
        
	// Recursion.
	if (t.status == int_true) {
	  t = PBEnumerate(d + 1, false);
	  nrback++;
	  image.depth = d;
	}
                    
	// If   <t>   =   fail, either   the   recursive   call  was
	//   unsuccessful,  or all new  elements   have been added  to
	//   levels  below  the current one   (this happens if  base =
	//   image up to current level).
	if (t.status != int_fail) {
	  // Representative case, element found: Return it.
	  // Subgroup case, base <> image  before current level:  We
	  //   need  only find  a representative  because we already
	  //   know the stabilizer of <L> at an earlier level.
	  if (repr || !wasTriv) {
	    std::cerr << "PBEnumerate, case 8\n";
	    return t;
	  }
	  else {
	    // Subgroup case, base  <> image at current level: Enlarge
	    //   <L>    with  <t>. Decrease <max>     according to the
	    //   enlarged <L>. Reset <R> to the enlarged <L>.
	    //	    for (int dd=0; dd<d; dd++)
	    //	      AddGeneratorsExtendSchreierTree( L[ dd ], {t});
	    AddGeneratorsExtendSchreierTree(L, 0, {t.val});
	    if (m < int(L.stabilizer[d].orbit.size())) {
	      std::cerr << "PBEnumerate, case 9\n";
	      return {int_fail,{}};
	    }
	    max = PositionNthTrueBlist( orB_sing, m - L.stabilizer[d].orbit.size());
	    SetStabChainFromLevel(R, L, d);
	  }
	}
        
	// Now  we can remove the   entire <R>-orbit  of <b> from  the
	// candidate list.
	if (R.stabilizer[d].transversal[b] != -1)
	  SubtractBlist(orb[d], BlistList(range, R.stabilizer[d].orbit));
	else
	  SubtractBlistOrbitStabChain(orb[d], StrongGeneratorsStabChain(R, d), b_int);
	b = orb[d].find_next(b);
      }
      
    }
    std::cerr << "PBEnumerate, case 10\n";
    return {int_fail, {}};
  };
  std::cerr << "PartitionBacktrack step 4\n";

  nrback=0; // count the number of times we jumped up

  // Trivial cases first.
  if (IsTrivial(G)) {
    if (!repr)
      return {int_group, G, {}};
    if (Pr(G.identity))
      return {int_perm,{},G.identity};
    else
      return {int_fail,{},{}};
  }
  std::cerr << "PartitionBacktrack step 5\n";
    
  // Construct the <image>.
  image.data=data;
  std::cerr << "PartitionBacktrack step 5.1\n";
  image.depth=1;
  if (repr) {
    image.partition = data.P;
  }
  else {
    image.partition = rbase.partition;
  }
  if (IsBool(rbase.level2)) {
    std::cerr << "PartitionBacktrack step 5.2\n";
    image.level2 = {int_false,-444,{},0};
    std::cerr << "PartitionBacktrack step 5.3\n";
  }
  else {
    std::cerr << "PartitionBacktrack step 5.4\n";
    std::cerr << "bool=" << rbase.level2.Stot.UseCycle << "\n";
    std::cerr << "PartitionBacktrack step 5.4.1\n";
    image.level2 = rbase.level2;
    std::cerr << "PartitionBacktrack step 5.5\n";
    image.perm2  = {int_perm, G.identity};
    std::cerr << "PartitionBacktrack step 5.6\n";
  }
  std::cerr << "PartitionBacktrack step 6\n";
    
  // If  <Pr> is  function,   multiply  permutations. Otherwise, keep   them
  // factorized.
  image.perm = {int_perm, G.identity};
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
      for (int i=0; i<int(fix.size()); i++)
	ProcessFixpoint_image(image, fix[i], fixP[i], -1);
    }
    // In   the representative case,   assign  to <L>  and <R>  stabilizer
    // chains.
    //    L := ListStabChain( CopyStabChain( StabChainMutable( L ) ) );
    //    R := ListStabChain( CopyStabChain( StabChainMutable( R ) ) );
  }
  std::cerr << "PartitionBacktrack step 7\n";
    
  int lenD=rbase.domain[rbase.domain.size()-1];
  for (int i=0; i<lenD; i++)
    range.push_back(i);
  permPlusBool<Telt> rep = PBEnumerate(0, !repr);
  std::cerr << "PartitionBacktrack step 8\n";
  if (!repr) {
    return {int_group, L, {}};
  }
  else {
    if (rep.status == int_perm)
      return {int_perm, {}, rep.val};
    else
      return {int_fail, {}, {}};
  }
  std::cerr << "PartitionBacktrack step 9\n";
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
  std::cerr << "Beginning of RepOpSetsPermGroup\n";
  std::cerr << "UseCycle=" << G.UseCycle << "\n";
  std::cerr << "After bool print\n";
  int n=G.n;
  std::vector<int> Omega = MovedPoints(G);
  if (repr && Phi.size() != Psi.size())
    return {int_fail, {}, {}};
  if (IsSubset(Phi, Omega) || ForAll(Omega, [&](int const &p) -> bool {return !Phi[p];})) {
    if (repr) {
      if (Difference_face(Phi, Omega) != Difference_face(Psi, Omega))
	return {int_fail, {}, {}};
      else
	return {int_perm,{},G.identity};
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
    for (auto & eVal : Omega)
      if (Ph[eVal] == 1)
	IntVect.push_back(eVal);
    //    std::vector<int> UnionVect;
    Face Ph_copy = Ph;
    for (auto & eVal : Omega)
      Ph_copy[eVal]=1;
    return GetPartition({IntVect, FaceToVector(Ph_copy)});
  };



  
  Partition P = GetPartitionFromPair(Phi);
  Partition Q = GetPartitionFromPair(Psi);


  auto GetSubgroup=[&](Face const& Ph) -> StabChain<Telt> {
    std::vector<Telt> LGen = StrongGeneratorsStabChain(G, 0);
    std::cerr << "GetSubgroup, |LGen|=" << LGen.size() << "\n";
    std::cerr << "GetSubgroup, LGen=";
    for (auto & eGen : LGen)
      std::cerr << " " << eGen;
    std::cerr << "\n";
    std::vector<Telt> sgs=Filtered(StrongGeneratorsStabChain(G, 0), [&](Telt const& g)->bool{return OnSets(Ph, g) == Ph;});
    return MinimalStabChain<Telt,Tint>(sgs, n);
  };
  
  
  StabChain<Telt> L = GetSubgroup(Phi);
  StabChain<Telt> R;
  if (repr)
    R = GetSubgroup(Psi);
  else
    R = L;
  std::cerr << "Orders: |R|=" << SizeStabChain<Telt,Tint>(R) << " |L|=" << SizeStabChain<Telt,Tint>(L) << "\n";
  rbaseType<Telt> rbase = EmptyRBase<Telt>({G, G}, true, Omega, P);
  std::cerr << "RepOpSetsPermGroup rbase.level2.status=" << GetIntTypeNature(rbase.level2.status) << "\n";
  std::vector<int> Phi_vect = FaceToVector(Phi);
  std::function<bool(Telt const&)> Pr=[&](Telt const& gen) -> bool {
    for (auto & i : Phi_vect) {
      int iImg=gen.at(i);
      if (Psi[iImg] == 0)
	return false;
    }
    return true;
  };
  std::cerr << "Before call to PartitionBacktrack\n";
  std::cerr << "bool=" << rbase.level2.Stot.UseCycle << "\n";
  std::cerr << "After bool print\n";
  return PartitionBacktrack<Telt,Tint>( G, Pr, repr, rbase, {Q}, L, R );
}

template<typename Telt,typename Tint>
StabChain<Telt> Stabilizer_OnSets(StabChain<Telt> const& G, Face const& Phi)
{
  std::cerr << "Beginning of Stabilizer_OnSets\n";
  bool repr=false;
  return RepOpSetsPermGroup<Telt,Tint>(G, repr, Phi, Phi).stab;
}


}

#endif
