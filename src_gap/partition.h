#ifndef DEFINE_PERMUTALIB_PARTITION_H
#define DEFINE_PERMUTALIB_PARTITION_H

#include "PermGroup.h"
#include "GapPrint.h"

#ifdef SYNCHRONIZED_DEBUG_GAP478
# define DEBUG_PARTITION
#endif

#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
# define CHECK_PARTITION
#endif



namespace permutalib {


template<typename Tidx>
struct Partition {
  std::vector<Tidx> points;
  std::vector<Tidx> firsts;
  std::vector<Tidx> lengths;
  std::vector<Tidx> cellno;
};


template<typename Tidx>
void NicePrintPartition(std::string const& str, Partition<Tidx> const& P)
{
  Tidx nbPart=P.firsts.size();
  std::cerr << str << " = [ ";
  for (Tidx iPart=0; iPart<nbPart; iPart++) {
    if (iPart > 0)
      std::cerr << ", ";
    Tidx eFirst=P.firsts[iPart];
    Tidx len=P.lengths[iPart];
    std::cerr << "[ ";
    for (Tidx i=0; i<len; i++) {
      if (i >0)
	std::cerr << ", ";
      Tidx ePt = P.points[eFirst + i];
      std::cerr << int(ePt+1);
    }
    std::cerr << " ]";
  }
  std::cerr << " ]\n";
}


template<typename Tidx>
void RawPrintPartition(Partition<Tidx> const& P)
{
  Tidx nbPart=P.lengths.size();
  std::vector<std::string> LPart(nbPart);
  for (Tidx iPart=0; iPart<nbPart; iPart++) {
    Tidx len=P.lengths[iPart];
    Tidx eFirst=P.firsts[iPart];
    std::vector<Tidx> eList(len);
    for (Tidx u=0; u<len; u++) {
      Tidx ePt = P.points[eFirst + u];
      eList[u] = ePt;
    }
    LPart[iPart] = GapStringIntVector(eList);
  }
  std::cerr << "CPP PART Partition=" << GapStringTVector(LPart) << "\n";
  std::cerr << "CPP PART points=" << GapStringIntVector(P.points) << "\n";
  std::cerr << "CPP PART firsts=" << GapStringIntVector(P.firsts) << "\n";
  std::cerr << "CPP PART lengths=" << GapStringTVector(P.lengths) << "\n";
  std::cerr << "CPP PART cellno=" << GapStringIntVector(P.cellno) << "\n";
}


template<typename Tidx>
void CheckConsistencyPartition(std::string const& str, Partition<Tidx> const& P)
{
  Tidx nbError=0;
  Tidx Max_points=0;
  bool ChangeVal=false;
  for (auto & eVal : P.points)
    if (eVal > Max_points) {
      Max_points = eVal;
      ChangeVal=true;
    }
  if (ChangeVal)
    Max_points++;
  if (Tidx(P.cellno.size()) != Max_points) {
    std::cerr << "P.cellno.size()=" << P.cellno.size() << " Max_points=" << int(Max_points) << "\n";
    std::cerr << "1: We should have |P.cellno| = Maximum(P.points) (at least, that's true at init)\n";
    nbError++;
  }
  if (P.firsts.size() != P.lengths.size()) {
    std::cerr << "2: P.firsts and P.lengths have different lengths\n";
    nbError++;
  }
  Tidx nbPoint=Tidx(P.cellno.size());
  Tidx nbPart=Tidx(P.lengths.size());
  std::vector<Tidx> MeasuredLength(nbPart,0);
  for (Tidx iPoint=0; iPoint<nbPoint; iPoint++) {
    Tidx iPart=P.cellno[iPoint];
    if ((iPart >= nbPart) && (iPart != std::numeric_limits<Tidx>::max())) {
      std::cerr << "3: Error, iPart=" << int(iPart) << " but nbPart=" << int(nbPart) << "\n";
      nbError++;
    }
    if (iPart >= 0 && iPart < nbPart)
      MeasuredLength[iPart]++;
  }
  for (Tidx iPart=0; iPart<nbPart; iPart++) {
    if (MeasuredLength[iPart] != P.lengths[iPart]) {
      std::cerr << "4: At iPart=" << int(iPart) << " we have error in lengths: MeasuredLength[iPart]=" << int(MeasuredLength[iPart]) << " and P.lengths[iPart]=" << int(P.lengths[iPart]) << "\n";
      nbError++;
    }
    if (MeasuredLength[iPart] == 0) {
      std::cerr << "5: At iPart=" << int(iPart) << " the length is zero\n";
      nbError++;
    }
  }
  for (Tidx iPart=0; iPart<nbPart; iPart++) {
    Tidx len=P.lengths[iPart];
    Tidx eFirst=P.firsts[iPart];
    for (Tidx u=0; u<len; u++) {
      Tidx ePt = P.points[eFirst + u];
      if (ePt < 0 || ePt >= nbPoint) {
	std::cerr << "6: At iPart=" << int(iPart) << " u=" << int(u) << " ePt=" << int(ePt) << " point out of range\n";
	nbError++;
      }
      if (ePt >= 0 && ePt < nbPoint) {
	if (P.cellno[ePt] != iPart) {
	  std::cerr << "7: At iPart=" << int(iPart) << " u=" << int(u) << " ePt=" << int(ePt) << " has wrong cellno\n";
	  nbError++;
	}
      }
    }
  }
  if (nbError > 0) {
    std::cerr << "Error at " << str << ". We found nbError=" << nbError << "\n";
    throw PermutalibException{1};
  }
}


template<typename Tidx>
Partition<Tidx> GetPartition(std::vector<std::vector<Tidx>> const& list)
{
#ifdef DEBUG_PARTITION
  std::vector<std::string> LStr;
  for (auto & eList : list)
    LStr.push_back(GapStringIntVector(eList));
  std::cerr << "CPP list=" << GapStringTVector(LStr) << "\n";
#endif
  std::vector<Tidx> points;
  Tidx Max_NbPoint = 0;
  bool ChangeVal = false;
  for (auto & eList : list)
    for (auto & eVal : eList) {
      points.push_back(eVal);
      if (eVal > Max_NbPoint) {
        Max_NbPoint = eVal;
        ChangeVal = true;
      }
    }
  if (ChangeVal)
    Max_NbPoint++;
  Tidx nbPoint=Max_NbPoint;
  Tidx nbPart=Tidx(list.size());
  std::vector<Tidx> firsts(nbPart);
  std::vector<Tidx> lengths(nbPart);
  std::vector<Tidx> cellno(nbPoint, std::numeric_limits<Tidx>::max());
  Tidx i=0;
  for (Tidx iPart=0; iPart<nbPart; iPart++) {
    Tidx len=Tidx(list[iPart].size());
    firsts[iPart]=i;
    lengths[iPart]=len;
    i += len;
    for (auto & eVal : list[iPart])
      cellno[eVal] = iPart;
  }
  Partition<Tidx> P{points, firsts, lengths, cellno};
#ifdef DEBUG_PARTITION
  std::cerr << "CPP After GetPartition operation P=\n";
  RawPrintPartition(P);
#endif
#ifdef CHECK_PARTITION
  CheckConsistencyPartition("GetPartition", P);
#endif
  return P;
}


template<typename Tidx>
Tidx GetNumberPoint(Partition<Tidx> const& ePartition)
{
  Tidx nbPoint = Tidx(ePartition.points.size());
  return nbPoint;
}


template<typename Tidx>
Tidx NumberCells(Partition<Tidx> const& ePartition)
{
  Tidx nbPart=Tidx(ePartition.firsts.size());
  return nbPart;
}


template<typename Tidx>
std::vector<Tidx> Cell(Partition<Tidx> const& ePartition, Tidx const& iPart)
{
  Tidx len=ePartition.lengths[iPart];
  Tidx eFirst=ePartition.firsts[iPart];
  std::vector<Tidx> eList(len);
  for (Tidx i=0; i<len; i++)
    eList[i]=ePartition.points[eFirst + i];
  return eList;
}


template<typename Tidx>
std::vector<std::vector<Tidx>> Cells(Partition<Tidx> const& ePartition)
{
  Tidx nbPart=ePartition.firsts.size();
  std::vector<std::vector<Tidx>> eListList(nbPart);
  for (Tidx iPart=0; iPart<nbPart; iPart++)
    eListList[iPart] = Cell(ePartition, iPart);
  return eListList;
}


template<typename Tidx>
Tidx CellNoPoint(Partition<Tidx> const& ePartition, Tidx const& pt)
{
  return ePartition.cellno[pt];
}


template<typename Tidx>
std::vector<Tidx> CellNoPoints(Partition<Tidx> const& ePartition, std::vector<Tidx> const& pts)
{
  Tidx len=pts.size();
  std::vector<Tidx> ret(len);
  for (Tidx i=0; i<len; i++)
    ret[i] = ePartition.cellno[pts[i]];
  return ret;
}


template<typename Tidx>
bool PointInCellNo(Partition<Tidx> const& ePartition, Tidx const& pt, Tidx const& iPart)
{
  return ePartition.cellno[pt] == iPart;
}


template<typename Tidx>
std::vector<Tidx> Fixcells(Partition<Tidx> const& ePartition)
{
  std::vector<Tidx> fix;
  Tidx nbPart=Tidx(ePartition.firsts.size());
  for (Tidx iPart=0; iPart<nbPart; iPart++) {
    if (ePartition.lengths[iPart] == 1) {
      Tidx eFirst=ePartition.firsts[iPart];
      Tidx eVal=ePartition.points[eFirst];
      fix.push_back(eVal);
    }
  }
  return fix;
}


template<typename Tidx>
Tidx SplitCell_Kernel(Partition<Tidx> & P, Tidx const& i, std::function<bool(Tidx)> const& test, Tidx const& out)
{
#ifdef DEBUG_PARTITION
  //  std::cerr << "CPP i=" << int(i+1) << " out=" << out << "\n";
  std::cerr << "CPP Before SplitCell_Kernel operation P=\n";
  RawPrintPartition(P);
#endif
#ifdef CHECK_PARTITION
  CheckConsistencyPartition("Input SplitCell_Kernel", P);
#endif
  Tidx a=P.firsts[i];
  Tidx b=a + P.lengths[i];
  Tidx l=b - 1;

  Tidx maxmov;
  if (out != std::numeric_limits<Tidx>::max() )
    maxmov = out;
  else
    maxmov = P.lengths[i]-1;
  Tidx B = l - maxmov;
#ifdef DEBUG_PARTITION
  std::cerr << "CPP maxmov=" << maxmov << " B=" << int(B+1) << "\n";
#endif
  a--;
  while (a<b) {
#ifdef DEBUG_PARTITION
    std::cerr << "CPP     1 a=" << int(a+1) << " b=" << int(b+1) << "\n";
#endif
    while(true) {
#ifdef DEBUG_PARTITION
      std::cerr << "CPP B LOOP\n";
#endif
      b--;
      if (b < B) {
#ifdef DEBUG_PARTITION
        std::cerr << "CPP exit 1\n";
#endif
        return std::numeric_limits<Tidx>::max();
      }
      if (!test(b))
        break;
    }
#ifdef DEBUG_PARTITION
    std::cerr << "CPP     2 a=" << int(a+1) << " b=" << int(b+1) << "\n";
#endif
    while(true) {
#ifdef DEBUG_PARTITION
      std::cerr << "CPP A LOOP\n";
#endif
      a++;
      if (a>b || test(a))
        break;
    }
#ifdef DEBUG_PARTITION
    std::cerr << "CPP     3 a=" << int(a+1) << " b=" << int(b+1) << "\n";
#endif
    if (a<b) {
      std::swap(P.points[a], P.points[b]);
    }
  }
#ifdef DEBUG_PARTITION
  std::cerr << "CPP a=" << int(a+1) << " l=" << int(l+1) << "\n";
#endif
  if (a > l) {
#ifdef DEBUG_PARTITION
    std::cerr << "CPP exit 2\n";
#endif
    return std::numeric_limits<Tidx>::max();
  }
  Tidx m=Tidx(P.firsts.size());
  for (Tidx idx=a; idx<=l; idx++)
    P.cellno[P.points[idx]] = m;
  P.firsts.push_back(a);
  P.lengths.push_back(l - a + 1);
  P.lengths[i] = P.lengths[i] - P.lengths[m];

#ifdef DEBUG_PARTITION
  std::cerr << "CPP After SplitCell_Kernel operation P=\n";
  RawPrintPartition(P);
#endif
#ifdef CHECK_PARTITION
  CheckConsistencyPartition("Output SplitCell_Kernel", P);
#endif
#ifdef DEBUG_PARTITION
  std::cerr << "CPP exit 3\n";
#endif
  return P.lengths[m];
}

template<typename Telt>
typename Telt::Tidx SplitCell_Partition(Partition<typename Telt::Tidx> & P, typename Telt::Tidx const& i, Partition<typename Telt::Tidx> const& Q, typename Telt::Tidx const& j, Telt const& g, typename Telt::Tidx const& out)
{
  using Tidx = typename Telt::Tidx;
#ifdef DEBUG_PARTITION
  std::cerr << "CPP SplitCell g=" << g << "\n";
  std::cerr << "CPP Q=\n";
  RawPrintPartition(Q);
#endif
  std::function<bool(Tidx)> test=[&](Tidx ePt) -> bool {
    Tidx fPt=PowAct(P.points[ePt], g);
#ifdef DEBUG_PARTITION
    std::cerr << "CPP SplitCellTestfun1 fPt=" << int(fPt+1) << "\n";
#endif
    return PointInCellNo(Q, fPt, j);
  };
  return SplitCell_Kernel(P, i, test, out);
}


template<typename Telt>
typename Telt::Tidx SplitCell_Face(Partition<typename Telt::Tidx> & P, typename Telt::Tidx const& i, Face const& f, typename Telt::Tidx const& j, Telt const& g, typename Telt::Tidx const& out)
{
  using Tidx = typename Telt::Tidx;
  std::function<bool(Tidx)> test=[&](Tidx ePt) -> bool {
    Tidx fPt=PowAct(P.points[ePt], g);
    if (j == 1) {
      if (f[fPt] == 1)
	return true;
      return false;
    } else {
      if (f[fPt] == 0)
	return true;
      return false;
    }
  };
  return SplitCell_Kernel(P, i, test, out);
}


template<typename Tidx>
Tidx IsolatePoint(Partition<Tidx> & P, Tidx const& a)
{
#ifdef CHECK_PARTITION
  CheckConsistencyPartition("Input IsolatePoint", P);
#endif
#ifdef DEBUG_PARTITION
  std::cerr << "CPP Input Partition\n";
  RawPrintPartition(P);
#endif
  Tidx nbPart=Tidx(P.firsts.size());
  Tidx iPart=P.cellno[a];
  Tidx eFirst=P.firsts[iPart];
  Tidx len=P.lengths[iPart];
  //  std::cerr << "a=" << a << " iPart=" << iPart << " eFirst=" << eFirst << " len=" << len << "\n";
  if (len == 1)
    return std::numeric_limits<Tidx>::max();
  Tidx pos=0;
  for (Tidx j=eFirst; j<eFirst + len; j++) {
    Tidx ePt=P.points[j];
    if (ePt == a)
      pos=j;
  }
  Tidx l=(eFirst + len) - 1;
  //  std::cerr << "pos=" << pos << " l=" << l << "\n";
  P.points[pos] = P.points[l];
  P.points[l]=a;
  //  int m=P.firsts.size() + 1;
  P.cellno[a]=nbPart;
  P.lengths[iPart]--;
  P.firsts.push_back(l);
  P.lengths.push_back(1);

#ifdef DEBUG_PARTITION
  std::cerr << "CPP After IsolatePoint operation P=\n";
  RawPrintPartition(P);
#endif
#ifdef CHECK_PARTITION
  CheckConsistencyPartition("Output IsolatePoint", P);
#endif
  return iPart;
}



template<typename Tidx>
Tidx UndoRefinement(Partition<Tidx> & P)
{
#ifdef CHECK_PARTITION
  CheckConsistencyPartition("Input UndoRefinement", P);
#endif
  Tidx nbPart=Tidx(P.firsts.size());
  Tidx pfm=P.firsts[nbPart-1];
  if (pfm == 0)
    return std::numeric_limits<Tidx>::max();
  Tidx plm=P.lengths[nbPart-1];
  Tidx m=P.cellno[P.points[pfm-1]];
  P.lengths[m] += plm;
  for (Tidx j=pfm; j<pfm + plm; j++) {
    Tidx ePt=P.points[j];
    P.cellno[ePt]=m;
  }
  P.firsts.pop_back();
  P.lengths.pop_back();
#ifdef DEBUG_PARTITION
  std::cerr << "CPP After UndoRefinement operation P=\n";
  RawPrintPartition(P);
#endif
#ifdef CHECK_PARTITION
  CheckConsistencyPartition("Output UndoRefinement", P);
#endif
  return m;
}


template<typename Tidx>
Tidx FixpointCellNo(Partition<Tidx> const& P, Tidx const& i)
{
  return P.points[P.firsts[i]];
}


template<typename Tidx>
Tidx FixcellPoint(Partition<Tidx> const& P, std::set<Tidx> & old)
{
  Tidx nbPart=P.lengths.size();
  std::vector<Tidx> poss;
  for (Tidx iPart=0; iPart<nbPart; iPart++) {
    if (P.lengths[iPart] == 1 && old.find(iPart) == old.end())
      poss.push_back(iPart);
  }
  Tidx nbPoss=poss.size();
  if (nbPoss == 0)
    return std::numeric_limits<Tidx>::max();
  Tidx idx=rand() % nbPoss;
  Tidx p=poss[idx];
  old.insert(p);
  return p;
}


template<typename Tidx>
struct typeFixcellsCell {
  bool res;
  std::vector<Tidx> K;
  std::vector<Tidx> I;
};


template<typename Tidx>
typeFixcellsCell<Tidx> FixcellsCell(Partition<Tidx> const& P, Partition<Tidx> const& Q, std::set<Tidx> & old)
{
  std::vector<Tidx> K, I;
  Tidx nbPart=P.firsts.size();
  for (Tidx iPart=0; iPart<nbPart; iPart++) {
    Tidx start=P.firsts[iPart];
    Tidx kPart=CellNoPoint(Q, P.points[start]);
    if (old.find(kPart) == old.end()) {
      std::function<bool()> eval=[&]() -> bool {
	for (int j=1; j<P.lengths[iPart]; j++) {
	  if (CellNoPoint(Q, P.points[j]) != kPart)
	    return false;
	}
	return true;
      };
      if (eval()) {
	old.insert(kPart);
	K.push_back(kPart);
	I.push_back(iPart);
      }
    }
  }
  if (K.size() == 0) {
    return {false,{},{}};
  }
  else {
    return {true, K, I};
  }
}


template<typename Tidx>
Partition<Tidx> TrivialPartition(std::vector<Tidx> const& Omega)
{
  return GetPartition({Omega});
}


template<typename Telt>
std::vector<std::vector<typename Telt::Tidx>> OrbitsPermsB(std::vector<Telt> const& gens, std::vector<typename Telt::Tidx> const& Omega)
{
  using Tidx = typename Telt::Tidx;
#ifdef DEBUG_PARTITION
  std::cerr << "DEBUG OrbitsPermB beginning\n";
#endif
  Tidx max=LargestMovedPoint(gens);
  Tidx maxP1 = 0;
  if (max != std::numeric_limits<Tidx>::max())
    maxP1 = max + 1;
  Face dom(maxP1);
  for (auto & eVal : Omega)
    if (eVal < maxP1)
      dom[eVal] = 1;
#ifdef DEBUG_PARTITION
  std::cerr << "DEBUG dom built\n";
#endif
  Face newF(maxP1);
  for (Tidx i=0; i<maxP1; i++)
    newF[i] = 1;
#ifdef DEBUG_PARTITION
  std::cerr << "DEBUG newF built\n";
#endif
  std::vector<std::vector<Tidx>> orbs;
  boost::dynamic_bitset<>::size_type fst=dom.find_first();
  while (fst != boost::dynamic_bitset<>::npos) {
    Tidx fst_i = Tidx(fst);
    std::vector<Tidx> orb{fst_i};
    newF[fst_i] = 0;
    dom [fst_i] = 0;
#ifdef DEBUG_PARTITION
    std::cerr << "DEBUG Before beginning of pnt loop\n";
#endif
    size_t posOrb=0;
    while(true) {
      Tidx pnt = orb[posOrb];
      for (auto & gen : gens) {
        Tidx img = PowAct(pnt, gen);
        if (newF[img]) {
          orb.push_back(img);
          newF[img] = 0;
          dom [img] = 0;
        }
      }
      posOrb++;
      if (posOrb >= orb.size())
        break;
    }
#ifdef DEBUG_PARTITION
    std::cerr << "DEBUG After the pnt loop\n";
#endif
    orbs.push_back(orb);
    fst=dom.find_first();
  }
  for (auto & pnt : Omega)
    if (pnt >= maxP1)
      orbs.push_back({pnt});
  return orbs;
}



template<typename Telt>
Partition<typename Telt::Tidx> OrbitsPartition(std::vector<Telt> const& gens, typename Telt::Tidx const&n, std::vector<typename Telt::Tidx> const& Omega)
{
#ifdef DEBUG_PARTITION
  std::cerr << "CPP OrbitsPartition, using OrbitsPerms\n";
  std::cerr << "CPP generators=" << GapStringTVector(gens) << "\n";
#endif
  return GetPartition(OrbitsPermsB(gens, Omega));
}


template<typename Tarith>
int SmallestPrimeDivisor(Tarith const& size)
{
  if (size == 1)
    return 1;
  int i=2;
  while(true) {
    int res=size % i;
    if (res == 0)
      return i;
    i++;
  }
}


template<typename Tarith, typename Tidx>
Partition<Tidx> CollectedPartition(Partition<Tidx> const& P, Tarith const& size)
{
  Partition<Tidx> C = P;
  Tidx div=SmallestPrimeDivisor(size);
  Tidx nbPart=P.firsts.size();
  Tidx nbPartTot=nbPart;
  for (Tidx iPart=0; iPart<nbPart; iPart++) {
    Tidx sizPart=P.lengths[iPart];
    if (sizPart < div) {
      Tidx eFirst=P.firsts[iPart];
      for (Tidx i=0; i<sizPart; i++) {
	Tidx ePt=P.points[eFirst + i];
	if (i == 0) {
	  C.lengths[iPart]=1;
	} else {
	  C.cellno[ePt] = nbPartTot;
	  C.firsts.push_back(eFirst + i);
	  C.lengths.push_back(1);
	}
      }
    }
  }
  return C;
}


}

#endif
