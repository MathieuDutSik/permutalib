#ifndef DEFINE_PARTITION_CODE
#define DEFINE_PARTITION_CODE

#include "PermGroup.h"

#define DEBUG_PARTITION
#define CHECK_PARTITION

namespace permutalib {

struct Partition {
  std::vector<int> points;
  std::vector<int> firsts;
  std::vector<int> lengths;
  std::vector<int> cellno;
};

void NicePrintPartition(std::string const& str, Partition const& P)
{
  int nbPart=P.firsts.size();
  std::cerr << str << " = [ ";
  for (int iPart=0; iPart<nbPart; iPart++) {
    if (iPart > 0)
      std::cerr << ", ";
    int eFirst=P.firsts[iPart];
    int len=P.lengths[iPart];
    std::cerr << "[ ";
    for (int i=0; i<len; i++) {
      if (i >0)
	std::cerr << ", ";
      int ePt = P.points[eFirst + i];
      std::cerr << (ePt+1);
    }
    std::cerr << " ]";
  }
  std::cerr << " ]\n";
}

 
void RawPrintPartition(Partition const& P)
{
  int nbPart=P.lengths.size();
  std::vector<std::string> LPart(nbPart);
  for (int iPart=0; iPart<nbPart; iPart++) {
    int len=P.lengths[iPart];
    int eFirst=P.firsts[iPart];
    std::vector<int> eList(len);
    for (int u=0; u<len; u++) {
      int ePt = P.points[eFirst + u];
      eList[u] = ePt;
    }
    LPart[iPart] = GapStringIntVector(eList);
  }
  std::cerr << "CPP Partition=" << GapStringTVector(LPart) << "\n";
  std::cerr << "CPP points=" << GapStringIntVector(P.points) << "\n";
  std::cerr << "CPP firsts=" << GapStringIntVector(P.firsts) << "\n";
  std::cerr << "CPP lengths=" << GapStringTVector(P.lengths) << "\n";
  std::cerr << "CPP cellno=" << GapStringIntVector(P.cellno) << "\n";
}


void CheckConsistencyPartition(std::string const& str, Partition const& P)
{
  int nbError=0;
  if (P.cellno.size() != P.points.size()) {
    std::cerr << "1: P.cellno and P.points have different lengths\n";
    nbError++;
  }
  if (P.firsts.size() != P.lengths.size()) {
    std::cerr << "2: P.firsts and P.lengths have different lengths\n";
    nbError++;
  }
  int nbPoint=P.cellno.size();
  int nbPart=P.lengths.size();
  std::vector<int> MeasuredLength(nbPart,0);
  for (int iPoint=0; iPoint<nbPoint; iPoint++) {
    int iPart=P.cellno[iPoint];
    if (iPart >= nbPart) {
      std::cerr << "3: Error, iPart=" << iPart << " but nbPart=" << nbPart << "\n";
      nbError++;
    }
    if (iPart < nbPart)
      MeasuredLength[iPart]++;
  }
  for (int iPart=0; iPart<nbPart; iPart++) {
    if (MeasuredLength[iPart] != P.lengths[iPart]) {
      std::cerr << "4: At iPart=" << iPart << " we have error in lengths\n";
      nbError++;
    }
    if (MeasuredLength[iPart] == 0) {
      std::cerr << "5: At iPart=" << iPart << " the length is zero\n";
      nbError++;
    }
  }
  for (int iPart=0; iPart<nbPart; iPart++) {
    int len=P.lengths[iPart];
    int eFirst=P.firsts[iPart];
    for (int u=0; u<len; u++) {
      int ePt = P.points[eFirst + u];
      if (ePt < 0 || ePt >= nbPoint) {
	std::cerr << "6: At iPart=" << iPart << " u=" << u << " ePt=" << ePt << " point out of range\n";
	nbError++;
      }
      if (ePt >= 0 && ePt < nbPoint) {
	if (P.cellno[ePt] != iPart) {
	  std::cerr << "7: At iPart=" << iPart << " u=" << u << " ePt=" << ePt << " has wrong cellno\n";
	  nbError++;
	}
      }
    }
  }
  if (nbError > 0) {
    std::cerr << "Error at " << str << ". We found nbError=" << nbError << "\n";
    throw TerminalException{1};
  }
}

 
Partition GetPartition(std::vector<std::vector<int>> const& list)
{
  std::vector<int> points;
  for (auto & eList : list)
    for (auto & eVal : eList)
      points.push_back(eVal);
  int nbPoint=points.size();
  int nbPart=list.size();
  std::vector<int> firsts(nbPart);
  std::vector<int> lengths(nbPart);
  std::vector<int> cellno(nbPoint);
  int i=0;
  for (int iPart=0; iPart<nbPart; iPart++) {
    int len=list[iPart].size();
    firsts[iPart]=i;
    lengths[iPart]=len;
    i += len;
    for (auto & eVal : list[iPart])
      cellno[eVal] = iPart;
  }
  Partition P{points, firsts, lengths, cellno};
#ifdef DEBUG_PARTITION
  std::cerr << "CPP After GetPartition operation P=\n";
  RawPrintPartition(P);
#endif
#ifdef CHECK_PARTITION
  CheckConsistencyPartition("GetPartition", P);
#endif
  return P;
}


int GetNumberPoint(Partition const& ePartition)
{
  int nbPoint = ePartition.points.size();
  return nbPoint;
}


int NumberCells(Partition const& ePartition)
{
  int nbPart=ePartition.firsts.size();
  return nbPart;
}


std::vector<int> Cell(Partition const& ePartition, int const& iPart)
{
  int len=ePartition.lengths[iPart];
  int eFirst=ePartition.firsts[iPart];
  std::vector<int> eList(len);
  for (int i=0; i<len; i++) {
    int eVal=ePartition.points[eFirst + i];
    eList[i]=eVal;
  }
  return eList;
}

std::vector<std::vector<int>> Cells(Partition const& ePartition)
{
  int nbPart=ePartition.firsts.size();
  std::vector<std::vector<int>> eListList(nbPart);
  for (int iPart=0; iPart<nbPart; iPart++)
    eListList[iPart] = Cell(ePartition, iPart);
  return eListList;
}


int CellNoPoint(Partition const& ePartition, int const& pt)
{
  return ePartition.cellno[pt];
}


std::vector<int> CellNoPoints(Partition const& ePartition, std::vector<int> const& pts)
{
  int len=pts.size();
  std::vector<int> ret(len);
  for (int i=0; i<len; i++)
    ret[i] = ePartition.cellno[pts[i]];
  return ret;
}

bool PointInCellNo(Partition const& ePartition, int const& pt, int const& iPart)
{
  return ePartition.cellno[pt] == iPart;
}

std::vector<int> Fixcells(Partition const& ePartition)
{
  std::vector<int> fix;
  int nbPart=ePartition.firsts.size();
  for (int iPart=0; iPart<nbPart; iPart++) {
    if (ePartition.lengths[iPart] == 1) {
      int eFirst=ePartition.firsts[iPart];
      int eVal=ePartition.points[eFirst];
      fix.push_back(eVal);
    }
  }
  return fix;
}


int SplitCell_Kernel(Partition & P, int const& i, std::function<bool(int)> const& test, int const& out)
{
#ifdef DEBUG_PARTITION
  //  std::cerr << "CPP i=" << (i+1) << " out=" << out << "\n";
  std::cerr << "CPP Before SplitCell_Kernel operation P=\n";
  RawPrintPartition(P);
#endif
#ifdef CHECK_PARTITION
  CheckConsistencyPartition("Input SplitCell_Kernel", P);
#endif
  int a=P.firsts[i];
  int b=a + P.lengths[i];
  int l=b-1;

  int maxmov;
  if (out >= 0)
    maxmov = out;
  else
    maxmov = P.lengths[i]-1;
  int B = l - maxmov;
  std::cerr << "CPP maxmov=" << maxmov << " B=" << (B+1) << "\n";
  a--;
  while (a<b) {
    std::cerr << "CPP 1 a=" << (a+1) << " b=" << (b+1) << "\n";
    while(true) {
      std::cerr << "CPP repeat loop on b\n";
      b--;
      if (b < B) {
        std::cerr << "CPP exit 1\n";
        return -1;
      }
      if (!test(b))
        break;
    }
    std::cerr << "CPP 2 a=" << (a+1) << " b=" << (b+1) << "\n";
    while(true) {
      a++;
      if (a>b || test(a))
        break;
    }
    std::cerr << "CPP 3 a=" << (a+1) << " b=" << (b+1) << "\n";
    if (a<b) {
      std::swap(P.points[a], P.points[b]);
    }
  }
  std::cerr << "CPP a=" << (a+1) << " l=" << (l+1) << "\n";
  if (a > l) {
    std::cerr << "CPP exit 2\n";
    return -1;
  }
  int m=P.firsts.size();
  for (int idx=a; idx<=l; idx++) {
    P.cellno[P.points[idx]] = m;
  }
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
  std::cerr << "CPP exit 3\n";
  return P.lengths[m];
}

template<typename Telt>
int SplitCell_Partition(Partition & P, int const& i, Partition const& Q, int const& j, Telt const& g, int const& out)
{
  std::cerr << "CPP Q=\n";
  RawPrintPartition(Q);
  std::function<bool(int)> test=[&](int const& ePt) -> bool {
    int fPt=PowAct(P.points[ePt], g);
    std::cerr << "CPP SplitCellTestfun1\n";
    return PointInCellNo(Q, fPt, j);
  };
  return SplitCell_Kernel(P, i, test, out);
}


template<typename Telt>
int SplitCell_Face(Partition & P, int const& i, Face const& f, int const& j, Telt const& g, int const& out)
{
  std::function<bool(int)> test=[&](int const& ePt) -> bool {
    int fPt=PowAct(P.points[ePt], g);
    if (j == 1) {
      if (f[fPt] == 1)
	return true;
      return false;
    }
    else {
      if (f[fPt] == 0)
	return true;
      return false;
    }
  };
  return SplitCell_Kernel(P, i, test, out);
}


int IsolatePoint(Partition & P, int const& a)
{
#ifdef CHECK_PARTITION
  CheckConsistencyPartition("Input IsolatePoint", P);
#endif
#ifdef DEBUG_PARTITION
  std::cerr << "CPP Input Partition\n";
  RawPrintPartition(P);
#endif
  int nbPart=P.firsts.size();
  int iPart=P.cellno[a];
  int eFirst=P.firsts[iPart];
  int len=P.lengths[iPart];
  //  std::cerr << "a=" << a << " iPart=" << iPart << " eFirst=" << eFirst << " len=" << len << "\n";
  if (len == 1)
    return -1;
  int pos=-1;
  for (int j=eFirst; j<eFirst + len; j++) {
    int ePt=P.points[j];
    if (ePt == a)
      pos=j;
  }
  int l=eFirst + len-1;
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



int UndoRefinement(Partition & P)
{
#ifdef CHECK_PARTITION
  CheckConsistencyPartition("Input UndoRefinement", P);
#endif
  int nbPart=P.firsts.size();
  int pfm=P.firsts[nbPart-1];
  if (pfm == 0)
    return -1;
  int plm=P.lengths[nbPart-1];
  int m=P.cellno[P.points[pfm-1]];
  P.lengths[m] += plm;
  for (int j=pfm; j<pfm + plm; j++) {
    int ePt=P.points[j];
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


int FixpointCellNo(Partition const& P, int const& i)
{
  return P.points[P.firsts[i]];
}




int FixcellPoint(Partition const& P, std::set<int> & old)
{
  int nbPart=P.lengths.size();
  std::vector<int> poss;
  for (int iPart=0; iPart<nbPart; iPart++) {
    if (P.lengths[iPart] == 1 && old.find(iPart) == old.end())
      poss.push_back(iPart);
  }
  int nbPoss=poss.size();
  if (nbPoss == 0)
    return -1;
  int idx=rand() % nbPoss;
  int p=poss[idx];
  old.insert(p);
  return p;
}


struct typeFixcellsCell {
  bool res;
  std::vector<int> K;
  std::vector<int> I;
};
 
typeFixcellsCell FixcellsCell(Partition const& P, Partition const& Q, std::set<int> & old)
{
  std::vector<int> K, I;
  int nbPart=P.firsts.size();
  for (int iPart=0; iPart<nbPart; iPart++) {
    int start=P.firsts[iPart];
    int kPart=CellNoPoint(Q, P.points[start]);
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


Partition TrivialPartition(std::vector<int> const& Omega)
{
  return GetPartition({Omega});
}

template<typename Telt>
Partition OrbitsPartition(std::vector<Telt> const& gens, int const&n, std::vector<int> const& Omega)
{
  return GetPartition(OrbitsPerms(gens, n, Omega));
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
 
template<typename Tarith>
Partition CollectedPartition(Partition const& P, Tarith const& size)
{
  Partition C=P;
  int div=SmallestPrimeDivisor(size);
  int nbPart=P.firsts.size();
  int nbPartTot=nbPart;
  for (int iPart=0; iPart<nbPart; iPart++) {
    int sizPart=P.lengths[iPart];
    if (sizPart < div) {
      int eFirst=P.firsts[iPart];
      for (int i=0; i<sizPart; i++) {
	int ePt=P.points[eFirst + i];
	if (i == 0) {
	  C.lengths[iPart]=1;
	}
	else {
	  C.cellno[ePt]=nbPartTot;
	  C.firsts.push_back(eFirst+i);
	  C.lengths.push_back(1);
	}
      }
    }
  }
  return C;
}

 
 
}

#endif
