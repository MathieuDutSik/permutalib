// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GAP_PARTITION_H_
#define SRC_GAP_PARTITION_H_

#include "GapPrint.h"
#include "PermGroup.h"
#include <limits>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef SYNCHRONIZED_DEBUG_GAP478
#define DEBUG_PARTITION
#endif

#ifdef PERMUTALIB_BLOCKING_SANITY_CHECK
#define CHECK_PARTITION
#endif

namespace permutalib {

template <typename Tidx> struct Partition {
  std::vector<Tidx> points;
  std::vector<Tidx> firsts;
  std::vector<Tidx> lengths;
  std::vector<Tidx> cellno;
};

/*
When working with partitions, we have to deal with the possibility
that the set of active points of the group is different from the set
of nodes being partitionned:
---The classical scenario is a group acting on vertices 0....n-1
but the set of moved points being e.g. {2,3}.
---The reverse occurs for intersections. We take the set Omega of all
points moved by G and H. Therefore there is movement outside of
Omega.



*/

//
// Print functionality
//

template <typename Tidx>
void NicePrintPartition(std::string const &str, Partition<Tidx> const &P) {
  Tidx nbPart = P.firsts.size();
  std::cerr << str << " = [ ";
  for (Tidx iPart = 0; iPart < nbPart; iPart++) {
    if (iPart > 0)
      std::cerr << ", ";
    Tidx eFirst = P.firsts[iPart];
    Tidx len = P.lengths[iPart];
    std::cerr << "[ ";
    for (Tidx i = 0; i < len; i++) {
      if (i > 0)
        std::cerr << ", ";
      Tidx ePt = P.points[eFirst + i];
      std::cerr << int(ePt + 1);
    }
    std::cerr << " ]";
  }
  std::cerr << " ]\n";
}

template <typename Tidx> void RawPrintPartition(Partition<Tidx> const &P) {
  Tidx nbPart = P.lengths.size();
  std::vector<std::string> LPart(nbPart);
  for (Tidx iPart = 0; iPart < nbPart; iPart++) {
    Tidx len = P.lengths[iPart];
    Tidx eFirst = P.firsts[iPart];
    std::vector<Tidx> eList(len);
    for (Tidx u = 0; u < len; u++) {
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

template <typename Tidx>
void CheckConsistencyPartition(std::string const &str,
                               Partition<Tidx> const &P) {
  Tidx nbError = 0;
  Tidx Max_points = 0;
  bool ChangeVal = false;
  for (auto &eVal : P.points)
    if (eVal > Max_points) {
      Max_points = eVal;
      ChangeVal = true;
    }
  if (ChangeVal)
    Max_points++;
  if (Tidx(P.cellno.size()) != Max_points) {
    std::cerr << "P.cellno.size()=" << P.cellno.size()
              << " Max_points=" << int(Max_points) << "\n";
    std::cerr << "1: We should have |P.cellno| = Maximum(P.points) (at least, "
                 "that's true at init)\n";
    nbError++;
  }
  if (P.firsts.size() != P.lengths.size()) {
    std::cerr << "2: P.firsts and P.lengths have different lengths\n";
    nbError++;
  }
  Tidx nbPoint = Tidx(P.cellno.size());
  Tidx nbPart = Tidx(P.lengths.size());
  std::vector<Tidx> MeasuredLength(nbPart, 0);
  for (Tidx iPoint = 0; iPoint < nbPoint; iPoint++) {
    Tidx iPart = P.cellno[iPoint];
    if ((iPart >= nbPart) && (iPart != std::numeric_limits<Tidx>::max())) {
      std::cerr << "3: Error, iPart=" << int(iPart)
                << " but nbPart=" << int(nbPart) << "\n";
      nbError++;
    }
    if (iPart >= 0 && iPart < nbPart)
      MeasuredLength[iPart]++;
  }
  for (Tidx iPart = 0; iPart < nbPart; iPart++) {
    if (MeasuredLength[iPart] != P.lengths[iPart]) {
      std::cerr << "4: At iPart=" << int(iPart)
                << " we have error in lengths: MeasuredLength[iPart]="
                << int(MeasuredLength[iPart])
                << " and P.lengths[iPart]=" << int(P.lengths[iPart]) << "\n";
      nbError++;
    }
    if (MeasuredLength[iPart] == 0) {
      std::cerr << "5: At iPart=" << int(iPart) << " the length is zero\n";
      nbError++;
    }
  }
  for (Tidx iPart = 0; iPart < nbPart; iPart++) {
    Tidx len = P.lengths[iPart];
    Tidx eFirst = P.firsts[iPart];
    for (Tidx u = 0; u < len; u++) {
      Tidx ePt = P.points[eFirst + u];
      if (ePt < 0 || ePt >= nbPoint) {
        std::cerr << "6: At iPart=" << int(iPart) << " u=" << int(u)
                  << " ePt=" << int(ePt) << " point out of range\n";
        nbError++;
      }
      if (ePt >= 0 && ePt < nbPoint) {
        if (P.cellno[ePt] != iPart) {
          std::cerr << "7: At iPart=" << int(iPart) << " u=" << int(u)
                    << " ePt=" << int(ePt) << " has wrong cellno\n";
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

//
// The real functionality
//

template <typename Tidx>
Partition<Tidx> GetPartition(std::vector<std::vector<Tidx>> const &list) {
#ifdef DEBUG_PARTITION
  std::vector<std::string> LStr;
  for (auto &eList : list)
    LStr.push_back(GapStringIntVector(eList));
  std::cerr << "CPP list=" << GapStringTVector(LStr) << "\n";
#endif
  size_t tot_siz = 0;
  for (auto &eList : list)
    tot_siz += eList.size();
  size_t pos_ins = 0;
  std::vector<Tidx> points(tot_siz);
  Tidx Max_NbPoint = 0;
  bool ChangeVal = false;
  for (auto &eList : list)
    for (auto &eVal : eList) {
      points[pos_ins] = eVal;
      pos_ins++;
      if (eVal > Max_NbPoint) {
        Max_NbPoint = eVal;
        ChangeVal = true;
      }
    }
  if (ChangeVal)
    Max_NbPoint++;
  Tidx nbPoint = Max_NbPoint;
  Tidx nbPart = Tidx(list.size());
  std::vector<Tidx> firsts(nbPart);
  std::vector<Tidx> lengths(nbPart);
  std::vector<Tidx> cellno(nbPoint, std::numeric_limits<Tidx>::max());
  Tidx i = 0;
  for (Tidx iPart = 0; iPart < nbPart; iPart++) {
    const std::vector<Tidx> &eL = list[iPart];
    Tidx len = Tidx(eL.size());
    firsts[iPart] = i;
    lengths[iPart] = len;
    i += len;
    for (auto &eVal : eL)
      cellno[eVal] = iPart;
  }
  Partition<Tidx> P{std::move(points), std::move(firsts), std::move(lengths),
                    std::move(cellno)};
#ifdef DEBUG_PARTITION
  std::cerr << "CPP After GetPartition operation P=\n";
  RawPrintPartition(P);
#endif
#ifdef CHECK_PARTITION
  CheckConsistencyPartition("GetPartition", P);
#endif
  return P;
}

template <typename Tidx>
Tidx GetNumberPoint(Partition<Tidx> const &ePartition) {
  Tidx nbPoint = Tidx(ePartition.points.size());
  return nbPoint;
}

template <typename Tidx> Tidx NumberCells(Partition<Tidx> const &ePartition) {
  Tidx nbPart = Tidx(ePartition.firsts.size());
  return nbPart;
}

template <typename Tidx>
std::vector<Tidx> Cell(Partition<Tidx> const &ePartition, Tidx const &iPart) {
  Tidx len = ePartition.lengths[iPart];
  Tidx eFirst = ePartition.firsts[iPart];
  std::vector<Tidx> eList(len);
  for (Tidx i = 0; i < len; i++)
    eList[i] = ePartition.points[eFirst + i];
  return eList;
}

template <typename Tidx>
std::vector<std::vector<Tidx>> Cells(Partition<Tidx> const &ePartition) {
  Tidx nbPart = ePartition.firsts.size();
  std::vector<std::vector<Tidx>> eListList(nbPart);
  for (Tidx iPart = 0; iPart < nbPart; iPart++)
    eListList[iPart] = Cell(ePartition, iPart);
  return eListList;
}

template <typename Tidx>
Tidx CellNoPoint(Partition<Tidx> const &ePartition, Tidx const &pt) {
  return ePartition.cellno[pt];
}

template <typename Tidx>
std::vector<Tidx> CellNoPoints(Partition<Tidx> const &ePartition,
                               std::vector<Tidx> const &pts) {
  Tidx len = pts.size();
  std::vector<Tidx> ret(len);
  for (Tidx i = 0; i < len; i++)
    ret[i] = ePartition.cellno[pts[i]];
  return ret;
}

template <typename Tidx>
bool PointInCellNo(Partition<Tidx> const &ePartition, Tidx const &pt,
                   Tidx const &iPart) {
  return ePartition.cellno[pt] == iPart;
}

template <typename Tidx>
std::vector<Tidx> Fixcells(Partition<Tidx> const &ePartition) {
#ifdef DEBUG_PARTITION
  std::cerr << "CPP beginning of Fixcells\n";
  RawPrintPartition(ePartition);
#endif
  std::vector<Tidx> fix;
  Tidx nbPart = Tidx(ePartition.firsts.size());
  for (Tidx iPart = 0; iPart < nbPart; iPart++) {
    if (ePartition.lengths[iPart] == 1) {
      Tidx eFirst = ePartition.firsts[iPart];
      Tidx eVal = ePartition.points[eFirst];
      fix.push_back(eVal);
    }
  }
  return fix;
}

template <typename Tidx, typename Ftest>
Tidx SplitCell_Kernel(Partition<Tidx> &P, Tidx const &i, Ftest test,
                      Tidx const &out) {
#ifdef DEBUG_PARTITION
  //  std::cerr << "CPP i=" << int(i+1) << " out=" << out << "\n";
  std::cerr << "CPP Before SplitCell_Kernel operation P=\n";
  RawPrintPartition(P);
#endif
#ifdef CHECK_PARTITION
  CheckConsistencyPartition("Input SplitCell_Kernel", P);
#endif
  // We need to shift the index by 1 because we are decreasing the index by 1
  // with a--. Since we can have Tidx = uint32_t this gets us to a=65536
  // and thus breaks a < b comparisons.
  Tidx a = P.firsts[i] + 1;
#ifdef DEBUG_PARTITION
  std::cerr << "CPP SplitCell i=" << (i + 1) << " a=" << int(a) << "\n";
#endif
  Tidx b = a + P.lengths[i];
  Tidx l = b - 1;

  Tidx maxmov;
  if (out != std::numeric_limits<Tidx>::max())
    maxmov = out;
  else
    maxmov = P.lengths[i] - 1;
  Tidx B = l - maxmov;
#ifdef DEBUG_PARTITION
  std::cerr << "CPP maxmov=" << maxmov << " B=" << int(B) << "\n";
#endif
  a--;
#ifdef DEBUG_PARTITION
  std::cerr << "CPP Before loop a=" << int(a) << " b=" << int(b) << "\n";
#endif
  while (a < b) {
#ifdef DEBUG_PARTITION
    std::cerr << "CPP     1 a=" << int(a) << " b=" << int(b) << "\n";
#endif
    while (true) {
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
      if (!test(b - 1))
        break;
    }
#ifdef DEBUG_PARTITION
    std::cerr << "CPP     2 a=" << int(a) << " b=" << int(b) << "\n";
#endif
    while (true) {
#ifdef DEBUG_PARTITION
      std::cerr << "CPP A LOOP\n";
#endif
      a++;
      if (a > b || test(a - 1))
        break;
    }
#ifdef DEBUG_PARTITION
    std::cerr << "CPP     3 a=" << int(a) << " b=" << int(b) << "\n";
#endif
    if (a < b) {
      std::swap(P.points[a - 1], P.points[b - 1]);
    }
  }
#ifdef DEBUG_PARTITION
  std::cerr << "CPP a=" << int(a) << " l=" << int(l) << "\n";
#endif
  if (a > l) {
#ifdef DEBUG_PARTITION
    std::cerr << "CPP exit 2\n";
#endif
    return std::numeric_limits<Tidx>::max();
  }
  Tidx m = Tidx(P.firsts.size());
  for (Tidx idx = a; idx <= l; idx++)
    P.cellno[P.points[idx - 1]] = m;
  P.firsts.push_back(a - 1);
  P.lengths.push_back((1 + l) - a);
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

template <typename Telt>
typename Telt::Tidx SplitCell_Partition_e(
    Partition<typename Telt::Tidx> &P, typename Telt::Tidx const &i,
    Partition<typename Telt::Tidx> const &Q, typename Telt::Tidx const &j,
    Telt const &g, typename Telt::Tidx const &out) {
  using Tidx = typename Telt::Tidx;
#ifdef DEBUG_PARTITION
  std::cerr << "CPP SplitCell g=" << g << "\n";
  std::cerr << "CPP Q=\n";
  RawPrintPartition(Q);
#endif
  auto test = [&](const Tidx &ePt) -> bool {
    Tidx fPt = PowAct(P.points[ePt], g);
#ifdef DEBUG_PARTITION
    std::cerr << "CPP SplitCellTestfun1 fPt=" << int(fPt + 1) << "\n";
#endif
    return PointInCellNo(Q, fPt, j);
  };
  return SplitCell_Kernel(P, i, test, out);
}

template <typename Tidx>
Tidx SplitCell_Partition(Partition<Tidx> &P, Tidx const &i,
                         Partition<Tidx> const &Q, Tidx const &j,
                         Tidx const &out) {
#ifdef DEBUG_PARTITION
  std::cerr << "CPP SplitCell g=()\n";
  std::cerr << "CPP Q=\n";
  RawPrintPartition(Q);
#endif
  auto test = [&](const Tidx &ePt) -> bool {
    Tidx fPt = P.points[ePt];
#ifdef DEBUG_PARTITION
    std::cerr << "CPP SplitCellTestfun1 fPt=" << int(fPt + 1) << "\n";
#endif
    return PointInCellNo(Q, fPt, j);
  };
  return SplitCell_Kernel(P, i, test, out);
}

template <typename Telt>
typename Telt::Tidx SplitCell_Face(Partition<typename Telt::Tidx> &P,
                                   typename Telt::Tidx const &i, Face const &f,
                                   typename Telt::Tidx const &j, Telt const &g,
                                   typename Telt::Tidx const &out) {
  using Tidx = typename Telt::Tidx;
  auto test = [&](const Tidx &ePt) -> bool {
    Tidx fPt = PowAct(P.points[ePt], g);
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

template <typename Tidx> Tidx IsolatePoint(Partition<Tidx> &P, Tidx const &a) {
#ifdef CHECK_PARTITION
  CheckConsistencyPartition("Input IsolatePoint", P);
#endif
#ifdef DEBUG_PARTITION
  std::cerr << "CPP Input Partition\n";
  RawPrintPartition(P);
#endif
  Tidx nbPart = Tidx(P.firsts.size());
  Tidx iPart = P.cellno[a];
  Tidx eFirst = P.firsts[iPart];
  Tidx len = P.lengths[iPart];
  //  std::cerr << "a=" << a << " iPart=" << iPart << " eFirst=" << eFirst << "
  //  len=" << len << "\n";
  if (len == 1)
    return std::numeric_limits<Tidx>::max();
  Tidx pos = 0;
  for (Tidx j = eFirst; j < eFirst + len; j++) {
    Tidx ePt = P.points[j];
    if (ePt == a)
      pos = j;
  }
  Tidx l = (eFirst + len) - 1;
  //  std::cerr << "pos=" << pos << " l=" << l << "\n";
  P.points[pos] = P.points[l];
  P.points[l] = a;
  //  int m=P.firsts.size() + 1;
  P.cellno[a] = nbPart;
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

template <typename Tidx> Tidx UndoRefinement(Partition<Tidx> &P) {
#ifdef CHECK_PARTITION
  CheckConsistencyPartition("Input UndoRefinement", P);
#endif
  Tidx nbPart = Tidx(P.firsts.size());
  Tidx pfm = P.firsts[nbPart - 1];
  if (pfm == 0)
    return std::numeric_limits<Tidx>::max();
  Tidx plm = P.lengths[nbPart - 1];
  Tidx m = P.cellno[P.points[pfm - 1]];
  P.lengths[m] += plm;
  for (Tidx j = pfm; j < pfm + plm; j++) {
    Tidx ePt = P.points[j];
    P.cellno[ePt] = m;
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

template <typename Tidx>
Tidx FixpointCellNo(Partition<Tidx> const &P, Tidx const &i) {
  return P.points[P.firsts[i]];
}

/*
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
  Tidx idx=random() % nbPoss;
  Tidx p=poss[idx];
  old.insert(p);
  return p;
}
*/

template <typename Tidx> struct typeFixcellsCell {
  bool res;
  std::vector<Tidx> K;
  std::vector<Tidx> I;
};

template <typename Tidx>
typeFixcellsCell<Tidx> FixcellsCell(Partition<Tidx> const &P,
                                    Partition<Tidx> const &Q,
                                    std::set<Tidx> &old) {
  std::vector<Tidx> K, I;
  Tidx nbPart = P.firsts.size();
  for (Tidx iPart = 0; iPart < nbPart; iPart++) {
    Tidx start = P.firsts[iPart];
    Tidx kPart = CellNoPoint(Q, P.points[start]);
    if (old.find(kPart) == old.end()) {
      auto eval = [&]() -> bool {
        for (int j = 1; j < P.lengths[iPart]; j++) {
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
    return {false, {}, {}};
  } else {
    return {true, std::move(K), std::move(I)};
  }
}

template <typename Tidx>
Partition<Tidx> TrivialPartition(std::vector<Tidx> const &Omega) {
  return GetPartition(std::vector<std::vector<Tidx>>({Omega}));
}

template <typename Telt>
std::vector<std::vector<typename Telt::Tidx>>
OrbitsPermsB(std::vector<Telt> const &gens,
             std::vector<typename Telt::Tidx> const &Omega) {
  using Tidx = typename Telt::Tidx;
#ifdef DEBUG_PARTITION
  std::cerr << "DEBUG OrbitsPermB beginning\n";
#endif
  Tidx max = LargestMovedPoint(gens);
  Tidx maxP1 = 0;
  if (max != std::numeric_limits<Tidx>::max())
    maxP1 = max + 1;
  Face dom(maxP1);
  for (auto &eVal : Omega)
    if (eVal < maxP1)
      dom[eVal] = 1;
#ifdef DEBUG_PARTITION
  std::cerr << "DEBUG dom built\n";
#endif
  Face newF(maxP1);
  for (Tidx i = 0; i < maxP1; i++)
    newF[i] = 1;
#ifdef DEBUG_PARTITION
  std::cerr << "DEBUG newF built\n";
#endif
  std::vector<std::vector<Tidx>> orbs;
  boost::dynamic_bitset<>::size_type fst = dom.find_first();
  while (fst != boost::dynamic_bitset<>::npos) {
    Tidx fst_i = Tidx(fst);
    std::vector<Tidx> orb{fst_i};
    newF[fst_i] = 0;
    dom[fst_i] = 0;
#ifdef DEBUG_PARTITION
    std::cerr << "DEBUG Before beginning of pnt loop\n";
#endif
    size_t posOrb = 0;
    while (true) {
      Tidx pnt = orb[posOrb];
      for (auto &gen : gens) {
        Tidx img = PowAct(pnt, gen);
        if (newF[img]) {
          orb.push_back(img);
          newF[img] = 0;
          dom[img] = 0;
        }
      }
      posOrb++;
      if (posOrb >= orb.size())
        break;
    }
#ifdef DEBUG_PARTITION
    std::cerr << "DEBUG After the pnt loop\n";
#endif
    orbs.emplace_back(std::move(orb));
    fst = dom.find_first();
  }
  for (auto &pnt : Omega)
    if (pnt >= maxP1)
      orbs.push_back({pnt});
  return orbs;
}

template <typename Telt>
Partition<typename Telt::Tidx>
OrbitsPartition(std::vector<Telt> const &gens,
                std::vector<typename Telt::Tidx> const &Omega) {
#ifdef DEBUG_PARTITION
  std::cerr << "CPP OrbitsPartition, using OrbitsPerms\n";
  std::cerr << "CPP generators=" << GapStringTVector(gens) << "\n";
#endif
  return GetPartition(OrbitsPermsB(gens, Omega));
}

template <typename Tarith> int SmallestPrimeDivisor(Tarith const &size) {
  if (size == 1)
    return 1;
  int i = 2;
  while (true) {
    Tarith res = size % Tarith(i);
    if (res == 0)
      return i;
    i++;
  }
}

template <typename Tarith, typename Tidx>
Partition<Tidx> CollectedPartition(Partition<Tidx> const &P,
                                   Tarith const &size) {
  Tidx div = SmallestPrimeDivisor(size);
#ifdef DEBUG_PARTITION
  std::cerr << "CPP div=" << div << "\n";
#endif
  Tidx nbPart = P.firsts.size();
  std::unordered_map<Tidx, std::vector<Tidx>> map;
  for (Tidx iPart = 0; iPart < nbPart; iPart++) {
    Tidx sizPart = P.lengths[iPart];
    map[sizPart].push_back(iPart);
  }
  std::vector<std::vector<Tidx>> C;
  for (auto &kv : map) {
    Tidx n_block = kv.second.size();
    if (n_block < div) {
      for (auto &iPart : kv.second) {
        C.push_back(Cell(P, iPart));
      }
    } else {
      std::vector<Tidx> NewConn;
      for (auto &iPart : kv.second) {
        for (auto &val : Cell(P, iPart))
          NewConn.push_back(val);
      }
      C.push_back(NewConn);
    }
  }
  return GetPartition(C);
}

// clang-format off
}  // namespace permutalib
#endif  // SRC_GAP_PARTITION_H_
// clang-format on
