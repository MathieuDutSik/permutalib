// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GAP_RESOLUTION_H_
#define SRC_GAP_RESOLUTION_H_

// clang-format off
#include "Group.h"
#include <algorithm>
#include <cstdint>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
// clang-format on

/*
  Free ZG-resolutions of the trivial module Z for a finite group G.

  This is a C++ transcription of the function ResolutionFiniteGroup of the
  GAP package HAP (Graham Ellis). The algorithm is the one described in
    G. Ellis, Computing group resolutions,
    Journal of Symbolic Computation 38 (2004) 1077-1118.
  The construction works in the universal cover X~ of a classifying space X
  of G, dimension by dimension:
  --- The (i-1)-cells of X~ are the pairs (j, g) with j a free generator of
      the (i-1)-th term of the resolution and g an element of G.
  --- A "maximal contractible subcomplex" Y(i-1) is stored as the set of
      (i-1)-cells marked in max_complex[i]. It is grown by pairing each
      (i-1)-cell with an i-cell whose boundary contains it exactly once and
      otherwise only cells that are already in the subcomplex (the pairing is
      recorded in contraction_matrix[i], a discrete vector field).
  --- When no (i-1)-cell can be paired anymore, the (i-1)-cell whose
      "differential" is shortest is chosen as the boundary of a new free
      generator of the i-th term of the resolution.
  The names of the internal functions follow the GAP code so that both can be
  compared line by line: Contraction, Differential, CellValue, FirstZero,
  FindConsequences and NextResTerm.

  Conventions:
  --- Elements of G are enumerated in a sorted list, with the identity first
      (as in GAP, Elements(G) is sorted and the identity is the smallest
      permutation), so the output can be compared with HAP for the same
      generators.
  --- A chain in the i-th term of the resolution is a list of terms
      coeff * (elt . e^i_cell) with coeff a nonzero integer, cell the 0-based
      index of the free generator and elt the 0-based index of the group
      element in elements(). Chains are always kept reduced: terms are sorted
      by (cell, elt) and distinct.
  --- The product of two elements is the product of permutations of GAP and
      permutalib: (g*h)(x) = h(g(x)).

  What is provided:
  --- dimension(i): number of free generators of the i-th term (rank of R_i).
  --- boundary(i, j): the boundary of the j-th generator of R_i, a chain of
      R_{i-1}.
  --- homotopy(i, t): the contracting homotopy of the resolution, mapping an
      element of R_i to R_{i+1}. It satisfies d h + h d = id on R_i for i>=1
      and d h(t) = t - eps(t) e^0_0 on R_0, with eps the augmentation.
  --- integral_homology(n): the abelian invariants of H_n(G, Z), obtained by
      tensoring with Z and Smith normal form. The invariants are returned in
      the divisibility form d_1 | d_2 | ... which is the output format of
      Homology(TensorWithIntegers(R), n) in HAP. The form of AbelianInvariants
      of GAP (prime powers) is given by integral_homology_prime_power(n).
  --- homology_dimension_mod_p(n, p): the dimension of H_n(G, F_p), computed
      independently by rank computations modulo p.
  --- check_homology_consistency(n): consistency checks of the homology (see
      below). They are run automatically in integral_homology when
      DEBUG_RESOLUTION is defined.
  --- sizes(): the sum of the lengths of the boundaries in each degree, that
      is Size(R) in HAP.
  --- extend(): computes one more term of the resolution.

  The consistency checks of check_homology_consistency(n) are:
  --- H_0(G, Z) = Z and, for n >= 1, H_n(G, Z) is finite with every
      invariant dividing the order of G.
  --- The universal coefficient theorem: for every prime p dividing |G|,
      dim H_n(G, F_p) = t_p(H_n) + t_p(H_{n-1}) with t_p the number of
      torsion invariants divisible by p (the free factor Z of H_0 counts
      for the tensor term but not for the Tor term).
      For a prime not dividing |G|, dim H_n(G, F_p) = 0 for n >= 1.
  --- H_1(G, Z) is the abelianization: its order is |G| / |[G,G]| with the
      derived subgroup computed by the stabilizer chain code.

  The resolution of HAP depends on whether it is "extendible": in the last
  term, the non extendible version (the default of HAP) scans the generators
  in increasing order only when pairing cells, while the extendible version
  scans them in increasing then decreasing order. The constructor has the
  same option so that both outputs of HAP can be reproduced; extend() is
  available in both cases.

  What is not transcribed from HAP:
  --- The Tietze reduction option (arg[3] of ResolutionFiniteGroup).
  --- The characteristic p reductions (arg[4]), i.e. resolutions over F_p G.
  --- The matrix group case (handled in HAP by isomorphism to a permutation
      group). Any permutation type Telt of permutalib can be used here.
 */

#ifdef DEBUG
#define DEBUG_RESOLUTION
#endif

namespace permutalib {

// One term coeff * (elt . e_cell) of a free ZG-module.
struct ResolutionTerm {
  int64_t coeff;
  size_t cell;
  size_t elt;
};

using ResolutionChain = std::vector<ResolutionTerm>;

// Sort the terms by (cell, elt), add the coefficients of identical
// (cell, elt) and remove the terms with zero coefficient.
inline void ReduceResolutionChain(ResolutionChain &chain) {
  std::sort(chain.begin(), chain.end(),
            [](ResolutionTerm const &a, ResolutionTerm const &b) {
              if (a.cell != b.cell)
                return a.cell < b.cell;
              return a.elt < b.elt;
            });
  size_t len = chain.size();
  size_t pos_write = 0;
  size_t pos_read = 0;
  while (pos_read < len) {
    ResolutionTerm term = chain[pos_read];
    pos_read++;
    while (pos_read < len && chain[pos_read].cell == term.cell &&
           chain[pos_read].elt == term.elt) {
      term.coeff += chain[pos_read].coeff;
      pos_read++;
    }
    if (term.coeff != 0) {
      chain[pos_write] = term;
      pos_write++;
    }
  }
  chain.resize(pos_write);
}

inline ResolutionChain NegateResolutionChain(ResolutionChain const &chain) {
  ResolutionChain ret;
  ret.reserve(chain.size());
  for (auto const &term : chain) {
    ret.push_back({-term.coeff, term.cell, term.elt});
  }
  return ret;
}

// Number of terms of the chain counted with multiplicity, that is the length
// of the corresponding HAP word.
inline size_t LengthResolutionChain(ResolutionChain const &chain) {
  size_t len = 0;
  for (auto const &term : chain) {
    if (term.coeff > 0)
      len += size_t(term.coeff);
    else
      len += size_t(-term.coeff);
  }
  return len;
}

// The string of the chain in the format of the words of HAP, that is a list
// of pairs [k, g] with k the signed 1-based index of the generator and g the
// 1-based index of the group element. A term of coefficient c is repeated
// |c| times, which is exactly the output of AlgebraicReduction in HAP.
inline std::string GapStringResolutionChain(ResolutionChain const &chain) {
  std::string str = "[ ";
  bool IsFirst = true;
  for (auto const &term : chain) {
    int64_t mult = term.coeff > 0 ? term.coeff : -term.coeff;
    int64_t k = term.coeff > 0 ? int64_t(term.cell + 1) : -int64_t(term.cell + 1);
    for (int64_t u = 0; u < mult; u++) {
      if (!IsFirst)
        str += ", ";
      IsFirst = false;
      str += "[ " + std::to_string(k) + ", " + std::to_string(term.elt + 1) + " ]";
    }
  }
  str += " ]";
  return str;
}

// The invariants of the Smith normal form of an integer matrix, that is the
// positive diagonal entries d_1 | d_2 | ... | d_r with r the rank.
// The matrix is given as a list of rows.
template <typename Tint>
std::vector<Tint> SmithNormalFormInvariants(std::vector<std::vector<Tint>> M) {
  size_t n_row = M.size();
  if (n_row == 0)
    return {};
  size_t n_col = M[0].size();
  if (n_col == 0)
    return {};
  auto abs_val = [](Tint const &x) -> Tint {
    if (x < 0)
      return -x;
    return x;
  };
  auto swap_rows = [&](size_t i1, size_t i2) {
    if (i1 != i2)
      std::swap(M[i1], M[i2]);
  };
  auto swap_cols = [&](size_t j1, size_t j2) {
    if (j1 != j2)
      for (size_t i = 0; i < n_row; i++)
        std::swap(M[i][j1], M[i][j2]);
  };
  std::vector<Tint> invariants;
  size_t t = 0;
  size_t n_min = std::min(n_row, n_col);
  while (t < n_min) {
    // Find the nonzero entry of the remaining submatrix of minimal absolute
    // value and move it to position (t, t).
    bool has_pivot = false;
    size_t i_piv = t;
    size_t j_piv = t;
    Tint val_piv = 0;
    for (size_t i = t; i < n_row; i++) {
      for (size_t j = t; j < n_col; j++) {
        if (M[i][j] != 0) {
          Tint a = abs_val(M[i][j]);
          if (!has_pivot || a < val_piv) {
            has_pivot = true;
            val_piv = a;
            i_piv = i;
            j_piv = j;
          }
        }
      }
    }
    if (!has_pivot)
      break;
    swap_rows(t, i_piv);
    swap_cols(t, j_piv);
    // Clear the row and the column of the pivot. Whenever a remainder is
    // nonzero it becomes the new (smaller) pivot, so this terminates.
    while (true) {
      bool clean = true;
      for (size_t i = t + 1; i < n_row; i++) {
        if (M[i][t] != 0) {
          Tint q = M[i][t] / M[t][t];
          if (q != 0)
            for (size_t j = t; j < n_col; j++)
              M[i][j] -= q * M[t][j];
          if (M[i][t] != 0) {
            swap_rows(t, i);
            clean = false;
            break;
          }
        }
      }
      if (!clean)
        continue;
      for (size_t j = t + 1; j < n_col; j++) {
        if (M[t][j] != 0) {
          Tint q = M[t][j] / M[t][t];
          if (q != 0)
            for (size_t i = t; i < n_row; i++)
              M[i][j] -= q * M[i][t];
          if (M[t][j] != 0) {
            swap_cols(t, j);
            clean = false;
            break;
          }
        }
      }
      if (!clean)
        continue;
      // The pivot must divide all the remaining entries. If some entry is
      // not divisible, add its row to the pivot row and clean again.
      bool divisible = true;
      for (size_t i = t + 1; i < n_row && divisible; i++) {
        for (size_t j = t + 1; j < n_col; j++) {
          if (M[i][j] % M[t][t] != 0) {
            for (size_t jj = t; jj < n_col; jj++)
              M[t][jj] += M[i][jj];
            divisible = false;
            break;
          }
        }
      }
      if (divisible)
        break;
    }
    if (M[t][t] < 0)
      M[t][t] = -M[t][t];
    invariants.push_back(M[t][t]);
    t++;
  }
  return invariants;
}

// The rank modulo a prime p of an integer matrix given as a list of rows.
template <typename Tint>
size_t RankModP(std::vector<std::vector<Tint>> M, Tint const &p) {
  size_t n_row = M.size();
  if (n_row == 0)
    return 0;
  size_t n_col = M[0].size();
  auto reduce = [&](Tint const &x) -> Tint {
    Tint r = x % p;
    if (r < 0)
      r += p;
    return r;
  };
  auto inverse = [&](Tint const &a) -> Tint {
    // Extended Euclid on (a, p) with 0 < a < p.
    Tint old_r = a, r = p, old_s = 1, s = 0;
    while (r != 0) {
      Tint q = old_r / r;
      Tint tmp = old_r - q * r;
      old_r = r;
      r = tmp;
      tmp = old_s - q * s;
      old_s = s;
      s = tmp;
    }
    return reduce(old_s);
  };
  for (size_t i = 0; i < n_row; i++)
    for (size_t j = 0; j < n_col; j++)
      M[i][j] = reduce(M[i][j]);
  size_t rank = 0;
  for (size_t j = 0; j < n_col && rank < n_row; j++) {
    size_t i_piv = rank;
    while (i_piv < n_row && M[i_piv][j] == 0)
      i_piv++;
    if (i_piv == n_row)
      continue;
    std::swap(M[rank], M[i_piv]);
    Tint inv = inverse(M[rank][j]);
    for (size_t jj = j; jj < n_col; jj++)
      M[rank][jj] = reduce(M[rank][jj] * inv);
    for (size_t i = 0; i < n_row; i++) {
      if (i != rank && M[i][j] != 0) {
        Tint c = M[i][j];
        for (size_t jj = j; jj < n_col; jj++)
          M[i][jj] = reduce(M[i][jj] - c * M[rank][jj]);
      }
    }
    rank++;
  }
  return rank;
}

// The prime factors of a positive integer, by trial division.
template <typename Tint> std::vector<Tint> PrimeFactors(Tint n) {
  std::vector<Tint> ret;
  Tint p = 2;
  while (p * p <= n) {
    if (n % p == 0) {
      ret.push_back(p);
      while (n % p == 0)
        n /= p;
    }
    p += 1;
  }
  if (n > 1)
    ret.push_back(n);
  return ret;
}

// The abelian invariants in the form of AbelianInvariants of GAP: the
// prime power factors of the given invariants sorted increasingly, followed
// by the zeros.
template <typename Tint>
std::vector<Tint> PrimePowerAbelianInvariants(std::vector<Tint> const &inv) {
  std::vector<Tint> ret;
  size_t n_zero = 0;
  for (auto const &d : inv) {
    if (d == 0) {
      n_zero++;
      continue;
    }
    Tint n = d;
    for (auto const &p : PrimeFactors(d)) {
      Tint q = 1;
      while (n % p == 0) {
        n /= p;
        q *= p;
      }
      ret.push_back(q);
    }
  }
  std::sort(ret.begin(), ret.end());
  for (size_t u = 0; u < n_zero; u++)
    ret.push_back(0);
  return ret;
}

template <typename Telt, typename Tint> class FiniteGroupResolution {
public:
  using Tidx = typename Telt::Tidx;
  using Tgroup = Group<Telt, Tint>;

private:
  // The pairing of an (i-1)-cell in the discrete vector field. In HAP this
  // is ContractionMatrix[i][j][g]: 1 if the cell is in the initial maximal
  // complex (state 1), or a pair [sign*j', g'] describing the i-cell through
  // which the cell is contracted (state 2).
  struct Pairing {
    uint8_t state; // 0: unset, 1: in initial complex, 2: paired with a cell
    int sign;
    size_t cell;
    size_t elt;
  };
  static constexpr size_t max_order_multiplication_table = 2048;
  // The group and its elements
  size_t N;
  std::vector<Telt> elts;
  std::unordered_map<Telt, size_t> map_position;
  bool use_table;
  std::vector<uint32_t> mult_table;
  std::vector<size_t> extended_elts;
  // The number of computed terms and the length requested at construction
  size_t K;
  size_t K_target;
  bool extendible;
  // The generators and identity, kept for the consistency checks
  std::vector<Telt> l_gen_store;
  Telt id_store;
  // pseudo_boundary[i][j] is the boundary of the j-th generator of R_i.
  std::vector<std::vector<ResolutionChain>> pseudo_boundary;
  // max_complex[i][j][g] = 1 if the (i-1)-cell (j,g) belongs to the maximal
  // contractible subcomplex Y(i-1).
  std::vector<std::vector<std::vector<uint8_t>>> max_complex;
  // contraction_matrix[i][j][g] is the pairing of the (i-1)-cell (j,g).
  std::vector<std::vector<std::vector<Pairing>>> contraction_matrix;
  // computed_contractions[i][j][g] is the memoized value of
  // Contraction(i, +(j,g)), and contraction_length[i][j][g] its length as
  // an unreduced HAP word (used for the choice of the next generator).
  std::vector<std::vector<std::vector<std::optional<ResolutionChain>>>>
      computed_contractions;
  std::vector<std::vector<std::vector<size_t>>> contraction_length;

  void build_elements(std::vector<Telt> const &l_gen, Telt const &id) {
    Tgroup G(l_gen, id);
    elts = G.get_all_element();
    std::sort(elts.begin(), elts.end());
    N = elts.size();
    for (size_t i = 0; i < N; i++) {
      if (elts[i] == id) {
        if (i != 0)
          std::swap(elts[0], elts[i]);
        break;
      }
    }
    if (N == 0 || !(elts[0] == id)) {
      std::cerr << "RES: The identity is missing from the list of elements\n";
      throw PermutalibException{1};
    }
    for (size_t i = 0; i < N; i++) {
      map_position[elts[i]] = i;
    }
    use_table = N <= max_order_multiplication_table;
    if (use_table) {
      mult_table.resize(N * N);
      for (size_t g = 0; g < N; g++) {
        for (size_t h = 0; h < N; h++) {
          Telt prod = elts[g] * elts[h];
          mult_table[g * N + h] = uint32_t(map_position.at(prod));
        }
      }
    }
    // The generators are sorted and the identity is removed, as in HAP.
    std::vector<Telt> gens = l_gen;
    std::sort(gens.begin(), gens.end());
    gens.erase(std::unique(gens.begin(), gens.end()), gens.end());
    for (auto const &g : gens) {
      if (!(g == id))
        extended_elts.push_back(map_position.at(g));
    }
    for (size_t g = 0; g < N; g++)
      extended_elts.push_back(g);
  }

  void init_state() {
    K = 0;
    pseudo_boundary.resize(2);
    max_complex.resize(2);
    contraction_matrix.resize(2);
    computed_contractions.resize(2);
    contraction_length.resize(2);
    // The 0-cells are the group elements; the identity vertex is the initial
    // maximal complex Y(0).
    std::vector<uint8_t> row(N, 0);
    row[0] = 1;
    max_complex[1] = {row};
  }

  // Index of elts[g] * elts[h]
  size_t product_index(size_t g, size_t h) const {
    if (use_table)
      return mult_table[g * N + h];
    Telt prod = elts[g] * elts[h];
    return map_position.at(prod);
  }

  // The boundary of the i-cell sign*(j, g): the chain g.d(e^i_j) up to sign.
  // Terms are not reduced but are distinct since pseudo_boundary is reduced.
  ResolutionChain boundary_cell(size_t i, int sign, size_t j, size_t g) const {
    ResolutionChain ret;
    if (i == 0)
      return ret;
    ResolutionChain const &bnd = pseudo_boundary[i][j];
    ret.reserve(bnd.size());
    for (auto const &term : bnd) {
      ret.push_back({sign * term.coeff, term.cell, product_index(g, term.elt)});
    }
    return ret;
  }

  // HAP: CellValue. Number v of terms (with multiplicity) in the boundary of
  // the i-cell (j,g) that are not in MC and, when v = 1, that term.
  std::pair<size_t, ResolutionTerm>
  cell_value(size_t i, std::vector<std::vector<uint8_t>> const &MC, size_t j,
             size_t g) const {
    size_t v = 0;
    ResolutionTerm q{0, 0, 0};
    for (auto const &term : pseudo_boundary[i][j]) {
      size_t elt = product_index(g, term.elt);
      if (MC[term.cell][elt] == 0) {
        v += size_t(term.coeff > 0 ? term.coeff : -term.coeff);
        q = {term.coeff, term.cell, elt};
      }
    }
    return {v, q};
  }

  // HAP: Contraction(i, x) for the (i-1)-cell x = sign*(j, g). Returns a
  // chain c of i-cells with d(c) = -x modulo the (i-1)-cells of Y(i-1).
  ResolutionChain contraction(size_t i, int sign, size_t j, size_t g) {
    if (i < 1) {
      return {{-1, 0, 0}};
    }
    std::optional<ResolutionChain> &memo = computed_contractions[i][j][g];
    if (!memo) {
      Pairing const &pr = contraction_matrix[i][j][g];
      if (pr.state == 1) {
        memo = ResolutionChain();
        contraction_length[i][j][g] = 0;
      } else {
        if (pr.state != 2) {
          std::cerr << "RES: Contraction called on an unpaired cell\n";
          throw PermutalibException{1};
        }
        // For x positive: c = -m + sum over the terms y of -d(m), y != x, of
        // Contraction(i, y).
        ResolutionChain c;
        c.push_back({-pr.sign, pr.cell, pr.elt});
        size_t len = 1;
        for (auto const &term : pseudo_boundary[i][pr.cell]) {
          size_t elt = product_index(pr.elt, term.elt);
          if (term.cell == j && elt == g)
            continue;
          int64_t coeff = -pr.sign * term.coeff;
          int s = coeff > 0 ? 1 : -1;
          int64_t mult = coeff > 0 ? coeff : -coeff;
          ResolutionChain sub = contraction(i, s, term.cell, elt);
          for (auto const &sub_term : sub) {
            c.push_back({mult * sub_term.coeff, sub_term.cell, sub_term.elt});
          }
          len += size_t(mult) * contraction_length[i][term.cell][elt];
        }
        ReduceResolutionChain(c);
        computed_contractions[i][j][g] = std::move(c);
        contraction_length[i][j][g] = len;
      }
    }
    if (sign > 0)
      return *memo;
    return NegateResolutionChain(*memo);
  }

  // The HAP length of the word Differential(i, (j,g)) before reduction.
  size_t differential_length(size_t i, size_t j, size_t g) {
    size_t len = 1;
    if (i == 1) {
      len += 1;
    } else {
      for (auto const &term : pseudo_boundary[i - 1][j]) {
        size_t elt = product_index(g, term.elt);
        int s = term.coeff > 0 ? 1 : -1;
        int64_t mult = term.coeff > 0 ? term.coeff : -term.coeff;
        // Forces the computation of the memoized contraction
        (void)contraction(i - 1, s, term.cell, elt);
        len += size_t(mult) * contraction_length[i - 1][term.cell][elt];
      }
    }
    return len;
  }

  // HAP: Differential(i, p). The boundary of the new i-cell attached along
  // the (i-1)-cell p = (j, g): the chain p + sum over the terms x of d(p) of
  // Contraction(i-1, x). It is a cycle of R_{i-1}.
  ResolutionChain differential(size_t i, size_t j, size_t g) {
    ResolutionChain diff;
    diff.push_back({1, j, g});
    if (i == 1) {
      diff.push_back({-1, 0, 0});
    } else {
      for (auto const &term : pseudo_boundary[i - 1][j]) {
        size_t elt = product_index(g, term.elt);
        int s = term.coeff > 0 ? 1 : -1;
        int64_t mult = term.coeff > 0 ? term.coeff : -term.coeff;
        ResolutionChain sub = contraction(i - 1, s, term.cell, elt);
        for (auto const &sub_term : sub) {
          diff.push_back({mult * sub_term.coeff, sub_term.cell, sub_term.elt});
        }
      }
    }
    ReduceResolutionChain(diff);
    return diff;
  }

  // HAP: FirstZero. Among the (i-1)-cells not in MC, the first one whose
  // differential is shortest.
  std::pair<size_t, size_t>
  first_zero(size_t i, std::vector<std::vector<uint8_t>> const &MC,
             std::vector<std::vector<size_t>> &diff_lengths) {
    bool found = false;
    size_t len_min = 0;
    std::pair<size_t, size_t> ret{0, 0};
    size_t dim = MC.size();
    for (size_t j = 0; j < dim; j++) {
      for (size_t g = 0; g < N; g++) {
        if (MC[j][g] == 0) {
          if (diff_lengths[j][g] == 0)
            diff_lengths[j][g] = differential_length(i, j, g);
          size_t len = diff_lengths[j][g];
          if (!found || len < len_min) {
            found = true;
            len_min = len;
            ret = {j, g};
          }
        }
      }
    }
    if (!found) {
      std::cerr << "RES: FirstZero called with a complete complex\n";
      throw PermutalibException{1};
    }
    return ret;
  }

  // HAP: FindConsequences. Pair every (i-1)-cell that can be contracted
  // through an i-cell whose other boundary terms are already in MC.
  void find_consequences(size_t i, std::vector<std::vector<uint8_t>> &MC,
                         size_t &n_zero) {
    size_t dim = pseudo_boundary[i].size();
    std::vector<size_t> iterset;
    for (size_t j = 0; j < dim; j++)
      iterset.push_back(j);
    if (i < K_target || extendible) {
      for (size_t j = dim; j > 0; j--)
        iterset.push_back(j - 1);
    }
    bool toggle = true;
    while (toggle) {
      toggle = false;
      for (size_t g : extended_elts) {
        for (size_t j : iterset) {
          std::pair<size_t, ResolutionTerm> cv = cell_value(i, MC, j, g);
          if (cv.first == 1) {
            ResolutionTerm const &p = cv.second;
            MC[p.cell][p.elt] = 1;
            n_zero--;
            int sign = p.coeff > 0 ? 1 : -1;
            contraction_matrix[i][p.cell][p.elt] = {2, sign, j, g};
            max_complex[i + 1][j][g] = 1;
            toggle = true;
          }
        }
      }
    }
  }

  // HAP: NextResTerm. Computes the i-th term of the resolution.
  void next_res_term(size_t i) {
    if (pseudo_boundary.size() < i + 2) {
      pseudo_boundary.resize(i + 2);
      max_complex.resize(i + 2);
      contraction_matrix.resize(i + 2);
      computed_contractions.resize(i + 2);
      contraction_length.resize(i + 2);
    }
    size_t dim_prev = max_complex[i].size();
    pseudo_boundary[i].clear();
    max_complex[i + 1].clear();
    contraction_matrix[i].assign(dim_prev, std::vector<Pairing>(N, Pairing{0, 0, 0, 0}));
    computed_contractions[i].assign(dim_prev, std::vector<std::optional<ResolutionChain>>(N));
    contraction_length[i].assign(dim_prev, std::vector<size_t>(N, 0));
    std::vector<std::vector<uint8_t>> MC = max_complex[i];
    size_t n_zero = 0;
    for (size_t j = 0; j < dim_prev; j++) {
      for (size_t g = 0; g < N; g++) {
        if (MC[j][g] == 1)
          contraction_matrix[i][j][g].state = 1;
        else
          n_zero++;
      }
    }
    std::vector<std::vector<size_t>> diff_lengths(dim_prev, std::vector<size_t>(N, 0));
    std::vector<uint8_t> row(N, 0);
    row[0] = 1;
    while (n_zero > 0) {
      std::pair<size_t, size_t> p = first_zero(i, MC, diff_lengths);
      ResolutionChain diff = differential(i, p.first, p.second);
      pseudo_boundary[i].push_back(std::move(diff));
      max_complex[i + 1].push_back(row);
      MC[p.first][p.second] = 1;
      n_zero--;
      contraction_matrix[i][p.first][p.second] = {2, 1, pseudo_boundary[i].size() - 1, 0};
      find_consequences(i, MC, n_zero);
    }
  }

public:
  // Resolution of length K of the group generated by l_gen. With
  // _extendible = false the output is the one of ResolutionFiniteGroup(G, K)
  // in HAP, with _extendible = true the one of
  // ResolutionFiniteGroup(G, K, false, 0, "extendible").
  FiniteGroupResolution(std::vector<Telt> const &l_gen, Telt const &id, size_t _K,
                        bool _extendible = false)
      : K_target(_K), extendible(_extendible), l_gen_store(l_gen), id_store(id) {
    build_elements(l_gen, id);
    init_state();
    for (size_t i = 0; i < _K; i++)
      extend();
  }
  FiniteGroupResolution(Tgroup const &G, size_t _K, bool _extendible = false)
      : FiniteGroupResolution(G.GeneratorsOfGroup(), G.get_identity(), _K, _extendible) {}

  // Computes one more term of the resolution.
  void extend() {
    next_res_term(K + 1);
    K++;
  }

  // The number of computed terms: R_0, ..., R_K are available.
  size_t length() const { return K; }
  size_t group_order() const { return N; }
  // The sum of the lengths of the boundaries of the generators of R_i for
  // i = 1..K, as Size(R) in HAP.
  std::vector<size_t> sizes() const {
    std::vector<size_t> ret;
    for (size_t i = 1; i <= K; i++) {
      size_t s = 0;
      for (auto const &bnd : pseudo_boundary[i])
        s += LengthResolutionChain(bnd);
      ret.push_back(s);
    }
    return ret;
  }
  std::vector<Telt> const &elements() const { return elts; }
  size_t position(Telt const &g) const { return map_position.at(g); }
  // Index of elements()[g] * elements()[h]
  size_t product(size_t g, size_t h) const { return product_index(g, h); }

  // The rank of the free ZG-module R_i.
  size_t dimension(int i) const {
    if (i < 0)
      return 0;
    if (i == 0)
      return 1;
    if (size_t(i) > K) {
      std::cerr << "RES: dimension(" << i << ") requested but only " << K
                << " terms were computed\n";
      throw PermutalibException{1};
    }
    return pseudo_boundary[i].size();
  }

  // The boundary of the j-th free generator of R_i (1 <= i <= K), a chain
  // of R_{i-1}. For i = 1 this is the chain g_j - 1 with g_j the group
  // element chosen as generator.
  ResolutionChain const &boundary(size_t i, size_t j) const {
    if (i < 1 || i > K) {
      std::cerr << "RES: boundary(" << i << ", " << j << ") out of range\n";
      throw PermutalibException{1};
    }
    return pseudo_boundary[i][j];
  }

  // The boundary of a chain of R_i, extended ZG-linearly.
  ResolutionChain boundary(size_t i, ResolutionChain const &chain) const {
    ResolutionChain ret;
    if (i < 1)
      return ret;
    for (auto const &term : chain) {
      for (auto const &bnd : pseudo_boundary[i][term.cell]) {
        ret.push_back({term.coeff * bnd.coeff, bnd.cell, product_index(term.elt, bnd.elt)});
      }
    }
    ReduceResolutionChain(ret);
    return ret;
  }

  // The augmentation R_0 -> Z.
  int64_t augmentation(ResolutionChain const &chain) const {
    int64_t sum = 0;
    for (auto const &term : chain)
      sum += term.coeff;
    return sum;
  }

  // The contracting homotopy h_i : R_i -> R_{i+1} applied to a term of R_i,
  // for 0 <= i < K. This is HAP's homotopy(i, p) = -Contraction(i+1, p).
  ResolutionChain homotopy(size_t i, ResolutionTerm const &term) {
    if (i + 1 > K) {
      std::cerr << "RES: homotopy(" << i << ", .) needs " << i + 1
                << " terms but only " << K << " were computed\n";
      throw PermutalibException{1};
    }
    int s = term.coeff > 0 ? 1 : -1;
    int64_t mult = term.coeff > 0 ? term.coeff : -term.coeff;
    ResolutionChain c = contraction(i + 1, -s, term.cell, term.elt);
    if (mult != 1) {
      for (auto &e : c)
        e.coeff *= mult;
    }
    return c;
  }

  ResolutionChain homotopy(size_t i, ResolutionChain const &chain) {
    ResolutionChain ret;
    for (auto const &term : chain) {
      ResolutionChain c = homotopy(i, term);
      ret.insert(ret.end(), c.begin(), c.end());
    }
    ReduceResolutionChain(ret);
    return ret;
  }

  // The matrix of the boundary d_i tensored with Z over ZG, that is the
  // integer matrix with dimension(i-1) rows and dimension(i) columns
  // obtained by mapping every group element to 1.
  std::vector<std::vector<Tint>> boundary_matrix_over_Z(size_t i) const {
    size_t n_row = dimension(int(i) - 1);
    size_t n_col = dimension(int(i));
    std::vector<std::vector<Tint>> M(n_row, std::vector<Tint>(n_col, 0));
    if (i < 1)
      return M;
    for (size_t j = 0; j < n_col; j++) {
      for (auto const &term : pseudo_boundary[i][j]) {
        M[term.cell][j] += Tint(long(term.coeff));
      }
    }
    return M;
  }

  // The abelian invariants of the integral homology H_n(G, Z) for
  // 0 <= n < K: the torsion coefficients d_1 | d_2 | ... (each > 1)
  // followed by one 0 for each infinite cyclic factor.
  std::vector<Tint> integral_homology(size_t n) const {
    std::vector<Tint> hom = integral_homology_kernel(n);
#ifdef DEBUG_RESOLUTION
    check_homology_consistency(n);
#endif
    return hom;
  }

  // The abelian invariants of H_n(G, Z) in the form of AbelianInvariants of
  // GAP, i.e. as prime powers.
  std::vector<Tint> integral_homology_prime_power(size_t n) const {
    return PrimePowerAbelianInvariants(integral_homology(n));
  }

  // The dimension of H_n(G, F_p) for a prime p, obtained from the ranks
  // modulo p of the boundary matrices.
  size_t homology_dimension_mod_p(size_t n, Tint const &p) const {
    if (n + 1 > K) {
      std::cerr << "RES: homology_dimension_mod_p(" << n << ") needs " << n + 1
                << " terms but only " << K << " were computed\n";
      throw PermutalibException{1};
    }
    size_t rank_n = RankModP(boundary_matrix_over_Z(n), p);
    size_t rank_np1 = RankModP(boundary_matrix_over_Z(n + 1), p);
    return dimension(int(n)) - rank_n - rank_np1;
  }

  // Consistency checks of the homology in degree n. See the description at
  // the top of the file. Throws a PermutalibException on failure.
  void check_homology_consistency(size_t n) const {
    std::vector<Tint> hom = integral_homology_kernel(n);
    Tint order = UnsignedToTint<Tint>(N);
    auto StringTint = [](Tint const &x) -> std::string {
      std::ostringstream os;
      os << x;
      return os.str();
    };
    auto fail = [&](std::string const &msg) {
      std::cerr << "RES: homology consistency failure in degree " << n << ": " << msg << "\n";
      throw PermutalibException{1};
    };
    if (n == 0) {
      if (hom.size() != 1 || hom[0] != 0)
        fail("H_0 is not Z");
      return;
    }
    for (auto const &d : hom) {
      if (d == 0)
        fail("H_n has an infinite factor");
      if (d < 2)
        fail("an invariant is smaller than 2");
      if (order % d != 0)
        fail("an invariant does not divide the order of the group");
    }
    // Universal coefficient theorem for the primes dividing the order.
    std::vector<Tint> hom_prev = integral_homology_kernel(n - 1);
    // Number of invariants d with d tensor F_p nonzero (tensor = true, then
    // a zero invariant counts) or with Tor(Z/d, F_p) nonzero (only torsion).
    auto count_div = [](std::vector<Tint> const &inv, Tint const &p, bool tensor) -> size_t {
      size_t t = 0;
      for (auto const &d : inv) {
        if (d == 0) {
          if (tensor)
            t++;
        } else {
          if (d % p == 0)
            t++;
        }
      }
      return t;
    };
    std::vector<Tint> primes = PrimeFactors(order);
    for (auto const &p : primes) {
      size_t dim_p = homology_dimension_mod_p(n, p);
      size_t expected = count_div(hom, p, true) + count_div(hom_prev, p, false);
      if (dim_p != expected)
        fail("the universal coefficient theorem fails for p=" + StringTint(p) + ": dim H_n(G,F_p)=" +
             std::to_string(dim_p) + " but t_p(H_n)+t_p(H_{n-1})=" + std::to_string(expected));
    }
    // A prime not dividing the order gives no homology.
    Tint q = 2;
    while (order % q == 0 || PrimeFactors(q).size() != 1 || PrimeFactors(q)[0] != q)
      q += 1;
    if (homology_dimension_mod_p(n, q) != 0)
      fail("nonzero homology modulo a prime not dividing the order");
    // H_1 is the abelianization
    if (n == 1) {
      Tgroup G(l_gen_store, id_store);
      Tint order_derived = G.DerivedSubgroup().size();
      Tint order_h1 = 1;
      for (auto const &d : hom)
        order_h1 *= d;
      if (order_h1 * order_derived != order)
        fail("the order of H_1 differs from the index of the derived subgroup");
    }
  }

private:
  std::vector<Tint> integral_homology_kernel(size_t n) const {
    if (n + 1 > K) {
      std::cerr << "RES: integral_homology(" << n << ") needs " << n + 1
                << " terms but only " << K << " were computed\n";
      throw PermutalibException{1};
    }
    std::vector<Tint> snf_n = SmithNormalFormInvariants(boundary_matrix_over_Z(n));
    std::vector<Tint> snf_np1 = SmithNormalFormInvariants(boundary_matrix_over_Z(n + 1));
    size_t rank_n = snf_n.size();
    size_t rank_np1 = snf_np1.size();
    std::vector<Tint> ret;
    for (auto const &d : snf_np1) {
      if (d > 1)
        ret.push_back(d);
    }
    size_t free_rank = dimension(int(n)) - rank_n - rank_np1;
    for (size_t u = 0; u < free_rank; u++)
      ret.push_back(0);
    return ret;
  }

public:
  // A GAP record with the same content as the HapResolution of HAP:
  // elts, the list of dimensions and the boundaries as HAP words.
  std::string GapString() const {
    std::string str = "rec(elts:=" + GapStringTVector(elts) + ",\n";
    str += "  dimension:=[ ";
    for (size_t i = 0; i <= K; i++) {
      if (i > 0)
        str += ", ";
      str += std::to_string(dimension(int(i)));
    }
    str += " ],\n  boundary:=[ ";
    for (size_t i = 1; i <= K; i++) {
      if (i > 1)
        str += ",\n    ";
      str += "[ ";
      size_t dim = pseudo_boundary[i].size();
      for (size_t j = 0; j < dim; j++) {
        if (j > 0)
          str += ", ";
        str += GapStringResolutionChain(pseudo_boundary[i][j]);
      }
      str += " ]";
    }
    str += " ])";
    return str;
  }
};

// clang-format off
}  // namespace permutalib
#endif  // SRC_GAP_RESOLUTION_H_
// clang-format on
