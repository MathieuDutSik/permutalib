// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Permutation.h"
#include "gmpxx.h"
#include <chrono>
#include <fstream>
#include <sstream>

#include "Group.h"
#include "Resolution.h"

/*
  Checks of the resolutions computed in Resolution.h:
  --- d o d = 0 and the augmentation vanishes on the image of d_1.
  --- The contracting homotopy satisfies d h + h d = id.
  --- The partial contracting homotopy y of a cycle z satisfies d y = z.
  --- The consistency checks of check_homology_consistency (finiteness,
      divisibility by the group order, universal coefficient theorem with
      the homology modulo p, abelianization).
  --- The integral homology matches the expected values of the input file,
      which come from HAP or from the literature.
  --- The homology recomputed from a conjugate of the group (so with a
      different enumeration of the elements) and in the extendible mode is
      the same, since it is an invariant of the group.
  --- Optionally, the ranks, the sizes of the boundaries and the dimensions
      of the homology modulo p match the values recorded by HAP.

  Input file format:
    nGroup
    then for each group:
      n nbGen                     (as in ReadGroupFromStream)
      the nbGen permutations
      K                           (length of the resolution to compute)
      then keyword lines until "end":
        homology n len d_1 ... d_len   expected invariants of H_n(G,Z) in
                                       Smith normal form (d_1 | d_2 | ...)
        ranks r_0 ... r_K              expected ranks (from HAP)
        sizes s_1 ... s_K              expected Size(R) (from HAP)
        modp p n dim                   expected dimension of H_n(G, F_p)
        norecompute                    skip the recomputation from the
                                       conjugate group (large groups)
 */

template <typename Tres>
void CheckBoundarySquare(Tres const &R) {
  size_t K = R.length();
  for (size_t i = 1; i <= K; i++) {
    size_t dim = R.dimension(int(i));
    for (size_t j = 0; j < dim; j++) {
      permutalib::ResolutionChain const &bnd = R.boundary(i, j);
      if (i == 1) {
        if (R.augmentation(bnd) != 0) {
          std::cerr << "The augmentation of d_1(e_" << j << ") is nonzero\n";
          throw permutalib::PermutalibException{1};
        }
      } else {
        permutalib::ResolutionChain bnd2 = R.boundary(i - 1, bnd);
        if (bnd2.size() != 0) {
          std::cerr << "d_" << i - 1 << " d_" << i << "(e_" << j << ") is nonzero\n";
          throw permutalib::PermutalibException{1};
        }
      }
    }
  }
}

template <typename Tres>
void CheckHomotopy(Tres &R) {
  size_t K = R.length();
  size_t N = R.group_order();
  for (size_t i = 0; i + 1 < K; i++) {
    size_t dim = R.dimension(int(i));
    for (size_t j = 0; j < dim; j++) {
      for (size_t g = 0; g < N; g++) {
        permutalib::ResolutionTerm x{1, j, g};
        permutalib::ResolutionChain sum = R.boundary(i + 1, R.homotopy(i, x));
        if (i == 0) {
          // d h(x) = x - eps(x) e^0_0
          sum.push_back({-1, 0, g});
          sum.push_back({1, 0, 0});
        } else {
          permutalib::ResolutionChain hd = R.homotopy(i - 1, R.boundary(i, {x}));
          sum.insert(sum.end(), hd.begin(), hd.end());
          sum.push_back({-1, j, g});
        }
        permutalib::ReduceResolutionChain(sum);
        if (sum.size() != 0) {
          std::cerr << "The homotopy identity fails in dimension " << i
                    << " for the cell (" << j << ", " << g << ")\n";
          throw permutalib::PermutalibException{1};
        }
      }
    }
  }
}

// The partial contracting homotopy y of a cycle z satisfies d y = z. It is
// compared with the contracting homotopy h(z), which also satisfies
// d h(z) = z, by the total length of the solutions and the running time.
// The cycles are z = d(c) for small chains c of R_{i+1} (z = g - 1 for
// i = 0) and z = x - h(d x) for x a cell of R_i, at most n_sample of each
// kind chosen by a fixed pseudo-random sequence (fewer for large groups).
template <typename Tres>
void CheckPartialHomotopy(Tres &R) {
  size_t K = R.length();
  size_t N = R.group_order();
  size_t n_sample = N > 1000 ? 20 : 200;
  uint64_t seed = 1;
  auto next_rand = [&](size_t n) -> size_t {
    seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
    return size_t((seed >> 33) % n);
  };
  for (size_t i = 0; i < K; i++) {
    std::vector<permutalib::ResolutionChain> l_cycle;
    size_t n_cell = R.dimension(int(i)) * N;
    for (size_t u = 0; u < std::min(n_sample, n_cell); u++) {
      permutalib::ResolutionChain z;
      if (i == 0) {
        z = {{-1, 0, 0}, {1, 0, next_rand(N)}};
      } else {
        permutalib::ResolutionChain c;
        for (size_t v = 0; v < 3; v++)
          c.push_back({int64_t(next_rand(3)) - 1, next_rand(R.dimension(int(i) + 1)), next_rand(N)});
        permutalib::ReduceResolutionChain(c);
        z = R.boundary(i + 1, c);
      }
      l_cycle.push_back(z);
      if (i > 0) {
        permutalib::ResolutionTerm x{1, next_rand(R.dimension(int(i))), next_rand(N)};
        z = R.homotopy(i - 1, R.boundary(i, {x}));
        for (auto &term : z)
          term.coeff = -term.coeff;
        z.push_back(x);
        permutalib::ReduceResolutionChain(z);
        l_cycle.push_back(z);
      }
    }
    size_t len_h = 0;
    size_t len_p = 0;
    int64_t time_h = 0;
    int64_t time_p = 0;
    size_t n_shorter = 0;
    size_t n_longer = 0;
    for (auto const &z : l_cycle) {
      std::chrono::time_point<std::chrono::system_clock> time1 =
          std::chrono::system_clock::now();
      permutalib::ResolutionChain y_p = R.partial_contracting_homotopy(i, z);
      std::chrono::time_point<std::chrono::system_clock> time2 =
          std::chrono::system_clock::now();
      permutalib::ResolutionChain y_h = R.homotopy(i, z);
      std::chrono::time_point<std::chrono::system_clock> time3 =
          std::chrono::system_clock::now();
      time_p += std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count();
      time_h += std::chrono::duration_cast<std::chrono::microseconds>(time3 - time2).count();
      permutalib::ResolutionChain diff = R.boundary(i + 1, y_p);
      for (auto const &term : z)
        diff.push_back({-term.coeff, term.cell, term.elt});
      permutalib::ReduceResolutionChain(diff);
      if (diff.size() != 0) {
        std::cerr << "The partial contracting homotopy fails in dimension " << i << "\n";
        throw permutalib::PermutalibException{1};
      }
      size_t l_p = permutalib::LengthResolutionChain(y_p);
      len_p += l_p;
      size_t l_h = permutalib::LengthResolutionChain(y_h);
      len_h += l_h;
      if (l_p < l_h)
        n_shorter++;
      if (l_p > l_h)
        n_longer++;
    }
    std::cerr << "  partial homotopy i=" << i << " |cycles|=" << l_cycle.size()
              << " length partial=" << len_p << " (" << time_p << "us)"
              << " homotopy=" << len_h << " (" << time_h << "us) shorter=" << n_shorter
              << " longer=" << n_longer << "\n";
  }
}

template <typename T>
std::string StringVector(std::vector<T> const &v) {
  std::ostringstream os;
  os << "[";
  for (size_t u = 0; u < v.size(); u++) {
    if (u > 0)
      os << ",";
    os << v[u];
  }
  os << "]";
  return os.str();
}

template <typename Tres>
std::vector<size_t> GetRanks(Tres const &R) {
  std::vector<size_t> ranks;
  for (size_t i = 0; i <= R.length(); i++)
    ranks.push_back(R.dimension(int(i)));
  return ranks;
}

template <typename Telt>
std::vector<Telt> ConjugateGenerators(std::vector<Telt> const &l_gen, Telt const &id) {
  // Conjugate by the permutation i -> n-1-i, which changes the sorted
  // enumeration of the elements and so the whole run of the algorithm.
  using Tidx = typename Telt::Tidx;
  Tidx n = id.size();
  std::vector<Tidx> eList(n);
  for (Tidx i = 0; i < n; i++)
    eList[i] = Tidx(n - 1 - i);
  Telt eConj(eList);
  Telt eConjInv = permutalib::Inverse(eConj);
  std::vector<Telt> ret;
  for (auto const &g : l_gen)
    ret.push_back(eConjInv * g * eConj);
  return ret;
}

int main(int argc, char *argv[]) {
  try {
    using Tidx = uint8_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tres = permutalib::FiniteGroupResolution<Telt, Tint>;
    if (argc != 3) {
      std::cerr << "TestResolution [InputFile] [case]\n";
      std::cerr << "\n";
      std::cerr << "InputFile: The file containing the groups and expected values\n";
      std::cerr << "case: \"check\" for running all the checks\n";
      std::cerr << "      \"print\" for printing the ranks, sizes and homology\n";
      throw permutalib::PermutalibException{1};
    }
    std::string InputFile = argv[1];
    std::string ecase = argv[2];
    if (ecase != "check" && ecase != "print") {
      std::cerr << "Available options are check and print\n";
      throw permutalib::PermutalibException{1};
    }
    std::ifstream is(InputFile);
    int nGroup;
    is >> nGroup;
    for (int iGroup = 0; iGroup < nGroup; iGroup++) {
      std::cerr << "iGroup=" << iGroup << "/" << nGroup << "\n";
      // The original generators are used (and not those of the stabilizer
      // chain) since the resolution of HAP depends on them.
      std::pair<std::vector<Telt>, Telt> pair = permutalib::ReadListGenFromStream<Telt>(is);
      std::vector<Telt> const &l_gen = pair.first;
      Telt const &id = pair.second;
      size_t K;
      is >> K;
      std::vector<std::vector<Tint>> expected_hom(K);
      std::vector<bool> has_expected_hom(K, false);
      std::vector<size_t> expected_ranks;
      std::vector<size_t> expected_sizes;
      std::vector<std::tuple<long, size_t, size_t>> expected_modp;
      bool recompute = true;
      while (true) {
        std::string keyword;
        is >> keyword;
        if (keyword == "end")
          break;
        if (keyword == "homology") {
          size_t n, len;
          is >> n >> len;
          if (n >= K) {
            std::cerr << "homology degree " << n << " needs K > " << n << "\n";
            throw permutalib::PermutalibException{1};
          }
          has_expected_hom[n] = true;
          for (size_t u = 0; u < len; u++) {
            long val;
            is >> val;
            expected_hom[n].push_back(Tint(val));
          }
        } else if (keyword == "ranks") {
          expected_ranks.resize(K + 1);
          for (size_t i = 0; i <= K; i++)
            is >> expected_ranks[i];
        } else if (keyword == "sizes") {
          expected_sizes.resize(K);
          for (size_t i = 0; i < K; i++)
            is >> expected_sizes[i];
        } else if (keyword == "modp") {
          long p;
          size_t n, dim;
          is >> p >> n >> dim;
          expected_modp.push_back({p, n, dim});
        } else if (keyword == "norecompute") {
          recompute = false;
        } else {
          std::cerr << "Unknown keyword " << keyword << " in the input file\n";
          throw permutalib::PermutalibException{1};
        }
      }
      std::chrono::time_point<std::chrono::system_clock> time1 =
          std::chrono::system_clock::now();
      Tres R(l_gen, id, K);
      std::chrono::time_point<std::chrono::system_clock> time2 =
          std::chrono::system_clock::now();
      std::vector<size_t> ranks = GetRanks(R);
      std::vector<size_t> sizes = R.sizes();
      std::cerr << "|G|=" << R.group_order() << " K=" << K
                << " ranks=" << StringVector(ranks) << " sizes=" << StringVector(sizes)
                << " milliseconds="
                << std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count()
                << "\n";
      std::vector<std::vector<Tint>> hom(K);
      for (size_t n = 0; n < K; n++)
        hom[n] = R.integral_homology(n);
      for (size_t n = 1; n < K; n++) {
        std::cerr << "  H_" << n << " = " << StringVector(hom[n])
                  << " GAP form " << StringVector(R.integral_homology_prime_power(n)) << "\n";
      }
      if (ecase == "check") {
        for (size_t n = 0; n < K; n++)
          R.check_homology_consistency(n);
        for (size_t n = 1; n < K; n++) {
          if (has_expected_hom[n] && hom[n] != expected_hom[n]) {
            std::cerr << "H_" << n << " = " << StringVector(hom[n])
                      << " but expected " << StringVector(expected_hom[n]) << "\n";
            throw permutalib::PermutalibException{1};
          }
        }
        if (expected_ranks.size() > 0 && ranks != expected_ranks) {
          std::cerr << "ranks=" << StringVector(ranks) << " but expected "
                    << StringVector(expected_ranks) << "\n";
          throw permutalib::PermutalibException{1};
        }
        if (expected_sizes.size() > 0 && sizes != expected_sizes) {
          std::cerr << "sizes=" << StringVector(sizes) << " but expected "
                    << StringVector(expected_sizes) << "\n";
          throw permutalib::PermutalibException{1};
        }
        for (auto const &e : expected_modp) {
          long p = std::get<0>(e);
          size_t n = std::get<1>(e);
          size_t dim = R.homology_dimension_mod_p(n, Tint(p));
          if (dim != std::get<2>(e)) {
            std::cerr << "dim H_" << n << "(G, F_" << p << ") = " << dim
                      << " but expected " << std::get<2>(e) << "\n";
            throw permutalib::PermutalibException{1};
          }
        }
        CheckBoundarySquare(R);
        CheckHomotopy(R);
        CheckPartialHomotopy(R);
        std::cerr << "  d o d = 0, homotopy, consistency and expected values are correct\n";
        if (recompute) {
          Tres R2(ConjugateGenerators(l_gen, id), id, K, true);
          std::cerr << "  conjugate extendible run: ranks=" << StringVector(GetRanks(R2)) << "\n";
          for (size_t n = 0; n < K; n++) {
            std::vector<Tint> hom2 = R2.integral_homology(n);
            if (hom2 != hom[n]) {
              std::cerr << "H_" << n << " = " << StringVector(hom2)
                        << " for the conjugate group but " << StringVector(hom[n])
                        << " for the original one\n";
              throw permutalib::PermutalibException{1};
            }
          }
          std::cerr << "  the homology of the conjugate group is the same\n";
        }
      }
    }
    std::cerr << "Normal completion of the program\n";
  } catch (permutalib::PermutalibException const &e) {
    std::cerr << "Erroneous completion of the program\n";
    exit(e.eVal);
  }
  return 0;
}
