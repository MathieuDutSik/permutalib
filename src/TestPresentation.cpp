// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Checks the presentations of Presentation.h on a list of groups, for the
// deterministic and the randomized (with the exact order) stabilizer
// chains:
// --- The relators of the presentation on the strong generators hold, and
//     the coset enumeration of each level, the index of <A_{i+1}> in
//     <A_i | T_i, ..., T_k>, gives the basic orbit length. This proves that
//     the presented group has the order of the group.
// --- The simplified presentation on the top level generators: its
//     relators hold, it starts with the top level generators and, when the
//     order is at most max_order_enum, the coset enumeration of the trivial
//     subgroup gives the order.
// --- For the orders up to max_order_res, the presentations from the
//     resolution (Resolution.h), on the generators chosen by the
//     resolution, on both these and the input generators, and on the input
//     generators: the relators hold and, when the order is at most
//     max_order_enum, the coset enumeration of the trivial subgroup gives
//     the order (it may overflow for the input generators).
// Before that, the coset enumeration is checked on presentations of known
// groups. The list file has one line "name degree order base_length" per
// group, the group being in GroupDir/name.grp. Prints one line
// "RESULT name ok|FAIL ..." per group and chain; exit code 1 on a failure.
#include "Group.h"
#include "Permutation.h"
#include "Presentation.h"
#include "Resolution.h"
#include "StabChainMain.h"
#include "gmpxx.h"
#include <chrono>
#include <fstream>

using Tint = mpz_class;
using permutalib::PresentationWord;

size_t const max_cosets = 5000000;

// Repeats the word w, e times.
PresentationWord power_word(PresentationWord const &w, int e) {
  PresentationWord ret;
  for (int i = 0; i < e; i++)
    ret.insert(ret.end(), w.begin(), w.end());
  return ret;
}

int check_known_presentations() {
  struct Case {
    std::string name;
    size_t n_gen;
    std::vector<PresentationWord> rels;
    std::vector<PresentationWord> subgens;
    size_t index;
  };
  std::vector<Case> cases;
  // A_5 = <a, b | a^2, b^3, (ab)^5>
  cases.push_back({"A5", 2, {{1, 1}, {2, 2, 2}, power_word({1, 2}, 5)}, {}, 60});
  cases.push_back({"A5/<b>", 2, {{1, 1}, {2, 2, 2}, power_word({1, 2}, 5)},
                   {{2}}, 20});
  // The dihedral group of order 2n: <a, b | a^n, b^2, (ab)^2>
  for (int n : {1, 2, 7, 30})
    cases.push_back({"D" + std::to_string(2 * n), 2,
                     {power_word({1}, n), {2, 2}, {1, 2, 1, 2}}, {},
                     static_cast<size_t>(2 * n)});
  // S_5 as a Coxeter group
  {
    std::vector<PresentationWord> rels;
    for (int i = 1; i <= 4; i++) {
      rels.push_back({i, i});
      for (int j = i + 1; j <= 4; j++)
        rels.push_back(power_word({i, j}, j == i + 1 ? 3 : 2));
    }
    cases.push_back({"S5", 4, rels, {}, 120});
    cases.push_back({"S5/S4", 4, rels, {{1}, {2}, {3}}, 5});
  }
  // A collapse: <a, b | a^-1 b a b^-2, b^-1 a b a^-2> is trivial.
  cases.push_back({"trivial", 2, {{-1, 2, 1, -2, -2}, {-2, 1, 2, -1, -1}}, {}, 1});
  // The quaternion group <a, b | a^4, a^2 b^-2, b^-1 a b a>
  cases.push_back({"Q8", 2, {{1, 1, 1, 1}, {1, 1, -2, -2}, {-2, 1, 2, 1}}, {}, 8});
  int n_fail = 0;
  for (auto &c : cases) {
    std::optional<size_t> index =
        permutalib::CosetEnumeration(c.n_gen, c.rels, c.subgens, max_cosets);
    bool ok = index && *index == c.index;
    std::cout << "KNOWN " << c.name << " " << (ok ? "ok" : "FAIL") << "\n";
    if (!ok)
      n_fail++;
  }
  return n_fail;
}

template <typename Telt>
bool relators_hold(permutalib::Presentation<Telt> const &P, Telt const &id) {
  for (auto &r : P.relators)
    if (!permutalib::EvaluateWord(P.generators, id, r).is_identity())
      return false;
  return true;
}

template <typename Telt>
int run(std::string const &name, std::string const &eFile,
        Tint const &order_gap, Tint const &max_order_enum) {
  using Tidx_label = uint16_t;
  std::ifstream is(eFile);
  std::pair<std::vector<Telt>, Telt> pair =
      permutalib::ReadListGenFromStream<Telt>(is);
  std::vector<Telt> const &LGen = pair.first;
  Telt const &id = pair.second;
  int n_fail = 0;
  for (bool randomized : {false, true}) {
    std::vector<std::string> fails;
    auto fail = [&](std::string const &s) { fails.push_back(s); };
    permutalib::StabChainOptions<Tint, Telt> options =
        permutalib::GetStandardOptions<Tint, Telt>(id);
    options.randomized = randomized;
    if (randomized)
      options.size = order_gap;
    permutalib::StabChain<Telt, Tidx_label> S =
        permutalib::StabChainOp_listgen<Telt, Tidx_label, Tint>(LGen, options);
    permutalib::UnbindCycles(S);
    if (permutalib::Order<Telt, Tidx_label, Tint>(S) != order_gap)
      fail("wrong order of the chain");
    auto t0 = std::chrono::steady_clock::now();
    permutalib::StabChainPresentation<Telt> SP =
        permutalib::Kernel_StabChainPresentation(S);
    auto t1 = std::chrono::steady_clock::now();
    if (!permutalib::VerifyStabChainPresentation(SP, id, max_cosets))
      fail("the levels of the presentation are wrong");
    auto t2 = std::chrono::steady_clock::now();
    permutalib::Presentation<Telt> P_strong =
        permutalib::FlattenStabChainPresentation(SP);
    if (!relators_hold(P_strong, id))
      fail("a relator on the strong generators does not hold");
    auto t3 = std::chrono::steady_clock::now();
    permutalib::Presentation<Telt> P =
        permutalib::Kernel_Presentation(S, 10000000);
    auto t4 = std::chrono::steady_clock::now();
    if (!relators_hold(P, id))
      fail("a relator of the simplified presentation does not hold");
    std::vector<Telt> top;
    for (auto &g : permutalib::GetListGenerators(S))
      if (!g.is_identity())
        top.push_back(g);
    if (P.generators.size() < top.size() ||
        !std::equal(top.begin(), top.end(), P.generators.begin()))
      fail("the simplified presentation does not start with the top "
           "level generators");
    std::string enum_str = "skipped";
    if (order_gap <= max_order_enum) {
      std::optional<size_t> index = permutalib::CosetEnumeration(
          P.generators.size(), P.relators, {}, max_cosets);
      if (!index) {
        enum_str = "overflow";
        fail("the coset enumeration of the simplified presentation overflows");
      } else {
        enum_str = std::to_string(*index);
        if (Tint(static_cast<unsigned long>(*index)) != order_gap)
          fail("the simplified presentation has order " + enum_str);
      }
    }
    auto secs = [](auto a, auto b) {
      return std::chrono::duration<double>(b - a).count();
    };
    std::cout << "RESULT " << name << " " << (randomized ? "rnd" : "det")
              << " " << (fails.empty() ? "ok" : "FAIL")
              << " strong: gens=" << P_strong.generators.size()
              << " rels=" << P_strong.relators.size()
              << " len=" << permutalib::TotalLength(P_strong.relators)
              << " simplified: gens=" << P.generators.size()
              << " top=" << top.size() << " rels=" << P.relators.size()
              << " len=" << permutalib::TotalLength(P.relators)
              << " enum=" << enum_str << " t_pres=" << secs(t0, t1)
              << " t_verify=" << secs(t1, t2) << " t_simplify=" << secs(t3, t4);
    for (auto &s : fails)
      std::cout << " | " << s;
    std::cout << "\n";
    if (!fails.empty())
      n_fail++;
  }
  return n_fail;
}

// The order of the presented group by coset enumeration, as a string.
template <typename Telt>
std::string enumerated_order(permutalib::Presentation<Telt> const &P) {
  std::optional<size_t> index = permutalib::CosetEnumeration(
      P.generators.size(), P.relators, {}, max_cosets);
  if (!index)
    return "overflow";
  return std::to_string(*index);
}

template <typename Telt>
int run_resolution(std::string const &name, std::string const &eFile,
                   Tint const &order_gap, Tint const &max_order_enum) {
  std::ifstream is(eFile);
  std::pair<std::vector<Telt>, Telt> pair =
      permutalib::ReadListGenFromStream<Telt>(is);
  std::vector<Telt> const &LGen = pair.first;
  Telt const &id = pair.second;
  std::vector<std::string> fails;
  auto fail = [&](std::string const &s) { fails.push_back(s); };
  auto t0 = std::chrono::steady_clock::now();
  permutalib::FiniteGroupResolution<Telt, Tint> R(LGen, id, 2);
  permutalib::Presentation<Telt> P_y = R.presentation();
  auto t1 = std::chrono::steady_clock::now();
  permutalib::Presentation<Telt> P_xy =
      permutalib::PresentationByResolutionMixed<Telt, Tint>(LGen, id);
  permutalib::Presentation<Telt> P_x =
      permutalib::PresentationByResolution<Telt, Tint>(LGen, id);
  auto t2 = std::chrono::steady_clock::now();
  if (!relators_hold(P_y, id))
    fail("a relator on the generators of the resolution does not hold");
  if (!relators_hold(P_xy, id))
    fail("a relator on both generators does not hold");
  if (!relators_hold(P_x, id))
    fail("a relator on the input generators does not hold");
  std::vector<Telt> gens = permutalib::DistinctGenerators(LGen, id);
  if (P_x.generators != gens)
    fail("the generators are not the input generators");
  // The coset enumeration of the presentation on the input generators may
  // overflow since its relators are long; the one on both generators,
  // from which it is obtained by the Tietze transformations of
  // SimplifyPresentation, is checked instead.
  std::string enum_y = "skipped", enum_xy = "skipped", enum_x = "skipped";
  if (order_gap <= max_order_enum) {
    std::string ord = order_gap.get_str();
    enum_y = enumerated_order(P_y);
    enum_xy = enumerated_order(P_xy);
    enum_x = enumerated_order(P_x);
    if (enum_y != ord)
      fail("the presentation on the generators of the resolution has order " +
           enum_y);
    if (enum_xy != ord)
      fail("the presentation on both generators has order " + enum_xy);
    if (enum_x != ord && enum_x != "overflow")
      fail("the presentation on the input generators has order " + enum_x);
  }
  auto secs = [](auto a, auto b) {
    return std::chrono::duration<double>(b - a).count();
  };
  std::cout << "RESULT " << name << " res " << (fails.empty() ? "ok" : "FAIL")
            << " resolution: gens=" << P_y.generators.size()
            << " rels=" << P_y.relators.size()
            << " len=" << permutalib::TotalLength(P_y.relators)
            << " enum=" << enum_y << " both: enum=" << enum_xy
            << " input: gens=" << P_x.generators.size()
            << " rels=" << P_x.relators.size()
            << " len=" << permutalib::TotalLength(P_x.relators)
            << " enum=" << enum_x << " t_res=" << secs(t0, t1)
            << " t_input=" << secs(t1, t2);
  for (auto &s : fails)
    std::cout << " | " << s;
  std::cout << "\n";
  return fails.empty() ? 0 : 1;
}

int main(int argc, char *argv[]) {
  try {
    if (argc != 5) {
      std::cerr << "TestPresentation ListFile GroupDir max_order_enum "
                   "max_order_res\n";
      return 1;
    }
    std::ifstream is_list(argv[1]);
    std::string dir = argv[2];
    Tint max_order_enum(argv[3]);
    Tint max_order_res(argv[4]);
    int n_fail = check_known_presentations();
    int n_case = 0;
    std::string name, degree, order_str, base;
    while (is_list >> name >> degree >> order_str >> base) {
      Tint order_gap(order_str);
      std::string eFile = dir + "/" + name + ".grp";
      n_fail += run<permutalib::SingleSidedPerm<uint16_t>>(name, eFile,
                                                           order_gap,
                                                           max_order_enum);
      n_case += 2;
      if (order_gap <= max_order_res) {
        n_fail += run_resolution<permutalib::SingleSidedPerm<uint16_t>>(
            name, eFile, order_gap, max_order_enum);
        n_case++;
      }
    }
    std::cout << "n_case=" << n_case << " n_fail=" << n_fail << "\n";
    return n_fail == 0 ? 0 : 1;
  } catch (permutalib::PermutalibException const &e) {
    std::cerr << "PermutalibException\n";
    return 2;
  }
}
