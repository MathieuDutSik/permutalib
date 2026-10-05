// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Compares the deterministic and the randomized stabilizer chain
// constructions (StabChainRandomized, with no order, the exact order and a
// wrong order, twice the exact one) on a list of groups,
// against the orders computed by GAP. The list file has one line
// "name degree order base_length" per group, the group being in
// GroupDir/name.grp. Prints one line "RESULT name ok|FAIL ..." per group and
// permutation type; exit code 1 on a failure.
#include "Group.h"
#include "Permutation.h"
#include "StabChainMain.h"
#include "gmpxx.h"
#include <chrono>
#include <fstream>
#include <random>

using Tint = mpz_class;

template <typename Telt, typename Tidx_label>
bool check_chain_invariants(permutalib::StabChain<Telt, Tidx_label> const &S,
                            std::string &msg) {
  using Tidx = typename Telt::Tidx;
  permutalib::StabChain<Telt, Tidx_label> Sptr = S;
  std::vector<Tidx> base_before;
  int level = 0;
  while (Sptr != nullptr && Sptr->orbit.size() > 0) {
    Tidx bpt = Sptr->orbit[0];
    size_t n = Sptr->transversal.size();
    std::vector<bool> in_orbit(n, false);
    for (auto &p : Sptr->orbit) {
      if (in_orbit[p]) {
        msg = "level " + std::to_string(level) + ": repeated orbit point";
        return false;
      }
      in_orbit[p] = true;
    }
    for (size_t p = 0; p < n; p++) {
      bool has_label =
          Sptr->transversal[p] != std::numeric_limits<Tidx_label>::max();
      if (has_label != in_orbit[p]) {
        msg = "level " + std::to_string(level) +
              ": transversal does not match the orbit";
        return false;
      }
    }
    for (auto &lbl : Sptr->genlabels) {
      Telt const &g = Sptr->comm->labels[lbl];
      for (auto &b : base_before) {
        if (permutalib::PowAct(b, g) != b) {
          msg = "level " + std::to_string(level) +
                ": a generator moves an earlier base point";
          return false;
        }
      }
      for (auto &p : Sptr->orbit) {
        if (!in_orbit[permutalib::PowAct(p, g)]) {
          msg = "level " + std::to_string(level) +
                ": the orbit is not closed under a generator";
          return false;
        }
      }
    }
    for (auto &p : Sptr->orbit) {
      Telt x = permutalib::InverseRepresentative(Sptr, p);
      if (permutalib::PowAct(p, x) != bpt) {
        msg = "level " + std::to_string(level) +
              ": the inverse representative does not map to the base point";
        return false;
      }
      for (auto &b : base_before) {
        if (permutalib::PowAct(b, x) != b) {
          msg = "level " + std::to_string(level) +
                ": a representative moves an earlier base point";
          return false;
        }
      }
    }
    base_before.push_back(bpt);
    Sptr = Sptr->stabilizer;
    level++;
  }
  return true;
}

template <typename Telt> int run(std::string const &eFile, Tint const &order_gap, int n_set_tests) {
  using Tidx = typename Telt::Tidx;
  using Tidx_label = uint16_t;
  using Tgroup = permutalib::Group<Telt, Tint>;
  std::ifstream is(eFile);
  std::pair<std::vector<Telt>, Telt> pair =
      permutalib::ReadListGenFromStream<Telt>(is);
  std::vector<Telt> const &LGen = pair.first;
  Telt const &id = pair.second;
  Tidx n = id.size();
  std::vector<std::string> fails;
  auto fail = [&](std::string const &s) { fails.push_back(s); };
  auto build = [&](bool randomized, Tint const &size, double &secs)
      -> permutalib::StabChain<Telt, Tidx_label> {
    permutalib::StabChainOptions<Tint, Telt> options =
        permutalib::GetStandardOptions<Tint, Telt>(id);
    options.randomized = randomized;
    options.size = size;
    auto t0 = std::chrono::steady_clock::now();
    permutalib::StabChain<Telt, Tidx_label> S =
        permutalib::StabChainOp_listgen<Telt, Tidx_label, Tint>(LGen, options);
    permutalib::UnbindCycles(S);
    secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
               .count();
    return S;
  };
  double t_det, t_rnd, t_known;
  permutalib::StabChain<Telt, Tidx_label> S_det = build(false, Tint(0), t_det);
  permutalib::StabChain<Telt, Tidx_label> S_rnd = build(true, Tint(0), t_rnd);
  permutalib::StabChain<Telt, Tidx_label> S_known = build(true, order_gap, t_known);
  // A wrong order, too large: never reached, the completion runs.
  double t_wrong;
  permutalib::StabChain<Telt, Tidx_label> S_wrong =
      build(true, Tint(2) * order_gap, t_wrong);
  std::vector<std::pair<std::string, permutalib::StabChain<Telt, Tidx_label> *>>
      chains = {{"det", &S_det}, {"rnd", &S_rnd}, {"known", &S_known},
                {"wrong", &S_wrong}};
  for (auto &[name, Sp] : chains) {
    Tint ord = permutalib::Order<Telt, Tidx_label, Tint>(*Sp);
    if (ord != order_gap)
      fail(name + ": order " + ord.get_str() + " instead of " + order_gap.get_str());
    std::string msg;
    if (!check_chain_invariants(*Sp, msg))
      fail(name + ": " + msg);
    for (auto &g : LGen)
      if (!permutalib::IsElementInStabChain(*Sp, g))
        fail(name + ": an input generator is not a member");
  }
  if (!permutalib::EqualityTest(S_det, S_rnd))
    fail("det and rnd differ");
  if (!permutalib::EqualityTest(S_det, S_known))
    fail("det and known differ");
  if (!permutalib::EqualityTest(S_det, S_wrong))
    fail("det and wrong differ");
  // The Group constructors with an order. GeneratorsOfGroup is the strong
  // generating set, which depends on the construction: only its size is
  // reported.
  size_t n_strong_det = 0, n_strong_rnd = 0;
  {
    Tgroup G1(LGen, id, order_gap);
    if (G1.size() != order_gap)
      fail("the Group constructor with the order gives a wrong size");
    Tgroup G0(LGen, id);
    if (G1 != G0)
      fail("the Group constructor with the order gives another group");
    n_strong_det = G0.GeneratorsOfGroup().size();
    n_strong_rnd = G1.GeneratorsOfGroup().size();
  }
  // The top level keeps the input generators, as in the deterministic chain.
  {
    std::vector<Telt> top_det = permutalib::GetListGenerators(S_det);
    std::vector<Telt> top_rnd = permutalib::GetListGenerators(S_rnd);
    if (top_det != top_rnd)
      fail("the top level generators differ");
  }
  // Membership of uniformly random elements of each chain in the other,
  // and agreement on random permutations of the whole domain.
  std::mt19937_64 eng(12345);
  for (int iter = 0; iter < 30; iter++) {
    Telt x = permutalib::UniformRandomElement(S_det);
    Telt y = permutalib::UniformRandomElement(S_rnd);
    if (!permutalib::IsElementInStabChain(S_rnd, x) ||
        !permutalib::IsElementInStabChain(S_det, y))
      fail("a random element is not in the other chain");
    std::vector<Tidx> v(n);
    for (Tidx i = 0; i < n; i++)
      v[i] = i;
    std::shuffle(v.begin(), v.end(), eng);
    Telt z(v);
    if (permutalib::IsElementInStabChain(S_det, z) !=
        permutalib::IsElementInStabChain(S_rnd, z))
      fail("the membership of a random permutation differs");
  }
  // Group level: canonical images (which must not depend on the chain),
  // set stabilizers, representative actions.
  Tgroup G_det(S_det);
  Tgroup G_rnd(S_rnd);
  if (G_det.size() != G_rnd.size())
    fail("the Group sizes differ");
  for (int iter = 0; iter < n_set_tests; iter++) {
    permutalib::Face f(n);
    for (Tidx i = 0; i < n; i++)
      if (eng() % 3 == 0)
        f[i] = 1;
    if (G_det.CanonicalImage(f) != G_rnd.CanonicalImage(f))
      fail("CanonicalImage differs");
    if (G_det.OptCanonicalImage(f) != G_rnd.OptCanonicalImage(f))
      fail("OptCanonicalImage differs");
    Tgroup Stab_det = G_det.Stabilizer_OnSets(f);
    Tgroup Stab_rnd = G_rnd.Stabilizer_OnSets(f);
    if (Stab_det.size() != Stab_rnd.size() || Stab_det != Stab_rnd)
      fail("Stabilizer_OnSets differs");
    Telt u = permutalib::UniformRandomElement(S_rnd);
    permutalib::Face f2(n);
    for (Tidx i = 0; i < n; i++)
      if (f[i])
        f2[permutalib::PowAct(i, u)] = 1;
    std::optional<Telt> opt = G_rnd.RepresentativeAction_OnSets(f, f2);
    if (!opt) {
      fail("RepresentativeAction_OnSets found nothing");
    } else {
      for (Tidx i = 0; i < n; i++)
        if (f[i] && !f2[permutalib::PowAct(i, *opt)])
          fail("RepresentativeAction_OnSets is wrong");
    }
  }
  std::cout << "RESULT " << eFile.substr(eFile.find_last_of('/') + 1) << " "
            << (fails.empty() ? "ok" : "FAIL") << " t_det=" << t_det
            << " t_rnd=" << t_rnd << " t_known=" << t_known
            << " t_wrong=" << t_wrong
            << " base_det=" << permutalib::BaseStabChain(S_det).size()
            << " base_rnd=" << permutalib::BaseStabChain(S_rnd).size()
            << " strong_det=" << n_strong_det << " strong_rnd=" << n_strong_rnd;
  for (auto &s : fails)
    std::cout << " | " << s;
  std::cout << "\n";
  return fails.empty() ? 0 : 1;
}

int main(int argc, char *argv[]) {
  try {
    if (argc != 4) {
      std::cerr << "TestStabChainRandomized ListFile GroupDir n_set_tests\n";
      return 1;
    }
    std::ifstream is_list(argv[1]);
    std::string dir = argv[2];
    int n_set_tests = atoi(argv[3]);
    std::string name, degree, order_str, base;
    int n_fail = 0, n_case = 0;
    while (is_list >> name >> degree >> order_str >> base) {
      Tint order_gap(order_str);
      std::string eFile = dir + "/" + name + ".grp";
      n_fail += run<permutalib::SingleSidedPerm<uint16_t>>(eFile, order_gap,
                                                          n_set_tests);
      n_fail += run<permutalib::DoubleSidedPerm<uint16_t>>(eFile, order_gap,
                                                          n_set_tests);
      n_case += 2;
    }
    std::cout << "n_case=" << n_case << " n_fail=" << n_fail << "\n";
    return n_fail == 0 ? 0 : 1;
  } catch (permutalib::PermutalibException const &e) {
    std::cerr << "PermutalibException\n";
    return 2;
  }
}
