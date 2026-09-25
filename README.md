Permutation Group Library
=========================

This code contains permutation group code adapted from GAP.
The goal is to have permutation groups and the partition
backtrack.



Downloading code
----------------

Full source code to be downloaded is from

```sh
$ git clone git@github.com:MathieuDutSik/permutalib.git --recursive
```

or

```sh
$ git clone http://github.com/MathieuDutSik/permutalib.git --recursive
```






Usage
-----

* The first step is to choose the data type for the permutation. For example:

```cpp
using Tidx = int16_t;
using Telt = permutalib::SingleSidedPerm<Tidx>;
```

The number of elements on which the group can act is 2^16 - 1.

* The permutation element is build as

```cpp
std::vector<Tidx> eList(10);
Telt eElt(eList);
```

* The second step is to choose an integer type. For example:

```cpp
using Tint = mpz_class;
```

This integer type is used for computing order of groups. Since the size can grow pretty large (The symmetric group has n! elements) an arbitrary large integer type is needed.

* The permutation group is built as

```cpp
std::vector<Telt> ListGen;
permutalib::Group<Telt,Tint> eG(ListGen, n);
```

* The subset of 0..n-1 is built using the boost dynamic bitset (https://www.boost.org/doc/libs/1_36_0/libs/dynamic_bitset/dynamic_bitset.html). The corresponding code is:

```cpp
permutalib::Face subset1, subset1;
permutalib::Group<Telt,Tint> stab = eG.Stabilizer_OnSets(subset1);
std::pair<bool, Telt> test = eG.RepresentativeAction_OnSets(subset1, subset2);
Face subset1_can = eG.CanonicalImage(subset1);
```

* See Group.h for the full functionality and the examples.


Group resolutions
-----------------

`Resolution.h` computes free ZG-resolutions of Z for a finite permutation
group G, following the algorithm of G. Ellis, "Computing group resolutions",
J. Symbolic Computation 38 (2004), as implemented by `ResolutionFiniteGroup`
in the GAP package HAP. It is used as

```cpp
permutalib::FiniteGroupResolution<Telt, Tint> R(eG, 5);
size_t rank2 = R.dimension(2);                      // rank of R_2
permutalib::ResolutionChain const& d = R.boundary(2, 0); // boundary of the first generator of R_2
std::vector<Tint> H3 = R.integral_homology(3);      // abelian invariants of H_3(G, Z)
R.extend();                                          // one more term
```

A chain is a list of terms `coeff * (elt . e_cell)` with `elt` an index in
`R.elements()`. The contracting homotopy is available as `R.homotopy(i, term)`
and `R.GapString()` prints the resolution in the word format of HAP.

The homology is returned in the divisibility form `d_1 | d_2 | ...` of
`Homology(TensorWithIntegers(R), n)` in HAP; `R.integral_homology_prime_power(n)`
gives the prime power form of `AbelianInvariants` in GAP. The dimension of
`H_n(G, F_p)` is `R.homology_dimension_mod_p(n, p)` and
`R.check_homology_consistency(n)` verifies the finiteness of the homology, the
universal coefficient theorem and the abelianization (these checks run
automatically when compiling with `-DDEBUG_RESOLUTION`). With the same
generators, the resolution is the one of `ResolutionFiniteGroup` in HAP; for
instance `SymmetricGroup(5)` with length 5 gives the ranks `[1,4,10,20,35,56]`
and sizes `[8,38,100,204,340]` recorded in the tests of HAP.
`TestResolution` checks the resolutions of the groups in
`CI_tests/11_TestResolution/GroupsHomology` against values from HAP and the
literature, and `GapResolution` writes a resolution for comparison with HAP.




Rationale
---------

The code in GAP is a very good basis for computing with
permutation group and this author was very satisfied with its
speed and functionality (only one case related to shortest
vectors of Leech lattice created problems).

However, GAP itself had some problems:

  * While the permutation code is very good, GAP itself is slow
  * We want the code of permutation group accessible as a library.
  * We want parallel code using groups.
  * The benefits of C++ (speed, templates) are very attractive for this kind of code.

The permlib code (https://github.com/tremlin/PermLib) provided a solution to this problem.
It is a reimplementation of the algorithms and generally very satisfying.
However, in some cases it was very very slow compared to the GAP
code. This made it unusable in many contexts.

Thus the idea is to simply recode the GAP code into C++ in order
to achieve this. Idea is really not to try to be too clever and
adapt the code accordingly. We work with gap-4.7.8 as reference
gap source code implementation.

The library is single threaded but it can be used into multithreaded
code because the code is thread safe. There is no global variable
used.



Design choices
--------------

There are some differences that we have decided to do with
the existing GAP code:

  * The GAP uses a shared pointer (i.e. std::shared_ptr) semantic for storing the permutation. This can be seen by taking a very long permutation and taking 10000 copy of it: memory usage barely changes. We do not use shared pointer and instead use a std::vector and we store the common vectors.
  * The GAP code uses a linked list recursive data structure (i.e. struct GRP { ....., GRP* stab}). We follow this convention as well after trying to use a std::vector for storing.
  * The GAP code has deterministic random algorithms. It means that if you run again a GAP program you get exactly the same result. The C++ code uses simple rand() but we have the option of getting random number as in GAP in order to reproduce bugs.
  * The permutalib code uses just a single type for the permutation.



Licensing
---------

The GPLv2 licensing is available from the GAP.
The LGPL licensing is available so that the library can be used in other programs.



Contact information
-------------------

Contact Mathieu Dutour Sikiric at Mathieu.Dutour@gmail.com in case of any question.
