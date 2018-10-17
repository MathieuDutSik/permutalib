Permutation Group Library
=========================

This code contains permutation group code adapted from GAP.
The goal is to have permutation groups and the partition
backtrack.


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


Programming differences
-----------------------

There are some differences that we have decided to do with
the existing GAP code:

  * The GAP uses a shared pointer (i.e. std::shared_ptr) semantic for storing the permutation. This can be seen by taking a very long permutation and taking 10000 copy of it: memory usage barely changes. We do not use shared pointer and instead use a std::vector and we store the common vectors.
  * The GAP code uses a recursive data structure (i.e. struct GRP { ....., GRP* stab}). We decided to put this as a vector. This means that we have to put slices when doing the calls, that is put the object and the stabilizer level.
  * The GAP code has deterministic random algorithms. It means that if you run again a GAP program you get exactly the same result. The C++ code uses simple rand() but we have the option of getting random number as in GAP in order to reproduce bugs.



Downloading code
----------------

Full source code to be downloaded is from

git clone git@github.com:MathieuDutSik/permutalib.git --recursive

or

git clone http://github.com/MathieuDutSik/permutalib.git --recursive



Licensing
---------

The GPLv2 licensing is available from the GAP.
The LGPL licensing is available so that the library can be used in other programs.

Contact information
-------------------

Contact Mathieu Dutour Sikiric at Mathieu.Dutour@gmail.com in case of any question.
