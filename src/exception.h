// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GAP_EXCEPTION_H_
#define SRC_GAP_EXCEPTION_H_

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <type_traits>

/*
  The random numbers of permutalib, the same on every platform.

  When the program has the generator of basic_common_cpp (Basic_random.h on
  the include path, as in polyhedral_common), permutalib draws from it, so
  that one set_random_seed drives both libraries. Otherwise permutalib uses
  its own: a std::mt19937_64 per thread, whose output the C++ standard fixes,
  with a fixed seed, and the ranges formed here rather than by the std
  distributions, whose algorithms are left to each library. The bounded draws
  use the method of Basic_random.h; the two seed their engines differently,
  so their sequences differ.
 */
#if __has_include("Basic_random.h")
#include "Basic_random.h"
#define PERMUTALIB_BASIC_RANDOM
#else
#include <random>
#endif

namespace permutalib {

struct PermutalibException {
  int eVal;
};

#ifdef PERMUTALIB_BASIC_RANDOM

// Uniform in [0, n), n > 0.
inline size_t random_index(size_t n) { return ::random_index(n); }

// Uniform in the closed range [lo, hi], lo <= hi.
inline int random_int(int lo, int hi) { return ::random_int(lo, hi); }

#else

namespace random_detail {

inline uint64_t random_u64() {
  thread_local std::mt19937_64 engine(0x5eed5eed5eed5eedULL);
  return engine();
}

// Uniform in [0, n), n > 0, without bias: the draws below 2^64 mod n are
// rejected.
inline uint64_t random_below(uint64_t n) {
  uint64_t threshold = (0 - n) % n;
  while (true) {
    uint64_t x = random_u64();
    if (x >= threshold)
      return x % n;
  }
}

} // namespace random_detail

// Uniform in [0, n), n > 0.
inline size_t random_index(size_t n) {
  return static_cast<size_t>(
      random_detail::random_below(static_cast<uint64_t>(n)));
}

// Uniform in the closed range [lo, hi], lo <= hi.
inline int random_int(int lo, int hi) {
  int64_t lo64 = lo;
  int64_t hi64 = hi;
  uint64_t width = static_cast<uint64_t>(hi64 - lo64) + 1;
  int64_t offset = static_cast<int64_t>(random_detail::random_below(width));
  return static_cast<int>(lo64 + offset);
}

#endif

// Construct an integer type Tint from an unsigned (size / count) value. On LLP64
// (Windows) size_t and other 64-bit types are `long long`, for which gmpxx
// provides no constructor, making Tint(val) ambiguous when Tint is mpz_class /
// mpq_class. Build the value from 32-bit halves, which every Tint accepts;
// values that fit in 32 bits take the direct constructor. The source must be an
// unsigned integer type (the >> 32 decomposition assumes a non-negative value).
template <typename Tint, typename Tin>
inline Tint UnsignedToTint(Tin const &val) {
  static_assert(std::is_unsigned_v<Tin>,
                "UnsignedToTint expects an unsigned integer source type");
  if constexpr (sizeof(Tin) <= 4) {
    return Tint(val);
  } else {
    // 64-bit source: gmpxx has no long long constructor (ambiguous on LLP64),
    // so build the value from 32-bit halves. Multiply by 2^16 twice rather than
    // using <<= or a 2^32 literal: some Tint types (e.g. SafeInt64) provide no
    // operator<<=, and 2^32 is itself a long long literal that gmpxx would find
    // ambiguous.
    uint64_t v = static_cast<uint64_t>(val);
    Tint ret = Tint(static_cast<uint32_t>(v >> 32));
    ret = ret * Tint(65536) * Tint(65536);
    ret = ret + Tint(static_cast<uint32_t>(v & 0xFFFFFFFFu));
    return ret;
  }
}

// clang-format off
}  // namespace permutalib
#endif  // SRC_GAP_EXCEPTION_H_
// clang-format on
