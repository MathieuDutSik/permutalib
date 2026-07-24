// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GAP_EXCEPTION_H_
#define SRC_GAP_EXCEPTION_H_

#include <cstdint>
#include <cstdlib>
#include <type_traits>

// permutalib calls the unqualified global random(). On POSIX and on MinGW,
// <cstdlib> declares it; including it here (before the templates that use it)
// satisfies two-phase name lookup. Only non-MinGW Windows (e.g. MSVC) lacks
// POSIX random()/srandom(), so provide a shim there built on rand(). RAND_MAX
// can be as low as 32767, so combine two draws to widen the range.
#if defined(_WIN32) && !defined(__MINGW32__)
inline long random() {
  return (static_cast<long>(std::rand()) << 15) ^ static_cast<long>(std::rand());
}
inline void srandom(unsigned int seed) { std::srand(seed); }
#endif

namespace permutalib {

struct PermutalibException {
  int eVal;
};

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
    uint64_t v = static_cast<uint64_t>(val);
    Tint ret = Tint(static_cast<uint32_t>(v >> 32));
    ret <<= 32;
    ret += Tint(static_cast<uint32_t>(v & 0xFFFFFFFFu));
    return ret;
  }
}

// clang-format off
}  // namespace permutalib
#endif  // SRC_GAP_EXCEPTION_H_
// clang-format on
