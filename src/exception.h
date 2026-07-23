// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GAP_EXCEPTION_H_
#define SRC_GAP_EXCEPTION_H_

#include <cstdint>
#include <type_traits>

#if defined(_WIN32)
// MinGW / Windows provides no POSIX random() / srandom(). permutalib calls the
// unqualified global random(); supply portable equivalents built on rand().
// POSIX random() returns a value in [0, 2^31-1], but MinGW RAND_MAX is only
// 32767, so combine two rand() draws to widen the range.
#include <cstdlib>
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
