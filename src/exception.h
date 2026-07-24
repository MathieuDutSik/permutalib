// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GAP_EXCEPTION_H_
#define SRC_GAP_EXCEPTION_H_

#include <cstdint>
#include <cstdlib>
#include <type_traits>

namespace permutalib {

struct PermutalibException {
  int eVal;
};

// Portable pseudo-random source. POSIX random() is unavailable on Windows
// (MinGW / MSVC), and defining a global random() shim collides with other
// libraries that provide their own (e.g. basic_common_cpp), so use a namespaced
// helper instead. On platforms that have POSIX random() we keep it; on Windows
// we combine two std::rand() draws (RAND_MAX can be as low as 32767) to widen
// the range to the [0, 2^31 - 1] that callers expect.
inline long permutalib_random() {
#if defined(_WIN32)
  return (static_cast<long>(std::rand()) << 15) ^ static_cast<long>(std::rand());
#else
  return random();
#endif
}

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
