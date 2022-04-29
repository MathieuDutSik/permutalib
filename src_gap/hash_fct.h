#ifndef SRC_GAP_HASH_FCT_H_
#define SRC_GAP_HASH_FCT_H_

#include <cstring>

namespace permutalib {

static inline uint32_t murmur_32_scramble(uint32_t k) {
  k *= 0xcc9e2d51;
  k = (k << 15) | (k >> 17);
  k *= 0x1b873593;
  return k;
}

uint32_t murmur3_32(const uint8_t *key, size_t len, uint32_t seed) {
  uint32_t h = seed;
  uint32_t k;
  /* Read in groups of 4. */
  for (size_t i = len >> 2; i; i--) {
    // Here is a source of differing results across endiannesses.
    // A swap here has no effects on hash properties though.
    memcpy(&k, key, sizeof(uint32_t));
    key += sizeof(uint32_t);
    h ^= murmur_32_scramble(k);
    h = (h << 13) | (h >> 19);
    h = h * 5 + 0xe6546b64;
  }
  /* Read the rest. */
  k = 0;
  for (size_t i = len & 3; i; i--) {
    k <<= 8;
    k |= key[i - 1];
  }
  // A swap is *not* necessary here because the preceding loop already
  // places the low bytes in the low places according to whatever endianness
  // we use. Swaps only apply when the memory is copied in a chunk.
  h ^= murmur_32_scramble(k);
  /* Finalize. */
  h ^= uint32_t(len);
  h ^= h >> 16;
  h *= 0x85ebca6b;
  h ^= h >> 13;
  h *= 0xc2b2ae35;
  h ^= h >> 16;
  return h;
}

template <typename T> inline T unaligned_load(void const *ptr) noexcept {
  // using memcpy so we don't get into unaligned load problems.
  // compiler should optimize this very well anyways.
  T t;
  std::memcpy(&t, ptr, sizeof(T));
  return t;
}

inline size_t robin_hood_hash_bytes(void const *ptr, size_t len,
                                    const uint64_t &seed) noexcept {
  static constexpr uint64_t m = UINT64_C(0xc6a4a7935bd1e995);
  //    static constexpr uint64_t seed = UINT64_C(0xe17a1465);
  static constexpr unsigned int r = 47;

  auto const *const data64 = static_cast<uint64_t const *>(ptr);
  uint64_t h = seed ^ (len * m);

  size_t const n_blocks = len / 8;
  for (size_t i = 0; i < n_blocks; ++i) {
    auto k = unaligned_load<uint64_t>(data64 + i);

    k *= m;
    k ^= k >> r;
    k *= m;

    h ^= k;
    h *= m;
  }
  auto const *const data8 =
      reinterpret_cast<uint8_t const *>(data64 + n_blocks);
  switch (len & 7U) {
  case 7:
    h ^= static_cast<uint64_t>(data8[6]) << 48U;
    [[fallthrough]];
  case 6:
    h ^= static_cast<uint64_t>(data8[5]) << 40U;
    [[fallthrough]];
  case 5:
    h ^= static_cast<uint64_t>(data8[4]) << 32U;
    [[fallthrough]];
  case 4:
    h ^= static_cast<uint64_t>(data8[3]) << 24U;
    [[fallthrough]];
  case 3:
    h ^= static_cast<uint64_t>(data8[2]) << 16U;
    [[fallthrough]];
  case 2:
    h ^= static_cast<uint64_t>(data8[1]) << 8U;
    [[fallthrough]];
  case 1:
    h ^= static_cast<uint64_t>(data8[0]);
    h *= m;
    [[fallthrough]];
  default:
    break;
  }

  h ^= h >> r;

  // not doing the final step here, because this will be done by keyToIdx
  // anyways h *= m; h ^= h >> r;
  return static_cast<size_t>(h);
}

// clang-format off
}  // namespace permutalib
#endif  // SRC_GAP_HASH_FCT_H_
// clang-format on
