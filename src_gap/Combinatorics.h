#ifndef SRC_GAP_COMBINATORICS_H_
#define SRC_GAP_COMBINATORICS_H_

#include <algorithm>
#include <vector>

template <typename T1, typename T2>
void SortParallel_PairList(std::vector<T1> &list1, std::vector<T2> &list2) {
  struct PairData {
    T1 val1;
    T2 val2;
  };
  size_t len = list1.size();
  std::vector<PairData> ListPair(len);
  for (size_t i = 0; i < len; i++)
    ListPair[i] = {list1[i], list2[i]};
  sort(ListPair.begin(), ListPair.end(),
       [](PairData const &a, PairData const &b) -> bool {
         return a.val1 < b.val1;
       });
  for (size_t i = 0; i < len; i++) {
    list1[i] = ListPair[i].val1;
    list2[i] = ListPair[i].val2;
  }
}

// clang-format off
#endif  // SRC_GAP_COMBINATORICS_H_
// clang-format on
