// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GAP_TESTINGFCT_H_
#define SRC_GAP_TESTINGFCT_H_

namespace permutalib {

template<typename T>
struct TimeEval_perm {
  std::chrono::time_point<std::chrono::system_clock> time;
  TimeEval_perm() { time = std::chrono::system_clock::now(); }
  int64_t eval() {
    std::chrono::time_point<std::chrono::system_clock> timeNew =
        std::chrono::system_clock::now();
    int64_t delta = std::chrono::duration_cast<T>(timeNew - time).count();
    time = timeNew;
    return delta;
  }
};

template <typename T>
std::ostream &operator<<(std::ostream &os, TimeEval_perm<T> &x) {
  os << x.eval();
  return os;
}

using MicrosecondTime_perm = TimeEval_perm<std::chrono::microseconds>;
using NanosecondTime_perm = TimeEval_perm<std::chrono::nanoseconds>;


// clang-format off
}  // namespace permutalib
#endif  // SRC_GAP_GAPPRINT_H_
// clang-format on
