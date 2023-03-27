// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GAP_TESTINGFCT_H_
#define SRC_GAP_TESTINGFCT_H_

namespace permutalib {

template<typename T>
struct TimeEval {
  std::chrono::time_point<std::chrono::system_clock> time;
  TimeEval() { time = std::chrono::system_clock::now(); }
  int64_t eval() {
    std::chrono::time_point<std::chrono::system_clock> timeNew =
        std::chrono::system_clock::now();
    int64_t delta = std::chrono::duration_cast<T>(timeNew - time).count();
    time = timeNew;
    return delta;
  }
};

using MicrosecondTime = TimeEval<std::chrono::microseconds>;
using NanosecondTime = TimeEval<std::chrono::nanoseconds>;


// clang-format off
}  // namespace permutalib
#endif  // SRC_GAP_GAPPRINT_H_
// clang-format on
