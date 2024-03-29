// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GAP_GAPPRINT_H_
#define SRC_GAP_GAPPRINT_H_

#include "Face_basic.h"
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

namespace permutalib {

static const int int_reducedm1 = -1;
static const int int_false = 0;
static const int int_true = 1;
static const int int_fail = 2;
static const int int_int = 3;
static const int int_perm = 4;
static const int int_group = 5;
static const int int_stablev = 6;


template<typename T>
T ReadScalar(std::string const& estr) {
  T ret_val;
  std::istringstream is(estr);
  is >> ret_val;
  return ret_val;
}

template <typename F> std::string GapString_F(size_t const &len, F const &f) {
  std::string str = "[ ";
  for (size_t i = 0; i < len; i++) {
    if (i > 0)
      str += ", ";
    str += f(i);
  }
  str += " ]";
  return str;
}

template <typename T> std::string GapStringTVector(std::vector<T> const &f) {
  std::ostringstream os;
  os << "[ ";
  size_t len = f.size();
  for (size_t i = 0; i < len; i++) {
    if (i > 0)
      os << ", ";
    os << f[i];
  }
  os << " ]";
  return os.str();
}

template <typename T> std::string GapStringTVectorB(std::vector<T> const &f) {
  if (f.size() == 0)
    return std::string("[ ]");
  return GapStringTVector(f);
}

template <typename Tidx>
std::string GapStringIntVector(std::vector<Tidx> const &ev) {
  auto f = [&](size_t const &i) -> std::string {
    return std::to_string(ev[i] + 1);
  };
  return GapString_F(ev.size(), f);
}

std::string GapStringTrueFalseFail(int const &v) {
  if (v == int_true)
    return std::string("true");
  if (v == int_false)
    return std::string("false");
  if (v == int_fail)
    return std::string("fail");
  return std::string("Unsupported case in GapStringTrueFalseFail\n");
}

std::string GapStringBool(bool const &v) {
  if (v)
    return std::string("true");
  else
    return std::string("false");
}

std::string GapStringBoolVector(Face const &f) {
  auto f_str = [&](size_t const &i) -> std::string {
    return GapStringBool(f[i]);
  };
  return GapString_F(f.size(), f_str);
}

std::string GapStringSetBoolVector(Face const &f) {
  std::vector<int> eS;
  size_t len = f.size();
  for (size_t i = 0; i < len; i++)
    if (f[i])
      eS.push_back(static_cast<int>(i));
  return GapStringIntVector(eS);
}

std::string GapStringListBoolVector(std::vector<Face> const &f) {
  if (f.size() == 0)
    return std::string("[ ]");
  std::vector<std::string> Lstr(f.size());
  for (size_t i = 0; i < f.size(); i++)
    Lstr[i] = GapStringBoolVector(f[i]);
  return GapStringTVector(Lstr);
}

std::string GapStringBoolVectorB(std::vector<int8_t> const &f) {
  auto f_str = [&](size_t const &i) -> std::string {
    if (f[i] == 1)
      return "true";
    return "false";
  };
  return GapString_F(f.size(), f_str);
}

std::string ConstrainedIntInfinity_to_string(int const &val, int const &n) {
  if (val < n) {
    return std::to_string(val + 1);
  } else {
    return std::string("infinity");
  }
}

template <typename Tidx> std::string PosFail_to_string(Tidx const &pos) {
  if (pos == std::numeric_limits<Tidx>::max())
    return std::string("fail");
  return std::to_string(pos + 1);
}

template <typename Tidx> std::string PosFalse_to_string(Tidx const &pos) {
  if (pos == std::numeric_limits<Tidx>::max())
    return std::string("false");
  return std::to_string(pos + 1);
}

template <typename T>
std::string GapStringMissingTVector(std::vector<std::optional<T>> const &f) {
  size_t len = f.size();
  std::vector<bool> Vstatus(len);
  for (size_t i = 0; i < len; i++) {
    bool status = false;
    if (f[i])
      status = true;
    Vstatus[i] = status;
  }
  size_t last_nz = 0;
  size_t first_nz = std::numeric_limits<size_t>::max();
  for (size_t i = 0; i < len; i++) {
    if (Vstatus[i])
      last_nz = i + 1;
    if (first_nz == std::numeric_limits<size_t>::max()) {
      if (Vstatus[i]) {
        first_nz = i;
      }
    }
  }
  std::ostringstream os;
  if (first_nz == std::numeric_limits<size_t>::max()) {
    return "[  ]";
  }
  os << "[";
  if (first_nz > 0) {
    os << " ";
    for (size_t i = 0; i < first_nz; i++)
      os << ",";
  }
  for (size_t i = first_nz; i < last_nz; i++) {
    if (i > first_nz)
      os << ",";
    if (Vstatus[i])
      os << " " << *(f[i]);
  }
  os << " ]";
  return os.str();
}

// clang-format off
}  // namespace permutalib
#endif  // SRC_GAP_GAPPRINT_H_
// clang-format on
