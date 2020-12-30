#ifndef GAP_PRINT_INCLUDE
#define GAP_PRINT_INCLUDE


std::string GapStringIntVector(std::vector<int> const& f)
{
  std::string str;
  str += "[ ";
  int len=f.size();
  for (int i=0; i<len; i++) {
    if (i>0)
      str += ", ";
    str += std::to_string(f[i]+1);
  }
  str += " ]";
  return str;
}

std::string GapStringBool(bool const& v)
{
  if (v)
    return std::string("true");
  else
    return std::string("false");
}


std::string GapStringBoolVector(Face const& f)
{
  std::string str;
  str += "[ ";
  int len=f.size();
  for (int i=0; i<len; i++) {
    if (i>0)
      str += ", ";
    if (f[i])
      str += "true";
    else
      str += "false";
  }
  str += " ]";
  return str;
}

std::string ConstrainedIntInfinity_to_string(int const& val, int const& n)
{
  if (val < n) {
    return std::to_string(val+1);
  } else {
    return std::string("infinity");
  }
}


std::string PosFail_to_string(int const& pos)
{
  if (pos == -1)
    return std::string("fail");
  return std::to_string(pos + 1);
}


std::string PosFalse_to_string(int const& pos)
{
  if (pos == -1)
    return std::string("false");
  return std::to_string(pos + 1);
}




template<typename T>
std::string GapStringTVector(std::vector<T> const& f)
{
  std::ostringstream os;
  os << "[ ";
  int len=f.size();
  for (int i=0; i<len; i++) {
    if (i>0)
      os << ", ";
    os << f[i];
  }
  os << " ]";
  return os.str();
}


template<typename T>
std::string GapStringMissingTVector(std::vector<std::optional<T>> const& f)
{
  int len=f.size();
  std::vector<bool> Vstatus(len);
  for (int i=0; i<len; i++) {
    bool status=false;
    if (f[i])
      status=true;
    Vstatus[i] = status;
  }
  int last_nz = 0;
  int first_nz = -1;
  for (int i=0; i<len; i++) {
    if (Vstatus[i])
      last_nz = i+1;
    if (first_nz == -1) {
      if (Vstatus[i]) {
        first_nz = i;
      }
    }
  }
  std::ostringstream os;
  if (first_nz == -1) {
    return "[  ]";
  }
  os << "[";
  if (first_nz>0) {
    os << " ";
    for (int i=0; i<first_nz; i++)
      os << ",";
  }
  for (int i=first_nz; i<last_nz; i++) {
    if (i > first_nz)
      os << ",";
    if (Vstatus[i])
      os << " " << *(f[i]);
  }
  os << " ]";
  return os.str();
}


#endif
