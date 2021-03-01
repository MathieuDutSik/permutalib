#ifndef PLUS_INFINITY_INCLUDE
#define PLUS_INFINITY_INCLUDE

#include <iostream>



template<typename T>
struct Tplusinfinity {
  Tplusinfinity(bool const& val1, T const& val2)
  {
    IsInf = val1;
    value = val2;
  }
  Tplusinfinity(T const& val)
  {
    IsInf = false;
    value = val;
  }
  void SetToInfinity()
  {
    IsInf = true;
  }
  Tplusinfinity<T> operator=(Tplusinfinity<T> const& x)
  {
    IsInf = x.IsInf;
    value = x.value;
    return *this;
  }
  Tplusinfinity<T> operator=(T const& val)
  {
    IsInf = false;
    value = val;
    return *this;
  }
  bool IsInfinity() const
  {
    return IsInf;
  }
  T GetValue() const
  {
    return value;
  }
private:
  bool IsInf;
  T value;
};


template<typename T>
bool operator==(Tplusinfinity<T> const& x, Tplusinfinity<T> const& y)
{
  if (x.IsInfinity() && y.IsInfinity())
    return true;
  if (x.IsInfinity() != y.IsInfinity())
    return false;
  return x.GetValue() == y.GetValue();
}


template<typename T>
bool operator<(Tplusinfinity<T> const& x, Tplusinfinity<T> const& y)
{
  if (x.IsInfinity() && y.IsInfinity())
    return false;
  if (!x.IsInfinity() && y.IsInfinity())
    return true;
  if (x.IsInfinity() && !y.IsInfinity())
    return false;
  return x.GetValue() < y.GetValue();
}

template<typename T>
std::ostream& operator<<(std::ostream& os, Tplusinfinity<T> const& x)
{
  if (x.IsInfinity())
    os << "infinity";
  else
    os << x.GetValue();
  return os;
}


template<typename T>
bool operator<=(Tplusinfinity<T> const& x, Tplusinfinity<T> const& y)
{
  if (x.IsInfinity() && y.IsInfinity())
    return true;
  if (!x.IsInfinity() && y.IsInfinity())
    return true;
  if (x.IsInfinity() && !y.IsInfinity())
    return false;
  return x.GetValue() <= y.GetValue();
}

template<typename T>
bool operator<=(T const& x, Tplusinfinity<T> const& y)
{
  if (y.IsInfinity())
    return true;
  return x <= y.GetValue();
}


template<typename T>
bool operator>(Tplusinfinity<T> const& x, Tplusinfinity<T> const& y)
{
  if (x.IsInfinity() && y.IsInfinity())
    return false;
  if (!x.IsInfinity() && y.IsInfinity())
    return false;
  if (x.IsInfinity() && !y.IsInfinity())
    return true;
  return x.GetValue() > y.GetValue();
}


template<typename T>
Tplusinfinity<T> operator+(Tplusinfinity<T> const& x, Tplusinfinity<T> const& y)
{
  if (x.IsInfinity() || y.IsInfinity())
    return Tplusinfinity<T>(true, 0);
  return Tplusinfinity<T>(false, x.GetValue() + y.GetValue());
}

template<typename T>
Tplusinfinity<T> operator+(Tplusinfinity<T> const& x, T const& y)
{
  return Tplusinfinity<T>(x.IsInfinity(), x.GetValue() + y);
}

template<typename T>
Tplusinfinity<T> operator+(T const& x, Tplusinfinity<T> const& y)
{
  return Tplusinfinity<T>(y.IsInfinity(), x + y.GetValue());
}

template<typename T>
Tplusinfinity<T> operator*(Tplusinfinity<T> const& x, Tplusinfinity<T> const& y)
{
  if (x.IsInfinity() || y.IsInfinity())
    return Tplusinfinity<T>(true, 0);
  return Tplusinfinity<T>(false, x.GetValue() * y.GetValue());
}

template<typename T>
Tplusinfinity<T> operator*(Tplusinfinity<T> const& x, T const& y)
{
  return Tplusinfinity<T>(x.IsInfinity(), x.GetValue() * y);
}

template<typename T>
Tplusinfinity<T> operator*(T const& x, Tplusinfinity<T> const& y)
{
  return Tplusinfinity<T>(y.IsInfinity(), x * y.GetValue());
}




#endif
