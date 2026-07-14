/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   sfloat: project-wide scalar type for the physics core (AD seam).

   Stock build (default)      : sfloat == double, zero overhead,
                                behavior must be bit-identical to
                                upstream SPARTA.
   AD build (-DSPARTA_AD)     : sfloat == forward-mode dual number with
                                4 derivative components (the 4 Bezier
                                design parameters alpha = x1,y1,x2,y2).

   Design rules of the AD type:
   - implicit construction FROM double (safe: derivative = 0)
   - NO implicit conversion TO double (derivative would silently die);
     use spval(x) or x.val(); C-style casts like (int)x still work via
     the explicit operator
   - user-provided copy ctor makes the type non-trivially copyable, so
     passing it through printf-style varargs is a COMPILE ERROR: every
     print site must wrap with spval(), making derivative truncation
     visible and deliberate
   - default ctor zero-initializes (value and derivatives): malloc'd
     arrays are handled separately (memory.cpp uses calloc under
     SPARTA_AD so array derivatives start at zero, not garbage)
   - comparisons compare values only (branches follow the primal flow;
     this IS the frozen-stream/pathwise semantics)
   - layout is exactly 5 doubles, no padding: MPI stub datatype
     MPI_5DOUBLE (see STUBS/mpi.h) moves it by size; MPI_SFLOAT below
     selects the right datatype per build

   Passive by design (never sfloat): RNG internals/outputs, timers,
   binary file formats, ubuf bit-punning unions.
------------------------------------------------------------------------- */

#ifndef SPARTA_SFLOAT_H
#define SPARTA_SFLOAT_H

#include <cmath>

#ifndef SPARTA_AD

/* ---------------- stock build: plain double ---------------- */

typedef double sfloat;
template <typename T> inline T spval(T x) { return x; }   // identity
#define MPI_SFLOAT MPI_DOUBLE

#else

/* ---------------- AD build: forward dual, 4 directions ---------------- */

#define SFLOAT_NDIR 4

class sfloat {
 public:
  double v;                 // value
  double d[SFLOAT_NDIR];    // derivative components

  sfloat() : v(0.0) { zerod(); }
  sfloat(double x) : v(x) { zerod(); }              // implicit from double
  sfloat(const sfloat &o) : v(o.v) {                // user-provided on
    for (int i = 0; i < SFLOAT_NDIR; i++) d[i] = o.d[i];   // purpose:
  }                                                 // varargs guard
  sfloat &operator=(const sfloat &o) {
    v = o.v;
    for (int i = 0; i < SFLOAT_NDIR; i++) d[i] = o.d[i];
    return *this;
  }

  double val() const { return v; }
  explicit operator double() const { return v; }    // (double)x with spval semantics
  explicit operator int() const { return (int) v; } // static_cast<int>(x)
  explicit operator long() const { return (long) v; }
  explicit operator long long() const { return (long long) v; }
  explicit operator bool() const { return v != 0.0; }
  void zerod() { for (int i = 0; i < SFLOAT_NDIR; i++) d[i] = 0.0; }

  // compound assignment

  sfloat &operator+=(const sfloat &o) {
    v += o.v;
    for (int i = 0; i < SFLOAT_NDIR; i++) d[i] += o.d[i];
    return *this;
  }
  sfloat &operator-=(const sfloat &o) {
    v -= o.v;
    for (int i = 0; i < SFLOAT_NDIR; i++) d[i] -= o.d[i];
    return *this;
  }
  sfloat &operator*=(const sfloat &o) {
    for (int i = 0; i < SFLOAT_NDIR; i++) d[i] = d[i]*o.v + v*o.d[i];
    v *= o.v;
    return *this;
  }
  sfloat &operator/=(const sfloat &o) {
    for (int i = 0; i < SFLOAT_NDIR; i++)
      d[i] = (d[i]*o.v - v*o.d[i])/(o.v*o.v);
    v /= o.v;
    return *this;
  }
};

// increment/decrement (converted double counters): value +/- 1, deriv kept

inline sfloat &operator++(sfloat &a) { a.v += 1.0; return a; }
inline sfloat &operator--(sfloat &a) { a.v -= 1.0; return a; }
inline sfloat operator++(sfloat &a, int) { sfloat t(a); a.v += 1.0; return t; }
inline sfloat operator--(sfloat &a, int) { sfloat t(a); a.v -= 1.0; return t; }

// arithmetic (mixed double/sfloat handled by the implicit ctor)

inline sfloat operator-(const sfloat &a) {
  sfloat r(a);
  r.v = -r.v;
  for (int i = 0; i < SFLOAT_NDIR; i++) r.d[i] = -r.d[i];
  return r;
}
inline sfloat operator+(const sfloat &a) { return a; }

inline sfloat operator+(const sfloat &a, const sfloat &b) {
  sfloat r(a); r += b; return r;
}
inline sfloat operator-(const sfloat &a, const sfloat &b) {
  sfloat r(a); r -= b; return r;
}
inline sfloat operator*(const sfloat &a, const sfloat &b) {
  sfloat r(a); r *= b; return r;
}
inline sfloat operator/(const sfloat &a, const sfloat &b) {
  sfloat r(a); r /= b; return r;
}

// comparisons: value-based (primal control flow)

inline bool operator<(const sfloat &a, const sfloat &b)  { return a.v <  b.v; }
inline bool operator>(const sfloat &a, const sfloat &b)  { return a.v >  b.v; }
inline bool operator<=(const sfloat &a, const sfloat &b) { return a.v <= b.v; }
inline bool operator>=(const sfloat &a, const sfloat &b) { return a.v >= b.v; }
inline bool operator==(const sfloat &a, const sfloat &b) { return a.v == b.v; }
inline bool operator!=(const sfloat &a, const sfloat &b) { return a.v != b.v; }

// math functions (chain rule); cover what the SPARTA core uses

inline sfloat sqrt(const sfloat &a) {
  sfloat r;
  r.v = std::sqrt(a.v);
  double g = (r.v == 0.0) ? 0.0 : 0.5/r.v;
  for (int i = 0; i < SFLOAT_NDIR; i++) r.d[i] = g*a.d[i];
  return r;
}
inline sfloat fabs(const sfloat &a) {
  sfloat r(a);
  if (a.v < 0.0) { r.v = -r.v;
    for (int i = 0; i < SFLOAT_NDIR; i++) r.d[i] = -r.d[i]; }
  return r;
}
inline sfloat exp(const sfloat &a) {
  sfloat r;
  r.v = std::exp(a.v);
  for (int i = 0; i < SFLOAT_NDIR; i++) r.d[i] = r.v*a.d[i];
  return r;
}
inline sfloat log(const sfloat &a) {
  sfloat r;
  r.v = std::log(a.v);
  for (int i = 0; i < SFLOAT_NDIR; i++) r.d[i] = a.d[i]/a.v;
  return r;
}
inline sfloat log2(const sfloat &a) {
  return log(a)/std::log(2.0);
}
inline sfloat log10(const sfloat &a) {
  return log(a)/std::log(10.0);
}
inline sfloat sin(const sfloat &a) {
  sfloat r;
  r.v = std::sin(a.v);
  double c = std::cos(a.v);
  for (int i = 0; i < SFLOAT_NDIR; i++) r.d[i] = c*a.d[i];
  return r;
}
inline sfloat cos(const sfloat &a) {
  sfloat r;
  r.v = std::cos(a.v);
  double s = -std::sin(a.v);
  for (int i = 0; i < SFLOAT_NDIR; i++) r.d[i] = s*a.d[i];
  return r;
}
inline sfloat tan(const sfloat &a) {
  sfloat r;
  r.v = std::tan(a.v);
  double g = 1.0 + r.v*r.v;
  for (int i = 0; i < SFLOAT_NDIR; i++) r.d[i] = g*a.d[i];
  return r;
}
inline sfloat asin(const sfloat &a) {
  sfloat r;
  r.v = std::asin(a.v);
  double g = 1.0/std::sqrt(1.0 - a.v*a.v);
  for (int i = 0; i < SFLOAT_NDIR; i++) r.d[i] = g*a.d[i];
  return r;
}
inline sfloat acos(const sfloat &a) {
  sfloat r;
  r.v = std::acos(a.v);
  double g = -1.0/std::sqrt(1.0 - a.v*a.v);
  for (int i = 0; i < SFLOAT_NDIR; i++) r.d[i] = g*a.d[i];
  return r;
}
inline sfloat atan(const sfloat &a) {
  sfloat r;
  r.v = std::atan(a.v);
  double g = 1.0/(1.0 + a.v*a.v);
  for (int i = 0; i < SFLOAT_NDIR; i++) r.d[i] = g*a.d[i];
  return r;
}
inline sfloat atan2(const sfloat &y, const sfloat &x) {
  sfloat r;
  r.v = std::atan2(y.v,x.v);
  double den = x.v*x.v + y.v*y.v;
  for (int i = 0; i < SFLOAT_NDIR; i++)
    r.d[i] = (x.v*y.d[i] - y.v*x.d[i])/den;
  return r;
}
inline sfloat pow(const sfloat &a, const sfloat &b) {
  sfloat r;
  r.v = std::pow(a.v,b.v);
  // d(a^b) = a^b * (b'*ln a + b*a'/a); guard a<=0 (integer-like use)
  double la = (a.v > 0.0) ? std::log(a.v) : 0.0;
  double g  = (a.v != 0.0) ? b.v*r.v/a.v : 0.0;
  for (int i = 0; i < SFLOAT_NDIR; i++)
    r.d[i] = g*a.d[i] + r.v*la*b.d[i];
  return r;
}
inline sfloat pow(const sfloat &a, double b) { return pow(a,sfloat(b)); }
inline sfloat pow(const sfloat &a, int b)    { return pow(a,sfloat((double)b)); }
inline sfloat pow(double a, const sfloat &b) { return pow(sfloat(a),b); }

// derivative-zero a.e. functions: derivative deliberately dropped

inline sfloat floor(const sfloat &a) { return sfloat(std::floor(a.v)); }
inline sfloat ceil(const sfloat &a)  { return sfloat(std::ceil(a.v)); }
inline sfloat round(const sfloat &a) { return sfloat(std::round(a.v)); }
inline sfloat trunc(const sfloat &a) { return sfloat(std::trunc(a.v)); }
inline sfloat fmod(const sfloat &a, const sfloat &b) {
  sfloat r(a);                     // d/da fmod = 1 a.e.
  r.v = std::fmod(a.v,b.v);
  return r;
}

inline sfloat erf(const sfloat &a) {
  sfloat r;
  r.v = std::erf(a.v);
  double g = 2.0/std::sqrt(M_PI)*std::exp(-a.v*a.v);
  for (int i = 0; i < SFLOAT_NDIR; i++) r.d[i] = g*a.d[i];
  return r;
}
inline sfloat erfc(const sfloat &a) {
  sfloat r;
  r.v = std::erfc(a.v);
  double g = -2.0/std::sqrt(M_PI)*std::exp(-a.v*a.v);
  for (int i = 0; i < SFLOAT_NDIR; i++) r.d[i] = g*a.d[i];
  return r;
}

// digamma via recurrence + asymptotic series (for tgamma's derivative)

inline double sfloat_digamma(double x) {
  double r = 0.0;
  while (x < 6.0) { r -= 1.0/x; x += 1.0; }
  double f = 1.0/(x*x);
  return r + std::log(x) - 0.5/x
         - f*(1.0/12.0 - f*(1.0/120.0 - f/252.0));
}
inline sfloat tgamma(const sfloat &a) {
  sfloat r;
  r.v = std::tgamma(a.v);
  double g = r.v*sfloat_digamma(a.v);
  for (int i = 0; i < SFLOAT_NDIR; i++) r.d[i] = g*a.d[i];
  return r;
}
inline sfloat lgamma(const sfloat &a) {
  sfloat r;
  r.v = std::lgamma(a.v);
  double g = sfloat_digamma(a.v);
  for (int i = 0; i < SFLOAT_NDIR; i++) r.d[i] = g*a.d[i];
  return r;
}

inline sfloat fmax(const sfloat &a, const sfloat &b) { return a.v >= b.v ? a : b; }
inline sfloat fmin(const sfloat &a, const sfloat &b) { return a.v <= b.v ? a : b; }

inline bool isnan(const sfloat &a) { return std::isnan(a.v); }
inline bool isinf(const sfloat &a) { return std::isinf(a.v); }

inline sfloat abs(const sfloat &a) { return fabs(a); }

template <typename T> inline T spval(T x)   { return x; }    // identity
inline double spval(const sfloat &x)        { return x.v; }  // deriv ends here

#define MPI_SFLOAT MPI_5DOUBLE

#endif  /* SPARTA_AD */

#endif
