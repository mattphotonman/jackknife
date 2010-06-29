#ifndef INCLUDE_DOUBLEJACKKNIFE
#define INCLUDE_DOUBLEJACKKNIFE

#include <complex>
#include "jackknife.h"

class DoubleJackknifeCmplx;
//class DoubleJackknifeTimeSeries;
//class DoubleJackknifeCmplxTimeSeries;
class WMEDataReadIn;
class WMECmplxDataReadIn;

class DoubleJackknife
{
 private:
  Jackknife *jk;
  int N;
  double ave;
  double err;
 public:
  friend class DoubleJackknifeCmplx;
  explicit DoubleJackknife(int n=3);
  DoubleJackknife(const DoubleJackknife & J);
  DoubleJackknife(int n, double d);
  virtual ~DoubleJackknife();
  DoubleJackknife & operator=(const DoubleJackknife & J);
  DoubleJackknife & CalcErr();
  DoubleJackknife & CalcAll();
  double ReturnAve() const;
  double ReturnErr() const;
  double ReturnJk(int i) const;
  double ReturnJkErr(int i) const;
  double ReturnDoubleJk(int i, int j) const;
  int ReturnN() const;
  int IsComplex() const { return 0; }
  double Unjackknife(int i) const;
  DoubleJackknife & FromSingleJk(const Jackknife & J);
  friend DoubleJackknife FromSingleJk(const Jackknife & J);
  DoubleJackknife operator+(const DoubleJackknife & J) const;
  DoubleJackknife & operator+=(const DoubleJackknife & J);
  DoubleJackknife operator-(const DoubleJackknife & J) const;
  DoubleJackknife operator-() const;
  DoubleJackknife & operator-=(const DoubleJackknife & J);
  DoubleJackknife operator*(const DoubleJackknife & J) const;
  DoubleJackknife operator*(double d) const;
  friend DoubleJackknife operator*(double d, const DoubleJackknife & J);
  DoubleJackknife & operator*=(const DoubleJackknife & J);
  DoubleJackknife & operator*=(double d);
  DoubleJackknife operator/(const DoubleJackknife & J) const;
  DoubleJackknife operator/(double d) const;
  friend DoubleJackknife operator/(double d, const DoubleJackknife & J);
  DoubleJackknife & operator/=(const DoubleJackknife & J);
  DoubleJackknife & operator/=(double d);
  DoubleJackknife & Log();
  friend DoubleJackknife Log(const DoubleJackknife & J);
  DoubleJackknife & Abs();
  friend DoubleJackknife Abs(const DoubleJackknife & J);
  friend DoubleJackknife Abs(const DoubleJackknifeCmplx & J);
  friend DoubleJackknife Real(const DoubleJackknifeCmplx & J);
  friend DoubleJackknife Imag(const DoubleJackknifeCmplx & J);
  friend DoubleJackknife Combine(const DoubleJackknife & J1, const DoubleJackknife & J2);
  //friend DoubleJackknifeTimeSeries ReadWMEData(const WMEDataReadIn & D);
  //friend void OutputJk(string filename, const DoubleJackknifeTimeSeries & JT);
};

class DoubleJackknifeCmplx
{
 private:
  JackknifeCmplx *jk;
  int N;
  complex<double> ave;
  double err;
 public:
  friend class DoubleJackknife;
  explicit DoubleJackknifeCmplx(int n=3);
  DoubleJackknifeCmplx(const DoubleJackknifeCmplx & J);
  DoubleJackknifeCmplx(int n, complex<double> c);
  virtual ~DoubleJackknifeCmplx();
  DoubleJackknifeCmplx & operator=(const DoubleJackknifeCmplx & J);
  DoubleJackknifeCmplx & CalcErr();
  DoubleJackknifeCmplx & CalcAll();
  complex<double> ReturnAve() const;
  double ReturnErr() const;
  complex<double> ReturnJk(int i) const;
  double ReturnJkErr(int i) const;
  complex<double> ReturnDoubleJk(int i, int j) const;
  int ReturnN() const;
  int IsComplex() const { return 1; }
  complex<double> Unjackknife(int i) const;
  DoubleJackknifeCmplx & FromSingleJk(const JackknifeCmplx & J);
  friend DoubleJackknifeCmplx FromSingleJk(const JackknifeCmplx & J);
  DoubleJackknifeCmplx operator+(const DoubleJackknifeCmplx & J) const;
  DoubleJackknifeCmplx & operator+=(const DoubleJackknifeCmplx & J);
  DoubleJackknifeCmplx operator-(const DoubleJackknifeCmplx & J) const;
  DoubleJackknifeCmplx operator-() const;
  DoubleJackknifeCmplx & operator-=(const DoubleJackknifeCmplx & J);
  DoubleJackknifeCmplx operator*(const DoubleJackknifeCmplx & J) const;
  DoubleJackknifeCmplx operator*(complex<double> c) const;
  friend DoubleJackknifeCmplx operator*(complex<double> c, const DoubleJackknifeCmplx & J);
  DoubleJackknifeCmplx & operator*=(const DoubleJackknifeCmplx & J);
  DoubleJackknifeCmplx & operator*=(complex<double> c);
  DoubleJackknifeCmplx operator/(const DoubleJackknifeCmplx & J) const;
  DoubleJackknifeCmplx operator/(complex<double> c) const;
  friend DoubleJackknifeCmplx operator/(complex<double> c, const DoubleJackknifeCmplx & J);
  DoubleJackknifeCmplx & operator/=(const DoubleJackknifeCmplx & J);
  DoubleJackknifeCmplx & operator/=(complex<double> c);
  DoubleJackknifeCmplx & Abs();
  friend DoubleJackknife Abs(const DoubleJackknifeCmplx & J);
  DoubleJackknifeCmplx & Real();
  friend DoubleJackknife Real(const DoubleJackknifeCmplx & J);
  DoubleJackknifeCmplx & Imag();
  friend DoubleJackknife Imag(const DoubleJackknifeCmplx & J);
  DoubleJackknifeCmplx & Conj();
  friend DoubleJackknifeCmplx Conj(const DoubleJackknifeCmplx & J);
  friend DoubleJackknifeCmplx Combine(const DoubleJackknifeCmplx & J1, const DoubleJackknifeCmplx & J2);
  //friend DoubleJackknifeCmplxTimeSeries ReadWMECmplxData(const WMECmplxDataReadIn & D);
  //friend void OutputJk(string filename, const DoubleJackknifeCmplxTimeSeries & JT);
};

#endif
