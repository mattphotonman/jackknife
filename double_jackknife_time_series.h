#ifndef INCLUDE_DOUBLE_JACKKNIFE_TIME_SERIES
#define INCLUDE_DOUBLE_JACKKNIFE_TIME_SERIES

#include <complex>
#include "jackknife.h"
#include "double_jackknife.h"
#include "jackknife_time_series.h"

class DoubleJackknifeTimeSeries
{
 private:
  int N;
  int Nt;
 protected:
  DoubleJackknife *Jk;
  int *is_defined;
 public:
  friend class DoubleJackknifeCmplxTimeSeries;
  explicit DoubleJackknifeTimeSeries(int nt=1, int n=3);
  DoubleJackknifeTimeSeries(const DoubleJackknifeTimeSeries & JT);
  DoubleJackknifeTimeSeries(int nt, const DoubleJackknife & J);
  virtual ~DoubleJackknifeTimeSeries();
  DoubleJackknifeTimeSeries & operator=(const DoubleJackknifeTimeSeries & JT);
  int ReturnIsDefined(int t) const;
  double ReturnAve(int t) const;
  double ReturnErr(int t) const;
  double ReturnJk(int t, int i) const;
  double ReturnJkErr(int t, int i) const;
  double ReturnDoubleJk(int t, int i, int j) const;
  int ReturnNt() const;
  int ReturnN() const;
  int IsComplex() const { return 0; }
  DoubleJackknifeTimeSeries & FromSingleJk(const JackknifeTimeSeries & JT);
  friend DoubleJackknifeTimeSeries FromSingleJk(const JackknifeTimeSeries & JT);
  DoubleJackknifeTimeSeries operator+(const DoubleJackknifeTimeSeries & JT) const;
  DoubleJackknifeTimeSeries & operator+=(const DoubleJackknifeTimeSeries & JT);
  DoubleJackknifeTimeSeries operator-(const DoubleJackknifeTimeSeries & JT) const;
  DoubleJackknifeTimeSeries operator-() const;
  DoubleJackknifeTimeSeries & operator-=(const DoubleJackknifeTimeSeries & JT);
  DoubleJackknifeTimeSeries operator*(const DoubleJackknifeTimeSeries & JT) const;
  DoubleJackknifeTimeSeries operator*(double d) const;
  friend DoubleJackknifeTimeSeries operator*(double d, const DoubleJackknifeTimeSeries & JT);
  DoubleJackknifeTimeSeries & operator*=(const DoubleJackknifeTimeSeries & JT);
  DoubleJackknifeTimeSeries & operator*=(double d);
  DoubleJackknifeTimeSeries operator/(const DoubleJackknifeTimeSeries & JT) const;
  DoubleJackknifeTimeSeries operator/(double d) const;
  friend DoubleJackknifeTimeSeries operator/(double d, const DoubleJackknifeTimeSeries & JT);
  DoubleJackknifeTimeSeries & operator/=(const DoubleJackknifeTimeSeries & JT);
  DoubleJackknifeTimeSeries & operator/=(double d);
  DoubleJackknifeTimeSeries & Reverse();
  friend DoubleJackknifeTimeSeries Reverse(const DoubleJackknifeTimeSeries & JT);
  DoubleJackknifeTimeSeries & Translate(int t);
  friend DoubleJackknifeTimeSeries Translate(int t, const DoubleJackknifeTimeSeries & JT);
  DoubleJackknifeTimeSeries & TranslateReverse(int t);
  friend DoubleJackknifeTimeSeries TranslateReverse(int t, const DoubleJackknifeTimeSeries & JT);
  DoubleJackknifeTimeSeries & M_eff();
  friend DoubleJackknifeTimeSeries M_eff(const DoubleJackknifeTimeSeries & JT);
  DoubleJackknifeTimeSeries & Abs();
  friend DoubleJackknifeTimeSeries Abs(const DoubleJackknifeTimeSeries & JT);
  friend DoubleJackknifeTimeSeries Abs(const DoubleJackknifeCmplxTimeSeries & JT);
  friend DoubleJackknifeTimeSeries Real(const DoubleJackknifeCmplxTimeSeries & JT);
  friend DoubleJackknifeTimeSeries Imag(const DoubleJackknifeCmplxTimeSeries & JT);
  double Normalize();
  friend DoubleJackknifeTimeSeries Combine(const DoubleJackknifeTimeSeries & JT1, const DoubleJackknifeTimeSeries & JT2);
  //friend DoubleJackknifeTimeSeries ReadWMEData(const WMEDataReadIn & D);
  friend void OutputJk(string filename, const DoubleJackknifeTimeSeries & JT);
};

class DoubleJackknifeCmplxTimeSeries
{
 private:
  int N;
  int Nt;
 protected:
  DoubleJackknifeCmplx *Jk;
  int *is_defined;
 public:
  friend class DoubleJackknifeTimeSeries;
  explicit DoubleJackknifeCmplxTimeSeries(int nt=1, int n=3);
  DoubleJackknifeCmplxTimeSeries(const DoubleJackknifeCmplxTimeSeries & JT);
  DoubleJackknifeCmplxTimeSeries(int nt, const DoubleJackknifeCmplx & J);
  virtual ~DoubleJackknifeCmplxTimeSeries();
  DoubleJackknifeCmplxTimeSeries & operator=(const DoubleJackknifeCmplxTimeSeries & JT);
  int ReturnIsDefined(int t) const;
  complex<double> ReturnAve(int t) const;
  double ReturnErr(int t) const;
  complex<double> ReturnJk(int t, int i) const;
  double ReturnJkErr(int t, int i) const;
  complex<double> ReturnDoubleJk(int t, int i, int j) const;
  int ReturnNt() const;
  int ReturnN() const;
  int IsComplex() const { return 1; }
  DoubleJackknifeCmplxTimeSeries & FromSingleJk(const JackknifeCmplxTimeSeries & JT);
  friend DoubleJackknifeCmplxTimeSeries FromSingleJk(const JackknifeCmplxTimeSeries & JT);
  DoubleJackknifeCmplxTimeSeries operator+(const DoubleJackknifeCmplxTimeSeries & JT) const;
  DoubleJackknifeCmplxTimeSeries & operator+=(const DoubleJackknifeCmplxTimeSeries & JT);
  DoubleJackknifeCmplxTimeSeries operator-(const DoubleJackknifeCmplxTimeSeries & JT) const;
  DoubleJackknifeCmplxTimeSeries operator-() const;
  DoubleJackknifeCmplxTimeSeries & operator-=(const DoubleJackknifeCmplxTimeSeries & JT);
  DoubleJackknifeCmplxTimeSeries operator*(const DoubleJackknifeCmplxTimeSeries & JT) const;
  DoubleJackknifeCmplxTimeSeries operator*(complex<double> c) const;
  friend DoubleJackknifeCmplxTimeSeries operator*(complex<double> c, const DoubleJackknifeCmplxTimeSeries & JT);
  DoubleJackknifeCmplxTimeSeries & operator*=(const DoubleJackknifeCmplxTimeSeries & JT);
  DoubleJackknifeCmplxTimeSeries & operator*=(complex<double> c);
  DoubleJackknifeCmplxTimeSeries operator/(const DoubleJackknifeCmplxTimeSeries & JT) const;
  DoubleJackknifeCmplxTimeSeries operator/(complex<double> c) const;
  friend DoubleJackknifeCmplxTimeSeries operator/(complex<double> c, const DoubleJackknifeCmplxTimeSeries & JT);
  DoubleJackknifeCmplxTimeSeries & operator/=(const DoubleJackknifeCmplxTimeSeries & JT);
  DoubleJackknifeCmplxTimeSeries & operator/=(complex<double> c);
  DoubleJackknifeCmplxTimeSeries & Reverse();
  friend DoubleJackknifeCmplxTimeSeries Reverse(const DoubleJackknifeCmplxTimeSeries & JT);
  DoubleJackknifeCmplxTimeSeries & Translate(int t);
  friend DoubleJackknifeCmplxTimeSeries Translate(int t, const DoubleJackknifeCmplxTimeSeries & JT);
  DoubleJackknifeCmplxTimeSeries & TranslateReverse(int t);
  friend DoubleJackknifeCmplxTimeSeries TranslateReverse(int t, const DoubleJackknifeCmplxTimeSeries & JT);
  DoubleJackknifeCmplxTimeSeries & Abs();
  friend DoubleJackknifeTimeSeries Abs(const DoubleJackknifeCmplxTimeSeries & JT);
  DoubleJackknifeCmplxTimeSeries & Real();
  friend DoubleJackknifeTimeSeries Real(const DoubleJackknifeCmplxTimeSeries & JT);
  DoubleJackknifeCmplxTimeSeries & Imag();
  friend DoubleJackknifeTimeSeries Imag(const DoubleJackknifeCmplxTimeSeries & JT);
  DoubleJackknifeCmplxTimeSeries & Conj();
  friend DoubleJackknifeCmplxTimeSeries Conj(const DoubleJackknifeCmplxTimeSeries & JT);
  double Normalize();
  friend DoubleJackknifeCmplxTimeSeries Combine(const DoubleJackknifeCmplxTimeSeries & JT1, const DoubleJackknifeCmplxTimeSeries & JT2);
  //friend DoubleJackknifeCmplxTimeSeries ReadWMECmplxData(const WMECmplxDataReadIn & D);
  friend void OutputJk(string filename, const DoubleJackknifeCmplxTimeSeries & JT);
};

void OutputAve(string filename, const DoubleJackknifeTimeSeries & JT);
void OutputAve(string filename, const DoubleJackknifeCmplxTimeSeries & JT);



#endif
