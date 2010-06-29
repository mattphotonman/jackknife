#ifndef INCLUDE_JACKKNIFE_TIME_SERIES
#define INCLUDE_JACKKNIFE_TIME_SERIES

#include <complex>
#include "jackknife.h"
#include "data_read_in.h"

class DoubleJackknifeTimeSeries;
class DoubleJackknifeCmplxTimeSeries;

class JackknifeTimeSeries
{
 private:
  int N;
  int Nt;
 protected:
  Jackknife *Jk;
  int *is_defined;
 public:
  friend class JackknifeCmplxTimeSeries;
  friend class DoubleJackknifeTimeSeries;
  explicit JackknifeTimeSeries(int nt=1, int n=2);
  JackknifeTimeSeries(const JackknifeTimeSeries & JT);
  JackknifeTimeSeries(int nt, const Jackknife & J);
  virtual ~JackknifeTimeSeries();
  JackknifeTimeSeries & operator=(const JackknifeTimeSeries & JT);
  int ReturnIsDefined(int t) const;
  double ReturnAve(int t) const;
  double ReturnErr(int t) const;
  double ReturnJk(int t, int i) const;
  int ReturnNt() const;
  int ReturnN() const;
  int IsComplex() const { return 0; }
  JackknifeTimeSeries operator+(const JackknifeTimeSeries & JT) const;
  JackknifeTimeSeries & operator+=(const JackknifeTimeSeries & JT);
  JackknifeTimeSeries operator-(const JackknifeTimeSeries & JT) const;
  JackknifeTimeSeries operator-() const;
  JackknifeTimeSeries & operator-=(const JackknifeTimeSeries & JT);
  JackknifeTimeSeries operator*(const JackknifeTimeSeries & JT) const;
  JackknifeTimeSeries operator*(double d) const;
  friend JackknifeTimeSeries operator*(double d, const JackknifeTimeSeries & JT);
  JackknifeTimeSeries & operator*=(const JackknifeTimeSeries & JT);
  JackknifeTimeSeries & operator*=(double d);
  JackknifeTimeSeries operator/(const JackknifeTimeSeries & JT) const;
  JackknifeTimeSeries operator/(double d) const;
  friend JackknifeTimeSeries operator/(double d, const JackknifeTimeSeries & JT);
  JackknifeTimeSeries & operator/=(const JackknifeTimeSeries & JT);
  JackknifeTimeSeries & operator/=(double d);
  JackknifeTimeSeries & Reverse();
  friend JackknifeTimeSeries Reverse(const JackknifeTimeSeries & JT);
  JackknifeTimeSeries & Translate(int t);
  friend JackknifeTimeSeries Translate(int t, const JackknifeTimeSeries & JT);
  JackknifeTimeSeries & TranslateReverse(int t);
  friend JackknifeTimeSeries TranslateReverse(int t, const JackknifeTimeSeries & JT);
  JackknifeTimeSeries & M_eff();
  friend JackknifeTimeSeries M_eff(const JackknifeTimeSeries & JT);
  JackknifeTimeSeries & Abs();
  friend JackknifeTimeSeries Abs(const JackknifeTimeSeries & JT);
  friend JackknifeTimeSeries Abs(const JackknifeCmplxTimeSeries & JT);
  friend JackknifeTimeSeries Real(const JackknifeCmplxTimeSeries & JT);
  friend JackknifeTimeSeries Imag(const JackknifeCmplxTimeSeries & JT);
  double Normalize();
  friend JackknifeTimeSeries Combine(const JackknifeTimeSeries & JT1, const JackknifeTimeSeries & JT2);
  friend JackknifeTimeSeries ReadWMEData(const WMEDataReadIn & D);
  friend void OutputJk(string filename, const JackknifeTimeSeries & JT);
  friend JackknifeTimeSeries ReadTimeSeriesFile(string ave_file, string jks_file, int Nt, int N);
};

class JackknifeCmplxTimeSeries
{
 private:
  int N;
  int Nt;
 protected:
  JackknifeCmplx *Jk;
  int *is_defined;
 public:
  friend class JackknifeTimeSeries;
  friend class DoubleJackknifeCmplxTimeSeries;
  explicit JackknifeCmplxTimeSeries(int nt=1, int n=2);
  JackknifeCmplxTimeSeries(const JackknifeCmplxTimeSeries & JT);
  JackknifeCmplxTimeSeries(int nt, const JackknifeCmplx & J);
  virtual ~JackknifeCmplxTimeSeries();
  JackknifeCmplxTimeSeries & operator=(const JackknifeCmplxTimeSeries & JT);
  int ReturnIsDefined(int t) const;
  complex<double> ReturnAve(int t) const;
  double ReturnErr(int t) const;
  complex<double> ReturnJk(int t, int i) const;
  int ReturnNt() const;
  int ReturnN() const;
  int IsComplex() const { return 1; }
  JackknifeCmplxTimeSeries operator+(const JackknifeCmplxTimeSeries & JT) const;
  JackknifeCmplxTimeSeries & operator+=(const JackknifeCmplxTimeSeries & JT);
  JackknifeCmplxTimeSeries operator-(const JackknifeCmplxTimeSeries & JT) const;
  JackknifeCmplxTimeSeries operator-() const;
  JackknifeCmplxTimeSeries & operator-=(const JackknifeCmplxTimeSeries & JT);
  JackknifeCmplxTimeSeries operator*(const JackknifeCmplxTimeSeries & JT) const;
  JackknifeCmplxTimeSeries operator*(complex<double> c) const;
  friend JackknifeCmplxTimeSeries operator*(complex<double> c, const JackknifeCmplxTimeSeries & JT);
  JackknifeCmplxTimeSeries & operator*=(const JackknifeCmplxTimeSeries & JT);
  JackknifeCmplxTimeSeries & operator*=(complex<double> c);
  JackknifeCmplxTimeSeries operator/(const JackknifeCmplxTimeSeries & JT) const;
  JackknifeCmplxTimeSeries operator/(complex<double> c) const;
  friend JackknifeCmplxTimeSeries operator/(complex<double> c, const JackknifeCmplxTimeSeries & JT);
  JackknifeCmplxTimeSeries & operator/=(const JackknifeCmplxTimeSeries & JT);
  JackknifeCmplxTimeSeries & operator/=(complex<double> c);
  JackknifeCmplxTimeSeries & Reverse();
  friend JackknifeCmplxTimeSeries Reverse(const JackknifeCmplxTimeSeries & JT);
  JackknifeCmplxTimeSeries & Translate(int t);
  friend JackknifeCmplxTimeSeries Translate(int t, const JackknifeCmplxTimeSeries & JT);
  JackknifeCmplxTimeSeries & TranslateReverse(int t);
  friend JackknifeCmplxTimeSeries TranslateReverse(int t, const JackknifeCmplxTimeSeries & JT);
  JackknifeCmplxTimeSeries & Abs();
  friend JackknifeTimeSeries Abs(const JackknifeCmplxTimeSeries & JT);
  JackknifeCmplxTimeSeries & Real();
  friend JackknifeTimeSeries Real(const JackknifeCmplxTimeSeries & JT);
  JackknifeCmplxTimeSeries & Imag();
  friend JackknifeTimeSeries Imag(const JackknifeCmplxTimeSeries & JT);
  JackknifeCmplxTimeSeries & Conj();
  friend JackknifeCmplxTimeSeries Conj(const JackknifeCmplxTimeSeries & JT);
  double Normalize();
  friend JackknifeCmplxTimeSeries Combine(const JackknifeCmplxTimeSeries & JT1, const JackknifeCmplxTimeSeries & JT2);
  friend JackknifeCmplxTimeSeries ReadWMECmplxData(const WMECmplxDataReadIn & D);
  friend void OutputJk(string filename, const JackknifeCmplxTimeSeries & JT);
};

void OutputAve(string filename, const JackknifeTimeSeries & JT);
void OutputAve(string filename, const JackknifeCmplxTimeSeries & JT);



#endif
