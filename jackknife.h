#ifndef INCLUDE_JACKKNIFE
#define INCLUDE_JACKKNIFE

#include <complex>

class JackknifeCmplx;
class JackknifeTimeSeries;
class JackknifeCmplxTimeSeries;
class WMEDataReadIn;
class WMECmplxDataReadIn;
class DoubleJackknife;
class DoubleJackknifeCmplx;

class Jackknife
{
 private:
  double *jk;
  int N;
  double ave;
  double err;
 public:
  friend class JackknifeCmplx;
  friend class DoubleJackknife;
  explicit Jackknife(int n=2);
  Jackknife(const Jackknife & J);
  Jackknife(int n, double d);
  virtual ~Jackknife();
  Jackknife & operator=(const Jackknife & J);
  Jackknife & CalcErr();
  Jackknife & CalcAll();
  double ReturnAve() const;
  double ReturnErr() const;
  double ReturnJk(int i) const;
  int ReturnN() const;
  int IsComplex() const { return 0; }
  double Unjackknife(int i) const;
  friend DoubleJackknife FromSingleJk(const Jackknife & J);
  Jackknife operator+(const Jackknife & J) const;
  Jackknife & operator+=(const Jackknife & J);
  Jackknife operator-(const Jackknife & J) const;
  Jackknife operator-() const;
  Jackknife & operator-=(const Jackknife & J);
  Jackknife operator*(const Jackknife & J) const;
  Jackknife operator*(double d) const;
  friend Jackknife operator*(double d, const Jackknife & J);
  Jackknife & operator*=(const Jackknife & J);
  Jackknife & operator*=(double d);
  Jackknife operator/(const Jackknife & J) const;
  Jackknife operator/(double d) const;
  friend Jackknife operator/(double d, const Jackknife & J);
  Jackknife & operator/=(const Jackknife & J);
  Jackknife & operator/=(double d);
  Jackknife & Log();
  friend Jackknife Log(const Jackknife & J);
  Jackknife & Sqrt();
  friend Jackknife Sqrt(const Jackknife & J);
  Jackknife & Exp();
  friend Jackknife Exp(const Jackknife & J);
  Jackknife & Abs();
  friend Jackknife Abs(const Jackknife & J);
  friend Jackknife Abs(const JackknifeCmplx & J);
  friend Jackknife Real(const JackknifeCmplx & J);
  friend Jackknife Imag(const JackknifeCmplx & J);
  friend Jackknife Combine(const Jackknife & J1, const Jackknife & J2);
  friend JackknifeTimeSeries ReadWMEData(const WMEDataReadIn & D);
  friend void OutputJk(string filename, const JackknifeTimeSeries & JT);
  friend DoubleJackknife Combine(const DoubleJackknife & J1, const DoubleJackknife & J2);
  friend JackknifeTimeSeries ReadTimeSeriesFile(string ave_file, string jks_file, int Nt, int N);
  friend Jackknife ReadTextFile(string file, int N);
};

class JackknifeCmplx
{
 private:
  complex<double> *jk;
  int N;
  complex<double> ave;
  double err;
 public:
  friend class Jackknife;
  friend class DoubleJackknifeCmplx;
  explicit JackknifeCmplx(int n=2);
  JackknifeCmplx(const JackknifeCmplx & J);
  JackknifeCmplx(int n, complex<double> c);
  virtual ~JackknifeCmplx();
  JackknifeCmplx & operator=(const JackknifeCmplx & J);
  JackknifeCmplx & CalcErr();
  JackknifeCmplx & CalcAll();
  complex<double> ReturnAve() const;
  double ReturnErr() const;
  complex<double> ReturnJk(int i) const;
  int ReturnN() const;
  int IsComplex() const { return 1; }
  complex<double> Unjackknife(int i) const;
  friend DoubleJackknifeCmplx FromSingleJk(const JackknifeCmplx & J);
  JackknifeCmplx operator+(const JackknifeCmplx & J) const;
  JackknifeCmplx & operator+=(const JackknifeCmplx & J);
  JackknifeCmplx operator-(const JackknifeCmplx & J) const;
  JackknifeCmplx operator-() const;
  JackknifeCmplx & operator-=(const JackknifeCmplx & J);
  JackknifeCmplx operator*(const JackknifeCmplx & J) const;
  JackknifeCmplx operator*(complex<double> c) const;
  friend JackknifeCmplx operator*(complex<double> c, const JackknifeCmplx & J);
  JackknifeCmplx & operator*=(const JackknifeCmplx & J);
  JackknifeCmplx & operator*=(complex<double> c);
  JackknifeCmplx operator/(const JackknifeCmplx & J) const;
  JackknifeCmplx operator/(complex<double> c) const;
  friend JackknifeCmplx operator/(complex<double> c, const JackknifeCmplx & J);
  JackknifeCmplx & operator/=(const JackknifeCmplx & J);
  JackknifeCmplx & operator/=(complex<double> c);
  JackknifeCmplx & Abs();
  friend Jackknife Abs(const JackknifeCmplx & J);
  JackknifeCmplx & Real();
  friend Jackknife Real(const JackknifeCmplx & J);
  JackknifeCmplx & Imag();
  friend Jackknife Imag(const JackknifeCmplx & J);
  JackknifeCmplx & Conj();
  friend JackknifeCmplx Conj(const JackknifeCmplx & J);
  friend JackknifeCmplx Combine(const JackknifeCmplx & J1, const JackknifeCmplx & J2);
  friend JackknifeCmplxTimeSeries ReadWMECmplxData(const WMECmplxDataReadIn & D);
  friend void OutputJk(string filename, const JackknifeCmplxTimeSeries & JT);
  friend DoubleJackknifeCmplx Combine(const DoubleJackknifeCmplx & J1, const DoubleJackknifeCmplx & J2);
};

#endif
