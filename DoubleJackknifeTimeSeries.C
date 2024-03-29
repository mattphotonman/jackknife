#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
using namespace std;
#include "jackknife.h"
#include "double_jackknife.h"
#include "double_jackknife_time_series.h"


//DoubleJackknifeTimeSeries class

//Standard and Default Constructor
DoubleJackknifeTimeSeries::DoubleJackknifeTimeSeries(int nt, int n)
{
  if (n<3) {
    //Error: I think I'd rather have it abort here.
    cout << "N must be >= 3 for a DoubleJackknifeTimeSeries object, but was given N=" << n << ".  Setting N=3.\n";
    N=2;
  } else {
    N=n;
  }
  if (nt<1) {
    //Error: I think I'd rather have it abort here.
    cout << "Nt must be >= 1 for a DoubleJackknifeTimeSeries object, but was given Nt=" << nt << ".  Setting Nt=1.\n";
    Nt=1;
  } else {
    Nt=nt;
  }
  Jk=new DoubleJackknife [Nt+1];
  is_defined=new int [Nt+1];
  for (int t=0; t<=Nt; t++) {
    is_defined[t]=0;
    //Have to initialize a DoubleJackknife object of size N for each Jk[t]
    DoubleJackknife tmp(N);
    Jk[t]=tmp;
  }
}

//Copy Constructor
DoubleJackknifeTimeSeries::DoubleJackknifeTimeSeries(const DoubleJackknifeTimeSeries & JT)
{
  N=JT.N;
  Nt=JT.Nt;
  Jk=new DoubleJackknife [Nt+1];
  is_defined=new int [Nt+1];
  for (int t=0; t<=Nt; t++) {
    Jk[t]=JT.Jk[t];
    is_defined[t]=JT.is_defined[t];
  }
}

//Constructor with a DoubleJackknife object
DoubleJackknifeTimeSeries::DoubleJackknifeTimeSeries(int nt, const DoubleJackknife & J)
{
  if (nt<1) {
    //Error: I think I'd rather have it abort here.
    cout << "Nt must be >= 1 for a DoubleJackknifeTimeSeries object, but was given Nt=" << nt << ".  Setting Nt=1.\n";
    Nt=1;
  } else {
    Nt=nt;
  }
  N=J.ReturnN();
  Jk=new DoubleJackknife [Nt+1];
  is_defined=new int [Nt+1];
  for (int t=0; t<=Nt; t++) {
    Jk[t]=J;
    is_defined[t]=1;
  }
}

//Destructor
DoubleJackknifeTimeSeries::~DoubleJackknifeTimeSeries() 
{
  delete [] Jk;
  delete [] is_defined;
}

//Assignment operator
//Deep copy
DoubleJackknifeTimeSeries & DoubleJackknifeTimeSeries::operator=(const DoubleJackknifeTimeSeries & JT)
{
  if (this == &JT)
    return *this;
  delete [] Jk;
  delete [] is_defined;
  N=JT.N;
  Nt=JT.Nt;
  Jk=new DoubleJackknife [Nt+1];
  is_defined=new int [Nt+1];
  for (int t=0; t<=Nt; t++) {
    Jk[t]=JT.Jk[t];
    is_defined[t]=JT.is_defined[t];
  }
  return *this;
}

//Say whether or not there is data at t
int DoubleJackknifeTimeSeries::ReturnIsDefined(int t) const
{
  if (t<0 || t>Nt)
    return 0;
  else
    return is_defined[t];
}

//Return the average at t
double DoubleJackknifeTimeSeries::ReturnAve(int t) const
{
  if (t<0 || t>Nt || !is_defined[t]) {
    //Error:  I think I'd rather have it abort here.
    cout << "No data at t=" << t << " in DoubleJackknifeTimeSeries.\n";
    return 0.0;
  } else {
    return Jk[t].ReturnAve();
  }
}

//Return the error at t
double DoubleJackknifeTimeSeries::ReturnErr(int t) const
{
  if (t<0 || t>Nt || !is_defined[t]) {
    //Error:  I think I'd rather have it abort here.
    cout << "No data at t=" << t << " in DoubleJackknifeTimeSeries.\n";
    return 0.0;
  } else {
    return Jk[t].ReturnErr();
  }
}

//Return the ith jackknife value at t
double DoubleJackknifeTimeSeries::ReturnJk(int t, int i) const
{
  if (t<0 || t>Nt || !is_defined[t]) {
    //Error:  I think I'd rather have it abort here.
    cout << "No data at t=" << t << " in DoubleJackknifeTimeSeries.\n";
    return 0.0;
  }
  //Use DoubleJackknife's own bounds checking for i.
  return Jk[t].ReturnJk(i);
}

//Return the error on the ith jackknife value at t
double DoubleJackknifeTimeSeries::ReturnJkErr(int t, int i) const
{
  if (t<0 || t>Nt || !is_defined[t]) {
    //Error:  I think I'd rather have it abort here.
    cout << "No data at t=" << t << " in DoubleJackknifeTimeSeries.\n";
    return 0.0;
  }
  //Use DoubleJackknife's own bounds checking for i.
  return Jk[t].ReturnJkErr(i);
}

//Return double jackknife value i, j at t
double DoubleJackknifeTimeSeries::ReturnDoubleJk(int t, int i, int j) const
{
  if (t<0 || t>Nt || !is_defined[t]) {
    //Error:  I think I'd rather have it abort here.
    cout << "No data at t=" << t << " in DoubleJackknifeTimeSeries.\n";
    return 0.0;
  }
  //Use DoubleJackknife's own bounds checking for i, j.
  return Jk[t].ReturnDoubleJk(i,j);
}

//Return Nt
int DoubleJackknifeTimeSeries::ReturnNt() const
{
  return Nt;
}

//Return N
int DoubleJackknifeTimeSeries::ReturnN() const
{
  return N;
}

//Takes a JackknifeTimeSeries object, undoes it, then calculates 
//double jackknife values to put into the current DoubleJackknifeTimeSeries
//object.
//NOTE: ONLY WORKS for quantities that depend linearly on original data.
DoubleJackknifeTimeSeries & DoubleJackknifeTimeSeries::FromSingleJk(const JackknifeTimeSeries & JT)
{
  delete [] Jk;
  delete [] is_defined;
  N=JT.N;
  Nt=JT.Nt;
  Jk=new DoubleJackknife[Nt+1];
  is_defined=new int [Nt+1];
  for (int t=0; t<=Nt; t++) {
    Jk[t].FromSingleJk(JT.Jk[t]);
    if (JT.ReturnIsDefined(t)) {
      is_defined[t]=1;
    } else {
      is_defined[t]=0;
    }
  }
  return *this;
}

//Friend function version of the above.
DoubleJackknifeTimeSeries FromSingleJk(const JackknifeTimeSeries & JT)
{
  DoubleJackknifeTimeSeries tmp;
  tmp.FromSingleJk(JT);
  return tmp;
}

//Addition
DoubleJackknifeTimeSeries DoubleJackknifeTimeSeries::operator+(const DoubleJackknifeTimeSeries & JT) const
{
  DoubleJackknifeTimeSeries sum(*this);
  sum+=JT;
  return sum;
}

//+=
DoubleJackknifeTimeSeries & DoubleJackknifeTimeSeries::operator+=(const DoubleJackknifeTimeSeries & JT)
{
  if (N!=JT.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't add DoubleJackknifeTimeSeries objects with different number of jackknife values.  Object 1 N=" << N << " Object 2 N=" << JT.N << "\n";
  } else if (Nt!=JT.Nt) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't add DoubleJackknifeTimeSeries objects with different Nt.  Object 1 Nt=" << Nt << " Object 2 Nt=" << JT.Nt << "\n";
  } else {
    for (int t=0; t<=Nt; t++) {
      if (is_defined[t] && JT.is_defined[t])
	Jk[t]+=JT.Jk[t];
      else
	is_defined[t]=0;
    }
  }
  return *this;
}

//Subtraction
DoubleJackknifeTimeSeries DoubleJackknifeTimeSeries::operator-(const DoubleJackknifeTimeSeries & JT) const
{
  DoubleJackknifeTimeSeries difference(*this);
  difference-=JT;
  return difference;
}

//Negative (Unary operator)
DoubleJackknifeTimeSeries DoubleJackknifeTimeSeries::operator-() const
{
  DoubleJackknife tmp(N,0.0);
  DoubleJackknifeTimeSeries neg(Nt,tmp);
  neg-=*this;
  return neg;
}

//-=
DoubleJackknifeTimeSeries & DoubleJackknifeTimeSeries::operator-=(const DoubleJackknifeTimeSeries & JT)
{
  if (N!=JT.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't subtract DoubleJackknifeTimeSeries objects with different number of jackknife values.  Object 1 N=" << N << " Object 2 N=" << JT.N << "\n";
  } else if (Nt!=JT.Nt) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't subtract DoubleJackknifeTimeSeries objects with different Nt.  Object 1 Nt=" << Nt << " Object 2 Nt=" << JT.Nt << "\n";
  } else {
    for (int t=0; t<=Nt; t++) {
      if (is_defined[t] && JT.is_defined[t])
	Jk[t]-=JT.Jk[t];
      else
	is_defined[t]=0;
    }
  }
  return *this;
}

//Multiplication of two DoubleJackknifeTimeSeries objects
DoubleJackknifeTimeSeries DoubleJackknifeTimeSeries::operator*(const DoubleJackknifeTimeSeries & JT) const
{
  DoubleJackknifeTimeSeries prod(*this);
  prod*=JT;
  return prod;
}

//Multiplication of a DoubleJackknifeTimeSeries object by a scalar on the right
DoubleJackknifeTimeSeries DoubleJackknifeTimeSeries::operator*(double d) const
{
  DoubleJackknifeTimeSeries prod(*this);
  prod*=d;
  return prod;
}

//Multiplication by a scalar on the left (friend function)
DoubleJackknifeTimeSeries operator*(double d, const DoubleJackknifeTimeSeries & JT)
{
  return JT*d;
}

//*=DoubleJackknifeTimeSeries object
DoubleJackknifeTimeSeries & DoubleJackknifeTimeSeries::operator*=(const DoubleJackknifeTimeSeries & JT)
{
  if (N!=JT.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't multiply DoubleJackknifeTimeSeries objects with different number of jackknife values.  Object 1 N=" << N << " Object 2 N=" << JT.N << "\n";
  } else if (Nt!=JT.Nt) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't multiply DoubleJackknifeTimeSeries objects with different Nt.  Object 1 Nt=" << Nt << " Object 2 Nt=" << JT.Nt << "\n";
  } else {
    for (int t=0; t<=Nt; t++) {
      if (is_defined[t] && JT.is_defined[t])
	Jk[t]*=JT.Jk[t];
      else
	is_defined[t]=0;
    }
  }
  return *this;
}

//*=scalar
DoubleJackknifeTimeSeries & DoubleJackknifeTimeSeries::operator*=(double d)
{
  for (int t=0; t<=Nt; t++)
    if (is_defined[t])
      Jk[t]*=d;
  return *this;
}

//Division of two JackknifeTimeSeries objects
DoubleJackknifeTimeSeries DoubleJackknifeTimeSeries::operator/(const DoubleJackknifeTimeSeries & JT) const
{
  DoubleJackknifeTimeSeries q(*this);
  q/=JT;
  return q;
}

//Division by a scalar
DoubleJackknifeTimeSeries DoubleJackknifeTimeSeries::operator/(double d) const
{
  DoubleJackknifeTimeSeries q(*this);
  q/=d;
  return q;
}

//Division of a scalar by a DoubleJackknifeTimeSeries object (friend function)
DoubleJackknifeTimeSeries operator/(double d, const DoubleJackknifeTimeSeries & JT)
{
  int Nt=JT.Nt;
  int N=JT.N;
  DoubleJackknife tmp(N,d);
  DoubleJackknifeTimeSeries q(Nt,tmp);
  q/=JT;
  return q;
}

//slash equals DoubleJackknifeTimeSeries object
DoubleJackknifeTimeSeries & DoubleJackknifeTimeSeries::operator/=(const DoubleJackknifeTimeSeries & JT)
{
  if (N!=JT.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't divide DoubleJackknifeTimeSeries objects with different number of jackknife values.  Object 1 N=" << N << " Object 2 N=" << JT.N << "\n";
  } else if (Nt!=JT.Nt) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't divide DoubleJackknifeTimeSeries objects with different Nt.  Object 1 Nt=" << Nt << " Object 2 Nt=" << JT.Nt << "\n";
  } else {
    for (int t=0; t<=Nt; t++) {
      if (is_defined[t] && JT.is_defined[t])
	Jk[t]/=JT.Jk[t];
      else
	is_defined[t]=0;
    }
  }
  return *this;
}

//slash equals scalar
DoubleJackknifeTimeSeries & DoubleJackknifeTimeSeries::operator/=(double d)
{
  for (int t=0; t<=Nt; t++)
    if (is_defined[t])
      Jk[t]/=d;
  return *this;
}

DoubleJackknifeTimeSeries & DoubleJackknifeTimeSeries::Reverse()
{
  DoubleJackknifeTimeSeries tmp(*this);
  for (int t=0; t<=Nt; t++) {
    Jk[Nt-t]=tmp.Jk[t];
    is_defined[Nt-t]=tmp.is_defined[t];
  }
  return *this;
}

DoubleJackknifeTimeSeries Reverse(const DoubleJackknifeTimeSeries & JT)
{
  DoubleJackknifeTimeSeries tmp(JT);
  tmp.Reverse();
  return tmp;
}

DoubleJackknifeTimeSeries & DoubleJackknifeTimeSeries::Translate(int t)
{
  if (t>0 && t<=Nt) {
    for (int tt=Nt; tt>=0; tt--) {
      if (ReturnIsDefined(tt-t)) {
	Jk[tt]=Jk[tt-t];
	is_defined[tt]=1;
      } else {
	is_defined[tt]=0;
      }
    }
  } else if (t<0 && t>=-Nt) {
    for (int tt=0; tt<=Nt; tt++) {
      if (ReturnIsDefined(tt-t)) {
	Jk[tt]=Jk[tt-t];
	is_defined[tt]=1;
      } else {
	is_defined[tt]=0;
      }
    }
  } else if (t<-Nt || t>Nt) {
    cout << "Can only translate by amount -Nt<=delta_t<=Nt, but received delta_t=" << t << ".  Do nothing.\n";
  }
  //If t==0 then do nothing.
  return *this;
}

DoubleJackknifeTimeSeries Translate(int t, const DoubleJackknifeTimeSeries & JT)
{
  DoubleJackknifeTimeSeries tmp(JT);
  tmp.Translate(t);
  return tmp;
}

DoubleJackknifeTimeSeries & DoubleJackknifeTimeSeries::TranslateReverse(int t)
{
  if (t>=0 && t<=Nt) {
    DoubleJackknifeTimeSeries tmp(*this);
    for (int tt=0; tt<=t; tt++) {
      Jk[t-tt]=tmp.Jk[tt];
      is_defined[t-tt]=tmp.is_defined[tt];
    }
    for (int tt=t+1; tt<=Nt; tt++)
      is_defined[tt]=0;
  } else if (t>Nt && t<=2*Nt) {
    DoubleJackknifeTimeSeries tmp(*this);
    for (int tt=0; tt<=Nt; tt++)
      if (t-tt <= Nt) {
	Jk[t-tt]=tmp.Jk[tt];
	is_defined[t-tt]=tmp.is_defined[tt];
      }
    for (int tt=t-Nt-1; tt>=0; tt--)
      is_defined[tt]=0;
  } else {
    cout << "Can only translate and reverse for translations 0<=delta_t<=2*Nt, but received delta_t=" << t << ".  Do nothing.\n";
  }
  return *this;
}

DoubleJackknifeTimeSeries TranslateReverse(int t, const DoubleJackknifeTimeSeries & JT)
{
  DoubleJackknifeTimeSeries tmp(JT);
  tmp.TranslateReverse(t);
  return tmp;
}

//Note that there will be a problem if the argument of the
//log is negative, which happens when Jk[t] and Jk[t-1]
//have different signs.  To avoid this it's best to only
//do M_eff to something that's positive definite (e.g.
//something you've already taken the Abs of).
DoubleJackknifeTimeSeries & DoubleJackknifeTimeSeries::M_eff()
{
  DoubleJackknifeTimeSeries tmp(*this);
  is_defined[0]=0;
  for (int t=1; t<=Nt; t++) {
    if (tmp.is_defined[t] && tmp.is_defined[t-1]) {
      Jk[t]=-Log(tmp.Jk[t]/tmp.Jk[t-1]);
      is_defined[t]=1;
    } else
      is_defined[t]=0;
  }
  return *this;
}

DoubleJackknifeTimeSeries M_eff(const DoubleJackknifeTimeSeries & JT)
{
  DoubleJackknifeTimeSeries tmp(JT);
  tmp.M_eff();
  return tmp;
}

//Abs (member function)
DoubleJackknifeTimeSeries & DoubleJackknifeTimeSeries::Abs()
{
  for (int t=0; t<=Nt; t++)
    if (is_defined[t]) {
      DoubleJackknife tmp(Jk[t]);
      tmp.Abs(); //Have to do it this way or else it gets
                 //confused with DoubleJackknifeTimeSeries::Abs
      Jk[t]=tmp;
    }
  return *this;
}

//Abs (friend function)
DoubleJackknifeTimeSeries Abs(const DoubleJackknifeTimeSeries & JT)
{
  DoubleJackknifeTimeSeries tmp(JT);
  tmp.Abs();
  return tmp;
}

//Abs (friend function with a DoubleJackknifeCmplxTimeSeries argument)
DoubleJackknifeTimeSeries Abs(const DoubleJackknifeCmplxTimeSeries & JT)
{
  DoubleJackknifeTimeSeries tmp(JT.Nt,JT.N);
  for (int t=0; t<=tmp.Nt; t++) {
    if (JT.is_defined[t]) {
      tmp.Jk[t]=Abs(JT.Jk[t]);
      tmp.is_defined[t]=1;
    } else {
      tmp.is_defined[t]=0;
    }
  }
  return tmp;
}

//Real (friend function)
DoubleJackknifeTimeSeries Real(const DoubleJackknifeCmplxTimeSeries & JT)
{
  DoubleJackknifeTimeSeries tmp(JT.Nt,JT.N);
  for (int t=0; t<=tmp.Nt; t++) {
    if (JT.is_defined[t]) {
      tmp.Jk[t]=Real(JT.Jk[t]);
      tmp.is_defined[t]=1;
    } else {
      tmp.is_defined[t]=0;
    }
  }
  return tmp;
}

//Imag (friend function)
DoubleJackknifeTimeSeries Imag(const DoubleJackknifeCmplxTimeSeries & JT)
{
  DoubleJackknifeTimeSeries tmp(JT.Nt,JT.N);
  for (int t=0; t<=tmp.Nt; t++) {
    if (JT.is_defined[t]) {
      tmp.Jk[t]=Imag(JT.Jk[t]);
      tmp.is_defined[t]=1;
    } else {
      tmp.is_defined[t]=0;
    }
  }
  return tmp;
}

//Normalize (member function)
//Divides everything by the largest average value on the time slices.
//Returns the factor by which the object was divided.
//This is needed because Sam's fitter has problems with large numbers.
double DoubleJackknifeTimeSeries::Normalize()
{
  double fac=1.0;
  for (int t=0; t<=Nt; t++)
    if (is_defined[t])
      if (abs(Jk[t].ReturnAve())>fac)
	fac=abs(Jk[t].ReturnAve());
  *this/=fac;
  return fac;
}

//Combine (friend function)
//Takes two DoubleJackknifeTimeSeries objects with N jaccknife values and puts
//them together to make another DoubleJackknifeTimeSeries object with 2*N
//jackknife values.
//IMPORTANT: NOTE THAT it can only be used on quantities that depend linearly
//on the data
DoubleJackknifeTimeSeries Combine(const DoubleJackknifeTimeSeries & JT1, const DoubleJackknifeTimeSeries & JT2)
{
  if (JT1.N != JT2.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can only combine DoubleJackknifeTimeSeries objects with same number of jackknife values.  Object 1 N=" << JT1.N << " Object 2 N=" << JT2.N << "\n";
    return JT1;
  }
  if (JT1.Nt != JT2.Nt) {
    //Error: I think I'd rather have it abort here.
    cout << "Can only combine DoubleJackknifeTimeSeries objects with same Nt.  Object 1 Nt=" << JT1.Nt << " Object 2 Nt=" << JT2.Nt << "\n";
    return JT1;
  }
  int N=JT1.N;
  int Nt=JT1.Nt;
  DoubleJackknifeTimeSeries tmp(Nt,2*N);
  for (int t=0; t<=Nt; t++) {
    if (JT1.is_defined[t] && JT2.is_defined[t]) {
      tmp.Jk[t]=Combine(JT1.Jk[t],JT2.Jk[t]);
      tmp.is_defined[t]=1;
    } else {
      tmp.is_defined[t]=0;
    }
  }
  return tmp;
}

//Output jackknife values and errors to a file (friend function)
void OutputJk(string filename, const DoubleJackknifeTimeSeries & JT)
{
  ofstream fout(filename.c_str());
  fout.precision(16);
  int Nt=JT.ReturnNt();
  int N=JT.ReturnN();
  for (int i=0; i<N; i++)
    for (int t=0; t<=Nt; t++)
      if (JT.ReturnIsDefined(t))
	fout << t << " " << JT.ReturnJk(t,i) << " " << JT.ReturnJkErr(t,i) << "\n";
  fout.close();
}

//Output average and errors to a file
//not a member or friend function
void OutputAve(string filename, const DoubleJackknifeTimeSeries & JT)
{
  ofstream fout(filename.c_str());
  fout.precision(16);
  int Nt=JT.ReturnNt();
  for (int t=0; t<=Nt; t++)
    if (JT.ReturnIsDefined(t))
      fout << t << " " << JT.ReturnAve(t) << " " << JT.ReturnErr(t) << "\n";
  fout.close();
}




//DoubleJackknifeCmplxTimeSeries class

//Standard and Default Constructor
DoubleJackknifeCmplxTimeSeries::DoubleJackknifeCmplxTimeSeries(int nt, int n)
{
  if (n<3) {
    //Error: I think I'd rather have it abort here.
    cout << "N must be >= 3 for a DoubleJackknifeCmplxTimeSeries object, but was given N=" << n << ".  Setting N=3.\n";
    N=3;
  } else {
    N=n;
  }
  if (nt<1) {
    //Error: I think I'd rather have it abort here.
    cout << "Nt must be >= 1 for a DoubleJackknifeCmplxTimeSeries object, but was given Nt=" << nt << ".  Setting Nt=1.\n";
    Nt=1;
  } else {
    Nt=nt;
  }
  Jk=new DoubleJackknifeCmplx [Nt+1];
  is_defined=new int [Nt+1];
  for (int t=0; t<=Nt; t++) {
    is_defined[t]=0;
    //Have to initialize a DoubleJackknifeCmplx object of size N for each Jk[t]
    DoubleJackknifeCmplx tmp(N);
    Jk[t]=tmp;
  }
}

//Copy Constructor
DoubleJackknifeCmplxTimeSeries::DoubleJackknifeCmplxTimeSeries(const DoubleJackknifeCmplxTimeSeries & JT)
{
  N=JT.N;
  Nt=JT.Nt;
  Jk=new DoubleJackknifeCmplx [Nt+1];
  is_defined=new int [Nt+1];
  for (int t=0; t<=Nt; t++) {
    Jk[t]=JT.Jk[t];
    is_defined[t]=JT.is_defined[t];
  }
}

//Constructor with a DoubleJackknifeCmplx object
DoubleJackknifeCmplxTimeSeries::DoubleJackknifeCmplxTimeSeries(int nt, const DoubleJackknifeCmplx & J)
{
  if (nt<1) {
    //Error: I think I'd rather have it abort here.
    cout << "Nt must be >= 1 for a DoubleJackknifeCmplxTimeSeries object, but was given Nt=" << nt << ".  Setting Nt=1.\n";
    Nt=1;
  } else {
    Nt=nt;
  }
  N=J.ReturnN();
  Jk=new DoubleJackknifeCmplx [Nt+1];
  is_defined=new int [Nt+1];
  for (int t=0; t<=Nt; t++) {
    Jk[t]=J;
    is_defined[t]=1;
  }
}

//Destructor
DoubleJackknifeCmplxTimeSeries::~DoubleJackknifeCmplxTimeSeries() 
{
  delete [] Jk;
  delete [] is_defined;
}

//Assignment operator
//Deep copy
DoubleJackknifeCmplxTimeSeries & DoubleJackknifeCmplxTimeSeries::operator=(const DoubleJackknifeCmplxTimeSeries & JT)
{
  if (this == &JT)
    return *this;
  delete [] Jk;
  delete [] is_defined;
  N=JT.N;
  Nt=JT.Nt;
  Jk=new DoubleJackknifeCmplx [Nt+1];
  is_defined=new int [Nt+1];
  for (int t=0; t<=Nt; t++) {
    Jk[t]=JT.Jk[t];
    is_defined[t]=JT.is_defined[t];
  }
  return *this;
}

//Say whether or not there is data at t
int DoubleJackknifeCmplxTimeSeries::ReturnIsDefined(int t) const
{
  if (t<0 || t>Nt)
    return 0;
  else
    return is_defined[t];
}

//Return the average at t
complex<double> DoubleJackknifeCmplxTimeSeries::ReturnAve(int t) const
{
  if (t<0 || t>Nt || !is_defined[t]) {
    //Error:  I think I'd rather have it abort here.
    cout << "No data at t=" << t << " in DoubleJackknifeCmplxTimeSeries.\n";
    return 0.0;
  } else {
    return Jk[t].ReturnAve();
  }
}

//Return the error at t
double DoubleJackknifeCmplxTimeSeries::ReturnErr(int t) const
{
  if (t<0 || t>Nt || !is_defined[t]) {
    //Error:  I think I'd rather have it abort here.
    cout << "No data at t=" << t << " in DoubleJackknifeCmplxTimeSeries.\n";
    return 0.0;
  } else {
    return Jk[t].ReturnErr();
  }
}

//Return the ith jackknife value at t
complex<double> DoubleJackknifeCmplxTimeSeries::ReturnJk(int t, int i) const
{
  if (t<0 || t>Nt || !is_defined[t]) {
    //Error:  I think I'd rather have it abort here.
    cout << "No data at t=" << t << " in DoubleJackknifeCmplxTimeSeries.\n";
    return 0.0;
  }
  //Use DoubleJackknifeCmplx's own bounds checking for i.
  return Jk[t].ReturnJk(i);
}

//Return the error on the ith jackknife value at t
double DoubleJackknifeCmplxTimeSeries::ReturnJkErr(int t, int i) const
{
  if (t<0 || t>Nt || !is_defined[t]) {
    //Error:  I think I'd rather have it abort here.
    cout << "No data at t=" << t << " in DoubleJackknifeCmplxTimeSeries.\n";
    return 0.0;
  }
  //Use DoubleJackknifeCmplx's own bounds checking for i.
  return Jk[t].ReturnJkErr(i);
}

//Return double jackknife value i, j at t
complex<double> DoubleJackknifeCmplxTimeSeries::ReturnDoubleJk(int t, int i, int j) const
{
  if (t<0 || t>Nt || !is_defined[t]) {
    //Error:  I think I'd rather have it abort here.
    cout << "No data at t=" << t << " in DoubleJackknifeCmplxTimeSeries.\n";
    return 0.0;
  }
  //Use DoubleJackknifeCmplx's own bounds checking for i, j.
  return Jk[t].ReturnDoubleJk(i,j);
}


//Return Nt
int DoubleJackknifeCmplxTimeSeries::ReturnNt() const
{
  return Nt;
}

//Return N
int DoubleJackknifeCmplxTimeSeries::ReturnN() const
{
  return N;
}

//Takes a JackknifeCmplxTimeSeries object, undoes it, then calculates 
//double jackknife values to put into the current
//DoubleJackknifeCmplxTimeSeries object.
//NOTE: ONLY WORKS for quantities that depend linearly on original data.
DoubleJackknifeCmplxTimeSeries & DoubleJackknifeCmplxTimeSeries::FromSingleJk(const JackknifeCmplxTimeSeries & JT)
{
  delete [] Jk;
  delete [] is_defined;
  N=JT.N;
  Nt=JT.Nt;
  Jk=new DoubleJackknifeCmplx[Nt+1];
  is_defined=new int [Nt+1];
  for (int t=0; t<=Nt; t++) {
    Jk[t].FromSingleJk(JT.Jk[t]);
    if (JT.ReturnIsDefined(t)) {
      is_defined[t]=1;
    } else {
      is_defined[t]=0;
    }
  }
  return *this;
}

//Friend function version of the above.
DoubleJackknifeCmplxTimeSeries FromSingleJk(const JackknifeCmplxTimeSeries & JT)
{
  DoubleJackknifeCmplxTimeSeries tmp;
  tmp.FromSingleJk(JT);
  return tmp;
}

//Addition
DoubleJackknifeCmplxTimeSeries DoubleJackknifeCmplxTimeSeries::operator+(const DoubleJackknifeCmplxTimeSeries & JT) const
{
  DoubleJackknifeCmplxTimeSeries sum(*this);
  sum+=JT;
  return sum;
}

//+=
DoubleJackknifeCmplxTimeSeries & DoubleJackknifeCmplxTimeSeries::operator+=(const DoubleJackknifeCmplxTimeSeries & JT)
{
  if (N!=JT.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't add DoubleJackknifeCmplxTimeSeries objects with different number of jackknife values.  Object 1 N=" << N << " Object 2 N=" << JT.N << "\n";
  } else if (Nt!=JT.Nt) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't add DoubleJackknifeCmplxTimeSeries objects with different Nt.  Object 1 Nt=" << Nt << " Object 2 Nt=" << JT.Nt << "\n";
  } else {
    for (int t=0; t<=Nt; t++) {
      if (is_defined[t] && JT.is_defined[t])
	Jk[t]+=JT.Jk[t];
      else
	is_defined[t]=0;
    }
  }
  return *this;
}

//Subtraction
DoubleJackknifeCmplxTimeSeries DoubleJackknifeCmplxTimeSeries::operator-(const DoubleJackknifeCmplxTimeSeries & JT) const
{
  DoubleJackknifeCmplxTimeSeries difference(*this);
  difference-=JT;
  return difference;
}

//Negative (Unary operator)
DoubleJackknifeCmplxTimeSeries DoubleJackknifeCmplxTimeSeries::operator-() const
{
  DoubleJackknifeCmplx tmp(N,0.0);
  DoubleJackknifeCmplxTimeSeries neg(Nt,tmp);
  neg-=*this;
  return neg;
}

//-=
DoubleJackknifeCmplxTimeSeries & DoubleJackknifeCmplxTimeSeries::operator-=(const DoubleJackknifeCmplxTimeSeries & JT)
{
  if (N!=JT.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't subtract DoubleJackknifeCmplxTimeSeries objects with different number of jackknife values.  Object 1 N=" << N << " Object 2 N=" << JT.N << "\n";
  } else if (Nt!=JT.Nt) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't subtract DoubleJackknifeCmplxTimeSeries objects with different Nt.  Object 1 Nt=" << Nt << " Object 2 Nt=" << JT.Nt << "\n";
  } else {
    for (int t=0; t<=Nt; t++) {
      if (is_defined[t] && JT.is_defined[t])
	Jk[t]-=JT.Jk[t];
      else
	is_defined[t]=0;
    }
  }
  return *this;
}

//Multiplication of two DoubleJackknifeCmplxTimeSeries objects
DoubleJackknifeCmplxTimeSeries DoubleJackknifeCmplxTimeSeries::operator*(const DoubleJackknifeCmplxTimeSeries & JT) const
{
  DoubleJackknifeCmplxTimeSeries prod(*this);
  prod*=JT;
  return prod;
}

//Multiplication of a DoubleJackknifeCmplxTimeSeries object by a complex scalar on the right
DoubleJackknifeCmplxTimeSeries DoubleJackknifeCmplxTimeSeries::operator*(complex<double> c) const
{
  DoubleJackknifeCmplxTimeSeries prod(*this);
  prod*=c;
  return prod;
}

//Multiplication by a complex scalar on the left (friend function)
DoubleJackknifeCmplxTimeSeries operator*(complex<double> c, const DoubleJackknifeCmplxTimeSeries & JT)
{
  return JT*c;
}

//*=DoubleJackknifeCmplxTimeSeries object
DoubleJackknifeCmplxTimeSeries & DoubleJackknifeCmplxTimeSeries::operator*=(const DoubleJackknifeCmplxTimeSeries & JT)
{
  if (N!=JT.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't multiply DoubleJackknifeCmplxTimeSeries objects with different number of jackknife values.  Object 1 N=" << N << " Object 2 N=" << JT.N << "\n";
  } else if (Nt!=JT.Nt) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't multiply DoubleJackknifeTimeCmplxSeries objects with different Nt.  Object 1 Nt=" << Nt << " Object 2 Nt=" << JT.Nt << "\n";
  } else {
    for (int t=0; t<=Nt; t++) {
      if (is_defined[t] && JT.is_defined[t])
	Jk[t]*=JT.Jk[t];
      else
	is_defined[t]=0;
    }
  }
  return *this;
}

//*=scalar
DoubleJackknifeCmplxTimeSeries & DoubleJackknifeCmplxTimeSeries::operator*=(complex<double> c)
{
  for (int t=0; t<=Nt; t++)
    if (is_defined[t])
      Jk[t]*=c;
  return *this;
}

//Division of two DoubleJackknifeCmplxTimeSeries objects
DoubleJackknifeCmplxTimeSeries DoubleJackknifeCmplxTimeSeries::operator/(const DoubleJackknifeCmplxTimeSeries & JT) const
{
  DoubleJackknifeCmplxTimeSeries q(*this);
  q/=JT;
  return q;
}

//Division by a complex scalar
DoubleJackknifeCmplxTimeSeries DoubleJackknifeCmplxTimeSeries::operator/(complex<double> c) const
{
  DoubleJackknifeCmplxTimeSeries q(*this);
  q/=c;
  return q;
}

//Division of a complex scalar by a DoubleJackknifeCmplxTimeSeries object (friend function)
DoubleJackknifeCmplxTimeSeries operator/(complex<double> c, const DoubleJackknifeCmplxTimeSeries & JT)
{
  int Nt=JT.Nt;
  int N=JT.N;
  DoubleJackknifeCmplx tmp(N,c);
  DoubleJackknifeCmplxTimeSeries q(Nt,tmp);
  q/=JT;
  return q;
}

//slash equals DoubleJackknifeCmplxTimeSeries object
DoubleJackknifeCmplxTimeSeries & DoubleJackknifeCmplxTimeSeries::operator/=(const DoubleJackknifeCmplxTimeSeries & JT)
{
  if (N!=JT.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't divide DoubleJackknifeCmplxTimeSeries objects with different number of jackknife values.  Object 1 N=" << N << " Object 2 N=" << JT.N << "\n";
  } else if (Nt!=JT.Nt) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't divide DoubleJackknifeCmplxTimeSeries objects with different Nt.  Object 1 Nt=" << Nt << " Object 2 Nt=" << JT.Nt << "\n";
  } else {
    for (int t=0; t<=Nt; t++) {
      if (is_defined[t] && JT.is_defined[t])
	Jk[t]/=JT.Jk[t];
      else
	is_defined[t]=0;
    }
  }
  return *this;
}

//slash equals complex scalar
DoubleJackknifeCmplxTimeSeries & DoubleJackknifeCmplxTimeSeries::operator/=(complex<double> c)
{
  for (int t=0; t<=Nt; t++)
    if (is_defined[t])
      Jk[t]/=c;
  return *this;
}

DoubleJackknifeCmplxTimeSeries & DoubleJackknifeCmplxTimeSeries::Reverse()
{
  DoubleJackknifeCmplxTimeSeries tmp(*this);
  for (int t=0; t<=Nt; t++) {
    Jk[Nt-t]=tmp.Jk[t];
    is_defined[Nt-t]=tmp.is_defined[t];
  }
  return *this;
}

DoubleJackknifeCmplxTimeSeries Reverse(const DoubleJackknifeCmplxTimeSeries & JT)
{
  DoubleJackknifeCmplxTimeSeries tmp(JT);
  tmp.Reverse();
  return tmp;
}

DoubleJackknifeCmplxTimeSeries & DoubleJackknifeCmplxTimeSeries::Translate(int t)
{
  if (t>0 && t<=Nt) {
    for (int tt=Nt; tt>=0; tt--) {
      if (ReturnIsDefined(tt-t)) {
	Jk[tt]=Jk[tt-t];
	is_defined[tt]=1;
      } else {
	is_defined[tt]=0;
      }
    }
  } else if (t<0 && t>=-Nt) {
    for (int tt=0; tt<=Nt; tt++) {
      if (ReturnIsDefined(tt-t)) {
	Jk[tt]=Jk[tt-t];
	is_defined[tt]=1;
      } else {
	is_defined[tt]=0;
      }
    }
  } else if (t<-Nt || t>Nt) {
    cout << "Can only translate by amount -Nt<=delta_t<=Nt, but received delta_t=" << t << ".  Do nothing.\n";
  }
  //If t==0 then do nothing.
  return *this;
}

DoubleJackknifeCmplxTimeSeries Translate(int t, const DoubleJackknifeCmplxTimeSeries & JT)
{
  DoubleJackknifeCmplxTimeSeries tmp(JT);
  tmp.Translate(t);
  return tmp;
}

DoubleJackknifeCmplxTimeSeries & DoubleJackknifeCmplxTimeSeries::TranslateReverse(int t)
{
  if (t>=0 && t<=Nt) {
    DoubleJackknifeCmplxTimeSeries tmp(*this);
    for (int tt=0; tt<=t; tt++) {
      Jk[t-tt]=tmp.Jk[tt];
      is_defined[t-tt]=tmp.is_defined[tt];
    }
    for (int tt=t+1; tt<=Nt; tt++)
      is_defined[tt]=0;
  } else if (t>Nt && t<=2*Nt) {
    DoubleJackknifeCmplxTimeSeries tmp(*this);
    for (int tt=0; tt<=Nt; tt++)
      if (t-tt <= Nt) {
	Jk[t-tt]=tmp.Jk[tt];
	is_defined[t-tt]=tmp.is_defined[tt];
      }
    for (int tt=t-Nt-1; tt>=0; tt--)
      is_defined[tt]=0;
  } else {
    cout << "Can only translate and reverse for translations 0<=delta_t<=2*Nt, but received delta_t=" << t << ".  Do nothing.\n";
  }
  return *this;
}

DoubleJackknifeCmplxTimeSeries TranslateReverse(int t, const DoubleJackknifeCmplxTimeSeries & JT)
{
  DoubleJackknifeCmplxTimeSeries tmp(JT);
  tmp.TranslateReverse(t);
  return tmp;
}

//Abs (member function)
DoubleJackknifeCmplxTimeSeries & DoubleJackknifeCmplxTimeSeries::Abs()
{
  for (int t=0; t<=Nt; t++)
    if (is_defined[t]) {
      DoubleJackknifeCmplx tmp(Jk[t]);
      tmp.Abs(); //Have to do it this way or else it gets confused
                 //with DoubleJackknifeCmplxTimeSeries::Abs
      Jk[t]=tmp;
    }
  return *this;
}

//Real (member function)
DoubleJackknifeCmplxTimeSeries & DoubleJackknifeCmplxTimeSeries::Real()
{
  for (int t=0; t<=Nt; t++)
    if (is_defined[t]) {
      DoubleJackknifeCmplx tmp(Jk[t]);
      tmp.Real(); //Have to do it this way or else it gets confused
                  //with DoubleJackknifeCmplxTimeSeries::Real
      Jk[t]=tmp;
    }
  return *this;
}

//Imag (member function)
DoubleJackknifeCmplxTimeSeries & DoubleJackknifeCmplxTimeSeries::Imag()
{
  for (int t=0; t<=Nt; t++)
    if (is_defined[t]) {
      DoubleJackknifeCmplx tmp(Jk[t]);
      tmp.Imag(); //Have to do it this way or else it gets confused
                  //with DoubleJackknifeCmplxTimeSeries::Real
      Jk[t]=tmp;
    }
  return *this;
}

//Conj (member function)
DoubleJackknifeCmplxTimeSeries & DoubleJackknifeCmplxTimeSeries::Conj()
{
  for (int t=0; t<=Nt; t++)
    if (is_defined[t]) {
      DoubleJackknifeCmplx tmp(Jk[t]);
      tmp.Conj(); //Have to do it this way or else it gets confused
                  //with DoubleJackknifeCmplxTimeSeries::Abs
      Jk[t]=tmp;
    }
  return *this;
}

//Conj (friend function)
DoubleJackknifeCmplxTimeSeries Conj(const DoubleJackknifeCmplxTimeSeries & JT)
{
  DoubleJackknifeCmplxTimeSeries tmpc(JT);
  tmpc.Conj();
  return tmpc;
}

//Normalize (member function)
//Divides everything by the largest modulus of the average values on
//the time slices.
//Returns the factor by which the object was divided.
//This is needed because Sam's fitter has problems with large numbers.
double DoubleJackknifeCmplxTimeSeries::Normalize()
{
  double fac=1.0;
  for (int t=0; t<=Nt; t++)
    if (is_defined[t])
      if (abs(Jk[t].ReturnAve())>fac)
	fac=abs(Jk[t].ReturnAve());
  *this/=fac;
  return fac;
}

//Combine (friend function)
//Takes two DoubleJackknifeCmplxTimeSeries objects with N jaccknife values and
//puts them together to make another JackknifeCmplxTimeSeries object with 2*N
//jackknife values.
//IMPORTANT: NOTE THAT it can only be used on quantities that depend linearly
//on the data
DoubleJackknifeCmplxTimeSeries Combine(const DoubleJackknifeCmplxTimeSeries & JT1, const DoubleJackknifeCmplxTimeSeries & JT2)
{
  if (JT1.N != JT2.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can only combine DoubleJackknifeCmplxTimeSeries objects with same number of jackknife values.  Object 1 N=" << JT1.N << " Object 2 N=" << JT2.N << "\n";
    return JT1;
  }
  if (JT1.Nt != JT2.Nt) {
    //Error: I think I'd rather have it abort here.
    cout << "Can only combine DoubleJackknifeCmplxTimeSeries objects with same Nt.  Object 1 Nt=" << JT1.Nt << " Object 2 Nt=" << JT2.Nt << "\n";
    return JT1;
  }
  int N=JT1.N;
  int Nt=JT1.Nt;
  DoubleJackknifeCmplxTimeSeries tmp(Nt,2*N);
  for (int t=0; t<=Nt; t++) {
    if (JT1.is_defined[t] && JT2.is_defined[t]) {
      tmp.Jk[t]=Combine(JT1.Jk[t],JT2.Jk[t]);
      tmp.is_defined[t]=1;
    } else {
      tmp.is_defined[t]=0;
    }
  }
  return tmp;
}

//Output jackknife values and errors to a file (friend function)
void OutputJk(string filename, const DoubleJackknifeCmplxTimeSeries & JT)
{
  ofstream fout(filename.c_str());
  fout.precision(16);
  int Nt=JT.ReturnNt();
  int N=JT.ReturnN();
  for (int i=0; i<N; i++)
    for (int t=0; t<=Nt; t++)
      if (JT.ReturnIsDefined(t))
	fout << t << " " << real(JT.ReturnJk(t,i)) << " " << imag(JT.ReturnJk(t,i)) << " " << JT.ReturnJkErr(t,i) << "\n";
  fout.close();
}

//Output average and errors to a file
//not a member or friend function
void OutputAve(string filename, const DoubleJackknifeCmplxTimeSeries & JT)
{
  ofstream fout(filename.c_str());
  fout.precision(16);
  int Nt=JT.ReturnNt();
  for (int t=0; t<=Nt; t++)
    if (JT.ReturnIsDefined(t))
      fout << t << " " << real(JT.ReturnAve(t)) << " " << imag(JT.ReturnAve(t)) << " " << JT.ReturnErr(t) << "\n";
  fout.close();
}

