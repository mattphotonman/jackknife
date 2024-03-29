#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <cstdlib>
using namespace std;
#include "jackknife.h"
#include "jackknife_time_series.h"


//JackknifeTimeSeries class

//Standard and Default Constructor
JackknifeTimeSeries::JackknifeTimeSeries(int nt, int n)
{
  if (n<2) {
    //Error: I think I'd rather have it abort here.
    cout << "N must be >= 2 for a JackknifeTimeSeries object, but was given N=" << n << ".  Setting N=2.\n";
    N=2;
  } else {
    N=n;
  }
  if (nt<1) {
    //Error: I think I'd rather have it abort here.
    cout << "Nt must be >= 1 for a JackknifeTimeSeries object, but was given Nt=" << nt << ".  Setting Nt=1.\n";
    Nt=1;
  } else {
    Nt=nt;
  }
  Jk=new Jackknife [Nt+1];
  is_defined=new int [Nt+1];
  for (int t=0; t<=Nt; t++) {
    is_defined[t]=0;
    //Have to initialize a Jackknife object of size N for each Jk[t]
    Jackknife tmp(N);
    Jk[t]=tmp;
  }
}

//Copy Constructor
JackknifeTimeSeries::JackknifeTimeSeries(const JackknifeTimeSeries & JT)
{
  N=JT.N;
  Nt=JT.Nt;
  Jk=new Jackknife [Nt+1];
  is_defined=new int [Nt+1];
  for (int t=0; t<=Nt; t++) {
    Jk[t]=JT.Jk[t];
    is_defined[t]=JT.is_defined[t];
  }
}

//Constructor with a Jackknife object
JackknifeTimeSeries::JackknifeTimeSeries(int nt, const Jackknife & J)
{
  if (nt<1) {
    //Error: I think I'd rather have it abort here.
    cout << "Nt must be >= 1 for a JackknifeTimeSeries object, but was given Nt=" << nt << ".  Setting Nt=1.\n";
    Nt=1;
  } else {
    Nt=nt;
  }
  N=J.ReturnN();
  Jk=new Jackknife [Nt+1];
  is_defined=new int [Nt+1];
  for (int t=0; t<=Nt; t++) {
    Jk[t]=J;
    is_defined[t]=1;
  }
}

//Destructor
JackknifeTimeSeries::~JackknifeTimeSeries() 
{
  delete [] Jk;
  delete [] is_defined;
}

//Assignment operator
//Deep copy
JackknifeTimeSeries & JackknifeTimeSeries::operator=(const JackknifeTimeSeries & JT)
{
  if (this == &JT)
    return *this;
  delete [] Jk;
  delete [] is_defined;
  N=JT.N;
  Nt=JT.Nt;
  Jk=new Jackknife [Nt+1];
  is_defined=new int [Nt+1];
  for (int t=0; t<=Nt; t++) {
    Jk[t]=JT.Jk[t];
    is_defined[t]=JT.is_defined[t];
  }
  return *this;
}

//Say whether or not there is data at t
int JackknifeTimeSeries::ReturnIsDefined(int t) const
{
  if (t<0 || t>Nt)
    return 0;
  else
    return is_defined[t];
}

//Return the average at t
double JackknifeTimeSeries::ReturnAve(int t) const
{
  if (t<0 || t>Nt || !is_defined[t]) {
    //Error:  I think I'd rather have it abort here.
    cout << "No data at t=" << t << " in JackknifeTimeSeries.\n";
    return 0.0;
  } else {
    return Jk[t].ReturnAve();
  }
}

//Return the error at t
double JackknifeTimeSeries::ReturnErr(int t) const
{
  if (t<0 || t>Nt || !is_defined[t]) {
    //Error:  I think I'd rather have it abort here.
    cout << "No data at t=" << t << " in JackknifeTimeSeries.\n";
    return 0.0;
  } else {
    return Jk[t].ReturnErr();
  }
}

//Return the ith jackknife value at t
double JackknifeTimeSeries::ReturnJk(int t, int i) const
{
  if (t<0 || t>Nt || !is_defined[t]) {
    //Error:  I think I'd rather have it abort here.
    cout << "No data at t=" << t << " in JackknifeTimeSeries.\n";
    return 0.0;
  }
  if (i<0 || i>=N) {
    //Error:  I think I'd rather have it abort here.
    cout << "i must be 0<=i<N, but received i=" << i << ".\n";
    return 0.0;
  }
  return Jk[t].ReturnJk(i);
}

//Return Nt
int JackknifeTimeSeries::ReturnNt() const
{
  return Nt;
}

//Return N
int JackknifeTimeSeries::ReturnN() const
{
  return N;
}

//Uses an array of averages (ave_in) at each time slice, an array 
//of jackknife values (jk_in) at each time slice for each configuration,
//and an array (is_defined_in) specifying which time slices actually
//contain data, to set the value of the JackknifeTimeSeries object.
JackknifeTimeSeries & JackknifeTimeSeries::FromArray(double ave_in [], double** jk_in, int is_defined_in [], int Nt_in, int N_in)
{
  //Assumed sizes of arrays:
  //double ave_in[Nt_in+1]
  //double jk_in[Nt_in+1][N_in]
  //int is_defined_in[Nt_in+1]
  
  delete [] Jk;
  delete [] is_defined;
  N=N_in;
  Nt=Nt_in;
  Jk=new Jackknife [Nt+1];
  is_defined=new int [Nt+1];
  for (int t=0; t<=Nt; t++) {
    if (is_defined_in[t])
      Jk[t].FromArray(ave_in[t],jk_in[t],N);
    else {
      //still fill it with a blank Jackknife object
      Jackknife tmp(N);
      Jk[t]=tmp;
    }
    is_defined[t]=is_defined_in[t];
  }
  
  return *this;

}

//Addition
JackknifeTimeSeries JackknifeTimeSeries::operator+(const JackknifeTimeSeries & JT) const
{
  JackknifeTimeSeries sum(*this);
  sum+=JT;
  return sum;
}

//+=
JackknifeTimeSeries & JackknifeTimeSeries::operator+=(const JackknifeTimeSeries & JT)
{
  if (N!=JT.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't add JackknifeTimeSeries objects with different number of jackknife values.  Object 1 N=" << N << " Object 2 N=" << JT.N << "\n";
  } else if (Nt!=JT.Nt) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't add JackknifeTimeSeries objects with different Nt.  Object 1 Nt=" << Nt << " Object 2 Nt=" << JT.Nt << "\n";
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
JackknifeTimeSeries JackknifeTimeSeries::operator-(const JackknifeTimeSeries & JT) const
{
  JackknifeTimeSeries difference(*this);
  difference-=JT;
  return difference;
}

//Negative (Unary operator)
JackknifeTimeSeries JackknifeTimeSeries::operator-() const
{
  Jackknife tmp(N,0.0);
  JackknifeTimeSeries neg(Nt,tmp);
  neg-=*this;
  return neg;
}

//-=
JackknifeTimeSeries & JackknifeTimeSeries::operator-=(const JackknifeTimeSeries & JT)
{
  if (N!=JT.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't subtract JackknifeTimeSeries objects with different number of jackknife values.  Object 1 N=" << N << " Object 2 N=" << JT.N << "\n";
  } else if (Nt!=JT.Nt) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't subtract JackknifeTimeSeries objects with different Nt.  Object 1 Nt=" << Nt << " Object 2 Nt=" << JT.Nt << "\n";
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

//Multiplication of two JackknifeTimeSeries objects
JackknifeTimeSeries JackknifeTimeSeries::operator*(const JackknifeTimeSeries & JT) const
{
  JackknifeTimeSeries prod(*this);
  prod*=JT;
  return prod;
}

//Multiplication of a JackknifeTimeSeries object by a scalar on the right
JackknifeTimeSeries JackknifeTimeSeries::operator*(double d) const
{
  JackknifeTimeSeries prod(*this);
  prod*=d;
  return prod;
}

//Multiplication by a scalar on the left (friend function)
JackknifeTimeSeries operator*(double d, const JackknifeTimeSeries & JT)
{
  return JT*d;
}

//*=JackknifeTimeSeries object
JackknifeTimeSeries & JackknifeTimeSeries::operator*=(const JackknifeTimeSeries & JT)
{
  if (N!=JT.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't multiply JackknifeTimeSeries objects with different number of jackknife values.  Object 1 N=" << N << " Object 2 N=" << JT.N << "\n";
  } else if (Nt!=JT.Nt) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't multiply JackknifeTimeSeries objects with different Nt.  Object 1 Nt=" << Nt << " Object 2 Nt=" << JT.Nt << "\n";
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
JackknifeTimeSeries & JackknifeTimeSeries::operator*=(double d)
{
  for (int t=0; t<=Nt; t++)
    if (is_defined[t])
      Jk[t]*=d;
  return *this;
}

//Division of two JackknifeTimeSeries objects
JackknifeTimeSeries JackknifeTimeSeries::operator/(const JackknifeTimeSeries & JT) const
{
  JackknifeTimeSeries q(*this);
  q/=JT;
  return q;
}

//Division by a scalar
JackknifeTimeSeries JackknifeTimeSeries::operator/(double d) const
{
  JackknifeTimeSeries q(*this);
  q/=d;
  return q;
}

//Division of a scalar by a JackknifeTimeSeries object (friend function)
JackknifeTimeSeries operator/(double d, const JackknifeTimeSeries & JT)
{
  int Nt=JT.Nt;
  int N=JT.N;
  Jackknife tmp(N,d);
  JackknifeTimeSeries q(Nt,tmp);
  q/=JT;
  return q;
}

//slash equals JackknifeTimeSeries object
JackknifeTimeSeries & JackknifeTimeSeries::operator/=(const JackknifeTimeSeries & JT)
{
  if (N!=JT.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't divide JackknifeTimeSeries objects with different number of jackknife values.  Object 1 N=" << N << " Object 2 N=" << JT.N << "\n";
  } else if (Nt!=JT.Nt) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't divide JackknifeTimeSeries objects with different Nt.  Object 1 Nt=" << Nt << " Object 2 Nt=" << JT.Nt << "\n";
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
JackknifeTimeSeries & JackknifeTimeSeries::operator/=(double d)
{
  for (int t=0; t<=Nt; t++)
    if (is_defined[t])
      Jk[t]/=d;
  return *this;
}

JackknifeTimeSeries & JackknifeTimeSeries::Reverse()
{
  JackknifeTimeSeries tmp(*this);
  for (int t=0; t<=Nt; t++) {
    Jk[Nt-t]=tmp.Jk[t];
    is_defined[Nt-t]=tmp.is_defined[t];
  }
  return *this;
}

JackknifeTimeSeries Reverse(const JackknifeTimeSeries & JT)
{
  JackknifeTimeSeries tmp(JT);
  tmp.Reverse();
  return tmp;
}

JackknifeTimeSeries & JackknifeTimeSeries::Translate(int t)
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

JackknifeTimeSeries Translate(int t, const JackknifeTimeSeries & JT)
{
  JackknifeTimeSeries tmp(JT);
  tmp.Translate(t);
  return tmp;
}

JackknifeTimeSeries & JackknifeTimeSeries::TranslateReverse(int t)
{
  if (t>=0 && t<=Nt) {
    JackknifeTimeSeries tmp(*this);
    for (int tt=0; tt<=t; tt++) {
      Jk[t-tt]=tmp.Jk[tt];
      is_defined[t-tt]=tmp.is_defined[tt];
    }
    for (int tt=t+1; tt<=Nt; tt++)
      is_defined[tt]=0;
  } else if (t>Nt && t<=2*Nt) {
    JackknifeTimeSeries tmp(*this);
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

JackknifeTimeSeries TranslateReverse(int t, const JackknifeTimeSeries & JT)
{
  JackknifeTimeSeries tmp(JT);
  tmp.TranslateReverse(t);
  return tmp;
}

//Note that there will be a problem if the argument of the
//log is negative, which happens when Jk[t] and Jk[t-1]
//have different signs.  To avoid this it's best to only
//do M_eff to something that's positive definite (e.g.
//something you've already taken the Abs of).
JackknifeTimeSeries & JackknifeTimeSeries::M_eff()
{
  JackknifeTimeSeries tmp(*this);
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

JackknifeTimeSeries M_eff(const JackknifeTimeSeries & JT)
{
  JackknifeTimeSeries tmp(JT);
  tmp.M_eff();
  return tmp;
}

JackknifeTimeSeries & JackknifeTimeSeries::M_eff_two_point()
{
  JackknifeTimeSeries tmp(*this);
  is_defined[0]=0;
  for (int t=1; t<=Nt; t++) {
    if (tmp.is_defined[t] && tmp.is_defined[t-1]) {
      Jk[t]=Meff_two_point(tmp.Jk[t-1],tmp.Jk[t],t,Nt);
      is_defined[t]=1;
    } else
      is_defined[t]=0;
  }
  return *this;
}

JackknifeTimeSeries M_eff_two_point(const JackknifeTimeSeries & JT)
{
  JackknifeTimeSeries tmp(JT);
  tmp.M_eff_two_point();
  return tmp;
}

JackknifeTimeSeries & JackknifeTimeSeries::M_eff_three_point()
{
  JackknifeTimeSeries tmp(*this);
  is_defined[0]=0;
  is_defined[Nt]=0;
  for (int t=1; t<Nt; t++) {
    if (tmp.is_defined[t+1] && tmp.is_defined[t] && tmp.is_defined[t-1]) {
      Jk[t]=Meff_three_point(tmp.Jk[t-1],tmp.Jk[t],tmp.Jk[t+1],t,Nt);
      is_defined[t]=1;
    } else
      is_defined[t]=0;
  }
  return *this;
}

JackknifeTimeSeries M_eff_three_point(const JackknifeTimeSeries & JT)
{
  JackknifeTimeSeries tmp(JT);
  tmp.M_eff_three_point();
  return tmp;
}

//Bin (member function)
//Note: Only valid if it depends linearly on original data.
JackknifeTimeSeries & JackknifeTimeSeries::Bin(int bin_size)
{
  for (int t=0; t<=Nt; t++)
    if (is_defined[t]) {
      Jackknife tmp(Jk[t]);
      tmp.Bin(bin_size); //Have to do it this way or else it gets
                         //confused with JackknifeTimeSeries::Bin
      N=tmp.ReturnN();
      Jk[t]=tmp;
    }
  return *this;
}

//Bin (friend function)
//Note: Only valid if it depends linearly on original data.
JackknifeTimeSeries Bin(const JackknifeTimeSeries & JT, int bin_size)
{
  JackknifeTimeSeries tmp(JT);
  tmp.Bin(bin_size);
  return tmp;
}

//Abs (member function)
JackknifeTimeSeries & JackknifeTimeSeries::Abs()
{
  for (int t=0; t<=Nt; t++)
    if (is_defined[t]) {
      Jackknife tmp(Jk[t]);
      tmp.Abs(); //Have to do it this way or else it gets
                 //confused with JackknifeTimeSeries::Abs
      Jk[t]=tmp;
    }
  return *this;
}

//Abs (friend function)
JackknifeTimeSeries Abs(const JackknifeTimeSeries & JT)
{
  JackknifeTimeSeries tmp(JT);
  tmp.Abs();
  return tmp;
}

//Abs (friend function with a JackknifeCmplxTimeSeries argument)
JackknifeTimeSeries Abs(const JackknifeCmplxTimeSeries & JT)
{
  JackknifeTimeSeries tmp(JT.Nt,JT.N);
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
JackknifeTimeSeries Real(const JackknifeCmplxTimeSeries & JT)
{
  JackknifeTimeSeries tmp(JT.Nt,JT.N);
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
JackknifeTimeSeries Imag(const JackknifeCmplxTimeSeries & JT)
{
  JackknifeTimeSeries tmp(JT.Nt,JT.N);
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
double JackknifeTimeSeries::Normalize()
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
//Takes two JackknifeTimeSeries objects with N jaccknife values and puts them
//together to make another JackknifeTimeSeries object with 2*N jackknife
//values.
//IMPORTANT: NOTE THAT it can only be used on quantities that depend linearly
//on the data
JackknifeTimeSeries Combine(const JackknifeTimeSeries & JT1, const JackknifeTimeSeries & JT2)
{
  if (JT1.N != JT2.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can only combine JackknifeTimeSeries objects with same number of jackknife values.  Object 1 N=" << JT1.N << " Object 2 N=" << JT2.N << "\n";
    return JT1;
  }
  if (JT1.Nt != JT2.Nt) {
    //Error: I think I'd rather have it abort here.
    cout << "Can only combine JackknifeTimeSeries objects with same Nt.  Object 1 Nt=" << JT1.Nt << " Object 2 Nt=" << JT2.Nt << "\n";
    return JT1;
  }
  int N=JT1.N;
  int Nt=JT1.Nt;
  JackknifeTimeSeries tmp(Nt,2*N);
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

//Output jackknife values to a file in the appropriate format for Sam's
//fitter (friend function)
void OutputJk(string filename, const JackknifeTimeSeries & JT)
{
  ofstream fout(filename.c_str());
  fout.precision(16);
  int Nt=JT.ReturnNt();
  int N=JT.ReturnN();
  for (int i=0; i<N; i++)
    for (int t=0; t<=Nt; t++)
      if (JT.ReturnIsDefined(t))
	fout << t << " " << JT.Jk[t].jk[i] << "\n";
  fout.close();
}

//Output average and errors to a file
//not a member or friend function
void OutputAve(string filename, const JackknifeTimeSeries & JT)
{
  ofstream fout(filename.c_str());
  fout.precision(16);
  int Nt=JT.ReturnNt();
  for (int t=0; t<=Nt; t++)
    if (JT.ReturnIsDefined(t))
      fout << t << " " << JT.ReturnAve(t) << " " << JT.ReturnErr(t) << "\n";
  fout.close();
}

//Reads in a JackknifeTimeSeries from two text files, one containing
//the average values, and one containing the jackknife values
JackknifeTimeSeries & JackknifeTimeSeries::ReadTimeSeriesFile(string ave_file, string jks_file, int Ntfile, int Nfile)
{
  double ave_tmp[Ntfile+1];
  //double jk_tmp[Ntfile+1][Nfile];
  double** jk_tmp=new double* [Ntfile+1];
  for (int t=0; t<=Ntfile; t++)
    jk_tmp[t]=new double [Nfile];
  int is_defined_tmp[Ntfile+1];
  double errors[Ntfile+1];

  int success=ReadTimeSeriesFiletoArray(ave_file,jks_file,Ntfile,Nfile,ave_tmp,jk_tmp,is_defined_tmp,errors);

  if (!success)
    return *this;
    
  //Read it in correctly, now use arrays to build the object.
  FromArray(ave_tmp,jk_tmp,is_defined_tmp,Ntfile,Nfile);

  //Have to delete jk_tmp now that we're done with it since it was
  //declared locally.
  for (int t=0; t<=Nfile; t++)
    delete [] jk_tmp[t];
  delete [] jk_tmp;
  
  //Compare errors read to errors calculated.
  //Won't abort if they don't agree, just give a warning.
  double tol=1E-6;
  int any_disagree=0;
  int disagree1[Nt+1], disagree2[Nt+1];
  for (int t=0; t<=Nt; t++) {
    disagree1[t]=0;
    disagree2[t]=0;
    if (ReturnIsDefined(t))
      if (errors[t] != 0) {
	double reldiff=abs((ReturnErr(t)-errors[t])/errors[t]);
	if (reldiff>tol) {
	  //Error read and calculated disagree at this time slice.
	  any_disagree=1;
	  disagree1[t]=1;
	} 
      } else if ( abs(ReturnErr(t))>tol ) {
	//Error read and calculated disagree at this time slice because
	//the error read was 0 and the calculated error is not 0.
	any_disagree=1;
	disagree2[t]=1;
      }
  }

  if (any_disagree) {
    cout << "ReadTimeSeriesFile:  Warning -- Error read and calculated disagree at the following time slices:";
    for (int t=0; t<=Nt; t++)
      if (disagree1[t])
	cout << " " << t;
    cout << "\n";
    cout << "ReadTimeSeriesFile:  Warning -- Error read and calculated disagree because error read was 0 at the following time slices:";
    for (int t=0; t<=Nt; t++)
      if (disagree2[t])
	cout << " " << t;
    cout << "\n";
  }
  
  return *this;
}

//Read in a Time Series from ave_file and its jackknife values from jks_file,
//and store the information in arrays ave, jk, is_defined, and errors.
//The arrays ave, jk, and is_defined can be put in JackknifeTimeSeries::FromArray
//to build an actual JackknifeTimeSeries object.  The errors array is
//just for checking.  The expected number of time slices (Nt) and jackknife values
//(N) in the files must be specified (they are inputs to the function).
//Not a member or friend function.
int ReadTimeSeriesFiletoArray(string ave_file, string jks_file, int Nt, int N, double ave [], double** jk, int is_defined [], double* errors)
{
  if (N<2) {
    //Error: I think I'd rather have it abort here.
    cout << "ReadTimeSeriesFiletoArray: N must be >= 2 for a JackknifeTimeSeries object, but was given N=" << N << ".\n";
    return 0;
  }
  if (Nt<1) {
    //Error: I think I'd rather have it abort here.
    cout << "ReadTimeSeriesFiletoArray: Nt must be >= 1 for a JackknifeTimeSeries object, but was given Nt=" << Nt << ".\n";
    return 0;
  }

  //Sizes of arrays assumed:
  //double ave[Nt+1]
  //double jk[Nt+1][N]
  //int is_defined[Nt+1]
  //double errors[Nt+1] (if it is input at all)

  for (int t=0; t<=Nt; t++)
    is_defined[t]=0;
  ifstream fin(ave_file.c_str());
  int t_last=-1, t, n_defined=0;
  string line;
  while (getline(fin,line)) {
    string vals[20];
    int nvals=split(line,vals,20);
    if (nvals!=3) {
      //Error: I think I'd rather have it abort here.
      cout << "ReadTimeSeriesFiletoArray: Should be three numbers in each column but read " << nvals << ".\n";
      return 0;
    }
    t=atoi(vals[0].c_str());
    if (t<=t_last) {
      //Error: I think I'd rather have it abort here.
      cout << "ReadTimeSeriesFiletoArray: Times should be in increasing order.\n";
      return 0;
    }
    if (t>Nt) {
      //Error: I think I'd rather have it abort here.
      cout << "ReadTimeSeriesFiletoArray: t must be <= Nt, but read t=" << t << ".\n";
      return 0;
    }
    ave[t]=atof(vals[1].c_str());
    is_defined[t]=1;
    if (errors!=0)
      errors[t]=atof(vals[2].c_str());
    n_defined++;
    t_last=t;
  }
  fin.close();
  ifstream fin2(jks_file.c_str());
  for (int config_no=0; config_no<N; config_no++) {
    for (t=0; t<=Nt; t++) {
      if (is_defined[t]) {
	if (!getline(fin2,line)) {
	  //Error: I think I'd rather have it abort here.
	  cout << "ReadTimeSeriesFiletoArray: Too few lines in jks file.\n";
	  return 0;
	}
	string vals[20];
	int nvals=split(line,vals,20);
	if (nvals!=2) {
	  //Error: I think I'd rather have it abort here.
	  cout << "ReadTimeSeriesFiletoArray: Wrong number of columns in jks file.\n";
	  return 0;
	}
	int t_read=atoi(vals[0].c_str());
	if (t_read != t) {
	  //Error: I think I'd rather have it abort here.
	  cout << "ReadTimeSeriesFiletoArray: Data for configuration " << config_no << " in wrong order.\n";
	  return 0;
	}
	jk[t][config_no]=atof(vals[1].c_str());
      }
    }
  }
  if (getline(fin2,line)) {
    //Error: I think I'd rather have it abort here.
    cout << "ReadTimeSeriesFiletoArray: Too many lines in jks file.\n";
    return 0;
  }
  
  return 1;
}




//JackknifeCmplxTimeSeries class

//Standard and Default Constructor
JackknifeCmplxTimeSeries::JackknifeCmplxTimeSeries(int nt, int n)
{
  if (n<2) {
    //Error: I think I'd rather have it abort here.
    cout << "N must be >= 2 for a JackknifeCmplxTimeSeries object, but was given N=" << n << ".  Setting N=2.\n";
    N=2;
  } else {
    N=n;
  }
  if (nt<1) {
    //Error: I think I'd rather have it abort here.
    cout << "Nt must be >= 1 for a JackknifeCmplxTimeSeries object, but was given Nt=" << nt << ".  Setting Nt=1.\n";
    Nt=1;
  } else {
    Nt=nt;
  }
  Jk=new JackknifeCmplx [Nt+1];
  is_defined=new int [Nt+1];
  for (int t=0; t<=Nt; t++) {
    is_defined[t]=0;
    //Have to initialize a JackknifeCmplx object of size N for each Jk[t]
    JackknifeCmplx tmp(N);
    Jk[t]=tmp;
  }
}

//Copy Constructor
JackknifeCmplxTimeSeries::JackknifeCmplxTimeSeries(const JackknifeCmplxTimeSeries & JT)
{
  N=JT.N;
  Nt=JT.Nt;
  Jk=new JackknifeCmplx [Nt+1];
  is_defined=new int [Nt+1];
  for (int t=0; t<=Nt; t++) {
    Jk[t]=JT.Jk[t];
    is_defined[t]=JT.is_defined[t];
  }
}

//Constructor with a JackknifeCmplx object
JackknifeCmplxTimeSeries::JackknifeCmplxTimeSeries(int nt, const JackknifeCmplx & J)
{
  if (nt<1) {
    //Error: I think I'd rather have it abort here.
    cout << "Nt must be >= 1 for a JackknifeCmplxTimeSeries object, but was given Nt=" << nt << ".  Setting Nt=1.\n";
    Nt=1;
  } else {
    Nt=nt;
  }
  N=J.ReturnN();
  Jk=new JackknifeCmplx [Nt+1];
  is_defined=new int [Nt+1];
  for (int t=0; t<=Nt; t++) {
    Jk[t]=J;
    is_defined[t]=1;
  }
}

//Destructor
JackknifeCmplxTimeSeries::~JackknifeCmplxTimeSeries() 
{
  delete [] Jk;
  delete [] is_defined;
}

//Assignment operator
//Deep copy
JackknifeCmplxTimeSeries & JackknifeCmplxTimeSeries::operator=(const JackknifeCmplxTimeSeries & JT)
{
  if (this == &JT)
    return *this;
  delete [] Jk;
  delete [] is_defined;
  N=JT.N;
  Nt=JT.Nt;
  Jk=new JackknifeCmplx [Nt+1];
  is_defined=new int [Nt+1];
  for (int t=0; t<=Nt; t++) {
    Jk[t]=JT.Jk[t];
    is_defined[t]=JT.is_defined[t];
  }
  return *this;
}

//Say whether or not there is data at t
int JackknifeCmplxTimeSeries::ReturnIsDefined(int t) const
{
  if (t<0 || t>Nt)
    return 0;
  else
    return is_defined[t];
}

//Return the average at t
complex<double> JackknifeCmplxTimeSeries::ReturnAve(int t) const
{
  if (t<0 || t>Nt || !is_defined[t]) {
    //Error:  I think I'd rather have it abort here.
    cout << "No data at t=" << t << " in JackknifeCmplxTimeSeries.\n";
    return 0.0;
  } else {
    return Jk[t].ReturnAve();
  }
}

//Return the error at t
double JackknifeCmplxTimeSeries::ReturnErr(int t) const
{
  if (t<0 || t>Nt || !is_defined[t]) {
    //Error:  I think I'd rather have it abort here.
    cout << "No data at t=" << t << " in JackknifeCmplxTimeSeries.\n";
    return 0.0;
  } else {
    return Jk[t].ReturnErr();
  }
}

//Return the ith jackknife value at t
complex<double> JackknifeCmplxTimeSeries::ReturnJk(int t, int i) const
{
  if (t<0 || t>Nt || !is_defined[t]) {
    //Error:  I think I'd rather have it abort here.
    cout << "No data at t=" << t << " in JackknifeCmplxTimeSeries.\n";
    return 0.0;
  }
  if (i<0 || i>=N) {
    //Error:  I think I'd rather have it abort here.
    cout << "i must be 0<=i<N, but received i=" << i << ".\n";
    return 0.0;
  }
  return Jk[t].ReturnJk(i);
}

//Return Nt
int JackknifeCmplxTimeSeries::ReturnNt() const
{
  return Nt;
}

//Return N
int JackknifeCmplxTimeSeries::ReturnN() const
{
  return N;
}

//Addition
JackknifeCmplxTimeSeries JackknifeCmplxTimeSeries::operator+(const JackknifeCmplxTimeSeries & JT) const
{
  JackknifeCmplxTimeSeries sum(*this);
  sum+=JT;
  return sum;
}

//+=
JackknifeCmplxTimeSeries & JackknifeCmplxTimeSeries::operator+=(const JackknifeCmplxTimeSeries & JT)
{
  if (N!=JT.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't add JackknifeCmplxTimeSeries objects with different number of jackknife values.  Object 1 N=" << N << " Object 2 N=" << JT.N << "\n";
  } else if (Nt!=JT.Nt) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't add JackknifeCmplxTimeSeries objects with different Nt.  Object 1 Nt=" << Nt << " Object 2 Nt=" << JT.Nt << "\n";
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
JackknifeCmplxTimeSeries JackknifeCmplxTimeSeries::operator-(const JackknifeCmplxTimeSeries & JT) const
{
  JackknifeCmplxTimeSeries difference(*this);
  difference-=JT;
  return difference;
}

//Negative (Unary operator)
JackknifeCmplxTimeSeries JackknifeCmplxTimeSeries::operator-() const
{
  JackknifeCmplx tmp(N,0.0);
  JackknifeCmplxTimeSeries neg(Nt,tmp);
  neg-=*this;
  return neg;
}

//-=
JackknifeCmplxTimeSeries & JackknifeCmplxTimeSeries::operator-=(const JackknifeCmplxTimeSeries & JT)
{
  if (N!=JT.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't subtract JackknifeCmplxTimeSeries objects with different number of jackknife values.  Object 1 N=" << N << " Object 2 N=" << JT.N << "\n";
  } else if (Nt!=JT.Nt) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't subtract JackknifeCmplxTimeSeries objects with different Nt.  Object 1 Nt=" << Nt << " Object 2 Nt=" << JT.Nt << "\n";
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

//Multiplication of two JackknifeCmplxTimeSeries objects
JackknifeCmplxTimeSeries JackknifeCmplxTimeSeries::operator*(const JackknifeCmplxTimeSeries & JT) const
{
  JackknifeCmplxTimeSeries prod(*this);
  prod*=JT;
  return prod;
}

//Multiplication of a JackknifeCmplxTimeSeries object by a complex scalar on the right
JackknifeCmplxTimeSeries JackknifeCmplxTimeSeries::operator*(complex<double> c) const
{
  JackknifeCmplxTimeSeries prod(*this);
  prod*=c;
  return prod;
}

//Multiplication by a complex scalar on the left (friend function)
JackknifeCmplxTimeSeries operator*(complex<double> c, const JackknifeCmplxTimeSeries & JT)
{
  return JT*c;
}

//*=JackknifeCmplxTimeSeries object
JackknifeCmplxTimeSeries & JackknifeCmplxTimeSeries::operator*=(const JackknifeCmplxTimeSeries & JT)
{
  if (N!=JT.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't multiply JackknifeCmplxTimeSeries objects with different number of jackknife values.  Object 1 N=" << N << " Object 2 N=" << JT.N << "\n";
  } else if (Nt!=JT.Nt) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't multiply JackknifeTimeCmplxSeries objects with different Nt.  Object 1 Nt=" << Nt << " Object 2 Nt=" << JT.Nt << "\n";
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
JackknifeCmplxTimeSeries & JackknifeCmplxTimeSeries::operator*=(complex<double> c)
{
  for (int t=0; t<=Nt; t++)
    if (is_defined[t])
      Jk[t]*=c;
  return *this;
}

//Division of two JackknifeCmplxTimeSeries objects
JackknifeCmplxTimeSeries JackknifeCmplxTimeSeries::operator/(const JackknifeCmplxTimeSeries & JT) const
{
  JackknifeCmplxTimeSeries q(*this);
  q/=JT;
  return q;
}

//Division by a complex scalar
JackknifeCmplxTimeSeries JackknifeCmplxTimeSeries::operator/(complex<double> c) const
{
  JackknifeCmplxTimeSeries q(*this);
  q/=c;
  return q;
}

//Division of a complex scalar by a JackknifeCmplxTimeSeries object (friend function)
JackknifeCmplxTimeSeries operator/(complex<double> c, const JackknifeCmplxTimeSeries & JT)
{
  int Nt=JT.Nt;
  int N=JT.N;
  JackknifeCmplx tmp(N,c);
  JackknifeCmplxTimeSeries q(Nt,tmp);
  q/=JT;
  return q;
}

//slash equals JackknifeCmplxTimeSeries object
JackknifeCmplxTimeSeries & JackknifeCmplxTimeSeries::operator/=(const JackknifeCmplxTimeSeries & JT)
{
  if (N!=JT.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't divide JackknifeCmplxTimeSeries objects with different number of jackknife values.  Object 1 N=" << N << " Object 2 N=" << JT.N << "\n";
  } else if (Nt!=JT.Nt) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't divide JackknifeCmplxTimeSeries objects with different Nt.  Object 1 Nt=" << Nt << " Object 2 Nt=" << JT.Nt << "\n";
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
JackknifeCmplxTimeSeries & JackknifeCmplxTimeSeries::operator/=(complex<double> c)
{
  for (int t=0; t<=Nt; t++)
    if (is_defined[t])
      Jk[t]/=c;
  return *this;
}

JackknifeCmplxTimeSeries & JackknifeCmplxTimeSeries::Reverse()
{
  JackknifeCmplxTimeSeries tmp(*this);
  for (int t=0; t<=Nt; t++) {
    Jk[Nt-t]=tmp.Jk[t];
    is_defined[Nt-t]=tmp.is_defined[t];
  }
  return *this;
}

JackknifeCmplxTimeSeries Reverse(const JackknifeCmplxTimeSeries & JT)
{
  JackknifeCmplxTimeSeries tmp(JT);
  tmp.Reverse();
  return tmp;
}

JackknifeCmplxTimeSeries & JackknifeCmplxTimeSeries::Translate(int t)
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

JackknifeCmplxTimeSeries Translate(int t, const JackknifeCmplxTimeSeries & JT)
{
  JackknifeCmplxTimeSeries tmp(JT);
  tmp.Translate(t);
  return tmp;
}

JackknifeCmplxTimeSeries & JackknifeCmplxTimeSeries::TranslateReverse(int t)
{
  if (t>=0 && t<=Nt) {
    JackknifeCmplxTimeSeries tmp(*this);
    for (int tt=0; tt<=t; tt++) {
      Jk[t-tt]=tmp.Jk[tt];
      is_defined[t-tt]=tmp.is_defined[tt];
    }
    for (int tt=t+1; tt<=Nt; tt++)
      is_defined[tt]=0;
  } else if (t>Nt && t<=2*Nt) {
    JackknifeCmplxTimeSeries tmp(*this);
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

JackknifeCmplxTimeSeries TranslateReverse(int t, const JackknifeCmplxTimeSeries & JT)
{
  JackknifeCmplxTimeSeries tmp(JT);
  tmp.TranslateReverse(t);
  return tmp;
}

//Abs (member function)
JackknifeCmplxTimeSeries & JackknifeCmplxTimeSeries::Abs()
{
  for (int t=0; t<=Nt; t++)
    if (is_defined[t]) {
      JackknifeCmplx tmp(Jk[t]);
      tmp.Abs(); //Have to do it this way or else it gets confused
                 //with JackknifeCmplxTimeSeries::Abs
      Jk[t]=tmp;
    }
  return *this;
}

//Real (member function)
JackknifeCmplxTimeSeries & JackknifeCmplxTimeSeries::Real()
{
  for (int t=0; t<=Nt; t++)
    if (is_defined[t]) {
      JackknifeCmplx tmp(Jk[t]);
      tmp.Real(); //Have to do it this way or else it gets confused
                  //with JackknifeCmplxTimeSeries::Real
      Jk[t]=tmp;
    }
  return *this;
}

//Imag (member function)
JackknifeCmplxTimeSeries & JackknifeCmplxTimeSeries::Imag()
{
  for (int t=0; t<=Nt; t++)
    if (is_defined[t]) {
      JackknifeCmplx tmp(Jk[t]);
      tmp.Imag(); //Have to do it this way or else it gets confused
                  //with JackknifeCmplxTimeSeries::Real
      Jk[t]=tmp;
    }
  return *this;
}

//Conj (member function)
JackknifeCmplxTimeSeries & JackknifeCmplxTimeSeries::Conj()
{
  for (int t=0; t<=Nt; t++)
    if (is_defined[t]) {
      JackknifeCmplx tmp(Jk[t]);
      tmp.Conj(); //Have to do it this way or else it gets confused
                  //with JackknifeCmplxTimeSeries::Abs
      Jk[t]=tmp;
    }
  return *this;
}

//Conj (friend function)
JackknifeCmplxTimeSeries Conj(const JackknifeCmplxTimeSeries & JT)
{
  JackknifeCmplxTimeSeries tmpc(JT);
  tmpc.Conj();
  return tmpc;
}

//Normalize (member function)
//Divides everything by the largest modulus of the average values on
//the time slices.
//Returns the factor by which the object was divided.
//This is needed because Sam's fitter has problems with large numbers.
double JackknifeCmplxTimeSeries::Normalize()
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
//Takes two JackknifeCmplxTimeSeries objects with N jaccknife values and puts
//them together to make another JackknifeCmplxTimeSeries object with 2*N
//jackknife values.
//IMPORTANT: NOTE THAT it can only be used on quantities that depend linearly
//on the data
JackknifeCmplxTimeSeries Combine(const JackknifeCmplxTimeSeries & JT1, const JackknifeCmplxTimeSeries & JT2)
{
  if (JT1.N != JT2.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can only combine JackknifeCmplxTimeSeries objects with same number of jackknife values.  Object 1 N=" << JT1.N << " Object 2 N=" << JT2.N << "\n";
    return JT1;
  }
  if (JT1.Nt != JT2.Nt) {
    //Error: I think I'd rather have it abort here.
    cout << "Can only combine JackknifeCmplxTimeSeries objects with same Nt.  Object 1 Nt=" << JT1.Nt << " Object 2 Nt=" << JT2.Nt << "\n";
    return JT1;
  }
  int N=JT1.N;
  int Nt=JT1.Nt;
  JackknifeCmplxTimeSeries tmp(Nt,2*N);
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

//Output jackknife values to a file in the appropriate format for Sam's
//fitter (friend function)
void OutputJk(string filename, const JackknifeCmplxTimeSeries & JT)
{
  ofstream fout(filename.c_str());
  fout.precision(16);
  int Nt=JT.ReturnNt();
  int N=JT.ReturnN();
  for (int i=0; i<N; i++)
    for (int t=0; t<=Nt; t++)
      if (JT.ReturnIsDefined(t))
	fout << t << " " << real(JT.Jk[t].jk[i]) << " " << imag(JT.Jk[t].jk[i]) << "\n";
  fout.close();
}

//Output average and errors to a file
//not a member or friend function
void OutputAve(string filename, const JackknifeCmplxTimeSeries & JT)
{
  ofstream fout(filename.c_str());
  fout.precision(16);
  int Nt=JT.ReturnNt();
  for (int t=0; t<=Nt; t++)
    if (JT.ReturnIsDefined(t))
      fout << t << " " << real(JT.ReturnAve(t)) << " " << imag(JT.ReturnAve(t)) << " " << JT.ReturnErr(t) << "\n";
  fout.close();
}

