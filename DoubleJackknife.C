#include <iostream>
#include <cmath>
#include <complex>
using namespace std;
#include "jackknife.h"
#include "double_jackknife.h"


//DoubleJackknife class

//Standard and Default Constructor
DoubleJackknife::DoubleJackknife(int n)
{
  if (n<3) {
    //Error: I think I'd rather have it abort here.
    cout << "N must be >= 3 for a DoubleJackknife object, but was given N=" << n << ".  Setting N=3.\n";
    N=3;
  } else {
    N=n;
  }
  jk=new Jackknife [N];
  for (int i=0; i<N; i++) {
    Jackknife tmp(N-1);
    jk[i]=tmp;
  }
}

//Copy Constructor
DoubleJackknife::DoubleJackknife(const DoubleJackknife & J)
{
  N=J.N;
  jk=new Jackknife [N];
  for (int i=0; i<N; i++)
    jk[i]=J.jk[i];
  ave=J.ave;
  err=J.err;
}

//Constructor with a scalar
//Sets all jackknife values equal to the input scalar d.
DoubleJackknife::DoubleJackknife(int n, double d)
{
  if (n<3) {
    //Error: I think I'd rather have it abort here.
    cout << "N must be >= 3 for a DoubleJackknife object, but was given N=" << n << ".  Setting N=3.\n";
    N=3;
  } else {
    N=n;
  }
  jk=new Jackknife [N];
  for (int i=0; i<N; i++) {
    Jackknife tmp(N-1,d);
    jk[i]=tmp;
  }
  ave=d;
  CalcAll();
}

//Destructor
DoubleJackknife::~DoubleJackknife()
{
  delete [] jk;
}

//Assignment operator
//Copies the jk array from object on the right hand side into a new block of
//memory which the object on the left hand side points to.  Also, the
//memory previously used by the jk array of the object on the left hand
//side is deleted.
DoubleJackknife & DoubleJackknife::operator=(const DoubleJackknife & J)
{
  if (this == &J)
    return *this;
  delete [] jk;
  N=J.N;
  jk=new Jackknife [N];
  for (int i=0; i<N; i++)
    jk[i]=J.jk[i];
  ave=J.ave;
  err=J.err;
  return *this;
}

//Calculate the jackknife error and store in err
DoubleJackknife & DoubleJackknife::CalcErr()
{
  err=0;
  for (int i=0; i<N; i++)
    err+=(jk[i].ReturnAve()-ave)*(jk[i].ReturnAve()-ave);
  err*=((double) N)/((double) (N-1));
  err=sqrt(err);
  return *this;
}

//Calculate everything that needs to be recalculated when a function is applied
//to a Jackknife object, or when two Jackknife objects are combined via a
//function.  For now this is just the error, but if in the future I add in
//things like the average of the jackknife values then this would also have
//to be recalculated.
DoubleJackknife & DoubleJackknife::CalcAll()
{
  CalcErr();
  return *this;
}

//Return the average
double DoubleJackknife::ReturnAve() const
{
  return ave;
}

//Return the error
double DoubleJackknife::ReturnErr() const
{
  return err;
}

//Return the ith jackknife value
double DoubleJackknife::ReturnJk(int i) const
{
  if (i<0 || i>=N) {
    //Error: I think I'd rather have it abort here.
    cout << "i must be 0<=i<N, but received i=" << i << ".  Returning first jackknife value instead.\n";
    i=0;
  }
  return jk[i].ReturnAve();
}

//Return the error on the ith jackknife value
double DoubleJackknife::ReturnJkErr(int i) const
{
  if (i<0 || i>=N) {
    //Error: I think I'd rather have it abort here.
    cout << "i must be 0<=i<N, but received i=" << i << ".  Returning error on first jackknife value instead.\n";
    i=0;
  }
  return jk[i].ReturnErr();
}

//Return double jackknife value i, j.
double DoubleJackknife::ReturnDoubleJk(int i, int j) const
{
  if (i<0 || i>=N) {
    //Error: I think I'd rather have it abort here.
    cout << "i must be 0<=i<N, but received i=" << i << ". Setting i=0 instead.\n";
    i=0;
  }
  if (j<0 || j>=N-1) {
    //Error: I think I'd rather have it abort here.
    cout << "j must be 0<=j<N-1, but received j=" << j << ".  Setting j=0 instead.\n";
    j=0;
  }
  return jk[i].ReturnJk(j);
}

//Return N
int DoubleJackknife::ReturnN() const
{
  return N;
}

//Unjackknife undoes the jackknife and returns the original data point
//Only valid for quantities that depend linearly on the original data
double DoubleJackknife::Unjackknife(int i) const
{
  if ( i<0 || i>=N ) {
    //Error: I think I'd rather have it abort here.
    cout << "Data only exists from DoubleJackknife object for 0<=i<" << N << " but was given i=" << i << ".\n";
    return 0.0;
  }
  double unjk=0.0;
  for (int j=0; j<N; j++)
    if (j!=i)
      unjk+=jk[j].ReturnAve();
  unjk-=((double) N-2)*jk[i].ReturnAve();
  return unjk;
}

//Takes a Jackknife object, undoes it, then calculates double jackknife
//values to put into the current DoubleJackknife object.
//NOTE: ONLY WORKS for quantities that depend linearly on original data.
DoubleJackknife & DoubleJackknife::FromSingleJk(const Jackknife & J)
{
  delete [] jk;
  N=J.N;
  jk=new Jackknife [N];
  double unjks[N];
  double sum=0.0;
  for (int i=0; i<N; i++) {
    unjks[i]=J.Unjackknife(i);
    sum+=unjks[i];
  }
  for (int i=0; i<N; i++) {
    Jackknife tmp(N-1);
    jk[i]=tmp;
    jk[i].ave=J.ReturnJk(i);
    for (int j=0; j<i; j++)
      jk[i].jk[j]=(sum-unjks[i]-unjks[j])/double(N-2);
    for (int j=i+1; j<N; j++)
      jk[i].jk[j-1]=(sum-unjks[i]-unjks[j])/double(N-2);
    jk[i].CalcAll();
  }
  ave=J.ave;
  err=J.err;
  return *this;
}

//Friend function version of the above.
DoubleJackknife FromSingleJk(const Jackknife & J)
{
  DoubleJackknife tmp;
  tmp.FromSingleJk(J);
  return tmp;
}

//Addition
DoubleJackknife DoubleJackknife::operator+(const DoubleJackknife & J) const
{
  DoubleJackknife sum(*this);
  sum+=J;
  return sum;
}

//+=
DoubleJackknife & DoubleJackknife::operator+=(const DoubleJackknife & J)
{
  if (N!=J.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't add DoubleJackknife objects with different number of jackknife values.  Object 1 N=" << N << " Object 2 N=" << J.N << "\n";
  } else {
    for (int i=0; i<N; i++)
      jk[i]+=J.jk[i];
    ave+=J.ave;
    CalcAll();
  }
  return *this;
}

//Subtraction
DoubleJackknife DoubleJackknife::operator-(const DoubleJackknife & J) const
{
  DoubleJackknife difference(*this);
  difference-=J;
  return difference;
}

//Negative (Unary operator)
DoubleJackknife DoubleJackknife::operator-() const
{
  DoubleJackknife neg(N,0.0);
  neg-=*this;
  return neg;
}

//-=
DoubleJackknife & DoubleJackknife::operator-=(const DoubleJackknife & J)
{
  if (N!=J.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't subtract DoubleJackknife objects with different number of jackknife values.  Object 1 N=" << N << " Object 2 N=" << J.N << "\n";
  } else {
    for (int i=0; i<N; i++)
      jk[i]-=J.jk[i];
    ave-=J.ave;
    CalcAll();
  }
  return *this;
}

//Multiplication of two DoubleJackknife objects
DoubleJackknife DoubleJackknife::operator*(const DoubleJackknife & J) const
{
  DoubleJackknife prod(*this);
  prod*=J;
  return prod;
}

//Multiplication by a scalar on the right
DoubleJackknife DoubleJackknife::operator*(double d) const
{
  DoubleJackknife prod(*this);
  prod*=d;
  return prod;
}

//Multiplication by a scalar on the left (friend function)
DoubleJackknife operator*(double d, const DoubleJackknife & J)
{
  return J*d;
}

//*=DoubleJackknife object
DoubleJackknife & DoubleJackknife::operator*=(const DoubleJackknife & J)
{
  if (N!=J.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't multiply DoubleJackknife objects with different number of jackknife values.  Object 1 N=" << N << " Object 2 N=" << J.N << "\n";
  } else {
    for (int i=0; i<N; i++)
      jk[i]*=J.jk[i];
    ave*=J.ave;
    CalcAll();
  }
  return *this;
}

//*=scalar
DoubleJackknife & DoubleJackknife::operator*=(double d)
{
  for (int i=0; i<N; i++)
    jk[i]*=d;
  ave*=d;
  CalcAll();
  return *this;
}

//Division of two DoubleJackknife objects
DoubleJackknife DoubleJackknife::operator/(const DoubleJackknife & J) const
{
  DoubleJackknife q(*this);
  q/=J;
  return q;
}

//Division by a scalar
DoubleJackknife DoubleJackknife::operator/(double d) const
{
  DoubleJackknife q(*this);
  q/=d;
  return q;
}

//Division of scalar by a DoubleJackknife object (friend function)
DoubleJackknife operator/(double d, const DoubleJackknife & J)
{
  int N=J.N;
  DoubleJackknife q(N,d);
  q/=J;
  return q;
}

//slash equals DoubleJackknife object
DoubleJackknife & DoubleJackknife::operator/=(const DoubleJackknife & J)
{
  if (N!=J.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't divide DoubleJackknife objects with different number of jackknife values.  Object 1 N=" << N << " Object 2 N=" << J.N << "\n";
  } else {
    for (int i=0; i<N; i++)
      jk[i]/=J.jk[i];
    ave/=J.ave;
    CalcAll();
  }
  return *this;
}

//slash equals scalar
DoubleJackknife & DoubleJackknife::operator/=(double d)
{
  for (int i=0; i<N; i++)
    jk[i]/=d;
  ave/=d;
  CalcAll();
  return *this;
}

//Log (member function)
DoubleJackknife & DoubleJackknife::Log()
{
  for (int i=0; i<N; i++) {
    jk[i].Log();
  }
  ave=log(ave);
  CalcAll();
  return *this;
}

//Log (friend function)
DoubleJackknife Log(const DoubleJackknife & J)
{
  DoubleJackknife tmp(J);
  tmp.Log();
  return tmp;
}

//Sqrt (member function)
DoubleJackknife & DoubleJackknife::Sqrt()
{
  for (int i=0; i<N; i++) {
    jk[i].Sqrt();
  }
  ave=sqrt(ave);
  CalcAll();
  return *this;
}

//Sqrt (friend function)
DoubleJackknife Sqrt(const DoubleJackknife & J)
{
  DoubleJackknife tmp(J);
  tmp.Sqrt();
  return tmp;
}

//Abs (member function)
DoubleJackknife & DoubleJackknife::Abs()
{
  for (int i=0; i<N; i++) {
    jk[i].Abs();
  }
  ave=abs(ave);
  CalcAll();
  return *this;
}

//Abs (friend function)
DoubleJackknife Abs(const DoubleJackknife & J)
{
  DoubleJackknife tmp(J);
  tmp.Abs();
  return tmp;
}

//Abs (friend function with a DoubleJackknifeCmplx argument)
DoubleJackknife Abs(const DoubleJackknifeCmplx & J)
{
  DoubleJackknife tmp(J.N);
  for (int i=0; i<tmp.N; i++)
    tmp.jk[i]=Abs(J.jk[i]);
  tmp.ave=abs(J.ave);
  tmp.CalcAll();
  return tmp;
}

//Real (friend function)
DoubleJackknife Real(const DoubleJackknifeCmplx & J)
{
  DoubleJackknife tmp(J.N);
  for (int i=0; i<tmp.N; i++)
    tmp.jk[i]=Real(J.jk[i]);
  tmp.ave=real(J.ave);
  tmp.CalcAll();
  return tmp;
}

//Imag (friend function)
DoubleJackknife Imag(const DoubleJackknifeCmplx & J)
{
  DoubleJackknife tmp(J.N);
  for (int i=0; i<tmp.N; i++)
    tmp.jk[i]=Imag(J.jk[i]);
  tmp.ave=imag(J.ave);
  tmp.CalcAll();
  return tmp;
}

//Combine (friend function)
//Takes two DoubleJackknife objects with N jaccknife values and puts them
//together to make another DoubleJackknife object with 2*N jackknife values.
//IMPORTANT: NOTE THAT it can only be used on quantities that depend linearly
//on the data
DoubleJackknife Combine(const DoubleJackknife & J1, const DoubleJackknife & J2)
{
  if (J1.N != J2.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can only combine DoubleJackknife objects with same number of jackknife values.  Object 1 N=" << J1.N << " Object 2 N=" << J2.N << "\n";
    return J1;
  }
  int N=J1.N;
  DoubleJackknife tmp(2*N);
  double unjk[2*N];
  double sum=0.0;
  for (int i=0; i<N; i++) {
    unjk[i]=J1.Unjackknife(i);
    unjk[i+N]=J2.Unjackknife(i);
    sum+=unjk[i]+unjk[i+N];
  }
  for (int i=0; i<2*N; i++) {
    tmp.jk[i].ave=(sum-unjk[i])/double(2*N-1);
    for (int j=0; j<i; j++)
      tmp.jk[i].jk[j]=(sum-unjk[i]-unjk[j])/double(2*N-2);
    for (int j=i+1; j<2*N; j++)
      tmp.jk[i].jk[j-1]=(sum-unjk[i]-unjk[j])/double(2*N-2);
    tmp.jk[i].CalcAll();
  }
  tmp.ave=0.5*(J1.ave+J2.ave);
  tmp.CalcAll();
  return tmp;
}



//DoubleJackknifeCmplx class

//Standard and Default Constructor
DoubleJackknifeCmplx::DoubleJackknifeCmplx(int n)
{
  if (n<3) {
    //Error: I think I'd rather have it abort here.
    cout << "N must be >= 3 for DoubleJackknifeCmplx object, but was given N=" << n << ".  Setting N=3.\n";
    N=3;
  } else {
    N=n;
  }
  jk=new JackknifeCmplx [N];
  for (int i=0; i<N; i++) {
    JackknifeCmplx tmp(N-1);
    jk[i]=tmp;
  }
}

//Copy Constructor
DoubleJackknifeCmplx::DoubleJackknifeCmplx(const DoubleJackknifeCmplx & J)
{
  N=J.N;
  jk=new JackknifeCmplx [N];
  for (int i=0; i<N; i++)
    jk[i]=J.jk[i];
  ave=J.ave;
  err=J.err;
}

//Constructor with a complex scalar
//Sets all jackknife values equal to the input scalar c.
//Do I need a separate one for a real scalar d?  Or will it convert a double
//automatically to a complex<double>?
DoubleJackknifeCmplx::DoubleJackknifeCmplx(int n, complex<double> c)
{
  if (n<3) {
    //Error: I think I'd rather have it abort here.
    cout << "N must be >= 3 for a DoubleJackknifeCmplx object, but was given N=" << n << ".  Setting N=3.\n";
    N=3;
  } else {
    N=n;
  }
  jk=new JackknifeCmplx [N];
  for (int i=0; i<N; i++) {
    JackknifeCmplx tmp(N-1,c);
    jk[i]=tmp;
  }
  ave=c;
  CalcAll();
}

//Destructor
DoubleJackknifeCmplx::~DoubleJackknifeCmplx()
{
  delete [] jk;
}

//Assignment operator
//Copies the jk array from object on the right hand side into a new block of
//memory which the object on the left hand side points to.  Also, the
//memory previously used by the jk array of the object on the left hand
//side is deleted.
DoubleJackknifeCmplx & DoubleJackknifeCmplx::operator=(const DoubleJackknifeCmplx & J)
{
  if (this == &J)
    return *this;
  delete [] jk;
  N=J.N;
  jk=new JackknifeCmplx [N];
  for (int i=0; i<N; i++)
    jk[i]=J.jk[i];
  ave=J.ave;
  err=J.err;
  return *this;
}

//Calculate the jackknife error and store in err
DoubleJackknifeCmplx & DoubleJackknifeCmplx::CalcErr()
{
  err=0;
  for (int i=0; i<N; i++)
    err+=real(conj(jk[i].ReturnAve()-ave)*(jk[i].ReturnAve()-ave)); //real just converts to a double,
                                                    //since this is the modulus squared
                                                    //it will already be a real number.
  err*=((double) N)/((double) (N-1));
  err=sqrt(err);
  return *this;
}

//Calculate everything that needs to be recalculated when a function is applied
//to a JackknifeCmplx object, or when two JackknifeCmplx objects are combined
//via a function.  For now this is just the error, but if in the future I add
//in things like the average of the jackknife values then this would also have
//to be recalculated.
DoubleJackknifeCmplx & DoubleJackknifeCmplx::CalcAll()
{
  CalcErr();
  return *this;
}

//Return the average
complex<double> DoubleJackknifeCmplx::ReturnAve() const
{
  return ave;
}

//Return the error
double DoubleJackknifeCmplx::ReturnErr() const
{
  return err;
}

//Return the ith jackknife value
complex<double> DoubleJackknifeCmplx::ReturnJk(int i) const
{
  if (i<0 || i>=N) {
    //Error: I think I'd rather have it abort here.
    cout << "i must be 0<=i<N, but received i=" << i << ".  Returning first jackknife value instead.\n";
    i=0;
  }
  return jk[i].ReturnAve();
}

//Return the error on the ith jackknife value
double DoubleJackknifeCmplx::ReturnJkErr(int i) const
{
  if (i<0 || i>=N) {
    //Error: I think I'd rather have it abort here.
    cout << "i must be 0<=i<N, but received i=" << i << ".  Returning first jackknife value instead.\n";
    i=0;
  }
  return jk[i].ReturnErr();
}

//Return double jackknife value i, j.
complex<double> DoubleJackknifeCmplx::ReturnDoubleJk(int i, int j) const
{
  if (i<0 || i>=N) {
    //Error: I think I'd rather have it abort here.
    cout << "i must be 0<=i<N, but received i=" << i << ".  Setting i=0 instead.\n";
    i=0;
  }
  if (j<0 || j>=N-1) {
    //Error: I think I'd rather have it abort here.
    cout << "j must be 0<=j<N-1, but received j=" << j << ".  Setting j=0 instead.\n";
    j=0;
  }
  return jk[i].ReturnJk(j);
}

//Return N
int DoubleJackknifeCmplx::ReturnN() const
{
  return N;
}

//Unjackknife undoes the jackknife and returns the original data point
//Only valid for quantities that depend linearly on the original data
complex<double> DoubleJackknifeCmplx::Unjackknife(int i) const
{
  if ( i<0 || i>=N ) {
    //Error: I think I'd rather have it abort here.
    cout << "Data only exists from DoubleJackknifeCmplx object for 0<=i<" << N << " but was given i=" << i << ".\n";
    return 0.0;
  }
  complex<double> unjk=0.0;
  for (int j=0; j<N; j++)
    if (j!=i)
      unjk+=jk[j].ReturnAve();
  unjk-=((double) N-2)*jk[i].ReturnAve();
  return unjk;
}

//Takes a JackknifeCmplx object, undoes it, then calculates double jackknife
//values to put into the current DoubleJackknifeCmplx object.
//NOTE: ONLY WORKS for quantities that depend linearly on original data.
DoubleJackknifeCmplx & DoubleJackknifeCmplx::FromSingleJk(const JackknifeCmplx & J)
{
  delete [] jk;
  N=J.N;
  jk=new JackknifeCmplx [N];
  complex<double> unjks[N];
  complex<double> sum=0.0;
  for (int i=0; i<N; i++) {
    unjks[i]=J.Unjackknife(i);
    sum+=unjks[i];
  }
  for (int i=0; i<N; i++) {
    JackknifeCmplx tmp(N-1);
    jk[i]=tmp;
    jk[i].ave=J.ReturnJk(i);
    for (int j=0; j<i; j++)
      jk[i].jk[j]=(sum-unjks[i]-unjks[j])/double(N-2);
    for (int j=i+1; j<N; j++)
      jk[i].jk[j-1]=(sum-unjks[i]-unjks[j])/double(N-2);
    jk[i].CalcAll();
  }
  ave=J.ave;
  err=J.err;
  return *this;
}

//Friend function version of the above.
DoubleJackknifeCmplx FromSingleJk(const JackknifeCmplx & J)
{
  DoubleJackknifeCmplx tmp;
  tmp.FromSingleJk(J);
  return tmp;
}

//Addition
DoubleJackknifeCmplx DoubleJackknifeCmplx::operator+(const DoubleJackknifeCmplx & J) const
{
  DoubleJackknifeCmplx sum(*this);
  sum+=J;
  return sum;
}

//+=
DoubleJackknifeCmplx & DoubleJackknifeCmplx::operator+=(const DoubleJackknifeCmplx & J)
{
  if (N!=J.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't add DoubleJackknifeCmplx objects with different number of jackknife values.  Object 1 N=" << N << " Object 2 N=" << J.N << "\n";
  } else {
    for (int i=0; i<N; i++)
      jk[i]+=J.jk[i];
    ave+=J.ave;
    CalcAll();
  }
  return *this;
}

//Subtraction
DoubleJackknifeCmplx DoubleJackknifeCmplx::operator-(const DoubleJackknifeCmplx & J) const
{
  DoubleJackknifeCmplx difference(*this);
  difference-=J;
  return difference;
}

//Negative (Unary operator)
DoubleJackknifeCmplx DoubleJackknifeCmplx::operator-() const
{
  DoubleJackknifeCmplx neg(N,0.0);
  neg-=*this;
  return neg;
}

//-=
DoubleJackknifeCmplx & DoubleJackknifeCmplx::operator-=(const DoubleJackknifeCmplx & J)
{
  if (N!=J.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't subtract DoubleJackknifeCmplx objects with different number of jackknife values.  Object 1 N=" << N << " Object 2 N=" << J.N << "\n";
  } else {
    for (int i=0; i<N; i++)
      jk[i]-=J.jk[i];
    ave-=J.ave;
    CalcAll();
  }
  return *this;
}

//Multiplication of two DoubleJackknifeCmplx objects
DoubleJackknifeCmplx DoubleJackknifeCmplx::operator*(const DoubleJackknifeCmplx & J) const
{
  DoubleJackknifeCmplx prod(*this);
  prod*=J;
  return prod;
}

//Multiplication by a complex scalar on the right
//Do I need a separate one for a real scalar d?  Or will it convert a double
//automatically to a complex<double>?
DoubleJackknifeCmplx DoubleJackknifeCmplx::operator*(complex<double> c) const
{
  DoubleJackknifeCmplx prod(*this);
  prod*=c;
  return prod;
}

//Multiplication by a complex scalar on the left (friend function)
//Do I need a separate one for a real scalar d?  Or will it convert a double
//automatically to a complex<double>?
DoubleJackknifeCmplx operator*(complex<double> c, const DoubleJackknifeCmplx & J)
{
  return J*c;
}

//*=DoubleJackknifeCmplx object
DoubleJackknifeCmplx & DoubleJackknifeCmplx::operator*=(const DoubleJackknifeCmplx & J)
{
  if (N!=J.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't multiply DoubleJackknifeCmplx objects with different number of jackknife values.  Object 1 N=" << N << " Object 2 N=" << J.N << "\n";
  } else {
    for (int i=0; i<N; i++)
      jk[i]*=J.jk[i];
    ave*=J.ave;
    CalcAll();
  }
  return *this;
}

//*=complex scalar
//Do I need a separate one for a real scalar d?  Or will it convert a double
//automatically to a complex<double>?
DoubleJackknifeCmplx & DoubleJackknifeCmplx::operator*=(complex<double> c)
{
  for (int i=0; i<N; i++)
    jk[i]*=c;
  ave*=c;
  CalcAll();
  return *this;
}

//Division of two DoubleJackknifeCmplx objects
DoubleJackknifeCmplx DoubleJackknifeCmplx::operator/(const DoubleJackknifeCmplx & J) const
{
  DoubleJackknifeCmplx q(*this);
  q/=J;
  return q;
}

//Division by a complex scalar
//Do I need a separate one for a real scalar d?  Or will it convert a double
//automatically to a complex<double>?
DoubleJackknifeCmplx DoubleJackknifeCmplx::operator/(complex<double> c) const
{
  DoubleJackknifeCmplx q(*this);
  q/=c;
  return q;
}

//Division of complex scalar by a DoubleJackknifeCmplx object (friend function)
//Do I need a separate one for a real scalar d?  Or will it convert a double
//automatically to a complex<double>?
DoubleJackknifeCmplx operator/(complex<double> c, const DoubleJackknifeCmplx & J)
{
  int N=J.N;
  DoubleJackknifeCmplx q(N,c);
  q/=J;
  return q;
}

//slash equals DoubleJackknifeCmplx object
DoubleJackknifeCmplx & DoubleJackknifeCmplx::operator/=(const DoubleJackknifeCmplx & J)
{
  if (N!=J.N) {
    //Error: I think I'd rather ahve it abort here.
    cout << "Can't divide DoubleJackknifeCmplx objects with different number of jackknife values.  Object 1 N=" << N << " Object 2 N=" << J.N << "\n";
  } else {
    for (int i=0; i<N; i++)
      jk[i]/=J.jk[i];
    ave/=J.ave;
    CalcAll();
  }
  return *this;
}

//slash equals complex scalar
//Do I need a separate one for a real scalar d?  Or will it convert a double
//automatically to a complex<double>?
DoubleJackknifeCmplx & DoubleJackknifeCmplx::operator/=(complex<double> c)
{
  for (int i=0; i<N; i++)
    jk[i]/=c;
  ave/=c;
  CalcAll();
  return *this;
}

//Abs (member function)
DoubleJackknifeCmplx & DoubleJackknifeCmplx::Abs()
{
  for (int i=0; i<N; i++)
    jk[i].Abs();
  ave=abs(ave);
  CalcAll();
  return *this;
}

//Real (member function)
DoubleJackknifeCmplx & DoubleJackknifeCmplx::Real()
{
  for (int i=0; i<N; i++)
    jk[i].Real();
  ave=real(ave);
  CalcAll();
  return *this;
}

//Imag (member function)
DoubleJackknifeCmplx & DoubleJackknifeCmplx::Imag()
{
  for (int i=0; i<N; i++)
    jk[i].Imag();
  ave=imag(ave);
  CalcAll();
  return *this;
}

//Conj (member function)
DoubleJackknifeCmplx & DoubleJackknifeCmplx::Conj()
{
  for (int i=0; i<N; i++)
    jk[i].Conj();
  ave=conj(ave);
  CalcAll();
  return *this;
}

//Conj (friend function)
DoubleJackknifeCmplx Conj(const DoubleJackknifeCmplx & J)
{
  DoubleJackknifeCmplx tmpc(J);
  tmpc.Conj();
  return tmpc;
}

//Combine (friend function)
//Takes two DoubleJackknifeCmplx objects with N jaccknife values and puts them
//together to make another DoubleJackknifeCmplx object with 2*N jackknife values.
//IMPORTANT: NOTE THAT it can only be used on quantities that depend linearly
//on the data
DoubleJackknifeCmplx Combine(const DoubleJackknifeCmplx & J1, const DoubleJackknifeCmplx & J2)
{
  if (J1.N != J2.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can only combine DoubleJackknifeCmplx objects with same number of jackknife values.  Object 1 N=" << J1.N << " Object 2 N=" << J2.N << "\n";
    return J1;
  }
  int N=J1.N;
  DoubleJackknifeCmplx tmp(2*N);
  complex<double> unjk[2*N];
  complex<double> sum(0.0,0.0);
  for (int i=0; i<N; i++) {
    unjk[i]=J1.Unjackknife(i);
    unjk[i+N]=J2.Unjackknife(i);
    sum+=unjk[i]+unjk[i+N];
  }
  for (int i=0; i<2*N; i++) {
    tmp.jk[i].ave=(sum-unjk[i])/double(2*N-1);
    for (int j=0; j<i; j++)
      tmp.jk[i].jk[j]=(sum-unjk[i]-unjk[j])/double(2*N-2);
    for (int j=i+1; j<2*N; j++)
      tmp.jk[i].jk[j-1]=(sum-unjk[i]-unjk[j])/double(2*N-2);
    tmp.jk[i].CalcAll();
  }
  tmp.ave=0.5*(J1.ave+J2.ave);
  tmp.CalcAll();
  return tmp;
}
