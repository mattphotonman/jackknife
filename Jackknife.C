#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <cstdlib>
using namespace std;
#include "jackknife.h"
#include "zeta_func.h"
#include "eff_masses.h"
#include "least_squares.h"
#include "global_const.h"
#include "random.h"


//Jackknife class

//Standard and Default Constructor
Jackknife::Jackknife(int n)
{
  if (n<2) {
    //Error: I think I'd rather have it abort here.
    cout << "N must be >= 2 for a Jackknife object, but was given N=" << n << ".  Setting N=2.\n";
    N=2;
  } else {
    N=n;
  }
  jk=new double [N];
}

//Copy Constructor
Jackknife::Jackknife(const Jackknife & J)
{
  N=J.N;
  jk=new double [N];
  for (int i=0; i<N; i++)
    jk[i]=J.jk[i];
  ave=J.ave;
  err=J.err;
}

//Constructor with a scalar
//Sets all jackknife values equal to the input scalar d.
Jackknife::Jackknife(int n, double d)
{
  if (n<2) {
    //Error: I think I'd rather have it abort here.
    cout << "N must be >= 2 for a Jackknife object, but was given N=" << n << ".  Setting N=2.\n";
    N=2;
  } else {
    N=n;
  }
  jk=new double [N];
  for (int i=0; i<N; i++)
    jk[i]=d;
  ave=d;
  CalcAll();
}

//Destructor
Jackknife::~Jackknife()
{
  delete [] jk;
}

//Assignment operator
//Copies the jk array from object on the right hand side into a new block of
//memory which the object on the left hand side points to.  Also, the
//memory previously used by the jk array of the object on the left hand
//side is deleted.
Jackknife & Jackknife::operator=(const Jackknife & J)
{
  if (this == &J)
    return *this;
  delete [] jk;
  N=J.N;
  jk=new double [N];
  for (int i=0; i<N; i++)
    jk[i]=J.jk[i];
  ave=J.ave;
  err=J.err;
  return *this;
}

//Calculate the jackknife error and store in err
Jackknife & Jackknife::CalcErr()
{
  err=0;

  if (jk_err_type<=0) {
    //Using error estimator
    //error^2 = N/(N-1) * sum( [ f(x^{jk}_i)-f(\bar{x}) ]^2 , i=0..N-1)
    //Now I can't figure why I used the factor N/(N-1) rather than (N-1)/N.

    for (int i=0; i<N; i++)
      err+=(jk[i]-ave)*(jk[i]-ave);
    err*=((double) N)/((double) (N-1));
    err=sqrt(err);

  } else if (jk_err_type>=1) {
    //Using error estimator
    //error^2 = (N-1)/N * sum( [ f(x^{jk}_i)-ave_jk ]^2 , i=0..N-1)
    //where ave_jk = 1/N * sum( f(x^{jk}_i) , i=0..N-1) is the average
    //of the jackknife values.
    //This seems to be the more standard (and possibly more accurate) method.
    
    double ave_jk=0.0;
    for (int i=0; i<N; i++)
      ave_jk+=jk[i];
    ave_jk/=double(N);
    for (int i=0; i<N; i++)
      err+=(jk[i]-ave_jk)*(jk[i]-ave_jk);
    err*=double(N-1)/double(N);
    err=sqrt(err);
    
  }
  
  return *this;
}

//Calculate everything that needs to be recalculated when a function is applied
//to a Jackknife object, or when two Jackknife objects are combined via a
//function.  For now this is just the error, but if in the future I add in
//things like the average of the jackknife values then this would also have
//to be recalculated.
Jackknife & Jackknife::CalcAll()
{
  CalcErr();
  return *this;
}

//Return the average
double Jackknife::ReturnAve() const
{
  return ave;
}

//Return the error
double Jackknife::ReturnErr() const
{
  return err;
}

//Return the ith jackknife value
double Jackknife::ReturnJk(int i) const
{
  if (i<0 || i>=N) {
    //Error: I think I'd rather have it abort here.
    cout << "i must be 0<=i<N, but received i=" << i << ".  Returning first jackknife value instead.\n";
    i=0;
  }
  return jk[i];
}

//Return N
int Jackknife::ReturnN() const
{
  return N;
}

//Unjackknife undoes the jackknife and returns the original data point
//Only valid for quantities that depend linearly on the original data
double Jackknife::Unjackknife(int i) const
{
  if ( i<0 || i>=N ) {
    //Error: I think I'd rather have it abort here.
    cout << "Data only exists from Jackknife object for 0<=i<" << N << " but was given i=" << i << ".\n";
    return 0.0;
  }
  double unjk=0.0;
  for (int j=0; j<N; j++)
    if (j!=i)
      unjk+=jk[j];
  unjk-=((double) N-2)*jk[i];
  return unjk;
}

//Takes an average and an array for the jackknife values and puts
//that into the Jackknife object.
Jackknife & Jackknife::FromArray(double ave_in, double jk_in[], int N_in)
{
  delete [] jk;
  N=N_in;
  jk=new double [N];
  for (int i=0; i<N; i++)
    jk[i]=jk_in[i];
  ave=ave_in;
  CalcAll();
  return *this;
}

//Addition
Jackknife Jackknife::operator+(const Jackknife & J) const
{
  Jackknife sum(*this);
  sum+=J;
  return sum;
}

//+=
Jackknife & Jackknife::operator+=(const Jackknife & J)
{
  if (N!=J.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't add Jackknife objects with different number of jackknife values.  Object 1 N=" << N << " Object 2 N=" << J.N << "\n";
  } else {
    for (int i=0; i<N; i++)
      jk[i]+=J.jk[i];
    ave+=J.ave;
    CalcAll();
  }
  return *this;
}

//Subtraction
Jackknife Jackknife::operator-(const Jackknife & J) const
{
  Jackknife difference(*this);
  difference-=J;
  return difference;
}

//Negative (Unary operator)
Jackknife Jackknife::operator-() const
{
  Jackknife neg(N,0.0);
  neg-=*this;
  return neg;
}

//-=
Jackknife & Jackknife::operator-=(const Jackknife & J)
{
  if (N!=J.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't subtract Jackknife objects with different number of jackknife values.  Object 1 N=" << N << " Object 2 N=" << J.N << "\n";
  } else {
    for (int i=0; i<N; i++)
      jk[i]-=J.jk[i];
    ave-=J.ave;
    CalcAll();
  }
  return *this;
}

//Multiplication of two Jackknife objects
Jackknife Jackknife::operator*(const Jackknife & J) const
{
  Jackknife prod(*this);
  prod*=J;
  return prod;
}

//Multiplication by a scalar on the right
Jackknife Jackknife::operator*(double d) const
{
  Jackknife prod(*this);
  prod*=d;
  return prod;
}

//Multiplication by a scalar on the left (friend function)
Jackknife operator*(double d, const Jackknife & J)
{
  return J*d;
}

//*=Jackknife object
Jackknife & Jackknife::operator*=(const Jackknife & J)
{
  if (N!=J.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't multiply Jackknife objects with different number of jackknife values.  Object 1 N=" << N << " Object 2 N=" << J.N << "\n";
  } else {
    for (int i=0; i<N; i++)
      jk[i]*=J.jk[i];
    ave*=J.ave;
    CalcAll();
  }
  return *this;
}

//*=scalar
Jackknife & Jackknife::operator*=(double d)
{
  for (int i=0; i<N; i++)
    jk[i]*=d;
  ave*=d;
  CalcAll();
  return *this;
}

//Division of two Jackknife objects
Jackknife Jackknife::operator/(const Jackknife & J) const
{
  Jackknife q(*this);
  q/=J;
  return q;
}

//Division by a scalar
Jackknife Jackknife::operator/(double d) const
{
  Jackknife q(*this);
  q/=d;
  return q;
}

//Division of scalar by a Jackknife object (friend function)
Jackknife operator/(double d, const Jackknife & J)
{
  int N=J.N;
  Jackknife q(N,d);
  q/=J;
  return q;
}

//slash equals Jackknife object
Jackknife & Jackknife::operator/=(const Jackknife & J)
{
  if (N!=J.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't divide Jackknife objects with different number of jackknife values.  Object 1 N=" << N << " Object 2 N=" << J.N << "\n";
  } else {
    for (int i=0; i<N; i++)
      jk[i]/=J.jk[i];
    ave/=J.ave;
    CalcAll();
  }
  return *this;
}

//slash equals scalar
Jackknife & Jackknife::operator/=(double d)
{
  for (int i=0; i<N; i++)
    jk[i]/=d;
  ave/=d;
  CalcAll();
  return *this;
}

//Log (member function)
Jackknife & Jackknife::Log()
{
  for (int i=0; i<N; i++) {
    jk[i]=log(jk[i]);
  }
  ave=log(ave);
  CalcAll();
  return *this;
}

//Log (friend function)
Jackknife Log(const Jackknife & J)
{
  Jackknife tmp(J);
  tmp.Log();
  return tmp;
}

//Sqrt (member function)
Jackknife & Jackknife::Sqrt()
{
  for (int i=0; i<N; i++) {
    jk[i]=sqrt(jk[i]);
  }
  ave=sqrt(ave);
  CalcAll();
  return *this;
}

//Sqrt (friend function)
Jackknife Sqrt(const Jackknife & J)
{
  Jackknife tmp(J);
  tmp.Sqrt();
  return tmp;
}

//Exp (member function)
Jackknife & Jackknife::Exp()
{
  for (int i=0; i<N; i++) {
    jk[i]=exp(jk[i]);
  }
  ave=exp(ave);
  CalcAll();
  return *this;
}

//Exp (friend function)
Jackknife Exp(const Jackknife & J)
{
  Jackknife tmp(J);
  tmp.Exp();
  return tmp;
}

//Meff_two_point (friend function)
Jackknife Meff_two_point(const Jackknife & C1, const Jackknife & C2, int t, int Nt)
{
  Jackknife result=C1;
  int N=result.ReturnN();
  result.ave=meff_two_point(C1.ReturnAve(),C2.ReturnAve(),t,Nt);
  for (int i=0; i<N; i++) {
    result.jk[i]=meff_two_point(C1.ReturnJk(i),C2.ReturnJk(i),t,Nt);
  }
  result.CalcAll();
  return result;
}

//Meff_three_point (friend function)
Jackknife Meff_three_point(const Jackknife & C1, const Jackknife & C2, const Jackknife & C3, int t, int Nt)
{
  Jackknife result=C1;
  int N=result.ReturnN();
  result.ave=meff_three_point(C1.ReturnAve(),C2.ReturnAve(),C3.ReturnAve(),t,Nt);
  for (int i=0; i<N; i++) {
    result.jk[i]=meff_three_point(C1.ReturnJk(i),C2.ReturnJk(i),C3.ReturnJk(i),t,Nt);
  }
  result.CalcAll();
  return result;
}

//Phi (member function)
Jackknife & Jackknife::Phi(int tau_x, int tau_y, int tau_z)
{
  ave=phi_func(ave,tau_x,tau_y,tau_z);
  for (int i=0; i<N; i++) {
    jk[i]=phi_func_phi0(jk[i],tau_x,tau_y,tau_z,ave);
  }
  CalcAll();
  return *this;
}

//Phi (friend function)
Jackknife Phi(const Jackknife & J, int tau_x, int tau_y, int tau_z)
{
  Jackknife tmp(J);
  tmp.Phi(tau_x,tau_y,tau_z);
  return tmp;
}

//Phi_Phi0 (member function)
Jackknife & Jackknife::Phi_Phi0(int tau_x, int tau_y, int tau_z, double phi0)
{
  ave=phi_func_phi0(ave,tau_x,tau_y,tau_z,phi0);
  for (int i=0; i<N; i++) {
    jk[i]=phi_func_phi0(jk[i],tau_x,tau_y,tau_z,ave);
  }
  CalcAll();
  return *this;
}

//Phi_Phi0 (friend function)
Jackknife Phi_Phi0(const Jackknife & J, int tau_x, int tau_y, int tau_z, double phi0)
{
  Jackknife tmp(J);
  tmp.Phi_Phi0(tau_x,tau_y,tau_z,phi0);
  return tmp;
}

//dPhi (member function)
Jackknife & Jackknife::dPhi(int tau_x, int tau_y, int tau_z)
{
  for (int i=0; i<N; i++) {
    jk[i]=dphi_func(jk[i],tau_x,tau_y,tau_z);
  }
  ave=dphi_func(ave,tau_x,tau_y,tau_z);
  CalcAll();
  return *this;
}

//dPhi (friend function)
Jackknife dPhi(const Jackknife & J, int tau_x, int tau_y, int tau_z)
{
  Jackknife tmp(J);
  tmp.dPhi(tau_x,tau_y,tau_z);
  return tmp;
}

//Bin (member function)
//Note: Only valid if quantity depends linearly on data
Jackknife & Jackknife::Bin(int bin_size)
{
  if (bin_size<1 || bin_size>=N) {
    //Error: I think I'd rather have it abort here.
    cout << "Bin size must be >=1 and <N.\n";
    return *this;
  } else if (bin_size==1)
    return *this;
  
  //Bin the unjackknifed data
  int Ntmp=N/bin_size;
  int Nbins=Ntmp;
  if (N%bin_size!=0) Nbins++;  //If it doesn't divide evenly then
                               //make an extra bin that's smaller.
  double binned[Nbins];
  for (int bin_num=0; bin_num<Ntmp; bin_num++) {
    binned[bin_num]=0.0;
    for (int i=bin_num*bin_size; i<(bin_num+1)*bin_size; i++)
      binned[bin_num]+=Unjackknife(i);
    binned[bin_num]/=double(bin_size);
  }
  if (N%bin_size!=0) {
    binned[Nbins-1]=0.0;
    for (int i=Ntmp*bin_size; i<N; i++)
      binned[Nbins-1]+=Unjackknife(i);
    binned[Nbins-1]/=double(N-Ntmp*bin_size);
  }

  //Replace the data for this object with binned data
  N=Nbins;
  double ave_bins=0.0;
  for (int i=0; i<N; i++)
    ave_bins+=binned[i];
  ave_bins/=double(N);
  if (jk_err_type<=0)
    ave=ave_bins;
  //if jk_err_type>=1 keep the old ave
  delete [] jk;
  jk=new double [N];
  for (int i=0; i<N; i++)
    jk[i]=(double(N)*ave_bins-binned[i])/double(N-1);

  CalcAll();
  return *this;
}

//Bin (friend function)
//Note: Only valid if quantity depends linearly on data
Jackknife Bin(const Jackknife & J, int bin_size)
{
  Jackknife tmp(J);
  tmp.Bin(bin_size);
  return tmp;
}

//Abs (member function)
Jackknife & Jackknife::Abs()
{
  for (int i=0; i<N; i++) {
    jk[i]=abs(jk[i]);
  }
  ave=abs(ave);
  CalcAll();
  return *this;
}

//Abs (friend function)
Jackknife Abs(const Jackknife & J)
{
  Jackknife tmp(J);
  tmp.Abs();
  return tmp;
}

//Abs (friend function with a JackknifeCmplx argument)
Jackknife Abs(const JackknifeCmplx & J)
{
  Jackknife tmp(J.N);
  for (int i=0; i<tmp.N; i++)
    tmp.jk[i]=abs(J.jk[i]);
  tmp.ave=abs(J.ave);
  tmp.CalcAll();
  return tmp;
}

//Real (friend function)
Jackknife Real(const JackknifeCmplx & J)
{
  Jackknife tmp(J.N);
  for (int i=0; i<tmp.N; i++)
    tmp.jk[i]=real(J.jk[i]);
  tmp.ave=real(J.ave);
  tmp.CalcAll();
  return tmp;
}

//Imag (friend function)
Jackknife Imag(const JackknifeCmplx & J)
{
  Jackknife tmp(J.N);
  for (int i=0; i<tmp.N; i++)
    tmp.jk[i]=imag(J.jk[i]);
  tmp.ave=imag(J.ave);
  tmp.CalcAll();
  return tmp;
}

//Combine (friend function)
//Takes two Jackknife objects with N jaccknife values and puts them together 
//to make another Jackknife object with 2*N jackknife values.
//IMPORTANT: NOTE THAT it can only be used on quantities that depend linearly
//on the data
Jackknife Combine(const Jackknife & J1, const Jackknife & J2)
{
  if (J1.N != J2.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can only combine Jackknife objects with same number of jackknife values.  Object 1 N=" << J1.N << " Object 2 N=" << J2.N << "\n";
    return J1;
  }
  int N=J1.N;
  Jackknife tmp(2*N);
  for (int i=0; i<N; i++) {
    tmp.jk[i]=0.0;
    for (int j=0; j<N; j++) {
      if (j!=i)
	tmp.jk[i]+=J1.Unjackknife(j);
      tmp.jk[i]+=J2.Unjackknife(j);
    }
    tmp.jk[i]/=((double) 2*N-1);
    tmp.jk[i+N]=0.0;
    for (int j=0; j<N; j++) {
      tmp.jk[i+N]+=J1.Unjackknife(j);
      if (j!=i)
	tmp.jk[i+N]+=J2.Unjackknife(j);
    }
    tmp.jk[i+N]/=((double) 2*N-1);
  }
  tmp.ave=0.5*(J1.ave+J2.ave);
  tmp.CalcAll();
  return tmp;
}

//Reads in a text file containing the average followed by the jackknife values
//on successive lines.
Jackknife & Jackknife::ReadTextFile(string file, int Nfile)
{
  double ave_tmp, jk_tmp[Nfile];
  int success=ReadTextFiletoArray(file,Nfile,ave_tmp,jk_tmp);

  if (!success)
    return *this;

  //Read it in properly, now put values into Jackknife object.
  FromArray(ave_tmp,jk_tmp,Nfile);

  return *this;
}

//Linear least squares fitting with one dependent variable.  This just
//calls Simul_Lin_Least_Squares in the appropriate way, and is here
//just for the convenience of not having to have an extra dimension
//in some arrays if you're only doing fitting with one dependent
//variable.
int Lin_Least_Squares(int Ncoeff, int Ndata, Jackknife* y, Jackknife** f, Jackknife C[])
{
  Jackknife** ytmp;
  Jackknife*** ftmp;
  
  ytmp=new Jackknife* [1];
  ytmp[0]=y;
  
  ftmp=new Jackknife** [1];
  ftmp[0]=f;

  int return_val=Simul_Lin_Least_Squares(1,Ncoeff,Ndata,ytmp,ftmp,C);
  
  //Free memory allocated in this function.
  delete [] ytmp;
  delete [] ftmp;

  return return_val;
}

//Linear least squares fitting with in general multiple dependent
//variables.  The format of the inputs and the outputs is exactly
//the same as the linear_least_square function (LeastSquares.C).
//The only difference between that and the current function is
//that the current one does the fit under the jackknife.  Since
//the objects are only single Jackknife objects, errors are "frozen".
int Simul_Lin_Least_Squares(int Nvar, int Ncoeff, int Ndata, Jackknife** y, Jackknife*** f, Jackknife C[])
{
  double** ytmp;
  double** sigma;
  double*** ftmp;
  double* Ctmp;

  ytmp=new double* [Nvar];
  for (int alpha=0; alpha<Nvar; alpha++)
    ytmp[alpha]=new double [Ndata];
  sigma=new double* [Nvar];
  for (int alpha=0; alpha<Nvar; alpha++)
    sigma[alpha]=new double [Ndata];
  ftmp=new double** [Nvar];
  for (int alpha=0; alpha<Nvar; alpha++) {
    ftmp[alpha]=new double* [Ncoeff];
    for (int beta=0; beta<Ncoeff; beta++)
      ftmp[alpha][beta]=new double [Ndata];
  }
  Ctmp=new double [Ncoeff];
  
  int N=y[0][0].ReturnN();
  for (int alpha=0; alpha<Nvar; alpha++)
    for (int i=0; i<Ndata; i++) {
      if (y[alpha][i].ReturnN() != N) {
	//Error: I think I'd rather have it abort here.
	cout << "All y[alpha][i] must have same number of jackknife values.\n";
	return 0;
      }
      for (int beta=0; beta<Ncoeff; beta++)
	if (f[alpha][beta][i].ReturnN() != N) {
	  //Error: I think I'd rather have it abort here.
	  cout << "All f[alpha][beta][i] must have same number of jackknife values.\n";
	  return 0;
	}
    }
  
  Jackknife tmp(N);
  for (int beta=0; beta<Ncoeff; beta++)
    C[beta]=tmp;

  //Frozen errors
  for (int alpha=0; alpha<Nvar; alpha++)
    for (int i=0; i<Ndata; i++)
      sigma[alpha][i]=y[alpha][i].ReturnErr();
  
  //Fit the average
  for (int alpha=0; alpha<Nvar; alpha++)
    for (int i=0; i<Ndata; i++) {
      ytmp[alpha][i]=y[alpha][i].ReturnAve();
      for (int beta=0; beta<Ncoeff; beta++)
	ftmp[alpha][beta][i]=f[alpha][beta][i].ReturnAve();
    }
  linear_least_square(Nvar,Ncoeff,Ndata,ytmp,sigma,ftmp,Ctmp);
  for (int beta=0; beta<Ncoeff; beta++)
    C[beta].ave=Ctmp[beta];
  
  //Fit the jackknife blocks
  for (int j=0; j<N; j++) {
    for (int alpha=0; alpha<Nvar; alpha++)
      for (int i=0; i<Ndata; i++) {
	ytmp[alpha][i]=y[alpha][i].ReturnJk(j);
	for (int beta=0; beta<Ncoeff; beta++)
	  ftmp[alpha][beta][i]=f[alpha][beta][i].ReturnJk(j);
      }
    linear_least_square(Nvar,Ncoeff,Ndata,ytmp,sigma,ftmp,Ctmp);
    for (int beta=0; beta<Ncoeff; beta++)
      C[beta].jk[j]=Ctmp[beta];
  }
  
  //Have to do CalcAll to all of the C[beta]
  for (int beta=0; beta<Ncoeff; beta++)
    C[beta].CalcAll();
  
  //Free memory allocated in this function
  for (int alpha=0; alpha<Nvar; alpha++) {
    delete [] ytmp[alpha];
    delete [] sigma[alpha];
    for (int beta=0; beta<Ncoeff; beta++)
      delete [] ftmp[alpha][beta];
    delete [] ftmp[alpha];
  }
  delete [] ytmp;
  delete [] sigma;
  delete [] ftmp;
  delete [] Ctmp;

  return 1;
}

//Same as Lin_Least_Squares, but fully covariant fit.
int Cov_Lin_Least_Squares(int Ncoeff, int Ndata, Jackknife* y, Jackknife** f, Jackknife C[])
{
  Jackknife** ytmp;
  Jackknife*** ftmp;
  
  ytmp=new Jackknife* [1];
  ytmp[0]=y;
  
  ftmp=new Jackknife** [1];
  ftmp[0]=f;

  int return_val=Simul_Cov_Lin_Least_Squares(1,Ncoeff,Ndata,ytmp,ftmp,C);
  
  //Free memory allocated in this function.
  delete [] ytmp;
  delete [] ftmp;

  return return_val;
}

//Same as Simul_Lin_Least_Squares, but doing a fully covariant fit.
//The covariance matrix can be determined from the Jackknife objects,
//and is frozen (i.e. same one used to fit all jackknife blocks and
//the central value), since these are only single Jackknife objects.
int Simul_Cov_Lin_Least_Squares(int Nvar, int Ncoeff, int Ndata, Jackknife** y, Jackknife*** f, Jackknife C[])
{
  double** ytmp;
  double**** cov; 
  double**** cov_inv;
  double*** ftmp;
  double* Ctmp;

  ytmp=new double* [Nvar];
  for (int alpha=0; alpha<Nvar; alpha++)
    ytmp[alpha]=new double [Ndata];
  
  cov=new double*** [Nvar];
  for (int alpha=0; alpha<Nvar; alpha++) {
    cov[alpha]=new double** [Ndata];
    for (int i=0; i<Ndata; i++) {
      cov[alpha][i]=new double* [Nvar];
      for (int alpha_p; alpha_p<Nvar; alpha_p++)
	cov[alpha][i][alpha_p]=new double [Ndata];
    }
  }

  cov_inv=new double*** [Nvar];
  for (int alpha=0; alpha<Nvar; alpha++) {
    cov_inv[alpha]=new double** [Ndata];
    for (int i=0; i<Ndata; i++) {
      cov_inv[alpha][i]=new double* [Nvar];
      for (int alpha_p; alpha_p<Nvar; alpha_p++)
	cov_inv[alpha][i][alpha_p]=new double [Ndata];
    }
  }

  ftmp=new double** [Nvar];
  for (int alpha=0; alpha<Nvar; alpha++) {
    ftmp[alpha]=new double* [Ncoeff];
    for (int beta=0; beta<Ncoeff; beta++)
      ftmp[alpha][beta]=new double [Ndata];
  }
  Ctmp=new double [Ncoeff];
  
  int N=y[0][0].ReturnN();
  for (int alpha=0; alpha<Nvar; alpha++)
    for (int i=0; i<Ndata; i++) {
      if (y[alpha][i].ReturnN() != N) {
	//Error: I think I'd rather have it abort here.
	cout << "All y[alpha][i] must have same number of jackknife values.\n";
	return 0;
      }
      for (int beta=0; beta<Ncoeff; beta++)
	if (f[alpha][beta][i].ReturnN() != N) {
	  //Error: I think I'd rather have it abort here.
	  cout << "All f[alpha][beta][i] must have same number of jackknife values.\n";
	  return 0;
	}
    }
  
  Jackknife tmp(N);
  for (int beta=0; beta<Ncoeff; beta++)
    C[beta]=tmp;

  //Frozen covariance matrix
  for (int alpha=0; alpha<Nvar; alpha++)
    for (int i=0; i<Ndata; i++)
      for (int alpha_p=0; alpha_p<Nvar; alpha_p++)
	for (int i_p=0; i_p<Ndata; i_p++)
	  cov[alpha][i][alpha_p][i_p]=Cov(y[alpha][i],y[alpha_p][i_p]);

  //Want to invert cov.  Need to map it to a two-dimensional
  //array.  Map the index pair (alpha,i) to alpha*Ndata+i+1, where
  //the +1 is so that numbering starts at 1 since I will use
  //Numerical Recipes functions to invert.
  double** tmp_mat=new double* [Nvar*Ndata+1];
  double** tmp_mat_inv=new double* [Nvar*Ndata+1];
  for (int ii=0; ii<Nvar*Ndata+1; ii++) {
    tmp_mat[ii]=new double [Nvar*Ndata+1];
    tmp_mat_inv[ii]=new double [Nvar*Ndata+1];
  }
  for (int alpha=0; alpha<Nvar; alpha++)
    for (int i=0; i<Ndata; i++)
      for (int alpha_p=0; alpha_p<Nvar; alpha_p++)
	for (int i_p=0; i_p<Ndata; i_p++)
	  tmp_mat[alpha*Ndata+i+1][alpha_p*Ndata+i_p+1]=cov[alpha][i][alpha_p][i_p];
  //Invert tmp_mat
  invert_mat(tmp_mat,tmp_mat_inv,Nvar*Ndata);
  //Put tmp_mat_inv into cov_inv.
  for (int alpha=0; alpha<Nvar; alpha++)
    for (int i=0; i<Ndata; i++)
      for (int alpha_p=0; alpha_p<Nvar; alpha_p++)
	for (int i_p=0; i_p<Ndata; i_p++)
	  cov_inv[alpha][i][alpha_p][i_p]=tmp_mat_inv[alpha*Ndata+i+1][alpha_p*Ndata+i_p+1];
  //Delete tmp_mat and tmp_mat_inv
  for (int ii=0; ii<Nvar*Ndata+1; ii++) {
    delete [] tmp_mat[ii];
    delete [] tmp_mat_inv[ii];
  }
  delete [] tmp_mat;
  delete [] tmp_mat_inv;
  
  
  //Fit the average
  for (int alpha=0; alpha<Nvar; alpha++)
    for (int i=0; i<Ndata; i++) {
      ytmp[alpha][i]=y[alpha][i].ReturnAve();
      for (int beta=0; beta<Ncoeff; beta++)
	ftmp[alpha][beta][i]=f[alpha][beta][i].ReturnAve();
    }
  cov_linear_least_square(Nvar,Ncoeff,Ndata,ytmp,cov_inv,ftmp,Ctmp);
  for (int beta=0; beta<Ncoeff; beta++)
    C[beta].ave=Ctmp[beta];
  
  //Fit the jackknife blocks
  for (int j=0; j<N; j++) {
    for (int alpha=0; alpha<Nvar; alpha++)
      for (int i=0; i<Ndata; i++) {
	ytmp[alpha][i]=y[alpha][i].ReturnJk(j);
	for (int beta=0; beta<Ncoeff; beta++)
	  ftmp[alpha][beta][i]=f[alpha][beta][i].ReturnJk(j);
      }
    cov_linear_least_square(Nvar,Ncoeff,Ndata,ytmp,cov_inv,ftmp,Ctmp);
    for (int beta=0; beta<Ncoeff; beta++)
      C[beta].jk[j]=Ctmp[beta];
  }
  
  //Have to do CalcAll to all of the C[beta]
  for (int beta=0; beta<Ncoeff; beta++)
    C[beta].CalcAll();
  
  //Free memory allocated in this function
  for (int alpha=0; alpha<Nvar; alpha++) {
    delete [] ytmp[alpha];
    
    for (int i=0; i<Ndata; i++) {
      for (int alpha_p=0; alpha_p<Nvar; alpha_p++)
	delete [] cov[alpha][i][alpha_p];
      delete [] cov[alpha][i];
    }
    delete [] cov[alpha];

    for (int i=0; i<Ndata; i++) {
      for (int alpha_p=0; alpha_p<Nvar; alpha_p++)
	delete [] cov_inv[alpha][i][alpha_p];
      delete [] cov_inv[alpha][i];
    }
    delete [] cov_inv[alpha];

    for (int beta=0; beta<Ncoeff; beta++)
      delete [] ftmp[alpha][beta];
    delete [] ftmp[alpha];
  }
  delete [] ytmp;
  delete [] cov;
  delete [] cov_inv;
  delete [] ftmp;
  delete [] Ctmp;

  return 1;
}

//Calculates chi^2 (under the jackknife) for the output of 
//Lin_Least_Squares.  Arguments are the same as for
//Lin_Least_Squares, but they must be given to this function
//AFTER first being given to Lin_Least_Squares.
Jackknife ChiSq(int Ncoeff, int Ndata, Jackknife* y, Jackknife** f, Jackknife C[])
{
  Jackknife** ytmp;
  Jackknife*** ftmp;
  
  ytmp=new Jackknife* [1];
  ytmp[0]=y;
  
  ftmp=new Jackknife** [1];
  ftmp[0]=f;

  Jackknife return_val=Simul_ChiSq(1,Ncoeff,Ndata,ytmp,ftmp,C);
  
  //Free memory allocated in this function.
  delete [] ytmp;
  delete [] ftmp;

  return return_val;
}

//chi^2 per d.o.f. = chi^2/(Ndata-Ncoeff) (note: Nvar=1)
Jackknife ChiSq_Per_Dof(int Ncoeff, int Ndata, Jackknife* y, Jackknife** f, Jackknife C[])
{
  return ChiSq(Ncoeff,Ndata,y,f,C)/double(Ndata-Ncoeff);
}

//Calculates chi^2 (under the jackknife) for the output of 
//Simul_Lin_Least_Squares.  Arguments are the same as for
//Simul_Lin_Least_Squares, but they must be given to this function
//AFTER first being given to Simul_Lin_Least_Squares.
Jackknife Simul_ChiSq(int Nvar, int Ncoeff, int Ndata, Jackknife** y, Jackknife*** f, Jackknife C[])
{
  Jackknife result=0.0*C[0];
  for (int alpha=0; alpha<Nvar; alpha++)
    for (int i=0; i<Ndata; i++) {
      Jackknife beta_sum=0.0*C[0];
      for (int beta=0; beta<Ncoeff; beta++)
	beta_sum+=C[beta]*f[alpha][beta][i];
      Jackknife tmp=(y[alpha][i]-beta_sum)/(y[alpha][i].ReturnErr());
      result+=tmp*tmp;
    }
  
  return result;
}

//chi^2 per d.o.f. = chi^2/(Nvar*Ndata-Ncoeff)
Jackknife Simul_ChiSq_Per_Dof(int Nvar, int Ncoeff, int Ndata, Jackknife** y, Jackknife*** f, Jackknife C[])
{
  return Simul_ChiSq(Nvar,Ncoeff,Ndata,y,f,C)/double(Nvar*Ndata-Ncoeff);
}

//Same as ChiSq, but for the output of Cov_Lin_Least_Squares.
Jackknife Cov_ChiSq(int Ncoeff, int Ndata, Jackknife* y, Jackknife** f, Jackknife C[])
{
  Jackknife** ytmp;
  Jackknife*** ftmp;
  
  ytmp=new Jackknife* [1];
  ytmp[0]=y;
  
  ftmp=new Jackknife** [1];
  ftmp[0]=f;

  Jackknife return_val=Simul_Cov_ChiSq(1,Ncoeff,Ndata,ytmp,ftmp,C);
  
  //Free memory allocated in this function.
  delete [] ytmp;
  delete [] ftmp;

  return return_val;
}

//chi^2 per d.o.f. = chi^2/(Ndata-Ncoeff) (note: Nvar=1)
Jackknife Cov_ChiSq_Per_Dof(int Ncoeff, int Ndata, Jackknife* y, Jackknife** f, Jackknife C[])
{
  return Cov_ChiSq(Ncoeff,Ndata,y,f,C)/double(Ndata-Ncoeff);
}

//Same as Simul_ChiSq but for output of Simul_Cov_Lin_Least_Squares
Jackknife Simul_Cov_ChiSq(int Nvar, int Ncoeff, int Ndata, Jackknife** y, Jackknife*** f, Jackknife C[])
{
  double**** cov;
  double**** cov_inv;
  
  cov=new double*** [Nvar];
  for (int alpha=0; alpha<Nvar; alpha++) {
    cov[alpha]=new double** [Ndata];
    for (int i=0; i<Ndata; i++) {
      cov[alpha][i]=new double* [Nvar];
      for (int alpha_p; alpha_p<Nvar; alpha_p++)
	cov[alpha][i][alpha_p]=new double [Ndata];
    }
  }
  
  cov_inv=new double*** [Nvar];
  for (int alpha=0; alpha<Nvar; alpha++) {
    cov_inv[alpha]=new double** [Ndata];
    for (int i=0; i<Ndata; i++) {
      cov_inv[alpha][i]=new double* [Nvar];
      for (int alpha_p; alpha_p<Nvar; alpha_p++)
	cov_inv[alpha][i][alpha_p]=new double [Ndata];
    }
  }


  //First we need to calculate the covariance matrix
  //and invert it.
  for (int alpha=0; alpha<Nvar; alpha++)
    for (int i=0; i<Ndata; i++)
      for (int alpha_p=0; alpha_p<Nvar; alpha_p++)
	for (int i_p=0; i_p<Ndata; i_p++)
	  cov[alpha][i][alpha_p][i_p]=Cov(y[alpha][i],y[alpha_p][i_p]);

  //Want to invert cov.  Need to map it to a two-dimensional
  //array.  Map the index pair (alpha,i) to alpha*Ndata+i+1, where
  //the +1 is so that numbering starts at 1 since I will use
  //Numerical Recipes functions to invert.
  double** tmp_mat=new double* [Nvar*Ndata+1];
  double** tmp_mat_inv=new double* [Nvar*Ndata+1];
  for (int ii=0; ii<Nvar*Ndata+1; ii++) {
    tmp_mat[ii]=new double [Nvar*Ndata+1];
    tmp_mat_inv[ii]=new double [Nvar*Ndata+1];
  }
  for (int alpha=0; alpha<Nvar; alpha++)
    for (int i=0; i<Ndata; i++)
      for (int alpha_p=0; alpha_p<Nvar; alpha_p++)
	for (int i_p=0; i_p<Ndata; i_p++)
	  tmp_mat[alpha*Ndata+i+1][alpha_p*Ndata+i_p+1]=cov[alpha][i][alpha_p][i_p];
  //Invert tmp_mat
  invert_mat(tmp_mat,tmp_mat_inv,Nvar*Ndata);
  //Put tmp_mat_inv into cov_inv.
  for (int alpha=0; alpha<Nvar; alpha++)
    for (int i=0; i<Ndata; i++)
      for (int alpha_p=0; alpha_p<Nvar; alpha_p++)
	for (int i_p=0; i_p<Ndata; i_p++)
	  cov_inv[alpha][i][alpha_p][i_p]=tmp_mat_inv[alpha*Ndata+i+1][alpha_p*Ndata+i_p+1];
  //Delete tmp_mat and tmp_mat_inv
  for (int ii=0; ii<Nvar*Ndata+1; ii++) {
    delete [] tmp_mat[ii];
    delete [] tmp_mat_inv[ii];
  }
  delete [] tmp_mat;
  delete [] tmp_mat_inv;

  //Now we can calculate chi^2

  Jackknife result=0.0*C[0];
  for (int alpha=0; alpha<Nvar; alpha++)
    for (int alpha_p=0; alpha_p<Nvar; alpha_p++)
      for (int i=0; i<Ndata; i++)
	for (int i_p=0; i_p<Ndata; i_p++) {
	  Jackknife beta_sum1=0.0*C[0];
	  Jackknife beta_sum2=beta_sum1;
	  for (int beta=0; beta<Ncoeff; beta++) {
	    beta_sum1+=C[beta]*f[alpha][beta][i];
	    beta_sum2+=C[beta]*f[alpha_p][beta][i_p];
	  }
	  Jackknife tmp1=y[alpha][i]-beta_sum1;
	  Jackknife tmp2=y[alpha_p][i_p]-beta_sum2;
	  result+=cov_inv[alpha][i][alpha_p][i_p]*tmp1*tmp2;
	}


  //Delete memory allocated
  for (int alpha=0; alpha<Nvar; alpha++) {
    
    for (int i=0; i<Ndata; i++) {
      for (int alpha_p=0; alpha_p<Nvar; alpha_p++)
	delete [] cov[alpha][i][alpha_p];
      delete [] cov[alpha][i];
    }
    delete [] cov[alpha];

    for (int i=0; i<Ndata; i++) {
      for (int alpha_p=0; alpha_p<Nvar; alpha_p++)
	delete [] cov_inv[alpha][i][alpha_p];
      delete [] cov_inv[alpha][i];
    }
    delete [] cov_inv[alpha];

  }
  delete [] cov;
  delete [] cov_inv;


  //Return chi^2  
  return result;
}

//chi^2 per d.o.f. = chi^2/(Nvar*Ndata-Ncoeff)
Jackknife Simul_Cov_ChiSq_Per_Dof(int Nvar, int Ncoeff, int Ndata, Jackknife** y, Jackknife*** f, Jackknife C[])
{
  return Simul_Cov_ChiSq(Nvar,Ncoeff,Ndata,y,f,C)/double(Nvar*Ndata-Ncoeff);
}

//Outputs the average followed by the jackknife values to a text file.
//Not a friend or a member function.
void OutputAveJk(string filename, const Jackknife & J)
{
  ofstream fout(filename.c_str());
  fout.precision(16);
  int N=J.ReturnN();
  fout << J.ReturnAve() << "\n";
  for (int i=0;i<N; i++)
    fout << J.ReturnJk(i) << "\n";
  fout.close();
}

//Reads in a text file containing the average followed by N jackknife values
//(where N is an input to the function), and puts them into the
//arrays ave and jk.
int ReadTextFiletoArray(string file, int N, double & ave, double jk[])
{
  if (N<2) {
    //Error: I think I'd rather have it abort here.
    cout << "ReadTextFiletoArray: N must be >=2 for a Jackknife object, but was given N=" << N << ".\n";
    return 0;
  }
  
  ifstream fin(file.c_str());
  string line;
  if (!getline(fin,line)) {
    //Error: I think I'd rather have it abort here.
    cout << "ReadTextFiletoArray: File empty.\n";
    return 0;
  }
  ave=atof(line.c_str());
  for (int i=0; i<N; i++) {
    if (!getline(fin,line)) {
      //Error: I think I'd rather have it abort here.
      cout << "ReadTextFiletoArray: File ended at jackknife value " << i << "\n";
      return 0;
    }
    jk[i]=atof(line.c_str());
  }
  if (getline(fin,line)) {
    //Error: I think I'd rather have it abort here.
    cout << "ReadTextFiletoArray: File didn't end after last jackknife value.\n";
    return 0;
  }
  fin.close();
  return 1;
}

//Takes a covariance matrix C[0..Nvar-1][0..Nvar-1] and a list of central
//values aves[0..Nvar-1], and returns an array of Jackknife objects
//Jacks[0..Nvar-1] that have the central values given by aves, are correlated
//as in the covariance matrix C, and that have Nsamp jackknife blocks (Nsamp
//is an input).  The Jackknife objects are constructed from fake data generated
//using Gaussian random numbers (which are generated in the function 
//gauss_corr_variables).
//(neither a member nor a friend function)
void JacksFromCovMat(double aves [], double** C, int Nvar, int Nsamp, Jackknife Jacks [])
{
  double** rands=new double* [Nvar];
  for (int i=0; i<Nvar; i++)
    rands[i]=new double [Nsamp];

  //NsmpC = Nsamp*C.
  double** NsmpC=new double* [Nvar];
  for (int i=0; i<Nvar; i++) {
    NsmpC[i]=new double [Nvar];
    for (int j=0; j<Nvar; j++)
      NsmpC[i][j]=Nsamp*C[i][j];
  }
  
  //Generate data points with covariance matrix Nsamp*C.
  //Since there will be Nsamp such data points, the averages
  //will have covariance matrix C, which is what we want.
  gauss_corr_variables(aves,NsmpC,Nvar,Nsamp,rands);

  //Now rands[i] contains an array of size Nsamp.  These
  //Nsamp values need to be jackknifed and fed into
  //Jacks[i], with central value coming from aves[i].
  double tot, jks[Nsamp];
  for (int i=0; i<Nvar; i++) {
    tot=0.0;
    for (int j=0; j<Nsamp; j++)
      tot+=rands[i][j];
    for (int j=0; j<Nsamp; j++)
      jks[j]=(tot-rands[i][j])/double(Nsamp-1);
    Jacks[i].FromArray(aves[i],jks,Nsamp);
  } 

  //Free memory allocated in this function by "new".
  for (int i=0; i<Nvar; i++) {
    delete [] rands[i];
    delete [] NsmpC[i];
  }
  delete [] rands;
  delete [] NsmpC;
}

//Calculates the covariance between two jackknife objects.
//The normalization is such that Cov(J,J)=J.err^2.
//(neither a member nor a friend function)
double Cov(const Jackknife & J1, const Jackknife & J2)
{
  int N=J1.ReturnN();
  if (J2.ReturnN()!=N) {
    //Error: I think I'd rather have it abort here.
    cout << "Cov: Can't take the covariance of Jackknife objects of different size.\n";
    return 0.0;
  }
  
  double cov=0.0;
  
  if (jk_err_type<=0) {
    //See CalcErr for explanation of jk_err_type.

    for (int i=0; i<N; i++)
      cov+=(J1.ReturnJk(i)-J1.ReturnAve())*(J2.ReturnJk(i)-J2.ReturnAve());
    cov*=double(N)/double(N-1);
    
  } else if (jk_err_type>=1) {
    //See CalcErr for explanation of jk_err_type.

    double ave_jk1=0.0, ave_jk2=0.0;
    for (int i=0; i<N; i++) {
      ave_jk1+=J1.ReturnJk(i);
      ave_jk2+=J2.ReturnJk(i);
    }
    ave_jk1/=double(N);
    ave_jk2/=double(N);
    for (int i=0; i<N; i++)
      cov+=(J1.ReturnJk(i)-ave_jk1)*(J2.ReturnJk(i)-ave_jk2);
    cov*=double(N-1)/double(N);

  }

  return cov;
  
}



//JackknifeCmplx class

//Standard and Default Constructor
JackknifeCmplx::JackknifeCmplx(int n)
{
  if (n<2) {
    //Error: I think I'd rather have it abort here.
    cout << "N must be >= 2 for JackknifeCmplx object, but was given N=" << n << ".  Setting N=2.\n";
    N=2;
  } else {
    N=n;
  }
  jk=new complex<double> [N];
}

//Copy Constructor
JackknifeCmplx::JackknifeCmplx(const JackknifeCmplx & J)
{
  N=J.N;
  jk=new complex<double> [N];
  for (int i=0; i<N; i++)
    jk[i]=J.jk[i];
  ave=J.ave;
  err=J.err;
}

//Constructor with a complex scalar
//Sets all jackknife values equal to the input scalar c.
//Do I need a separate one for a real scalar d?  Or will it convert a double
//automatically to a complex<double>?
JackknifeCmplx::JackknifeCmplx(int n, complex<double> c)
{
  if (n<2) {
    //Error: I think I'd rather have it abort here.
    cout << "N must be >= 2 for a JackknifeCmplx object, but was given N=" << n << ".  Setting N=2.\n";
    N=2;
  } else {
    N=n;
  }
  jk=new complex<double> [N];
  for (int i=0; i<N; i++)
    jk[i]=c;
  ave=c;
  CalcAll();
}

//Destructor
JackknifeCmplx::~JackknifeCmplx()
{
  delete [] jk;
}

//Assignment operator
//Copies the jk array from object on the right hand side into a new block of
//memory which the object on the left hand side points to.  Also, the
//memory previously used by the jk array of the object on the left hand
//side is deleted.
JackknifeCmplx & JackknifeCmplx::operator=(const JackknifeCmplx & J)
{
  if (this == &J)
    return *this;
  delete [] jk;
  N=J.N;
  jk=new complex<double> [N];
  for (int i=0; i<N; i++)
    jk[i]=J.jk[i];
  ave=J.ave;
  err=J.err;
  return *this;
}

//Calculate the jackknife error and store in err
JackknifeCmplx & JackknifeCmplx::CalcErr()
{
  err=0;

  if (jk_err_type<=0) {
    //Using error estimator
    //error^2 = N/(N-1) * sum( | f(x^{jk}_i)-f(\bar{x}) |^2 , i=0..N-1)
    //Now I can't figure why I used the factor N/(N-1) rather than (N-1)/N.
    
    for (int i=0; i<N; i++)
      err+=real(conj(jk[i]-ave)*(jk[i]-ave)); //real just converts to a double,
                                              //since this is the modulus squared
                                              //it will already be a real number.
    err*=((double) N)/((double) (N-1));
    err=sqrt(err);

  } else if (jk_err_type>=1) {
    //Using error estimator
    //error^2 = (N-1)/N * sum( | f(x^{jk}_i)-ave_jk |^2 , i=0..N-1)
    //where ave_jk = 1/N * sum( f(x^{jk}_i) , i=0..N-1) is the average
    //of the jackknife values.
    //This seems to be the more standard (and possibly more accurate) method.

    complex<double> ave_jk=0.0;
    for (int i=0; i<N; i++)
      ave_jk+=jk[i];
    ave_jk/=double(N);
    for (int i=0; i<N; i++)
      err+=real(conj(jk[i]-ave_jk)*(jk[i]-ave_jk));  //real just converts to a double,
                    //since this is the modulus squared
                    //it will already be a real number.
    err*=double(N-1)/double(N);
    err=sqrt(err);
    
  }

  return *this;
}

//Calculate everything that needs to be recalculated when a function is applied
//to a JackknifeCmplx object, or when two JackknifeCmplx objects are combined
//via a function.  For now this is just the error, but if in the future I add
//in things like the average of the jackknife values then this would also have
//to be recalculated.
JackknifeCmplx & JackknifeCmplx::CalcAll()
{
  CalcErr();
  return *this;
}

//Return the average
complex<double> JackknifeCmplx::ReturnAve() const
{
  return ave;
}

//Return the error
double JackknifeCmplx::ReturnErr() const
{
  return err;
}

//Return the ith jackknife value
complex<double> JackknifeCmplx::ReturnJk(int i) const
{
  if (i<0 || i>=N) {
    //Error: I think I'd rather have it abort here.
    cout << "i must be 0<=i<N, but received i=" << i << ".  Returning first jackknife value instead.\n";
    i=0;
  }
  return jk[i];
}

//Return N
int JackknifeCmplx::ReturnN() const
{
  return N;
}

//Unjackknife undoes the jackknife and returns the original data point
//Only valid for quantities that depend linearly on the original data
complex<double> JackknifeCmplx::Unjackknife(int i) const
{
  if ( i<0 || i>=N ) {
    //Error: I think I'd rather have it abort here.
    cout << "Data only exists from JackknifeCmplx object for 0<=i<" << N << " but was given i=" << i << ".\n";
    return 0.0;
  }
  complex<double> unjk=0.0;
  for (int j=0; j<N; j++)
    if (j!=i)
      unjk+=jk[j];
  unjk-=((double) N-2)*jk[i];
  return unjk;
}

//Addition
JackknifeCmplx JackknifeCmplx::operator+(const JackknifeCmplx & J) const
{
  JackknifeCmplx sum(*this);
  sum+=J;
  return sum;
}

//+=
JackknifeCmplx & JackknifeCmplx::operator+=(const JackknifeCmplx & J)
{
  if (N!=J.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't add JackknifeCmplx objects with different number of jackknife values.  Object 1 N=" << N << " Object 2 N=" << J.N << "\n";
  } else {
    for (int i=0; i<N; i++)
      jk[i]+=J.jk[i];
    ave+=J.ave;
    CalcAll();
  }
  return *this;
}

//Subtraction
JackknifeCmplx JackknifeCmplx::operator-(const JackknifeCmplx & J) const
{
  JackknifeCmplx difference(*this);
  difference-=J;
  return difference;
}

//Negative (Unary operator)
JackknifeCmplx JackknifeCmplx::operator-() const
{
  JackknifeCmplx neg(N,0.0);
  neg-=*this;
  return neg;
}

//-=
JackknifeCmplx & JackknifeCmplx::operator-=(const JackknifeCmplx & J)
{
  if (N!=J.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't subtract JackknifeCmplx objects with different number of jackknife values.  Object 1 N=" << N << " Object 2 N=" << J.N << "\n";
  } else {
    for (int i=0; i<N; i++)
      jk[i]-=J.jk[i];
    ave-=J.ave;
    CalcAll();
  }
  return *this;
}

//Multiplication of two JackknifeCmplx objects
JackknifeCmplx JackknifeCmplx::operator*(const JackknifeCmplx & J) const
{
  JackknifeCmplx prod(*this);
  prod*=J;
  return prod;
}

//Multiplication by a complex scalar on the right
//Do I need a separate one for a real scalar d?  Or will it convert a double
//automatically to a complex<double>?
JackknifeCmplx JackknifeCmplx::operator*(complex<double> c) const
{
  JackknifeCmplx prod(*this);
  prod*=c;
  return prod;
}

//Multiplication by a complex scalar on the left (friend function)
//Do I need a separate one for a real scalar d?  Or will it convert a double
//automatically to a complex<double>?
JackknifeCmplx operator*(complex<double> c, const JackknifeCmplx & J)
{
  return J*c;
}

//*=JackknifeCmplx object
JackknifeCmplx & JackknifeCmplx::operator*=(const JackknifeCmplx & J)
{
  if (N!=J.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can't multiply JackknifeCmplx objects with different number of jackknife values.  Object 1 N=" << N << " Object 2 N=" << J.N << "\n";
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
JackknifeCmplx & JackknifeCmplx::operator*=(complex<double> c)
{
  for (int i=0; i<N; i++)
    jk[i]*=c;
  ave*=c;
  CalcAll();
  return *this;
}

//Division of two JackknifeCmplx objects
JackknifeCmplx JackknifeCmplx::operator/(const JackknifeCmplx & J) const
{
  JackknifeCmplx q(*this);
  q/=J;
  return q;
}

//Division by a complex scalar
//Do I need a separate one for a real scalar d?  Or will it convert a double
//automatically to a complex<double>?
JackknifeCmplx JackknifeCmplx::operator/(complex<double> c) const
{
  JackknifeCmplx q(*this);
  q/=c;
  return q;
}

//Division of complex scalar by a JackknifeCmplx object (friend function)
//Do I need a separate one for a real scalar d?  Or will it convert a double
//automatically to a complex<double>?
JackknifeCmplx operator/(complex<double> c, const JackknifeCmplx & J)
{
  int N=J.N;
  JackknifeCmplx q(N,c);
  q/=J;
  return q;
}

//slash equal JackknifeCmplx object
JackknifeCmplx & JackknifeCmplx::operator/=(const JackknifeCmplx & J)
{
  if (N!=J.N) {
    //Error: I think I'd rather ahve it abort here.
    cout << "Can't divide JackknifeCmplx objects with different number of jackknife values.  Object 1 N=" << N << " Object 2 N=" << J.N << "\n";
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
JackknifeCmplx & JackknifeCmplx::operator/=(complex<double> c)
{
  for (int i=0; i<N; i++)
    jk[i]/=c;
  ave/=c;
  CalcAll();
  return *this;
}

//Abs (member function)
JackknifeCmplx & JackknifeCmplx::Abs()
{
  for (int i=0; i<N; i++)
    jk[i]=abs(jk[i]);
  ave=abs(ave);
  CalcAll();
  return *this;
}

//Real (member function)
JackknifeCmplx & JackknifeCmplx::Real()
{
  for (int i=0; i<N; i++)
    jk[i]=real(jk[i]);
  ave=real(ave);
  CalcAll();
  return *this;
}

//Imag (member function)
JackknifeCmplx & JackknifeCmplx::Imag()
{
  for (int i=0; i<N; i++)
    jk[i]=imag(jk[i]);
  ave=imag(ave);
  CalcAll();
  return *this;
}

//Conj (member function)
JackknifeCmplx & JackknifeCmplx::Conj()
{
  for (int i=0; i<N; i++)
    jk[i]=conj(jk[i]);
  ave=conj(ave);
  CalcAll();
  return *this;
}

//Conj (friend function)
JackknifeCmplx Conj(const JackknifeCmplx & J)
{
  JackknifeCmplx tmpc(J);
  tmpc.Conj();
  return tmpc;
}

//Combine (friend function)
//Takes two JackknifeCmplx objects with N jaccknife values and puts them
//together to make another JackknifeCmplx object with 2*N jackknife values.
//IMPORTANT: NOTE THAT it can only be used on quantities that depend linearly
//on the data
JackknifeCmplx Combine(const JackknifeCmplx & J1, const JackknifeCmplx & J2)
{
  if (J1.N != J2.N) {
    //Error: I think I'd rather have it abort here.
    cout << "Can only combine JackknifeCmplx objects with same number of jackknife values.  Object 1 N=" << J1.N << " Object 2 N=" << J2.N << "\n";
    return J1;
  }
  int N=J1.N;
  JackknifeCmplx tmp(2*N);
  for (int i=0; i<N; i++) {
    tmp.jk[i]=0.0;
    for (int j=0; j<N; j++) {
      if (j!=i)
	tmp.jk[i]+=J1.Unjackknife(j);
      tmp.jk[i]+=J2.Unjackknife(j);
    }
    tmp.jk[i]/=((double) 2*N-1);
    tmp.jk[i+N]=0.0;
    for (int j=0; j<N; j++) {
      tmp.jk[i+N]+=J1.Unjackknife(j);
      if (j!=i)
	tmp.jk[i+N]+=J2.Unjackknife(j);
    }
    tmp.jk[i+N]/=((double) 2*N-1);
  }
  tmp.ave=0.5*(J1.ave+J2.ave);
  tmp.CalcAll();
  return tmp;
}
