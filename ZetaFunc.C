#include <iostream>
#include <cmath>
using namespace std;
#include "zeta_func.h"
#include "global_const.h"

double phi_func_phi0(double q, int tau_x, int tau_y, int tau_z, double phi0)
{
  double phi=phi_func(q,tau_x,tau_y,tau_z);
  phi=phi-floor((phi-phi0+Pi/2.0)/Pi)*Pi;
  return phi;
}

double phi_func(double q, int tau_x, int tau_y, int tau_z)
{
  return Pi-atan2(pow(Pi,1.5)*q,zeta_func(q,tau_x,tau_y,tau_z));
}

double dphi_func(double q, int tau_x, int tau_y, int tau_z)
{
  double z=zeta_func(q,tau_x,tau_y,tau_z);
  double dz=dzeta_func(q,tau_x,tau_y,tau_z);
  return pow(Pi,1.5)*(q*dz-z)/(z*z+Pi*Pi*Pi*q*q);
}

double zeta_func(double q, int tau_x, int tau_y, int tau_z)
{
  int tau[3]={tau_x,tau_y,tau_z};
  return zeta(q,tau,1E-10,100);
}

double dzeta_func(double q, int tau_x, int tau_y, int tau_z)
{
  int tau[3]={tau_x,tau_y,tau_z};
  return dzeta(q,tau,1E-10,100);
}

double zeta(double q, int tau[], double epsilon, int N)
{
  return 1.0/sqrt(4.0*Pi)*(f1(q,epsilon)+f2(q,tau,epsilon,N)+f3(q,tau,epsilon));
}

double dzeta(double q, int tau[], double epsilon, int N)
{
  return 1.0/sqrt(4.0*Pi)*(df1(q,epsilon)+df2(q,tau,epsilon,N)+df3(q,tau,epsilon));
}

double f1(double q, double epsilon)
{
  //sum(Pi^(3/2)/(l-1/2)*(q^2)^l/l!,l=0..infinity)

  double sum=0.0;
  double powr=1.0;
  double fact=1.0;
  double term;
  double l=0.0;

  while (0==0) {
    term=powr/((l-0.5)*fact);
    sum+=term;
    if (fabs(term/sum)<epsilon/10.0) break;
    l++;
    powr*=q*q;
    fact*=l;
  }
  
  return pow(Pi,1.5)*sum;
}

double df1(double q, double epsilon)
{
  //2*q*sum(Pi^(3/2)/(l+1/2)*(q^2)^l/l!,l=0..infinity)

  double sum=0.0;
  double powr=1.0;
  double fact=1.0;
  double term;
  double l=0.0;

  while (0==0) {
    term=powr/((l+0.5)*fact);
    sum+=term;
    if (fabs(term/sum)<epsilon/10.0) break;
    l++;
    powr*=q*q;
    fact*=l;
  }
  
  return 2*q*pow(Pi,1.5)*sum;
}

double f2(double q, int tau[], double epsilon, int N)
{
  //int(exp(t*q^2)*(Pi/t)^(3/2)*sum^prime (-1)^(n dot tau)exp(-Pi^2*n^2/t),t=0..1)

  double lambda=0.1;
  double L2=lam2(lambda,epsilon);
  double out=0.0;

  for (int i=0;i<N;i++) {
    double t1=double(i)/double(N);
    double t2=double(i+1)/double(N);
    double f1, f2;
    if (i==0)
      f1=0.0;
    else
      f1=exp(t1*q*q)*pow(Pi/t1,1.5)*G_func(t1/(Pi*Pi),tau,lambda,L2);
    f2=exp(t2*q*q)*pow(Pi/t2,1.5)*G_func(t2/(Pi*Pi),tau,lambda,L2);
    out+=(f1+f2)/(2.0*double(N));
  }
  
  return out;
}

double df2(double q, int tau[], double epsilon, int N)
{
  //2*q*int(t*exp(t*q^2)*(Pi/t)^(3/2)*sum^prime (-1)^(n dot tau)exp(-Pi^2*n^2/t),t=0..1)

  double lambda=0.1;
  double L2=lam2(lambda,epsilon);
  double out=0.0;

  for (int i=0;i<N;i++) {
    double t1=double(i)/double(N);
    double t2=double(i+1)/double(N);
    double f1, f2;
    if (i==0)
      f1=0.0;
    else
      f1=t1*exp(t1*q*q)*pow(Pi/t1,1.5)*G_func(t1/(Pi*Pi),tau,lambda,L2);
    f2=t2*exp(t2*q*q)*pow(Pi/t2,1.5)*G_func(t2/(Pi*Pi),tau,lambda,L2);
    out+=(f1+f2)/(2.0*double(N));
  }
  
  return 2*q*out;
}

double f3(double q, int tau[], double epsilon)
{
  //sum(exp(-(n^2-q^2))/(n^2-q^2),n in Z^3 + tau/2)
  
  double out=0.0;
  double testval;
  int N=int(ceil(sqrt(1.0+q*q)));
  for (int m1=-N-tau[0]; m1<=N; m1++)
    for (int m2=-N-tau[1]; m2<=N; m2++)
      for (int m3=-N-tau[2]; m3<=N; m3++) {
	double n1=m1+0.5*tau[0];
	double n2=m2+0.5*tau[1];
	double n3=m3+0.5*tau[2];
	testval=n1*n1+n2*n2+n3*n3;
	if (testval<=1+q*q)
	  out+=exp(-(testval-q*q))/(testval-q*q);
      }
  testval=1.0/3.0*log(8.0/(fabs(out)*exp(-q*q)*my_pow(1.0-exp(-1.0),3)*epsilon));
  double lambda;
  if (testval>0.0)
    lambda=sqrt(testval)+1.0;
  else
    lambda=1.0;
  N=int(ceil(lambda));
  for (int m1=-N-tau[0]; m1<=N; m1++)
    for (int m2=-N-tau[1]; m2<=N; m2++)
      for (int m3=-N-tau[2]; m3<=N; m3++) {
	double n1=m1+0.5*tau[0];
	double n2=m2+0.5*tau[1];
	double n3=m3+0.5*tau[2];
	testval=n1*n1+n2*n2+n3*n3;
	if (testval>1+q*q)
	  out+=exp(-(testval-q*q))/(testval-q*q);
      }
  
  return out;
}

double df3(double q, int tau[], double epsilon)
{
  //2*q*sum(exp(-(n^2-q^2))*(1/(n^2-q^2)+1/(n^2-q^2)^2),n in Z^3 + tau/2)
  
  double out=0.0;
  double testval;
  int N=int(ceil(sqrt(1.0+q*q)));
  for (int m1=-N-tau[0]; m1<=N; m1++)
    for (int m2=-N-tau[1]; m2<=N; m2++)
      for (int m3=-N-tau[2]; m3<=N; m3++) {
	double n1=m1+0.5*tau[0];
	double n2=m2+0.5*tau[1];
	double n3=m3+0.5*tau[2];
	testval=n1*n1+n2*n2+n3*n3;
	if (testval<=1+q*q)
	  out+=exp(-(testval-q*q))*(1.0/(testval-q*q)+1.0/((testval-q*q)*(testval-q*q)));
      }
  testval=1.0/3.0*log(2.0/(fabs(out)*exp(-q*q)*my_pow(1.0-exp(-1.0),3)*epsilon));
  double lambda;
  if (testval>0.0)
    lambda=sqrt(testval)+1.0;
  else
    lambda=1.0;
  N=int(ceil(lambda));
  for (int m1=-N-tau[0]; m1<=N; m1++)
    for (int m2=-N-tau[1]; m2<=N; m2++)
      for (int m3=-N-tau[2]; m3<=N; m3++) {
	double n1=m1+0.5*tau[0];
	double n2=m2+0.5*tau[1];
	double n3=m3+0.5*tau[2];
	testval=n1*n1+n2*n2+n3*n3;
	if (testval>1+q*q)
	  out+=exp(-(testval-q*q))*(1.0/(testval-q*q)+1.0/((testval-q*q)*(testval-q*q)));
      }
  
  return 2*q*out;
}

double my_pow(double x, int n)
{
  //x^n

  int nn;
  double xx;
  if (n<0) {
    nn=-n;
    xx=1.0/x;
  } else {
    nn=n;
    xx=x;
  }
  
  double out=1.0;
  for (int i=0; i<nn; i++)
    out*=xx;
  
  return out;
}

int my_fac(int n)
{
  //n! , n>=0
  
  int out=1;
  for (int i=1; i<=n; i++)
    out*=i;
  
  return out;
}

double lam2(double lambda, double epsilon)
{
  //Upper cutoff in G_func for given lambda and tolerance
  
  int freq, min;
  Nlam(freq,min,lambda);
  double Nl=double(freq);
  double nl=sqrt(double(min));
  
  double out=sqrt(nl*nl-log((epsilon*Nl*my_pow(1.0-exp(-1.0),3))/(2*(13.0-12.0*exp(-1.0)+3.0*exp(-2.0)))));
  if (out>nl+1.0)
    return out;
  else
    return nl+1.0;
}


double G_func(double t, int tau[], double lambda, double lambda2)
{
  //sum((-1)^(n dot tau)*exp(-n^2/t),|n|>=lambda,|n_i|<lambda2)
  
  int N=int(floor(lambda2));
  if (N==lambda2) N--;
  double out=0.0;
  for (int n1=-N; n1<=N; n1++)
    for (int n2=-N; n2<=N; n2++)
      for (int n3=-N; n3<=N; n3++) {
	int testval=n1*n1+n2*n2+n3*n3;
	int xx=n1*tau[0]+n2*tau[1]+n3*tau[2];
	int sgn_mult=1-2*abs(xx%2);
	if (testval>=lambda*lambda)
	  out+=sgn_mult*exp(-double(testval)/t);
      }
  
  return out;
}

void Nlam(int& freq, int& min, double lambda)
{
  //Parameters for leading exponential approximation
  //freq is number of vectors in Z^3 with smallest norm >=lambda^2
  //min is that norm squared

  int N=int(ceil(fabs(lambda)));
  min=3*N*N;
  freq=0;
  
  for (int n1=N; n1>=0; n1--)
    for (int n2=N; n2>=0; n2--)
      for (int n3=N; n3>=0; n3--) {
	int testval=n1*n1+n2*n2+n3*n3;
	int mult=1;
	if (n1>0) mult*=2;
	if (n2>0) mult*=2;
	if (n3>0) mult*=2;
	if (testval==min)
	  freq+=mult;
	else if (testval<min && testval>=lambda*lambda) {
	  min=testval;
	  freq=mult;
	}
      }
}
