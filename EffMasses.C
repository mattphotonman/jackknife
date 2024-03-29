#include <iostream>
#include <cmath>
using namespace std;
#include "eff_masses.h"

double meff_two_point(double C1, double C2, int t, int Nt)
{
  int n=Nt-t;
  if (n<0) {
    cout << "Error: Must have t<=Nt in meff_two_point.\n";
    return 0.0;
  }
  double ratio=C2/C1;
  if (ratio>1 || ratio<=0) {
    cout << "Error: Invalid argument to meff_two_point.\n";
    return 0.0;
  }
  double a=1.0/ratio;
  double epsilon=1E-9;
  int max_itr=1000;
  int reached_max_itr=0;
  double f_val, df_val, change;
  int i;
  for (i=0; i<max_itr; i++) {
    f_val=f_func(a,n);
    df_val=df_func(a,n);
    change=(ratio-f_val)/df_val;
    if (fabs(change)<=epsilon) break;
    a+=change;
    if (i==max_itr-1) reached_max_itr=1;
  }
  if (reached_max_itr)
    cout << "Error: Exceeded max number of iterations meff_two_point.\n";
  return log(a);
}

double meff_three_point(double C1, double C2, double C3, int t, int Nt)
{
  int n=Nt-t;
  if (n<1) {
    cout << "Error: Must have t<=Nt-1 in meff_three_point.\n";
    return 0.0;
  }
  double ratio=(C3-C2)/(C2-C1);
  double max_ratio=(2.0*double(n)-1.0)/(2.0*double(n)+1.0);
  if (ratio>max_ratio || ratio<=0) {
    cout << "Error: Invalid argument to meff_three_point.\n";
    return 0.0;
  }
  double a;
  if (ratio/max_ratio>0.9)
    a=1.0+sqrt(3.0/double(n)*(1.0-ratio/max_ratio));
  else
    a=1.0/ratio;
  double epsilon=1E-9;
  int max_itr=1000;
  int reached_max_itr=0;
  double g_val, dg_val, change;
  int i;
  for (i=0; i<max_itr; i++) {
    g_val=g_func(a,n);
    dg_val=dg_func(a,n);
    change=(ratio-g_val)/dg_val;
    if (fabs(change)<=epsilon) break;
    a+=change;
    if (i==max_itr-1) reached_max_itr=1;
  }
  if (reached_max_itr)
    cout << "Error: Exceeded max number of iterations meff_three_point.\n";
  return log(a);
}  

double f_func(double a, int n)
{
  double pow1=1.0, pow2;
  for (int i=0; i<2*n+1; i++)
    pow1*=a;
  pow2=pow1*a;
  return (pow1+a)/(pow2+1.0);
}

double df_func(double a, int n)
{
  double pow1=1.0, pow2=1.0;
  for (int i=0; i<2*n+1; i++)
    pow1*=a;
  for (int i=0; i<4*n+3; i++)
    pow2*=a;
  double num=2.0*n*pow1-pow2+pow1-pow1*a*a+a-2.0*n*pow1*a*a;
  double denom=a*(pow1*a+1.0)*(pow1*a+1.0);
  return num/denom;
}

double g_func(double a, int n)
{
  double pow=1.0;
  for (int i=0; i<2*n; i++)
    pow*=a;
  return (pow-a)/(pow*a-1.0);
}

double dg_func(double a, int n)
{
  double pow1=1.0;
  double pow2=1.0;
  for (int i=0; i<2*n; i++)
    pow1*=a;
  for (int i=0; i<4*n+1; i++)
    pow2*=a;
  double num=-2.0*pow1*n+a+2.0*n*pow1*a*a-pow2;
  double denom=a*(pow1*a-1.0)*(pow1*a-1.0);
  return num/denom;
}
