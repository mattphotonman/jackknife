#ifndef INCLUDE_EFF_MASSES
#define INCLUDE_EFF_MASSES

double meff_two_point(double C1, double C2, int t, int Nt);
double meff_three_point(double C1, double C2, double C3, int t, int Nt);
double f_func(double a, int n);
double df_func(double a, int n);
double g_func(double a, int n);
double dg_func(double a, int n);

#endif
