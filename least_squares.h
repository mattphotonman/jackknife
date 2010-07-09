#ifndef INCLUDE_LEASTSQR
#define INCLUDE_LEASTSQR

void linear_least_square(int Nvar, int Ncoeff, int Ndata, double** y, double** sigma, double*** f, double C[]);
void ludcmp(double** a, int n, int *indx, double* d);
void lubksb(double** a, int n, int *indx, double b[]);
double chisq(int Nvar, int Ncoeff, int Ndata, double** y, double** sigma, double*** f, double C[]);
double chisq_per_dof(int Nvar, int Ncoeff, int Ndata, double** y, double** sigma, double*** f, double C[]);

#endif
