#ifndef INCLUDE_LEASTSQR
#define INCLUDE_LEASTSQR

void linear_least_square(int Nvar, int Ncoeff, int Ndata, double** y, double** sigma, double*** f, double C[]);
void cov_linear_least_square(int Nvar, int Ncoeff, int Ndata, double** y, double**** Cinv, double*** f, double Cf[]);
void ludcmp(double** a, int n, int *indx, double* d);
void lubksb(double** a, int n, int *indx, double b[]);
void invert_mat(double** mat, double** mat_inv, int n);
double chisq(int Nvar, int Ncoeff, int Ndata, double** y, double** sigma, double*** f, double C[]);
double chisq_per_dof(int Nvar, int Ncoeff, int Ndata, double** y, double** sigma, double*** f, double C[]);
double cov_chisq(int Nvar, int Ncoeff, int Ndata, double** y, double**** Cinv, double*** f, double Cf[]);
double cov_chisq_per_dof(int Nvar, int Ncoeff, int Ndata, double** y, double**** Cinv, double*** f, double Cf[]);

#endif
