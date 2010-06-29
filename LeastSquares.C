#include <iostream>
#include <cmath>
#include "nrutil.h"
#include "least_squares.h"
#define TINY 1.0e-20  //A small number.
using namespace std;

/*--------------------------------------------
Least squares fit to

y_alpha=sum_beta C_beta*f_{alpha beta}(x)

where x is the independent variable(s), C_beta are the coefficients to be
fit, y_alpha are the dependent variables, and f_{alpha beta} are the
functions we have chosen to use for the fitting form.  This is done by
minimizing the following function:

sum_{alpha,i}(y^i_alpha-sum_beta C_beta*f_{alpha beta}(x^i))^2/(sigma^i_alpha)^2

where i indexes the list of actual data points.

The linear_least_square requires the following inputs:

Nvar=number of independent variables (alpha=0,...,Nvar-1)
Ncoeff=number of coefficients to be fit (beta=0,...,Ncoeff-1)
Ndata=number of data points (i=0,...,Ndata-1)
y[alpha][i]=y^i_alpha
sigma[alpha][i]=sigma^i_alpha
f[alpha][beta][i]=f_{alpha beta}(x^i)

The return is

C[beta]=C_beta

which are the fitted coefficients.
--------------------------------------------*/
void linear_least_square(int Nvar, int Ncoeff, int Ndata, double** y, double** sigma, double*** f, double C[])
{
  int* indx=new int [Ncoeff+1];
  double d;
  double** a;
  a=new double* [Ncoeff+1];
  for (int gamma=0; gamma<Ncoeff; gamma++)
    a[gamma+1]=new double [Ncoeff+1];
  double* b=new double [Ncoeff+1];

  for (int gamma=0; gamma<Ncoeff; gamma++)
    b[gamma+1]=0.0;
  
  for (int gamma=0; gamma<Ncoeff; gamma++)
    for (int beta=0; beta<Ncoeff; beta++)
      a[gamma+1][beta+1]=0.0;
  
  for (int alpha=0; alpha<Nvar; alpha++)
    for (int gamma=0; gamma<Ncoeff; gamma++)
      for (int beta=0; beta<Ncoeff; beta++)
	for (int i=0; i<Ndata; i++)
	  a[gamma+1][beta+1]+=f[alpha][gamma][i]*f[alpha][beta][i]/(sigma[alpha][i]*sigma[alpha][i]);
  
  for (int alpha=0; alpha<Nvar; alpha++)
    for (int gamma=0; gamma<Ncoeff; gamma++)
      for (int i=0; i<Ndata; i++)
	b[gamma+1]+=y[alpha][i]*f[alpha][gamma][i]/(sigma[alpha][i]*sigma[alpha][i]);

  //ATTENTION: Notice that Numerical Recipes numbers vectors and matrices
  //from 1 to n rather than 0 to n-1, hence the indexing and sizes of
  //a, b, and indx above.

  ludcmp(a,Ncoeff,indx,&d);
  lubksb(a,Ncoeff,indx,b);
  
  for (int beta=0; beta<Ncoeff; beta++)
    C[beta]=b[beta+1];

  //Free memory allocated in this function.
  delete [] indx;
  for (int gamma=0; gamma<Ncoeff; gamma++)
    delete [] a[gamma+1];
  delete [] a;
  delete [] b;
}

/*--------------------------------------------
  Taken from section 2.3 for Numerical Recipes in C.

  ATTENTION:  Note that Numerical Recipes numbers matrices and vectors
  from 1 to n rather than 0 to n-1.

  Given a matrix a[1..n][1..n], this routine replaces it by the LU
  decomposition of a rowwise permutation of itself.  a and n are input.
  a is output, arranged as in equation (2.3.14) above; indx[1..n] is an
  output vector that records teh row permutation effected by the partial
  pivoting; d is output as +/-1 depending on whether the number of row
  interchanges was even or odd, respectively.  This routine is used in
  combination with lubksb to solve linear equations or invert a matrix.
--------------------------------------------*/
void ludcmp(double** a, int n, int *indx, double* d)
{
  int i, imax, j, k;
  double big, dum, sum, temp;
  double *vv;  //vv stores the implicit scaling of each row

  vv=vector(1,n);
  *d=1.0;  //No row interchanges yet.
  for (i=1; i<=n; i++) {  //Loop over rows to get the implicit
                          //scaling information.
    big=0.0;
    for (j=1;j<=n;j++)
      if ((temp=fabs(a[i][j])) > big) big=temp;
    if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
                           //No nonzero largest element.
    vv[i]=1.0/big;  //Save the scaling.
  }
  for (j=1; j<=n; j++) {  //This is the loop over columns of Crout's method.
    for (i=1;i<j;i++) {  //This is equation (2.3.12) except i=j.
      sum=a[i][j];
      for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;  //Initialize for the search for largest pivot element.
    for (i=j;i<=n;i++) {  //This is i=j of equation (2.3.12) and i=j+1...N
                          //of equation (2.3.13).
      sum=a[i][j];
      for (k=1;k<j;k++)
	sum-=a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
	//Is the figure of merit for the pivot better than the
	//best so far?
	big=dum;
	imax=i;
      }
    }
    if (j != imax) {  //Do we need to interchange rows?
      for (k=1;k<=n;k++) {  //Yes, do so...
	dum=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      *d = -(*d);  //...and change the parity of d.
      vv[imax]=vv[j];  //Also interchange the scale factor.
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=TINY;
    //If the pivot element is zero the matrix is singular (at least
    //to the precision of the algorithm).  For some applications on
    //singular matrix, it is desirable to substitute TINY for zero.
    if (j != n) {  //Now, finally, divide by the pivot element.
      dum=1.0/(a[j][j]);
      for (i=j+1;i<=n;i++) a[i][j]*=dum;
    }
  }  //Go back for the next column in the reduction.
  free_vector(vv,1,n);
}

/*--------------------------------------------
  Taken from section 2.3 of Numerical Recipes in C.

  ATTENTION:  Note that Numerical Recipes numbers matrices and vectors
  from 1 to n rather than 0 to n-1.

  Solves the set of n linear equations AX=B.  Here a[1..n][1..n] is input,
  not as the matrix A but rather as its LU decomposition, determined by the
  routine ludcmp.  indx[1..n] is input as the permutation vector returned by
  ludcmp.  b[1..n] is input as the right-hand side vector B, and returns
  with the solution vector X.  a, n, and indx are not modified by this
  routine and can be left in place for successive calls with different
  right-hand sides b.  This routine takes into account the possibility that
  b will begin with many zero elements, so it is efficient for use in
  matrix inversion.
--------------------------------------------*/
void lubksb(double** a, int n, int *indx, double b[])
{
  int i, ii=0, ip, j;
  double sum;
  
  for (i=1;i<=n;i++) {  //When ii is set to a positive value, it will become
                        //the index of the first nonvanishing element of b.
                        //We now do the forward substitution, equation (2.3.6).
                        //The only new wrinkle is to unscramble the permutation
                        //as we go.
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii)
      for (j=ii;j<=i-1;j++) sum-=a[i][j]*b[j];
    else if (sum) ii=i;  //A nonzero element was encountered, so from now on
                         //we will have to do the sums in the loop above.
    b[i]=sum;
  }
  for (i=n;i>=1;i--) {  //Now we do backsubstitution, equation (2.3.7).
    sum=b[i];
    for (j=i+1;j<=n;j++) sum-=a[i][j]*b[j];
    b[i]=sum/a[i][i];  //Store a component of the solution vector X.
  }  //All done!
}
