#include <iostream>
#include <cmath>
#include <cstdlib>
using namespace std;
#include "special_functions.h"


/* modified bessel function; needs also I1 */
double K1( double x) {
	double temp1,temp2,y;
	double a2,a3,a4,a5,a6,a8,a10,a12 ;
	if ( x > 2. ) {
           temp1 = 2./x;
           temp2 = 1.25331414 + .23498619*temp1;
           a2 = temp1*temp1 ;
           a3 = temp1*a2 ;
	   a4 = a2*a2 ;
           a5 = temp1*a4 ;
	   a6 = a2*a4 ;
           temp2 = temp2 -.03655620*a2 +.01504268*a3;
           temp2 = temp2 -.00780353*a4 +.00325614*a5;
           temp2 = temp2 -.00068245*a6;
           y = temp2/exp(x)/sqrt(x);
	}
	else if ( x > 0. ) {
           temp1 = x/2. ;
           a2 = temp1*temp1 ;
	   a4 = a2*a2 ;
	   a6 = a2*a4 ;
	   a8 = a2*a6 ;
	   a10 = a2*a8 ;
	   a12 = a2*a10 ;
           temp2 = x*log(temp1)*I1(x) + 1.0 + .15443144*a2 ;
           temp2 = temp2 - .67278579*a4 - .18156897*a6 ;
           temp2 = temp2 - .01919402*a8 - .00110404*a10 ;
           temp2 = temp2 - .00004686*a12 ;
           y = temp2/x ;
	}
        else {
	  cout << "sorry, don't have bessel function K1\n";
	  cout << " for argument x= " << x << "\n";
	  exit(1);
        }
	return(y);
}


double I1( double x) {
	double temp1,temp2,y;
	double a2,a4,a6,a8,a10,a12 ;
/* modified bessel function */
	if ( x <= 3.75 && -x <= 3.75 ) {
           temp1 = x/3.75 ;
           a2 = temp1*temp1 ;
	   a4 = a2*a2 ;
	   a6 = a2*a4 ;
	   a8 = a2*a6 ;
	   a10 = a2*a8 ;
	   a12 = a2*a10 ;
           temp2 = 0.5 +.87890594*a2 +.51498869*a4 ;
           temp2 = temp2 + .15084934*a6 + .02658733*a8 ;
           temp2 = temp2 + .00301532*a10 +.00032411*a12 ;
           y = x*temp2 ;
	}
        else {
	  cout << "sorry, don't have bessel function I1\n";
	  cout << " for argument x= " << x << "\n";
	  exit (1);
	}
      return(y) ;
}



/* modified bessel function; needs also I0 */
double K0( double x) {
	double temp1,temp2,y;
	double a2,a3,a4,a5,a6,a8,a10,a12 ;
	if ( x > 2. ) {
           temp1 = 2./x ;
           a2 = temp1*temp1 ;
           a3 = temp1*a2 ;
	   a4 = a2*a2 ;
           a5 = temp1*a4 ;
	   a6 = a2*a4 ;
           temp2 = 1.25331414 - .07832358*temp1 ;
           temp2 = temp2 +.02189568*a2 -.01062446*a3 ;
           temp2 = temp2 +.00587872*a4 -.00251540*a5 ;
           temp2 = temp2 +.00053208*a6  ;
           y = temp2/exp(x)/sqrt(x) ;
	}
	else if ( x > 0. ) {
           temp1 = x/2. ;
           a2 = temp1*temp1 ;
	   a4 = a2*a2 ;
	   a6 = a2*a4 ;
	   a8 = a2*a6 ;
	   a10 = a2*a8 ;
	   a12 = a2*a10 ;
           temp2 = -log(temp1)*I0(x)-.57721566 + .42278420*a2 ;
           temp2 = temp2 + .23069756*a4 + .03488590*a6 ;
           temp2 = temp2 + .00262698*a8 + .00010750*a10 ;
           temp2 = temp2 + .00000740*a12 ;
           y = temp2 ;
	}
	else {
	  cout << "sorry, don't have bessel function K0\n";
	  cout << " for argument x= " << x << "\n";
	  exit (1);
	}
      return (y);
}

/* modified bessel function;  */
double I0( double x) {
	double temp1,temp2,y;
	double a2,a4,a6,a8,a10,a12 ;
	if ( x <= 3.75 && -x <= 3.75 ) {
           temp1 = x/3.75 ;
           a2 = temp1*temp1 ;
	   a4 = a2*a2 ;
	   a6 = a2*a4 ;
	   a8 = a2*a6 ;
	   a10 = a2*a8 ;
	   a12 = a2*a10 ;
           temp2 = 1.0 +3.5156229*a2 +3.0899424*a4 ;
           temp2 = temp2 + 1.2067492*a6 + .2659732*a8 ;
           temp2 = temp2 + .0360768*a10 +.0045813*a12 ;
           y = temp2 ;
	}
	else {
	  cout << "sorry, don't have bessel function I0\n";
	  cout << " for argument x= " << x << "\n";
	  exit (1);
        }
        return (y);
       
}

//Bessel K_n function for n>=2.  Uses K1 and K0.
double besselK(int n, double x)
{
  int j;
  double bk,bkm,bkp,tox;

  if (n==0)
    return K0(x);
  else if (n==1)
    return K1(x);
  else if (n<0) {
    cout << "Index n less than 0 in besselK.";
    exit(1);
  }
  
  tox=2.0/x;
  bkm=K0(x); 
  bk=K1(x); 
  for (j=1;j<n;j++) {
    bkp=bkm+j*tox*bk; 
    bkm=bk; 
    bk=bkp;
  } 

  return bk;  
}
