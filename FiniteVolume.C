#include <iostream>
#include <cmath>
#include <cstdlib>
using namespace std;
#include "finite_volume.h"
#include "special_functions.h"

//Returns the multiplicity.  Use this function in other functions
//instead of calculating a mults array.  This way the mults array can
//be calculated internally to the return_mult function (using a
//static variable), and used anywhere, instead of each function
//calculating its own mults array.  The number of multiplicities to
//calculate initially is set inside the return_mult function by the
//variable nsq_last, and the function will throw an error if you ask
//it for a multiplicity for nsq higher than this.
//Choose nsq_last here, and if it is found that larger multiplicities
//are commonly needed, nsq_last should be increased here.
int return_mult(int nsq)
{
  const int nsq_last=100;
  if (nsq>nsq_last) {
    cout << "return_mult:  Error, only have mults up to nsq = " << nsq_last << " and you asked for nsq = " << nsq << ".  Abort.\n";
    exit(1);
  }

  if (nsq<0) {
    cout << "return_mult:  Error, nsq<0.  Abort.\n";
    exit(1);
  }
  
  static int mults[nsq_last+1];
  static int mults_calc_done=0;
  if (!mults_calc_done) {
    calc_mults(mults,nsq_last);
    mults_calc_done=1;
  }
  
  return mults[nsq];
}

//This function is called only by return_mult to avoid calculating
//the same multiplicities more than once in different places.
//NOTE that the second argument isn't the size of the array
//but is size(mults)-1 (the last index of the array).
void calc_mults(int mults[],int nsq_last)
{
  //Parameter check
  if (nsq_last<0) {
    cout << "calc_mults:  Error, nsq_last must be >=0.  Abort.\n";
    exit(1);
  }

  for (int nsq=0; nsq<=nsq_last; nsq++)
    mults[nsq]=0;
  
  int N=int(sqrt(double(nsq_last)));
  for (int nx=-N; nx<=N; nx++)
    for (int ny=-N; ny<=N; ny++)
      for (int nz=-N; nz<=N; nz++) {
	int nsq=nx*nx+ny*ny+nz*nz;
	if (nsq<=nsq_last)
	  mults[nsq]++;
      }
  
}

double delta1(double x)
{
  const int nsq_max=100;  //Should be <=nsq_last in return_mult function.
  const double epsilon=1.0E-9;
  
  double sum=0.0, term;
  for (int nsq=1; nsq<=nsq_max; nsq++) {
    term=return_mult(nsq)*K1(sqrt(double(nsq))*x)/sqrt(double(nsq));
    if (sum!=0.0)
      if (fabs(term/sum)<epsilon) break;
    sum+=term;
    
    if (nsq==nsq_max)
      cout << "delta1:  Warning, reached nsq_max and didn't converge.\n";
  }
  
  return 4.0/x*sum;
}
