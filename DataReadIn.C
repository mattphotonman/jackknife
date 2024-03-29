#include <iostream>
#include <stdio.h>
#include <complex>
using namespace std;
#include "data_read_in.h"


/*--------------------------------------------------------------------------
  WMEDataReadIn class
  
  The purpose of this class is just to read in and contain (not manipulate
  in any way) data from a WME or a WME.PIPI format file. N is the number 
  of configurations and Nt is the number of time slices, so "data" is a
  Nt x N array storing the data for a particular quantity at each time
  slice and for each configuration.  "was_read" is a Nt x N matrix
  containing either 0 or 1 in each element indicating whether data has
  been read for that configuration and that time slice.  When a line is
  read from a WME data file, the WMEDataReadIn object has to know if that
  line actually contains data for the particular quantity represented by
  the object.  To that end "col_nums" and "labels" are list of column
  numbers (integers) and strings respectively such that col_num[i] must
  contain the string labels[i] (for all i from 0 to n_labels-1) for the
  line to contain data for the quantity represented by the object.
  "t_col" is the column number that is read to find out for which time
  slice the line is giving data, and "data_col" is the column in which
  the actual data (a double) for the quantity at the given time slice and
  configuration number is found.  "info" is just a short string describing
  the quantity represented by the object for the user's information only
  and does not affect any class methods.

  Once data is read into an object of this class it is meant to be fed to
  a JackknifeTimeSeries object via the ReadData member function of the
  JackknifeTimeSeries class.
--------------------------------------------------------------------------*/

//Constructor
WMEDataReadIn::WMEDataReadIn(int nt, int n, int datcol, int tcol, string in)
{
  //Parameter bounds check
  if (n<2) {
    //Error: I think I'd rather have it abort here.
    cout << "N must be >= 2 for a WMEDataReadIn object, but was given N=" << n << ".  Setting N=2.\n";
    N=2;
  } else {
    N=n;
  }
  if (nt<1) {
    //Error: I think I'd rather have it abort here.
    cout << "Nt must be >= 1 for a WMEDataReadIn object, but was given Nt=" << nt << ".  Setting Nt=1.\n";
    Nt=1;
  } else {
    Nt=nt;
  }
  if (tcol<0) {
    //Error: I think I'd rather have it abort here.
    cout << "t_col must be >=0 for a WMEDataReadIn object, but was given t_col=" << tcol << ".  Setting t_col=2.\n";
    t_col=2;
  } else {
    t_col=tcol;
  }
  if (datcol<0) {
    //Error: I think I'd rather have it abort here.
    cout << "data_col must be >=0 for a WMEDataReadIn object, but was given data_col=" << datcol << ".  Setting data_col=" << t_col+1 << ".\n";
    data_col=t_col+1;
  } else if (datcol==t_col) {
    //Error: I think I'd rather have it abort here.
    cout << "data_col and t_col must be different in a WMEDataReadIn object, but was given " << datcol << " for both.  Setting data_col=" << t_col+1 << ".\n";
    data_col=t_col+1;
  } else {
    data_col=datcol;
  }

  //Set the remaining parameters that were passed to the constructor
  info=in;

  //Initialize everything else
  data=new double* [Nt];
  for (int t=0; t<Nt; t++)
    data[t]=new double [N];
  was_read=new int* [Nt];
  for (int t=0; t<Nt; t++) {
    was_read[t]=new int [N];
    for (int i=0; i<N; i++)
      was_read[t][i]=0;
  }
  col_nums=0; //null pointer
  labels=0; //null pointer
  n_labels=0;
}

//Destructor
WMEDataReadIn::~WMEDataReadIn()
{
  for (int t=0; t<Nt; t++)
    delete [] data[t];
  delete [] data;
  for (int t=0; t<Nt; t++)
    delete [] was_read[t];
  delete [] was_read;
  delete [] col_nums;
  delete [] labels;
}

//Assignment operator
//Deep copy
WMEDataReadIn & WMEDataReadIn::operator=(const WMEDataReadIn& D)
{
  if (this == &D)
    return *this;

  for (int t=0; t<Nt; t++)
    delete [] data[t];
  delete [] data;
  for (int t=0; t<Nt; t++)
    delete [] was_read[t];
  delete [] was_read;
  delete [] col_nums;
  delete [] labels;

  N=D.N;
  Nt=D.Nt;
  n_labels=D.n_labels;
  t_col=D.t_col;
  data_col=D.data_col;
  info=D.info;
  
  data=new double* [Nt];
  for (int t=0; t<Nt; t++) {
    data[t]=new double [N];
    for (int i=0; i<N; i++)
      data[t][i]=D.data[t][i];
  }
  was_read=new int* [Nt];
  for (int t=0; t<Nt; t++) {
    was_read[t]=new int [N];
    for (int i=0; i<N; i++)
      was_read[t][i]=D.was_read[t][i];
  }
  col_nums=new int [n_labels];
  labels=new string [n_labels];
  for (int j=0; j<n_labels; j++) {
    col_nums[j]=D.col_nums[j];
    labels[j]=D.labels[j];
  }
  return *this;
}

//Say if configuration i has been read for timeslice t
int WMEDataReadIn::ReturnIsDefined(int t, int i) const
{
  if (t<0 || t>=Nt || i<0 || i>=N) return 0;
      
  return was_read[t][i];
}
  
//Say if all configurations have been read for timeslice t
int WMEDataReadIn::ReturnIsDefined(int t) const
{
  if (t<0 || t>=Nt)
    return 0;
  for (int i=0; i<N; i++)
    if (!was_read[t][i])
      return 0;
  return 1;
}

//Say if all configurations have been read for all timeslices
int WMEDataReadIn::AllDataRead() const
{
  for (int t=0; t<Nt; t++)
    if (!ReturnIsDefined(t))
      return 0;
  return 1;
}

//Return Nt
int WMEDataReadIn::ReturnNt() const
{
  return Nt;
}

//Return N
int WMEDataReadIn::ReturnN() const
{
  return N;
}

//Return n_labels
int WMEDataReadIn::ReturnNlabels() const
{
  return n_labels;
}

//Return col_nums[i]
int WMEDataReadIn::ReturnColNum(int i) const
{
  if (i<0) {
    //Output to standard error instead of standard out?
    cout << "i must be >= 0 in ReturnColNum.\n";
    return -1;
  } else if (i>=n_labels) {
    //Output to standard error instead of standard out?
    cout << "Only " << n_labels << " col_nums in WMEDataReadIn object.\n";
    return -1;
  }
  return col_nums[i];
}

//Return labels[i]
string WMEDataReadIn::ReturnLabel(int i) const
{
  if (i<0) {
    //Output to standard error instead of standard out?
    cout << "i must be >= 0 in ReturnLabel.\n";
    return "";
  } else if (i>=n_labels) {
    //Output to standard error instead of standard out?
    cout << "Only " << n_labels << " labels in WMEDataReadIn object.\n";
    return "";
  }
  return labels[i];
}

//Return t_col
int WMEDataReadIn::ReturnTCol() const
{
  return t_col;
}

//Return data_col
int WMEDataReadIn::ReturnDataCol() const
{
  return data_col;
}

//Return info string
string WMEDataReadIn::Info() const
{
  return info;
}

//Add a column number to col_nums and corresponding string to labels
WMEDataReadIn & WMEDataReadIn::AddColInfo(int col, string lbl)
{
  if (col<0) {
    //Error: I think I'd rather have it abort here.
    cout << "Error in WMEDataReadIn::AddColInfo : tried to add column number " << col << " which is less than 0.\n";
    return *this;
  }
  int *tmp_col_nums=0;
  string *tmp_labels=0;
  if (n_labels>0) {
    tmp_col_nums=new int [n_labels];
    tmp_labels=new string [n_labels];
    for (int j=0; j<n_labels; j++) {
      tmp_col_nums[j]=col_nums[j];
      tmp_labels[j]=labels[j];
    }
  }
  
  delete [] col_nums;
  delete [] labels;
  
  n_labels++;
  col_nums=new int [n_labels];
  labels=new string [n_labels];
  
  for (int j=0; j<n_labels-1; j++) {
    col_nums[j]=tmp_col_nums[j];
    labels[j]=tmp_labels[j];
  }
  col_nums[n_labels-1]=col;
  labels[n_labels-1]=lbl;

  delete [] tmp_col_nums;
  delete [] tmp_labels;
  
  return *this;
}

//Read a line from a WME or WME.PIPI data file and add the info to the data array
//if this line contains data for the quantity represented by this object (check
//this by comparing with the col_nums and labels arrays).
//Return value is 0 if no data was read from the line and 1 if data was read.
//vrb is a flag for the verbosity (0=off 1=on, default off)
int WMEDataReadIn::ReadLine(int conf, const string & line, int vrb)
{
  if (conf>=N || conf<0) {
    //Error: I think I'd rather have it abort here.
    cout << "Error in WMEDataReadIn::ReadLine :  given configuration number " << conf << " which is outside acceptable range 0<=conf<" << N << ".\n";
    return 0;
  }
  string vals[20];
  int nvals=split(line,vals,20);
  if (nvals<=0) {
    if (vrb)
      cout << "nvals = " << nvals << " <= 0.\n";
    return 0;
  }
  if (n_labels==0) {
    //Error: I think I'd rather have it abort here.
    cout << "Error in WMEDataReadIn::ReadLine :  no labels information.\n";
    return 0;
  }
  //Check that column number col_nums[j] contains the string labels[j] for all j.
  for (int j=0; j<n_labels; j++) {
    if (col_nums[j]>=nvals) {
      if (vrb)
	cout << "col_nums[" << j << "] = " << col_nums[j] << " >= nvals = " << nvals << "\n";
      return 0;
    }
    if (labels[j] != vals[col_nums[j]]) {
      if (vrb)
	cout << "labels[" << j << "] = " << labels[j] << " != vals[" << col_nums[j] << "] = " << vals[col_nums[j]] << "\n";
      return 0;
    }
  }
  if (t_col >= nvals) {
    //Error: I think I'd rather have it abort here.
    cout << "Error in WMEDataReadIn::ReadLine:  t column supposed to be " << t_col << " but max column found was " << nvals-1 << ".\n";
    return 0;
  }
  int t;
  sscanf(vals[t_col].c_str(),"%d",&t);
  if (t<0 || t>=Nt) {
    //Error: I think I'd rather have it abort here.
    cout << "Error in WMEDataReadIn::ReadLine:  read t value " << t << " which is outside of the allowed range 0<=t<" << Nt << ".\n";
    return 0;
  }
  if (data_col >= nvals) {
    //Error: I think I'd rather have it abort here.
    cout << "Error in WMEDataReadIn::ReadLine:  data column supposed to be " << data_col << " but max column found was " << nvals-1 << ".\n";
    return 0;
  }
  double dat;
  sscanf(vals[data_col].c_str(),"%lf",&dat);
  if (was_read[t][conf]) {
    cout << "Already read data for t=" << t << " configuration " << conf << ".\n";
    return 0;
  }
  data[t][conf]=dat;
  was_read[t][conf]=1;
  return 1;
}

//Pass the data from a WMEDataReadIn object to a JackknifeTimeSeries object
//(friend function)
JackknifeTimeSeries ReadWMEData(const WMEDataReadIn & D)
{
  int Nt=D.Nt;
  int N=D.N;
  JackknifeTimeSeries JT(Nt,N);
  for (int t=0; t<Nt; t++) {
    if (D.ReturnIsDefined(t)) {
      for (int i=0; i<N; i++) {
	double jk_val=0.0;
	for (int ii=0; ii<N; ii++)
	  if (ii!=i)
	    jk_val+=D.data[t][ii];
	jk_val/=((double) N-1);
	JT.Jk[t].jk[i]=jk_val;
      }
      double ave_val=0.0;
      for (int i=0; i<N; i++)
	ave_val+=D.data[t][i];
      ave_val/=((double) N);
      JT.Jk[t].ave=ave_val;
      JT.Jk[t].CalcAll();
      JT.is_defined[t]=1;
    }
  }
  return JT;
}
      



//WMECmplxDataReadIn class
//Same as WMEDataReadIn class except it reads in complex numbers.


//Constructor
WMECmplxDataReadIn::WMECmplxDataReadIn(int nt, int n, int redatcol, int imdatcol, int tcol, string in)
{
  //Parameter bounds check
  if (n<2) {
    //Error: I think I'd rather have it abort here.
    cout << "N must be >= 2 for a WMECmplxDataReadIn object, but was given N=" << n << ".  Setting N=2.\n";
    N=2;
  } else {
    N=n;
  }
  if (nt<1) {
    //Error: I think I'd rather have it abort here.
    cout << "Nt must be >= 1 for a WMECmplxDataReadIn object, but was given Nt=" << nt << ".  Setting Nt=1.\n";
    Nt=1;
  } else {
    Nt=nt;
  }
  if (tcol<0) {
    //Error: I think I'd rather have it abort here.
    cout << "t_col must be >=0 for a WMECmplxDataReadIn object, but was given t_col=" << tcol << ".  Setting t_col=2.\n";
    t_col=2;
  } else {
    t_col=tcol;
  }
  if (redatcol<0) {
    //Error: I think I'd rather have it abort here.
    cout << "re_data_col must be >=0 for a WMECmplxDataReadIn object, but was given re_data_col=" << redatcol << ".  Setting re_data_col=" << t_col+1 << ".\n";
    re_data_col=t_col+1;
  } else if (redatcol==t_col) {
    //Error: I think I'd rather have it abort here.
    cout << "re_data_col and t_col must be different in a WMECmplxDataReadIn object, but was given " << redatcol << " for both.  Setting re_data_col=" << t_col+1 << ".\n";
    re_data_col=t_col+1;
  } else {
    re_data_col=redatcol;
  }
  if (imdatcol<0) {
    //Error: I think I'd rather have it abort here.
    int newval;
    if (re_data_col+1==t_col)
      newval=re_data_col+2;
    else
      newval=re_data_col+1;
    cout << "im_data_col must be >=0 for a WMECmplxDataReadIn object, but was given im_data_col=" << imdatcol << ".  Setting im_data_col=" << newval << ".\n";
    im_data_col=newval;
  } else if (imdatcol==t_col || imdatcol==re_data_col) {
    //Error: I think I'd rather have it abort here.
    int newval;
    if (re_data_col+1==t_col)
      newval=re_data_col+2;
    else
      newval=re_data_col+1;
    cout << "im_data_col must not be equal to re_data_col or t_col for a WMECmplxDataReadIn object, but was given im_data_col=" << imdatcol << ", re_data_col=" << re_data_col << ", t_col=" << t_col << ".  Setting im_data_col=" << newval << ".\n";
    im_data_col=newval;
  } else {
    im_data_col=imdatcol;
  } 

  //Set the remaining parameters that were passed to the constructor
  info=in;

  //Initialize everything else
  data=new complex<double>* [Nt];
  for (int t=0; t<Nt; t++)
    data[t]=new complex<double> [N];
  was_read=new int* [Nt];
  for (int t=0; t<Nt; t++) {
    was_read[t]=new int [N];
    for (int i=0; i<N; i++)
      was_read[t][i]=0;
  }
  col_nums=0; //null pointer
  labels=0; //null pointer
  n_labels=0;
}

//Destructor
WMECmplxDataReadIn::~WMECmplxDataReadIn()
{
  for (int t=0; t<Nt; t++)
    delete [] data[t];
  delete [] data;
  for (int t=0; t<Nt; t++)
    delete [] was_read[t];
  delete [] was_read;
  delete [] col_nums;
  delete [] labels;
}

//Assignment operator
//Deep copy
WMECmplxDataReadIn & WMECmplxDataReadIn::operator=(const WMECmplxDataReadIn& D)
{
  if (this == &D)
    return *this;

  for (int t=0; t<Nt; t++)
    delete [] data[t];
  delete [] data;
  for (int t=0; t<Nt; t++)
    delete [] was_read[t];
  delete [] was_read;
  delete [] col_nums;
  delete [] labels;

  N=D.N;
  Nt=D.Nt;
  n_labels=D.n_labels;
  t_col=D.t_col;
  re_data_col=D.re_data_col;
  im_data_col=D.im_data_col;
  info=D.info;
  
  data=new complex<double>* [Nt];
  for (int t=0; t<Nt; t++) {
    data[t]=new complex<double> [N];
    for (int i=0; i<N; i++)
      data[t][i]=D.data[t][i];
  }
  was_read=new int* [Nt];
  for (int t=0; t<Nt; t++) {
    was_read[t]=new int [N];
    for (int i=0; i<N; i++)
      was_read[t][i]=D.was_read[t][i];
  }
  col_nums=new int [n_labels];
  labels=new string [n_labels];
  for (int j=0; j<n_labels; j++) {
    col_nums[j]=D.col_nums[j];
    labels[j]=D.labels[j];
  }
  return *this;
}

//Say if configuration i has been read for timeslice t
int WMECmplxDataReadIn::ReturnIsDefined(int t, int i) const
{
  if (t<0 || t>=Nt || i<0 || i>=N) return 0;
      
  return was_read[t][i];
}
  
//Say if all configurations have been read for timeslice t
int WMECmplxDataReadIn::ReturnIsDefined(int t) const
{
  if (t<0 || t>=Nt)
    return 0;
  for (int i=0; i<N; i++)
    if (!was_read[t][i])
      return 0;
  return 1;
}

//Say if all configurations have been read for all timeslices
int WMECmplxDataReadIn::AllDataRead() const
{
  for (int t=0; t<Nt; t++)
    if (!ReturnIsDefined(t))
      return 0;
  return 1;
}

//Return Nt
int WMECmplxDataReadIn::ReturnNt() const
{
  return Nt;
}

//Return N
int WMECmplxDataReadIn::ReturnN() const
{
  return N;
}

//Return n_labels
int WMECmplxDataReadIn::ReturnNlabels() const
{
  return n_labels;
}

//Return col_nums[i]
int WMECmplxDataReadIn::ReturnColNum(int i) const
{
  if (i<0) {
    //Output to standard error instead of standard out?
    cout << "i must be >= 0 in ReturnColNum.\n";
    return -1;
  } else if (i>=n_labels) {
    //Output to standard error instead of standard out?
    cout << "Only " << n_labels << " col_nums in WMECmplxDataReadIn object.\n";
    return -1;
  }
  return col_nums[i];
}

//Return labels[i]
string WMECmplxDataReadIn::ReturnLabel(int i) const
{
  if (i<0) {
    //Output to standard error instead of standard out?
    cout << "i must be >= 0 in ReturnLabel.\n";
    return "";
  } else if (i>=n_labels) {
    //Output to standard error instead of standard out?
    cout << "Only " << n_labels << " labels in WMECmplxDataReadIn object.\n";
    return "";
  }
  return labels[i];
}

//Return t_col
int WMECmplxDataReadIn::ReturnTCol() const
{
  return t_col;
}

//Return re_data_col
int WMECmplxDataReadIn::ReturnReDataCol() const
{
  return re_data_col;
}

//Return im_data_col
int WMECmplxDataReadIn::ReturnImDataCol() const
{
  return im_data_col;
}

//Return info string
string WMECmplxDataReadIn::Info() const
{
  return info;
}

//Add a column number to col_nums and corresponding string to labels
WMECmplxDataReadIn & WMECmplxDataReadIn::AddColInfo(int col, string lbl)
{
  if (col<0) {
    //Error: I think I'd rather have it abort here.
    cout << "Error in WMECmplxDataReadIn::AddColInfo : tried to add column number " << col << " which is less than 0.\n";
    return *this;
  }
  int *tmp_col_nums=0;
  string *tmp_labels=0;
  if (n_labels>0) {
    tmp_col_nums=new int [n_labels];
    tmp_labels=new string [n_labels];
    for (int j=0; j<n_labels; j++) {
      tmp_col_nums[j]=col_nums[j];
      tmp_labels[j]=labels[j];
    }
  }
  
  delete [] col_nums;
  delete [] labels;
  
  n_labels++;
  col_nums=new int [n_labels];
  labels=new string [n_labels];
  
  for (int j=0; j<n_labels-1; j++) {
    col_nums[j]=tmp_col_nums[j];
    labels[j]=tmp_labels[j];
  }
  col_nums[n_labels-1]=col;
  labels[n_labels-1]=lbl;

  delete [] tmp_col_nums;
  delete [] tmp_labels;
  
  return *this;
}

//Read a line from a WME or WME.PIPI data file and add the info to the data array
//if this line contains data for the quantity represented by this object (check
//this by comparing with the col_nums and labels arrays).
//Return value is 0 if no data was read from the line and 1 if data was read.
//vrb is a flag for the verbosity (0=off 1=on, default off)
int WMECmplxDataReadIn::ReadLine(int conf, const string & line, int vrb)
{
  if (conf>=N || conf<0) {
    //Error: I think I'd rather have it abort here.
    cout << "Error in WMECmplxDataReadIn::ReadLine :  given configuration number " << conf << " which is outside acceptable range 0<=conf<" << N << ".\n";
    return 0;
  }
  string vals[20];
  int nvals=split(line,vals,20);
  if (nvals<=0) {
    if (vrb)
      cout << "nvals = " << nvals << " <= 0.\n";
    return 0;
  }
  if (n_labels==0) {
    //Error: I think I'd rather have it abort here.
    cout << "Error in WMECmplxDataReadIn::ReadLine :  no labels information.\n";
    return 0;
  }
  //Check that column number col_nums[j] contains the string labels[j] for all j.
  for (int j=0; j<n_labels; j++) {
    if (col_nums[j]>=nvals) {
      if (vrb)
	cout << "col_nums[" << j << "] = " << col_nums[j] << " >= nvals = " << nvals << "\n";
      return 0;
    }
    if (labels[j] != vals[col_nums[j]]) {
      if (vrb)
	cout << "labels[" << j << "] = " << labels[j] << " != vals[" << col_nums[j] << "] = " << vals[col_nums[j]] << "\n";
      return 0;
    }
  }
  if (t_col >= nvals) {
    //Error: I think I'd rather have it abort here.
    cout << "Error in WMECmplxDataReadIn::ReadLine:  t column supposed to be " << t_col << " but max column found was " << nvals-1 << ".\n";
    return 0;
  }
  int t;
  sscanf(vals[t_col].c_str(),"%d",&t);
  if (t<0 || t>=Nt) {
    //Error: I think I'd rather have it abort here.
    cout << "Error in WMECmplxDataReadIn::ReadLine:  read t value " << t << " which is outside of the allowed range 0<=t<" << Nt << ".\n";
    return 0;
  }
  if (re_data_col >= nvals) {
    //Error: I think I'd rather have it abort here.
    cout << "Error in WMECmplxDataReadIn::ReadLine:  real data column supposed to be " << re_data_col << " but max column found was " << nvals-1 << ".\n";
    return 0;
  }
  if (im_data_col >= nvals) {
    //Error: I think I'd rather have it abort here.
    cout << "Error in WMECmplxDataReadIn::ReadLine:  imaginary data column supposed to be " << im_data_col << " but max column found was " << nvals-1 << ".\n";
    return 0;
  }
  double redat, imdat;
  sscanf(vals[re_data_col].c_str(),"%lf",&redat);
  sscanf(vals[im_data_col].c_str(),"%lf",&imdat);
  if (was_read[t][conf]) {
    cout << "Already read data for t=" << t << " configuration " << conf << ".\n";
    return 0;
  }
  complex<double> dat(redat,imdat);
  data[t][conf]=dat;
  was_read[t][conf]=1;
  return 1;
}



//split function - splits a string by its white spaces
//not a member or friend function, but it is used by ReadLine.
int split(string line, string vals[], int max_vals)
{
  if (max_vals<1) {
    cout << "Error in split function: Too few elements in output array.\n";
    return -1;
  }
  string white_space(" \t");
  int pos1=line.find_first_not_of(white_space);
  if (pos1==string::npos)
    return 0;
  int pos2;
  int n_val=0;
  while (0==0) {
    pos2=line.find_first_of(white_space,pos1);
    if (pos2==string::npos) {
      vals[n_val].assign(line,pos1,line.length()-pos1);
      n_val++;
      break;
    }
    vals[n_val].assign(line,pos1,pos2-pos1);
    n_val++;
    pos1=line.find_first_not_of(white_space,pos2);
    if (pos1==string::npos)
      break;
    if (n_val>=max_vals) {
      cout << "Error in split function: Too many columns in input line.\n";
      return -1;
    }
  }
  return n_val;
}

//Pass the data from a WMECmplxDataReadIn object to a JackknifeCmplxTimeSeries
//object (friend function)
JackknifeCmplxTimeSeries ReadWMECmplxData(const WMECmplxDataReadIn & D)
{
  int Nt=D.Nt;
  int N=D.N;
  JackknifeCmplxTimeSeries JT(Nt,N);
  for (int t=0; t<Nt; t++) {
    if (D.ReturnIsDefined(t)) {
      for (int i=0; i<N; i++) {
	complex<double> jk_val=0.0;
	for (int ii=0; ii<N; ii++)
	  if (ii!=i)
	    jk_val+=D.data[t][ii];
	jk_val/=((double) N-1);
	JT.Jk[t].jk[i]=jk_val;
      }
      complex<double> ave_val=0.0;
      for (int i=0; i<N; i++)
	ave_val+=D.data[t][i];
      ave_val/=((double) N);
      JT.Jk[t].ave=ave_val;
      JT.Jk[t].CalcAll();
      JT.is_defined[t]=1;
    }
  }
  return JT;
}

