#include <iostream>
#include <fstream>
#include <stdio.h>
using namespace std;
#include "data_read_in_hisq.h"
#include "data_read_in.h"


/*--------------------------------------------------------------------------
  We provide a class Ensemb for entering the information about an
  ensemble, entering the data points we wish to take from that ensemble,
  and reading in data from a file, including the covariance matrix. Some
  non-member, non-friend functions take an array of Ensemb objects
  and help with the process of reading a data file into all of them,
  as well as outputting the result to a big covariance matrix.
--------------------------------------------------------------------------*/

//Constructor
Ensemb::Ensemb(double mL, double mS, double mC, double beta, double L)
{
  info=new double [5];
  info[0]=mL;
  info[1]=mS;
  info[2]=mC;
  info[3]=beta;
  info[4]=L;

  N_data=0;
  //Set arrays to null pointers until data points are added.
  data_x=0;
  data_y=0;
  cov_y=0;
  data_was_read=0;
  cov_was_read=0;
}

//Destructor
Ensemb::~Ensemb()
{
  //Sizes of arrays
  //double info[5];
  //double data_x[N_data][3]
  //double data_y[N_data]
  //double cov_y[N_data][N_data]
  //int data_was_read[N_data]
  //int cov_was_read[N_data][N_data]

  for (int i=0; i<N_data; i++) {
    delete [] data_x[i];
    delete [] cov_y[i];
    delete [] cov_was_read[i];
  }
  delete [] info;
  delete [] data_x;
  delete [] data_y;
  delete [] cov_y;
  delete [] data_was_read;
  delete [] cov_was_read;
}

//Say if y values and covariances have been read for all data points.
int Ensemb::AllDataCovRead() const
{
  return AllDataRead() && AllCovRead();
}

//Say if y values have been read for all data points.
int Ensemb::AllDataRead() const
{
  int allread=1;
  for (int i=0; i<N_data; i++)
    allread=allread && data_was_read[i];
  
  return allread;
}

//Say if covariances have been read for all data points.
int Ensemb::AllCovRead() const
{
  int allread=1;
  for (int i=0; i<N_data; i++)
    for (int j=0; j<N_data; j++)
      allread=allread && cov_was_read[i][j];

  return allread;
}

//Return N_data
int Ensemb::ReturnNdata() const
{
  return N_data;
}

//Return info[col]
double Ensemb::ReturnInfo(int col) const
{
  //Bounds checking
  if (col<0 || col>=5) {
    cout << "Ensemb::ReturnInfo:  Error, col must be >=0 and < 5.\n";
    return -1;
  }
  
  return info[col];
}

//Return data_x[i][col]
double Ensemb::ReturnDataX(int i, int col) const
{
  //Bounds checking
  if (i<0 || i>=N_data) {
    cout << "Ensemb::ReturnDataX:  Error, i must be >=0 and < N_data.\n";
    return -1;
  }
  if (col<0 || col>=3) {
    cout << "Ensemb::ReturnDataX:  Error, col must be >=0 and < 3.\n";
    return -1;
  }

  return data_x[i][col];
}

//Return data_y[i]
double Ensemb::ReturnDataY(int i) const
{
  //Bounds checking
  if (i<0 || i>=N_data) {
    cout << "Ensemb::ReturnDataY:  Error, i must be >=0 and < N_data.\n";
    return -1;
  }

  return data_y[i];
}

//Return cov_y[i][j] 
double Ensemb::ReturnCovY(int i, int j) const
{
  //Bounds checking
  if (i<0 || i>=N_data) {
    cout << "Ensemb::ReturnCovY:  Error, i must be >=0 and < N_data.\n";
    return -1;
  }
  //Bounds checking
  if (j<0 || j>=N_data) {
    cout << "Ensemb::ReturnCovY:  Error, j must be >=0 and < N_data.\n";
    return -1;
  }

  return cov_y[i][j];
}

//Return data_was_read[i]
int Ensemb::ReturnDataWasRead(int i) const
{
  //Bounds checking
  if (i<0 || i>=N_data) {
    cout << "Ensemb::ReturnDataWasRead:  Error, i must be >=0 and < N_data.\n";
    return -1;
  }
  
  return data_was_read[i];
}

//Return cov_was_read[i][j]
int Ensemb::ReturnCovWasRead(int i, int j) const
{
  //Bounds checking
  if (i<0 || i>=N_data) {
    cout << "Ensemb::ReturnCovWasRead:  Error, i must be >=0 and < N_data.\n";
    return -1;
  }
  //Bounds checking
  if (j<0 || j>=N_data) {
    cout << "Ensemb::ReturnCovWasRead:  Error, j must be >=0 and < N_data.\n";
    return -1;
  }

  return cov_was_read[i][j];
} 

//Add a new data point to be read for this ensemble by specifying its
//x values (flag, mA, mB).
void Ensemb::AddDataPt(double flag, double mA, double mB)
{
  //Sizes of arrays
  //double info[5];
  //double data_x[N_data][3]
  //double data_y[N_data]
  //double cov_y[N_data][N_data]
  //int data_was_read[N_data]
  //int cov_was_read[N_data][N_data]

  //Copy existing arrays
  double** old_data_x;
  double* old_data_y;
  double** old_cov_y;
  int* old_data_was_read;
  int** old_cov_was_read;

  if (N_data==0) {

    old_data_x=0;
    old_data_y=0;
    old_cov_y=0;
    old_data_was_read=0;
    old_cov_was_read=0;

  } else {

    old_data_x=new double* [N_data];
    for (int i=0; i<N_data; i++) {
      old_data_x[i]=new double [3];
      for (int col=0; col<3; col++)
	old_data_x[i][col]=data_x[i][col];
    }

    old_data_y=new double [N_data];
    for (int i=0; i<N_data; i++)
      old_data_y[i]=data_y[i];
    
    old_cov_y=new double* [N_data];
    for (int i=0; i<N_data; i++) {
      old_cov_y[i]=new double [N_data];
      for (int j=0; j<N_data; j++)
	old_cov_y[i][j]=cov_y[i][j];
    }
    
    old_data_was_read=new int [N_data];
    for (int i=0; i<N_data; i++)
      old_data_was_read[i]=data_was_read[i];
    
    old_cov_was_read=new int* [N_data];
    for (int i=0; i<N_data; i++) {
      old_cov_was_read[i]=new int [N_data];
      for (int j=0; j<N_data; j++)
	old_cov_was_read[i][j]=cov_was_read[i][j];
    }

  }

  //Delete existing arrays.
  for (int i=0; i<N_data; i++) {
    delete [] data_x[i];
    delete [] cov_y[i];
    delete [] cov_was_read[i];
  }
  delete [] data_x;
  delete [] data_y;
  delete [] cov_y;
  delete [] data_was_read;
  delete [] cov_was_read;


  //Update all of the arrays with the information of the new
  //point.

  //Updata N_data
  N_data++;

  //Update data_x
  data_x=new double* [N_data];
  //Copy over old data
  for (int i=0; i<N_data-1; i++) {
    data_x[i]=new double [3];
    for (int col=0; col<3; col++)
      data_x[i][col]=old_data_x[i][col];
  }
  //Enter new data
  data_x[N_data-1]=new double [3];
  data_x[N_data-1][0]=flag;
  data_x[N_data-1][1]=mA;
  data_x[N_data-1][2]=mB;
  
  //Update data_y
  data_y=new double [N_data];
  //Copy old data
  for (int i=0; i<N_data-1; i++)
    data_y[i]=old_data_y[i];
  //No new data to enter here
  
  //Update cov_y
  cov_y=new double* [N_data];
  //Copy old data
  for (int i=0; i<N_data-1; i++) {
    cov_y[i]=new double [N_data];
    for (int j=0; j<N_data-1; j++)
      cov_y[i][j]=old_cov_y[i][j];
  }
  //No new data to enter here
  cov_y[N_data-1]=new double [N_data];
  
  //Update data_was_read
  data_was_read=new int [N_data];
  //Copy old data
  for (int i=0; i<N_data-1; i++)
    data_was_read[i]=old_data_was_read[i];
  //Enter new data
  data_was_read[N_data-1]=0;
  
  //Update cov_was_read
  cov_was_read=new int* [N_data];
  //Copy old data
  for (int i=0; i<N_data-1; i++) {
    cov_was_read[i]=new int [N_data];
    for (int j=0; j<N_data-1; j++)
      cov_was_read[i][j]=old_cov_was_read[i][j];
    //Enter new data
    cov_was_read[i][N_data-1]=0;
  }  
  //Enter new data
  cov_was_read[N_data-1]=new int [N_data];
  for (int j=0; j<N_data; j++)
    cov_was_read[N_data-1][j]=0;

  //Delete memory allocated inside this function.
  for (int i=0; i<N_data-1; i++) {
    delete [] old_data_x[i];
    delete [] old_cov_y[i];
    delete [] old_cov_was_read[i];
  }
  delete [] old_data_x;
  delete [] old_data_y;
  delete [] old_cov_y;
  delete [] old_data_was_read;
  delete [] old_cov_was_read;
 
}

//Read in a line from a data file.  Have to check the number of columns
//to see if the line is giving you a data value, a covariance matrix
//element, or neither.  Then add that data in appropriately.
int Ensemb::ReadLine(const string & line)
{
  string vals[20];
  int nvals=split(line,vals,20);
  if (nvals<0) {
    return 0;
  }

  if (nvals!=9 && nvals!=10 && nvals!=17) {
    //Wrong number of columns for a data point or covariance
    //matrix element.  Just means that this line is a header
    //or something else, so skip it.
    return 1;
    //Note that this function only returns 0 if there is some 
    //sort of error.  Encountering a non-data line is not an error.
  }

  if (nvals==9 || nvals==10) {
    //Line containing a data point.
    //(Sometimes the error is included, in which case there are
    //10 columns.  The error is not actually read here since it
    //can be taken from the covariance matrix.)
    
    //Read ensemble info of this data point.
    double pt_info[5];
    for (int col=0; col<5; col++)
      pt_info[col]=atof(vals[col+3].c_str());

    //Check if it is in the ensemble of this object.
    int this_ensemb=1;
    for (int col=0; col<5; col++)
      this_ensemb=this_ensemb && (pt_info[col]==ReturnInfo(col));
    if (!this_ensemb)
      return 1;
    
    //Read x values of this data point.
    double pt_data_x[3];
    for (int col=0; col<3; col++)
      pt_data_x[col]=atof(vals[col].c_str());
    
    //Check which data point in this object has these x values,
    //if any.
    int i;
    for (i=0; i<N_data; i++) {
      int this_pt=1;
      for (int col=0; col<3; col++)
	this_pt=this_pt && (pt_data_x[col]==ReturnDataX(i,col));
      if (this_pt) break;
    }
    if (i==N_data) {
      //Data point not found in object, means we don't want to
      //read it.
      return 1;
    }

    //If we get here then the data point read is one for which we
    //want to collect info.
    
    //Check that we haven't already read something for this data
    //point.
    if (ReturnDataWasRead(i)) {
      cout << "Ensemb::ReadLine: Tried to read a data point for which data has already been read.\n";
      return 0;
    }
    
    //Read it in.
    data_y[i]=atof(vals[8].c_str());
    data_was_read[i]=1;
    
    return 1;
  }

  //If you get here then nvals=17 and this line contains a
  //covariance matrix element.
  
  //Read ensemble info of the two data points for which this line
  //gives the covariance.
  double pt_info[5], pt2_info;
  for (int col=0; col<5; col++) {
    pt_info[col]=atof(vals[col+3].c_str());
    pt2_info=atof(vals[col+11].c_str());
    //Only looking at covariance matrix elements within this ensemble,
    //so don't read this line if pt2_info!=pt_info[col].  (Covariance
    //matrix elements between ensembles should be 0 anyways if they're 
    //even written down.)
    if (pt2_info!=pt_info[col])
      return 1;
  }

  //Check if this line has a cov mat element for the ensemble of
  //this object.
  int this_ensemb=1;
  for (int col=0; col<5; col++)
    this_ensemb=this_ensemb && (pt_info[col]==ReturnInfo(col));
  if (!this_ensemb)
    return 1;
  
  //Read x values of the two data points.
  double pt1_data_x[3], pt2_data_x[3];
  for (int col=0; col<3; col++) {
    pt1_data_x[col]=atof(vals[col].c_str());
    pt2_data_x[col]=atof(vals[col+8].c_str());
  }

  //Check which of the data points have these x values, if any.
  int i, j;
  for (i=0; i<N_data; i++) {
    int this_pt=1;
    for (int col=0; col<3; col++)
      this_pt=this_pt && (pt1_data_x[col]==ReturnDataX(i,col));
    if (this_pt) break;
  }
  if (i==N_data) {
    //Data point not found in object, means we don't want to
    //read it.
    return 1;
  }
  for (j=0; j<N_data; j++) {
    int this_pt=1;
    for (int col=0; col<3; col++)
      this_pt=this_pt && (pt2_data_x[col]==ReturnDataX(j,col));
    if (this_pt) break;
  }
  if (j==N_data) {
    //Data point not found in object, means we don't want to
    //read it.
    return 1;
  }

  //If we get here then both the data points read are ones for which
  //we want to collect info, thus we should read this cov mat element.

  //Check that we haven't already read something for this cov mat
  //element.
  if (ReturnCovWasRead(i,j)) {
    cout << "Ensemb::ReadLine: Tried to read a covariance matrix element for which a value has already been read.\n";
    return 0;
  }

  //Read it in.
  cov_y[i][j]=atof(vals[16].c_str());
  cov_was_read[i][j]=1;
  
  return 1;
}
 

//Non-member, non-friend functions.

//Reads through a file and passes each line to each Ensemb object
//in an array E of Ensemb objects.  Array E has size N_ensemb.
void ReadFile(string file, Ensemb E[], int N_ensemb)
{
  ifstream fin(file.c_str());
  string line;
  
  while (getline(fin,line))
    for (int i=0; i<N_ensemb; i++) {
      int success=E[i].ReadLine(line);
      if (!success) {
	cout << "\nReadFile:  Ensemb::ReadLine failed. Abort.\n";
	exit(1);
      }
    }

  fin.close();
  
}

//Check if all Ensemb objects in an array E (of size N_ensemb) have
//read all of their data points and covariance matrices.
int AllRead(Ensemb E [], int N_ensemb)
{
  int all_read=1;
  for (int i=0; i<N_ensemb; i++)
    all_read=all_read && E[i].AllDataCovRead();
  
  return all_read;
}

//Take an array E (of size N_ensemb) of Ensemb objects, and put the
//data from all ensembles into aggregate arrays.  Array data_x
//contains the 8 x values for each point, array data_y contains
//the y values for each point, and cov_y a big covariance matrix
//for all points (covariances between ensembles are 0).  The indexing
//for data_x, data_y, and cov_y are in the same order and indices
//run over all points in all ensembles (in the order they were
//entered).
void EnsembsToArrays(Ensemb E [], int N_ensemb, double** data_x, double* data_y, double** cov_y)
{
  //Assumed array sizes
  //Ensemb E[N_ensemb]
  //double data_x[N][8]
  //double data_y[N]
  //double cov_y[N][N]
  //where N is the total number of data points in all ensembles.
  
  //Make an array with the number of data points contained in all
  //ensembles before a given one, and also count the total number
  //of data points N.
  int num_before[N_ensemb], N=0;
  for (int ii=0; ii<N_ensemb; ii++) {
    num_before[ii]=N;
    N+=E[ii].ReturnNdata();
  }

  //i_tot is an index from 0 to N-1, the overall index for the point.
  //ii is the ensemble number the point is in.
  //i is the index of the point within the ensemble.
  //Similarly for j_tot, jj, j.
  //Given i_tot, the other two are obtained via
  //inds(num_before,N_ensemb,i_tot,ii,i)
  //Use this notation to fill data_x, data_y, cov_y with the right
  //values.
  for (int i_tot=0; i_tot<N; i_tot++) {
    int ii, i;
    inds(num_before,N_ensemb,i_tot,ii,i);
    for (int col=0; col<3; col++)
      data_x[i_tot][col]=E[ii].ReturnDataX(i,col);
    for (int col=3; col<8; col++)
      data_x[i_tot][col]=E[ii].ReturnInfo(col-3);
    data_y[i_tot]=E[ii].ReturnDataY(i);

    for (int j_tot=0; j_tot<N; j_tot++) {
      int jj, j;
      inds(num_before,N_ensemb,j_tot,jj,j);
      if (ii!=jj)
	cov_y[i_tot][j_tot]=0.0;
      else
	cov_y[i_tot][j_tot]=E[ii].ReturnCovY(i,j);
    }

  }
  
}

//From an overall index i_tot, determines ensemble number ii and
//index within the ensemble i, (the latter two are outputs).
//Needs array num_before (of size N_ensemb) which gives the number
//of data points before a given ensemble in the overall indexing.
void inds(int num_before [], int N_ensemb, int i_tot, int & ii, int & i)
{
  int ii_tmp;

  for (ii_tmp=0; ii_tmp<N_ensemb-1; ii_tmp++)
    if ( i_tot>=num_before[ii_tmp] && i_tot<num_before[ii_tmp+1] ) {
      ii=ii_tmp;
      i=i_tot-num_before[ii];
      break;
    }
  
  if (ii_tmp==N_ensemb-1) {
    //In the last ensemble.
    ii=ii_tmp;
    i=i_tot-num_before[ii];
    //Since this function doesn't know how big the last ensemble is,
    //it could potentially return an i value out of the range of
    //the last ensemble if the i_tot given to this function is
    //outside the allowed range.  (Might have to worry about this
    //if I start using the inds function in more places.)
  }
    
}
