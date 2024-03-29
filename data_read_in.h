#ifndef INCLUDE_DATA_READ_IN
#define INCLUDE_DATA_READ_IN

#include <complex>
#include "jackknife.h"
#include "jackknife_time_series.h"

class WMEDataReadIn
{
 private:
  double **data;
  int N;
  int Nt;
  int **was_read;
  int *col_nums; //Column numbers (starting from 0 which must have certain 
                 //labels for the line to represent data that should be read.
  string *labels; //Labels that should be in the columns in col_nums.
  int n_labels; //Size of the col_nums and labels array.
  int t_col; //Column in which the time slice is found.
  int data_col; //Column in which the data is found.
  string info; //String describing the data being read in (for informational
               //purposes only).
 public:
  explicit WMEDataReadIn(int nt=1, int n=2, int datcol=3, int tcol=2, string in="");
  virtual ~WMEDataReadIn();
  WMEDataReadIn & operator=(const WMEDataReadIn & D);
  int ReturnIsDefined(int t, int i) const;
  int ReturnIsDefined(int t) const;
  int AllDataRead() const;
  int ReturnNt() const;
  int ReturnN() const;
  int ReturnNlabels() const;
  int ReturnColNum(int i) const;
  string ReturnLabel(int i) const;
  int ReturnTCol() const;
  int ReturnDataCol() const;
  string Info() const;
  int IsComplex() const { return 0; }
  WMEDataReadIn & AddColInfo(int col, string lbl);
  int ReadLine(int conf, const string & line, int vrb=0);
  friend JackknifeTimeSeries ReadWMEData(const WMEDataReadIn & D);
};

class WMECmplxDataReadIn
{
 private:
  complex<double> **data;
  int N;
  int Nt;
  int **was_read;
  int *col_nums; //Column numbers (starting from 0 which must have certain 
                 //labels for the line to represent data that should be read.
  string *labels; //Labels that should be in the columns in col_nums.
  int n_labels; //Size of the col_nums and labels array.
  int t_col; //Column which the time slice is found.
  int re_data_col; //Column in which the real part of the data value is found.
  int im_data_col; //Column in which the imaginary part of the data value is
                   //found.
  string info; //String describing the data being read in (for informational
               //purposes only).
 public:
  explicit WMECmplxDataReadIn(int nt=1, int n=2, int redatcol=3, int imdatcol=4, int tcol=2, string in="");
  virtual ~WMECmplxDataReadIn();
  WMECmplxDataReadIn & operator=(const WMECmplxDataReadIn & D);
  int ReturnIsDefined(int t, int i) const;
  int ReturnIsDefined(int t) const;
  int AllDataRead() const;
  int ReturnNt() const;
  int ReturnN() const;
  int ReturnNlabels() const;
  int ReturnColNum(int i) const;
  string ReturnLabel(int i) const;
  int ReturnTCol() const;
  int ReturnReDataCol() const;
  int ReturnImDataCol() const;
  string Info() const;
  int IsComplex() const { return 1; }
  WMECmplxDataReadIn & AddColInfo(int col, string lbl);
  int ReadLine(int conf, const string & line, int vrb=0);
  friend JackknifeCmplxTimeSeries ReadWMECmplxData(const WMECmplxDataReadIn & D);
};

int split(string line, string vals[], int max_vals);
    

#endif
