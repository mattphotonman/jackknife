#ifndef INCLUDE_DATA_READ_IN_HISQ
#define INCLUDE_DATA_READ_IN_HISQ

class Ensemb
{
 private:
  double *info;
  double **data_x;
  double *data_y;
  double **cov_y;
  int N_data;
  int *data_was_read;
  int **cov_was_read;
 public:
  explicit Ensemb(double mL=-1, double mS=-1, double mC=-1, double beta=-1, double L=-1);
  virtual ~Ensemb();
  Ensemb & operator=(const Ensemb & E);
  int AllDataCovRead() const;
  int AllDataRead() const;
  int AllCovRead() const;
  int ReturnNdata() const;
  double ReturnInfo(int col) const;
  double ReturnDataX(int i, int col) const;
  double ReturnDataY(int i) const;
  double ReturnCovY(int i, int j) const;
  int ReturnDataWasRead(int i) const;
  int ReturnCovWasRead(int i, int j) const;
  void AddDataPt(double flag, double mA, double mB);
  int ReadLine(const string & line);
};

void ReadFile(string file, Ensemb E [], int N_ensemb);
int AllRead(Ensemb E [], int N_ensemb);
void EnsembsToArrays(Ensemb E [], int N_ensemb, double** data_x, double* data_y, double** cov_y);
void inds(int num_before [], int N_ensemb, int i_tot, int & ii, int & i);

#endif
