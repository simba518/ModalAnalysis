#include <MatrixIO.h>
#include <SparseGenEigenSolver.h>
#include <SparseMatrixTools.h>
#include <Timer.h>
using namespace std;
using namespace Eigen;
using namespace UTILITY;
using namespace EIGEN3EXT;

#define ERROR_EXIT(msg, succ) { ERROR_LOG_COND(msg,succ); exit(0); };

// format of the matrix
// #Row #Column
// i a x  ------ For a line in the file, 'x' means the element value
// corresponding to the ith row and ath column of the whole matrix.
bool loadMat(SparseMatrix<double> &M, const string file){
    
  fstream inf;
  inf.open(file.c_str());
  if (!inf.is_open()){
	return false;
  }

  int row = 0, col = 0;
  inf >> row >> col;
  assert_gt(row,0);
  assert_gt(col,0);
  
  /// read the data
  M.resize(row, col);
      
  return true;
}

int main(int argc, char *argv[]){

  // check arg
  ERROR_EXIT("ussage: modal_analysis stiff_matrix_file mass_matrix_file number_of_eigs", (4 == argc));
  const string k_file = argv[1];
  const string m_file = argv[2];
  const int eig_num = atoi(argv[3]);
  ERROR_EXIT("eig number should > 6, while "<<eig_num<<" is given.", eig_num > 6);

  // load data
  INFO_LOG("loading data...");
  SparseMatrix<double> M, K;
  bool succ = loadMat(M, m_file);    ERROR_EXIT("failed to load: "<<m_file, succ);
  succ = loadMat(K, k_file);         ERROR_EXIT("failed to load: "<<k_file, succ);

  // convert
  INFO_LOG("convert to lower matrix...");
  const SparseMatrix<double> Klower = EIGEN3EXT::getLower(K);
  const SparseMatrix<double> Mlower = EIGEN3EXT::getLower(M);

  // solve
  INFO_LOG("solve the general eigenvalue problem...");
  MatrixXd W;
  VectorXd lambda;
  Timer timer;
  timer.start();
  succ = EigenSparseGenEigenSolver::solve(Klower,Mlower,W,lambda,12);
  timer.stop("time for solving: ");

  // save data
  const string W_f = k_file+".W";
  const string lambda_f = +".lambda";
  INFO_LOG("saving the data to ");
  succ = write(W_f, W);   ERROR_LOG_COND("failed to save W to "<< W_f, succ);
  succ = write(lambda_f, lambda);   ERROR_LOG_COND("failed to save lambda to "<< lambda_f, succ);
  
  return 0;
}
