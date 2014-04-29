#include <MatrixIO.h>
#include <SparseGenEigenSolver.h>
#include <SparseMatrixTools.h>
#include <AuxTools.h>
#include <Timer.h>
using namespace std;
using namespace Eigen;
using namespace UTILITY;
using namespace EIGEN3EXT;

#define ERROR_EXIT(msg, succ) { ERROR_LOG_COND(msg,succ); if(!succ) exit(0); };

// format of the matrix
// #Row #Column
// i a x  ------ For a line in the file, 'x' means the element value
// corresponding to the ith row and ath column of the whole matrix.
bool loadMat(SparseMatrix<double> &M, const string file, const bool check_sym=true){

  TRACE_FUN();
  INFO_LOG("loading file: "<< file);
    
  fstream inf;
  inf.open(file.c_str());
  if (!inf.is_open()){
	return false;
  }

  int row = 0, col = 0;
  inf >> row >> col;
  assert_gt(row,0);
  assert_gt(col,0);
  INFO_LOG("rows: "<< row);
  INFO_LOG("cols: "<< col);
  
  /// read the data
  typedef Eigen::Triplet<double> Tri;
  vector<Tri> data;

  Timer timer;
  timer.start();
  while (true) {
	int r,c;
	double v;
    inf >> r >> c >> v;
	data.push_back(Tri(r,c,v));
    if( inf.eof() )
	  break;
  }
  timer.stop("time for load data: ");
  inf.close();

  INFO_LOG("non zeros: "<< data.size());

  timer.start();
  M.resize(row, col);
  M.setFromTriplets(data.begin(), data.end());
  timer.stop("time for init M: ");

  if (check_sym){
	const SparseMatrix<double> Mt = M.transpose();
	assert_le((Mt-M).norm(),1e-12);
  }

  return true;
}

template <typename OS, typename INT>
void hex2vtk(OS &os,const double *node, size_t node_num,const INT *hex, size_t hex_num){

  os << "# vtk DataFile Version 2.0\nhex\nASCII\n\nDATASET UNSTRUCTURED_GRID\n";
  os << "POINTS " << node_num << " double\n";
  for(size_t i = 0; i < node_num; ++i)
	os << node[i*3+0] << " " << node[i*3+1] << " " << node[i*3+2] << "\n";

  os << "CELLS " << hex_num << " " << hex_num*9 << "\n";
  for(size_t i = 0; i < hex_num; ++i){
	os << 8 << "  "
	   << hex[i*8+0] << " " << hex[i*8+1] << " "
	   << hex[i*8+2] << " " << hex[i*8+3] << " "
	   << hex[i*8+4] << " " << hex[i*8+5] << " "
	   << hex[i*8+6] << " " << hex[i*8+7] << "\n";
  }

  os << "CELL_TYPES " << hex_num << "\n";
  for(size_t i = 0; i < hex_num; ++i)
	os << 12 << "\n";
}

void saveMeshes(const MatrixXd &U, const string input_mesh){
    
  TRACE_FUN();
  
  ifstream inf;
  inf.open(input_mesh);
  ERROR_EXIT("failed to open file: "<< input_mesh,inf.is_open());
  
  int num_nodes, num_hex;
  string tempt;
  inf >> tempt >> num_nodes >> num_hex >> tempt;
  assert_gt(num_nodes,0);
  assert_gt(num_hex,0);

  VectorXd nodes(num_nodes*3);
  vector<int> hex(num_hex*8);
  for (int i = 0; i < num_nodes; ++i){
    inf >> nodes[i*3+0];
    inf >> nodes[i*3+1];
    inf >> nodes[i*3+2];
	inf >> tempt;
  }
  for (int i = 0; i < num_hex; ++i){
	inf >> tempt;
    inf >> hex[i*8+0];
	inf >> hex[i*8+1];
	inf >> hex[i*8+2];
	inf >> hex[i*8+3];
	inf >> hex[i*8+4];
	inf >> hex[i*8+5];
	inf >> hex[i*8+6];
	inf >> hex[i*8+7];
	inf >> tempt;
	inf >> tempt;
	inf >> tempt;
  }
  inf.close();

  MatrixXd disp = U;
  if (disp.cols() <= 0){
	disp.resize(num_nodes*3,1);
	disp.setZero();
  }

  for (int i = 0; i < disp.cols(); ++i){

	const string fname = input_mesh+"_mode"+TOSTR(i)+".vtk";
	ofstream outf;
	outf.open(fname);
	ERROR_EXIT("failed to open file: "<< fname,outf.is_open());
	
	assert_eq(nodes.size(), disp.rows());
	// const VectorXd x = nodes+disp.col(i);
	VectorXd x = nodes;
	x[138] += -1941.05;
	x[139] += 1.84181e-05;
	x[140] += -5.51654e-05;

	hex2vtk(outf, &x[0], num_nodes, &hex[0], num_hex);
	outf.close();
  }
}

int main(int argc, char *argv[]){
  
  // check arg
  ERROR_EXIT("ussage: modal_analysis stiff_matrix_file mass_matrix_file number_of_eigs [input_mesh_file] [con_nodes_file]", (4 <= argc));
  const string k_file = argv[1];
  const string m_file = argv[2];
  const int eig_num = atoi(argv[3]);
  ERROR_EXIT("eig number should > 6, while "<<eig_num<<" is given.", eig_num > 6);

  // load data
  INFO_LOG("loading data...");
  SparseMatrix<double> M_all, K_all;
  bool succ = loadMat(M_all, m_file);    ERROR_EXIT("failed to load: "<<m_file, succ);
  succ = loadMat(K_all, k_file);         ERROR_EXIT("failed to load: "<<k_file, succ);

  // extract constrained nodes
  SparseMatrix<double> P = eye<double>(M_all.rows(), M_all.cols());
  if (argc > 5){
  	vector<int> con_nodes;
  	succ = loadVec(argv[5], con_nodes, UTILITY::TEXT);
  	ERROR_LOG_COND("failed to load the constrained nodes: "<<argv[5], succ);
  	set<int> con_set;
  	for (int i = 0; i < con_nodes.size(); ++i){
  	  assert_ge(con_nodes[i],0);
  	  con_set.insert(con_nodes[i]);
  	}
  	genReshapeMatrix(M_all.rows(), 3, con_set, P, false);
  }

  const SparseMatrix<double> M = P*M_all*P.transpose();
  const SparseMatrix<double> K = P*K_all*P.transpose();
  assert_eq(M.rows(), K.rows());
  assert_eq(M.cols(), K.cols());
  cout << "M.dim = " << M.rows() << endl;

  // convert
  INFO_LOG("convert to lower matrix...");
  const SparseMatrix<double> Klower = EIGEN3EXT::getLower(K);
  const SparseMatrix<double> Mlower = EIGEN3EXT::getLower(M);

  // solve
  INFO_LOG("solve the general eigenvalue problem...");
  MatrixXd W_sub;
  VectorXd lambda;
  Timer timer;
  timer.start();
  succ = EigenSparseGenEigenSolver::solve(Klower,Mlower,W_sub,lambda,eig_num);
  timer.stop("time for solving: ");
  const MatrixXd W = P.transpose()*W_sub;
  cout << "W.row = " << W.rows() << endl;
  cout << "W.col = " << W.cols() << endl;

  // save data
  INFO_LOG("saving the data...");
  const string W_f = k_file+".W";
  const string lambda_f = k_file +".lambda";
  ofstream outf;
  outf.open(W_f);
  ERROR_LOG_COND("failed to open file: "<< W_f<< "for saving W.", outf.is_open());
  outf << W.rows()<< "\t" << W.cols();
  for (int r = 0; r < W.rows(); ++r){
  	outf << endl;
  	for (int c = 0; c < W.cols(); ++c){
  	  outf << W(r,c);
  	  if (c != W.cols()-1)
  		outf << " ";
  	}
  }
  outf.close();

  outf.open(lambda_f);
  ERROR_LOG_COND("failed to open file: "<< lambda_f<< "for saving lambda.", outf.is_open());
  outf << lambda.size();
  for (int i = 0; i < lambda.size(); ++i){
    outf << endl << lambda[i];
  }
  outf.close();

  if (argc > 4){
  	saveMeshes(W, string(argv[4]));
  }
  
  return 0;
}
