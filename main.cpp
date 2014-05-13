#include <iomanip>
#include <MatrixIO.h>
#include <SparseGenEigenSolver.h>
#include <SparseMatrixTools.h>
#include <AuxTools.h>
#include <Timer.h>
#include <MassMatrix.h>
#include <ElasticForceTetFullStVK.h>
#include <EigenSolver.h>
using namespace std;
using namespace Eigen;
using namespace UTILITY;
using namespace EIGEN3EXT;

#define ERROR_EXIT(msg, succ) { ERROR_LOG_COND(msg,succ); if(!succ) exit(0); };
enum VOL_MESH_TYPE{TET_MESH, HEX_MESH};

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
	if (!inf.fail())
	  data.push_back(Tri(r,c,v));
	if( inf.eof() || inf.fail())
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

void loadVolMesh(const string vol_mesh, VectorXd &nodes, vector<int> &elem, const VOL_MESH_TYPE type){

  ifstream inf;
  inf.open(vol_mesh);
  ERROR_EXIT("failed to open file:"<<vol_mesh,inf.is_open());

  const int num_node_per_elem = type==TET_MESH?4:8;
  int num_nodes, num_elem;
  string tempt;
  inf >> tempt >> num_nodes >> num_elem >> tempt;
  assert_gt(num_nodes,0);
  assert_gt(num_elem,0);

  nodes.resize(num_nodes*3);
  elem.resize(num_elem*num_node_per_elem);
  int tempt_int;
  for (int i = 0; i < num_nodes; ++i){
	inf >> tempt_int;
	assert_eq(tempt_int,i+1);
    inf >> nodes[i*3+0];
    inf >> nodes[i*3+1];
    inf >> nodes[i*3+2];
  }
  for (int i = 0; i < num_elem; ++i){
	inf >> tempt_int;
	assert_eq(tempt_int,i+1);
	for (int j = 0; j < num_node_per_elem; ++j){
	  inf >> elem[i*num_node_per_elem+j]; 
	  elem[i*num_node_per_elem+j] -= 1;
	}
  }
  inf.close();
}

template <typename OS, typename INT>
void vol2vtk(OS &os,const double *node, size_t node_num,const INT *vol, size_t vol_num, const VOL_MESH_TYPE type){

  os << "# vtk DataFile Version 2.0\nvol\nASCII\n\nDATASET UNSTRUCTURED_GRID\n";
  os << "POINTS " << node_num << " double\n";
  for(size_t i = 0; i < node_num; ++i)
	os << node[i*3+0] << " " << node[i*3+1] << " " << node[i*3+2] << "\n";

  const int num_node_per_elem = type==TET_MESH?4:8;
  os << "CELLS " << vol_num << " " << vol_num*(num_node_per_elem+1) << "\n";
  for(size_t i = 0; i < vol_num; ++i){
	os << num_node_per_elem;
	for (int j = 0; j < num_node_per_elem; ++j)
	  os <<"  "<< vol[i*num_node_per_elem+j];
	os << "\n";
  }

  os << "CELL_TYPES " << vol_num << "\n";
  const int cell_type = type==TET_MESH?10:12;;
  for(size_t i = 0; i < vol_num; ++i)
	os << cell_type << "\n";
}

void saveVolMesh(const MatrixXd &U, const string vol_mesh, const VOL_MESH_TYPE type){
    
  TRACE_FUN();
  
  ifstream inf;
  inf.open(vol_mesh);
  ERROR_EXIT("failed to open file: "<<vol_mesh,inf.is_open());

  VectorXd nodes;
  vector<int> elem;
  loadVolMesh(vol_mesh,nodes,elem,type);
  
  const int num_node_per_elem = type==TET_MESH?4:8;
  MatrixXd disp = U;
  if (disp.cols() <= 0){
	disp.resize(nodes.size(),1);
	disp.setZero();
  }
  for (int i = 0; i < disp.cols(); ++i){

	const string fname = vol_mesh+"_mode"+TOSTR(i)+".vtk";
	ofstream outf;
	outf.open(fname);
	ERROR_EXIT("failed to open file: "<<fname,outf.is_open());
	
	assert_eq(nodes.size(), disp.rows());
	const VectorXd x = nodes+disp.col(i);
	vol2vtk(outf, &x[0], nodes.size()/3, &elem[0], elem.size()/num_node_per_elem, type);
	outf.close();
  }
}

void writeMatrix(const SparseMatrix<double> &M, const string filename){

  ofstream out;
  out.open(filename);
  out << M.rows() <<"\t"<< M.cols() <<endl;
  for (int k=0; k<M.outerSize(); ++k){
	for (SparseMatrix<double>::InnerIterator it(M,k); it; ++it){
	  out << it.row() << " "<< it.col() << " "<<it.value() << endl;
	}
  }
  out.close();
}

void checkMatrix(const SparseMatrix<double> &M, const SparseMatrix<double> &K, 
				 const string tet_file,const double density=1){

  VectorXd nodes;
  vector<int> tet;
  loadVolMesh(tet_file,nodes,tet,TET_MESH);
  pTetMesh tetmesh = pTetMesh(new (TetMesh));
  tetmesh->reset(nodes, tet);
  tetmesh->material().reset(1,1,0.45);
  
  MassMatrix mass;
  SparseMatrix<double> M_correct;
  mass.compute(M_correct, *tetmesh);
  assert_eq(M_correct.nonZeros(), M.nonZeros());
  cout<< "M_correct.norm() = " << M_correct.norm() << endl;
  cout<< "(M-M0).norm() = " << (M - M_correct).norm() << endl;

  ElasticForceTetFullStVK elas(tetmesh);
  elas.prepare();
  VectorXd x0;
  tetmesh->nodes(x0);
  const SparseMatrix<double> K_correct = elas.K(x0)*(-1.0f);
  writeMatrix(K_correct,tet_file+".stvk");

  assert_eq(K_correct.nonZeros(), K.nonZeros());
  cout<< "K_correct.norm() = " << K_correct.norm() << endl;
  cout<< "(K-K0).norm() = " << (K - K_correct).norm() << endl;

}

// solve the general eigenvalue problem using aparck
int generalMA(int argc, char *argv[]){

  // check arg
  ERROR_EXIT("ussage: modal_analysis stiff_matrix_file mass_matrix_file number_of_eigs [LARGEST/SMALLEST] [input_mesh_file] [4/8] [con_nodes_file]", (4 <= argc));
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
  if (argc > 7){
  	vector<int> con_nodes;
  	succ = loadVec(argv[7], con_nodes, UTILITY::TEXT);
  	ERROR_LOG_COND("failed to load the constrained nodes: "<<argv[7], succ);
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
  const char mode = argc>4&&string("LARGEST")==string(argv[4]) ? 'R':'S';
  succ = EigenSparseGenEigenSolver::solve(Klower,Mlower,W_sub,lambda,eig_num,mode,"LM");
  timer.stop("time for solving: ");

  // scale the basis W
  MatrixXd W = P.transpose()*W_sub;
  for (int i = 0; i < W.cols(); ++i)
	W.col(i).normalize();
  cout << "W.row = " << W.rows() << endl;
  cout << "W.col = " << W.cols() << endl;

  // save data
  INFO_LOG("saving the data...");
  const string W_f = k_file+".W"+mode;
  const string lambda_f = k_file +".lambda"+mode;
  ofstream outf;
  outf.open(W_f);
  ERROR_LOG_COND("failed to open file: "<< W_f<< "for saving W.", outf.is_open());
  outf << W.rows()<< "\t" << W.cols();
  for (int r = 0; r < W.rows(); ++r){
  	outf << endl;
  	for (int c = 0; c < W.cols(); ++c){
  	  outf<< setprecision(20) << W(r,c);
  	  if (c != W.cols()-1)
  		outf << " ";
  	}
  }
  outf.close();

  outf.open(lambda_f);
  ERROR_LOG_COND("failed to open file: "<< lambda_f<< "for saving lambda.", outf.is_open());
  outf << lambda.size();
  for (int i = 0; i < lambda.size(); ++i) outf<<setprecision(20)<< endl << lambda[i];
  outf.close();

  if (argc > 5){
	const VOL_MESH_TYPE type = atoi(argv[6])==4 ? TET_MESH:HEX_MESH;
	saveVolMesh(W, string(argv[5]), type); 
  }
  return 0;
}

// compute and save the largest eigenvalue and the corresponding
// eigvector using self-implemented power method.
int largestEigenValue(int argc, char *argv[]){

  // check arg
  ERROR_EXIT("ussage: modal_analysis stiff_matrix_file mass_matrix_file", (3 == argc));
  const string k_file = argv[1];
  const string m_file = argv[2];

  // load data
  INFO_LOG("loading data...");
  SparseMatrix<double> M_all, K_all;
  bool succ = loadMat(M_all, m_file);    ERROR_EXIT("failed to load: "<<m_file, succ);
  succ = loadMat(K_all, k_file);         ERROR_EXIT("failed to load: "<<k_file, succ);

  VectorXd eig_vec;
  double lambda = -1.0;
  INFO_LOG("solving for the largest eigenvalue...");
  Timer timer;
  timer.start();
  const int it = EIGEN3EXT::largestGenEigenSym(K_all,M_all,eig_vec,lambda,K_all.rows(),1e-8);
  timer.stop("time for solving: ");
  if (it <= 0){
	ERROR_LOG("failed to solve the problem.");
	return it;
  }

  // save data
  const string lambda_f = k_file +".largest_lambda";
  INFO_LOG("saving the data to: << " << lambda_f);
  ofstream outf;
  outf.open(lambda_f);
  ERROR_LOG_COND("failed to open file: "<< lambda_f<< "for saving lambda.", outf.is_open());
  outf << lambda;
  outf.close();
  
  return it;
}

int main(int argc, char *argv[]){
  
  return generalMA(argc, argv);
  // return largestEigenValue(argc, argv);
}
