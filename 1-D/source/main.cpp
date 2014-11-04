#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include "geometry.h"
//#include "hdf5.h"
#include <cmath>
#include <string>

namespace po = boost::program_options;
using namespace std;
using namespace Eigen;

// A helper function to simplify the main part.
template<class T>
ostream& operator<<(ostream& os, const vector<T>& v);
void computeNextTimeStep(Geometry *toGrid, string FDM, double deltaT, double h, double a); 
void stepInTime(SparseMatrix<double>* fd, SparseMatrix<double>* cd, Geometry* grid, double k, double h, double a);
double L1_norm(Geometry *toGrid, double *numerical_solution);
double L1_norm(Geometry *toGrid, double *numerical_solution);
double LInfinity_norm(Geometry *toGrid, double *numerical_solution);
void writeToFile(double *toData, double dim, ofstream &outputFile);
void writeToStdout(double *toData, double dim);
//void writeH5(std::ofstream xdm5, string outputDir, Geometry* toGrid, int timeStep);

//This program models the advection equation for various conservative finite difference methods
int main(int ac, char *av[]) {
    double cellCount = 0;
    double h = 0.0;
    double a = 0.0;
    double startT, endT, dt, ts;
    double CFL;
    string icType, bndryType;
    string FDM;
    string config_file;
    string outputDir;
    string data_fileName;

  try {
    //Set of options for Command Line
    po::options_description generic("Generic Options");
    generic.add_options() 
      ("help", "produce a help message")
      ("pF", po::value<string>(&config_file)->default_value("exampleParameters.cfg"), "name of a Parameter file")
    ;
    //Set of options for Config file
    po::options_description tokens("Parameters for stokes solver");
    tokens.add_options() 
      ("M", po::value<double>(&cellCount), "number of cells within the domain, must be a power of 2")
      ("startTime", po::value<double>(&startT), "Initial Time")
      ("endTime", po::value<double>(&endT), "End Time")
      ("finiteDifferenceMethod", po::value<string>(&FDM), "Possibilities are Upwind, Lax-Friedrichs, Lax-Wendroff")
      ("initialConditions", po::value<string>(&icType), "Type of initial condition, (Square Wave, Semicircle, Gaussian Pulse")
      ("boundaryConditions", po::value<string>(&bndryType), "Type of Boundary Condition, (Periodic")
      ("advectionConstant", po::value<double>(&a), "speed")
      ("outputDirectory", po::value<string>(&outputDir), "Output Directory")
      ("CFL", po::value<double>(&CFL), "CFL number")
      ("dataFilename", po::value<string>(&data_fileName), "Filename for initial data and final data to be written to")
    ;
    
    po::options_description cmdline_options;
    cmdline_options.add(generic).add(tokens);
    
    po::options_description config_file_options;
    config_file_options.add(tokens);
    
    po::positional_options_description pCF;
    pCF.add("pF", 1);
    
    po::variables_map vm;
    po::store(po::command_line_parser(ac, av).
              options(cmdline_options).positional(pCF).run(), vm);
    po::notify(vm);

    ifstream ifs(config_file.c_str());
    if(!ifs) {
      cout << "Could not open parameter file: " << config_file << "\n";
      return 0;
    }
    else
    {
      po::store(po::parse_config_file(ifs, config_file_options), vm);
      notify(vm);
    }

    if(vm.count("help")) {
      cout << cmdline_options << "\n";
      return 0;
    }
    else if(vm.count("pF")) {
      cout << "Parameter File has been set to: " << config_file << "\n";
    }
  }
  catch(exception& e) {
    cerr << "error: " << e.what() << "\n";
    return 1;
  }	
  catch(...) {
    cerr << "exception of unknown type";
    return 1;
  }
  
  //Opening output stream
  std::ofstream fs;
  //fs.open((outputDir + data_fileName).c_str(), ios::out);
  fs.open((outputDir + data_fileName).c_str());

  
  //Initializing the Problem
  h = 1/cellCount;
  dt = CFL*h/a;
  ts = (endT-startT)/dt;
  // Plus 2 since we need one extra cell on each far side due to periodicness. 
  Geometry grid(cellCount, h, icType); 
  Geometry* toGrid = &grid;
  int n = (int) toGrid->getM();
  double m = toGrid->getM();
  int i,j;
  cout << "Spacing of adjacent cells is: " << h << "\n";
  cout << "Number of Time Steps is set to: " << ts << "\n";
  double* cellSnapshot = toGrid->dispU();
  writeToStdout(cellSnapshot, m );
  cout << "At time step 0:" << endl;
  cout << "Initial conditions have been set to a " << icType << endl;
  cout << "Writing data to " << data_fileName << endl;
  for( int i=0; i<ts; i++) {
    computeNextTimeStep(toGrid, FDM, dt, h, a); 
    cellSnapshot = toGrid->dispU();
    //writeToStdout(cellSnapshot, m );
    writeToFile(cellSnapshot,m,fs);
  }

  cellSnapshot = toGrid->dispU();
  writeToStdout(cellSnapshot, m );
   
  fs.close();
  return 0;
}


double L1_norm(Geometry *toGrid, double *numerical_solution ) {
   double sum;
  double N = toGrid->getM();
  double* numerical_sol = toGrid->dispU();

  for(int i=0; i < N; i++) {
    sum += abs(numerical_sol[i]);
  }
  return sum;
}   

double L2_norm(Geometry *toGrid, double *numerical_solution) {
  double sum;
  double N = toGrid->getM();
  double* numerical_sol = toGrid->dispU();

  for(int i=0; i < N; i++) {
    sum += numerical_sol[i]*numerical_sol[i]; 
  }
  return sum;
}
  
double LInfinity_norm(Geometry *toGrid, double *numerical_solution) {
   double max=0;
  double N = toGrid->getM();
  double* numerical_sol = toGrid->dispU();

  for(int i=0; i < N; i++) {
    if ( abs(numerical_sol[i]) > max ) 
      max = numerical_sol[i];
  }
  return max;
}

void writeToFile(double *toData, double dim, ofstream &outputFile) {
  for(int i = 0; i < dim; i++) {
    outputFile << "x:" << i/dim << "Value:" << toData[i] << "| ";
  }
  outputFile << endl; 
}
void writeToStdout(double *toData, double dim){
  for(int i = 0; i < dim; i++) {
    cout << "x:" << i/dim << "Value:" << toData[i] << "| ";
  }
}
// ----------------------------
//First Order Difference
/*
  SparseMatrix<double,ColMajor> fda(grid.getM()+1, grid.getM()+1);  
//Second Order Difference
  SparseMatrix<double,ColMajor> identity(grid.getM()+1, grid.getM()+1);  
  //vector<Triplet<double> > fdTmp, cdTmp;
  //fdTmp.reserve(2);
  //cdTmp.reserve(2);
  //fdTmp.push_back(Triplet<double>(1,3,2));
  //cdTmp.push_back(Triplet<double>(1,3,2));
  //fd1.setFromTriplets(fdTmp.begin(),fdTmp.end());
  //cd1.setFromTriplets(cdTmp.begin(), cdTmp.end()); 
  
  fda.reserve(VectorXi::Constant(n+1,2));
  identity.reserve(VectorXi::Constant(n+1,2));
  
  for(j = 0; j < n+1; j++) {
    for(i= 0; i < n+1; i++) {  
      if(i==j) { 
        fda.insert(i,j) = -1;
        identity.insert(i,j) = 1;
      }
      else if((i-1)==j) {
        fda.insert(i,j) = 1;
      }
    }
  }
  
  fda.makeCompressed();
  identity.makeCompressed();
 //Initializing output
  //std::ofstream xdmfFile;
  DataWindow<double> uWindow (grid.getU(), grid.getM(), 1);
  cout << "original" << endl;
  cout << tWindow.displayMatrix() << endl; 

  cout << "Twindow" << endl;
  for( i=0; i<ts; i++) { 
    stepInTime(&fda, &identity, toGrid, dt, h, a_c);
    DataWindow<double> tWindow (grid.getT(), grid.getM(), 1);
    cout << tWindow.displayMatrix() << endl << endl; 
   // writeH5(xdmfFile, outputDir, toGrid, i);
  }
*/

template<class T>
ostream& operator<<(ostream& os, const vector<T>& v)
{
    copy(v.begin(), v.end(), ostream_iterator<T>(os, " "));
    return os;
}

void computeNextTimeStep(Geometry *toGrid, string FDM, double deltaT, double h, double a) {
  int N = (int) toGrid->getM();
  double *U = toGrid->getU();
  double *UNext = new double[N+2];
  double CFL = (deltaT*a)/h;

  for(int i=0; i < N+2; i++) 
    UNext[i] = 0;

  if ( FDM.compare("Upwind") == 0 ) {
    for (int i = 1; i < N+1; i++) 
      UNext[i] = U[i] + CFL * (U[i-1] - U[i]);  
  }
  else if ( FDM.compare("Lax-Friedrichs") == 0 ) {
    for (int i= 1; i < N+1; i++)
      UNext[i] = 0.5 * (U[i-1] + U[i+1]) + (CFL/2) * (U[i-1] - U[i+1]);
  }
  else if ( FDM.compare("Lax-Wendroff") == 0 ) {
    double uLeftHalf, uRightHalf;
    for (int i = 1; i < N+1; i++) {
      uRightHalf = 0.5 * (U[i] + U[i+1]) + (CFL/2) * (U[i] - U[i+1]);
      uLeftHalf = 0.5 * (U[i] + U[i-1]) + (CFL/2) * (U[i] - U[i-1]);
      UNext[i] = U[i] + CFL * (uLeftHalf - uRightHalf); 
    }
  }
  else
    cout << "Specified " << FDM <<" has not been implemented yet" << endl;
 
  for(int i=1; i < N+1; i++) 
    U[i] = UNext[i];

  U[0] = UNext[N];
  U[N+1] = UNext[1];
} 

void stepInTime(SparseMatrix<double>* fda, SparseMatrix<double>* identity, Geometry* grid, double k, double h, double a) {
  //cout << *lfda << endl;
  //cout << *lfds << endl;
  int n = (int) grid->getM();
  Eigen::Matrix<double, Dynamic, 1> u_j, u_j_1;
  u_j.resize(n+1,1);
  u_j_1.resize(n+1,1);
  for(int i = 0; i < n; i++) {
    u_j(i) = *(grid->getU()+i);
  }
  u_j(n) = 0;
  //Forward Difference Approximation Method
  double sigma = (k*a)/h;
  //u_j_1 = 1/2 * (*lfda) * u_j - c * (*lfds) * u_j;
  //u_j_1 = (*identity) * u_j + sigma * (*fda) * u_j;
  u_j_1 = (*identity) * u_j + (*fda) * u_j;
  if(u_j_1(n) != 0) {
    u_j_1(0) = u_j_1(n);
  }
  //Updating Temperature data 
  double* T = grid->getU();
  for(int i = 0; i < n; i++) {
    *(T+i) = u_j_1(i);
  }     
}

/*void writeH5(ofstream xdmfFile, string outputDir, Geometry* toGrid, int timeStep) {
  hid_t file_id;
  herr_t status;
  string fileName = "timestep";
  stringstream stmp;
  stmp << fileName << "-" << timeStep;  
  file_id = H5Fcreate((outputDir + "/" + stmp.str()  + ".h5").c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  cout << "<Outputting current data to \"" << outputDir << "/" << fileName  << "-" << timestep << ".h5\">" << endl;
  status = H5Fclose(file_id);
         
}*/

//void stepInTime(double* td, double n, double ts) {
/*void stepInTime(Geometry* grid) {
  int n = grid->getM();
  //MatrixXd A = MatrixXd::Constant(n,n,0);
  SparseMatrix<double, ColMajor> A(n,n);

  A.setFromTriplets(T.begin(),T.end());
  SparseLU<SparseMatrix<double, ColMajor> > solver;
  solver.analyzePattern(A);
  solver.factorize(A);

    
  cout << A << endl;
  return;
}*/

