#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include "geometry.h"
//#include "hdf5.h"
#include <cmath>
#include <string>
#include "math.h"

namespace po = boost::program_options;
using namespace std;
using namespace Eigen;

// A helper function to simplify the main part.
template<class T>
ostream& operator<<(ostream& os, const vector<T>& v);
void computeNextTimeStep(Geometry *toGrid, string FDM, double deltaT, double h, double a); 
void stepInTime(SparseMatrix<double>* fd, SparseMatrix<double>* cd, Geometry* grid, double k, double h, double a);
double L1_norm(double *exact_solution, double *numerical_solution, int m);
double L2_norm(double *exact_solution, double *numerical_solution, int m);
double LInfinity_norm(double *exact_solution, double *numerical_solution, int m);
void writeToFile(double *toData, double dim, double timeStp);
void writeToStdout(double *toData, double dim);
//void writeH5(std::ofstream xdm5, string outputDir, Geometry* toGrid, int timeStep);


//global variables
//Opening output stream
std::ofstream *ptrFs;

//This program models the advection equation for various conservative finite difference methods
int main(int ac, char *av[]) {
    int cellCount = 0;
    double h = 0.0;
    double a = 0.0;
    double startT, endT, dt, ts;
    double CFL;
    string icType, bndryType;
    string FDM;
    string config_file;
    string outputDir;
    string fileSuffix;
    string fileName;

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
      ("M", po::value<int>(&cellCount), "number of cells within the domain, must be a power of 2")
      ("startTime", po::value<double>(&startT), "Initial Time")
      ("endTime", po::value<double>(&endT), "End Time")
      ("finiteDifferenceMethod", po::value<string>(&FDM), "Possibilities are Upwind, Lax-Friedrichs, Lax-Wendroff")
      ("initialConditions", po::value<string>(&icType), "Type of initial condition, (Square Wave, Semicircle, Gaussian Pulse")
      ("boundaryConditions", po::value<string>(&bndryType), "Type of Boundary Condition, (Periodic")
      ("advectionConstant", po::value<double>(&a), "speed")
      ("outputDirectory", po::value<string>(&outputDir), "Output Directory")
      ("CFL", po::value<double>(&CFL), "CFL number")
      ("suffixFilename", po::value<string>(&fileSuffix), "Filename for initial data and final data to be written to")
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
  
  //DEfining the variables
  double conRateL1=0, conRateL2=0, conRateLInf=0;
  double lastL1=0,lastL2=0,lastLInf=0;
  bool firstItr = true;
  double *cellSnapshot;
  fileName = outputDir + '/' + FDM + "_" + icType + ".txt";
  ofstream stats(fileName.c_str());
  ostringstream sc;
  int n;
  double m;
  
  //Initializing the Problem
  //For assignment
  for( int k=6; k < 12; k++, cellCount *= 2) {
   //whats going on when one casts an int to a double. Are the terms distributed to the exponent, mantissa appropriately? 
    h = 1/cellCount;
    sc << cellCount;
    string count = sc.str();
    fileName = outputDir + '/' + FDM + "for" + icType + "_h_" + count + '.' + fileSuffix;
    ofstream fs(fileName.c_str());
    ptrFs = &fs;
    
    dt = CFL*h/a;
    ts = (endT-startT)/dt;
    // Plus 2 since we need one extra cell on each far side due to periodicness. 
    //Initializing the grid, and toGrid will then be refrence to the grid. This will allow me to pass the Grid's address location in memory to a function.
    Geometry grid(cellCount, icType); 
    Geometry* toGrid = &grid;
    n = (int) toGrid->getM();
    m = toGrid->getM();
    int i,j;
    
    cout << "Spacing of adjacent cells is: " << h << "\n";
    cout << "Number of Time Steps is set to: " << ts << "\n";
    cout << "Initial conditions have been set to a " << icType << endl;
    cout << "Writing data to " << fileName << endl;
 
    i = 0; 
    if(firstItr)
      writeToFile(toGrid->dispU(),n, i);
    //The main work:
    for( double i=0; i<(double)ts; i++) {
      //cout << "At time step " << i << ":" << endl;
      cellSnapshot = toGrid->dispU();
      computeNextTimeStep(toGrid, FDM, dt, h, a); 
    }

    Geometry grid_tmp(cellCount, h, icType); 
    writeToFile(cellSnapshot, n, (int) (ts-1));
    fs.close();

    double L1, L2, LInf;
    L1 = L1_norm(grid_tmp.dispU(), cellSnapshot, n);
    L2 = L2_norm(grid_tmp.dispU(), cellSnapshot, n);
    LInf = LInfinity_norm(grid_tmp.dispU(), cellSnapshot, n);
  
    if (!firstItr) {
      conRateL1 = log2 (lastL1/L1);
      conRateL2 = log2 (lastL2/L2);
      conRateLInf = log2 (lastLInf/LInf);
    }
    else 
      firstItr = false;

    stats << "Number of cells: " << cellCount  << endl;
    stats << "\t L1 Norm: " << L1 << endl;
    if(!firstItr)
      stats << "\t\t Convergence Rate: " << conRateL1 << endl;
    stats << "\t L2 Norm: " << L2 << endl;
    if(!firstItr)
      stats << "\t\t Convergence Rate: " << conRateL2 << endl;
    stats << "\t LInf Norm: " << LInf << endl << endl;
    if(!firstItr)
      stats << "\t\t Convergence Rate: " << conRateLInf << endl;

    lastL1 = L1;
    lastL2 = L2;
    lastLInf = LInf;
  }   
  
  stats.close();

  return 0;
}


double L1_norm(double *exact_solution, double *numerical_solution, int m ) {
  double sum;
  for(int i=0; i < m; i++) {
    sum += abs(exact_solution[i] - numerical_solution[i]);
  }
  return sum/m;
}   

double L2_norm(double *exact_sol, double *numerical_solution, int m) {
  double sum=0, error;
  double h = 1/(double) m;
  for(int i=0; i < m; i++) {
    error = abs(exact_sol[i] - numerical_solution[i]);
      sum = sum + (h*error*error);
//    sum += pow(error,2); 
  }
  
  //return pow(sum,0.5)/((double) m);
  return std::sqrt(sum);
}
  
double LInfinity_norm(double *exact_solution, double *numerical_solution, int m) {
  double max=0;

  for(int i=0; i < m; i++) {
    if ( abs(numerical_solution[i]-exact_solution[i]) > max ) 
      max = abs(numerical_solution[i]-exact_solution[i]);
  }
  return max;
}

void writeToFile(double *toData, double dim, double timeStp) {
  for(int i = 0; i < dim; i++) {
    //Format: <time step> <space coords> <u (numerical approximation)>
    *ptrFs << timeStp << ',' << i/dim << ',' << toData[i] << endl; 
  }
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
//Add tow ghost cells per side
void computeNextTimeStep(Geometry *toGrid, string FDM, double deltaT, double h, double a) {
  int N = (int) toGrid->getM();
  double *grid = toGrid->getGrid();
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
      uRightHalf = 0.5 * (U[i] + U[i+1]) - (CFL/2) * (U[i+1] - U[i]);
      uLeftHalf = 0.5 * (U[i] + U[i-1]) - (CFL/2) * (U[i] - U[i-1]);
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

