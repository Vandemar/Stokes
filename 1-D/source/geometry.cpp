#include "geometry.h"

#define PI 3.14159265

using namespace std;

Geometry::Geometry(double numberOfCells, double dx, string initialCond ) {
  M = numberOfCells;
  h = dx;
  double x;
  int N = (int) M;  
  U = new double[N+2];
  //initializing temperature data to a square wave
  if ( initialCond.compare("Square Wave") == 0) {
    for(int i =1; i <= N/2;  ++i) {
      U[i] = U[N+1-i] = 0;
      if ( (i-1)/M >= 0.25 && (i-1)/M <= 0.75 )
        U[i] = U[N+1-i] = 1;
    }
  }
  else if ( initialCond.compare("Semicircle") == 0 ) {
    double tmp;
    for(int i =1; i <= N/2;  ++i) {
      x = (i-1)/M;
      tmp = (x-0.5)*(x-0.5);
      U[i] = U[N+1-i] = sqrt(0.25 - tmp);
    }
  }
  else if ( initialCond.compare("Gaussian Pulse") == 0 ) {
     double tmp;
     for(int i = 1; i <= N/2; i++) { 
      x = (i-1)/M;
      tmp = -256 * (x-0.5)*(x-0.5);
      U[i] = U[N-i+1] = exp(tmp);
     }
  }
  else {
    U[1] = 1; 
    for( int i = 2; i <= N/2; i++) 
      U[i] = U[N+1-i] = 0;
  } 
  
  //Periodic boundary conditions
  U[0] = U[N];
  U[N+1] = U[1];
}

Geometry::~Geometry() {
  delete[] U;
}

double* Geometry::dispU() {
  double *disp_u = new double[M];
  double *full_u = getU();
  for ( int j=1; j<=M; j++) 
    disp_u[j-1] = full_u[j];
  return disp_u;
} 
