#include "geometry.h"

#define PI 3.14159265

using namespace std;

Geometry::Geometry(double numberOfCells, double dx, string initialCond ) {
  M = numberOfCells;
  h = dx;
  double x;
  int N = (int) M;  
  U = new double[N];
  //initializing temperature data to a square wave
  if ( initialCond.compare("Square Wave") == 0) {
    for(int i =0; i < N/2;  ++i) {
      U[i] = U[N-i-1] = 0;
      if ( i/M >= 0.25 && i/M <= 0.75 )
        U[i] = U[N-1-i] = 1;
    }
  }
  else if ( initialCond.compare("Semicircle") == 0 ) {
    double tmp;
    for(int i =0; i < N/2;  ++i) {
      x = i/M;
      tmp = (x-0.5)*(x-0.5);
      U[i] = U[N-1-i] = sqrt(0.25 - tmp);
    }
  }
  else if ( initialCond.compare("Gaussian Pulse") == 0 ) {
     double tmp;
     for(int i = 0; i < N/2; i++) { 
      x = i/M;
      tmp = -256 * (x-0.5)*(x-0.5);
      U[i] = U[N-i-1] = exp(tmp);
     }
  }
  else {
    U[0] = 1; 
    for( int i = 1; i < N; i++) 
      U[i] = 0;
  } 
}

Geometry::~Geometry() {
  delete[] U;
}


