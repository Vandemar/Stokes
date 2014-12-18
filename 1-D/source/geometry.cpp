#include "geometry.h"

#define PI 3.14159265

using namespace std;

Geometry::Geometry(int nCells, string initialCond ) {
  M = nCells;
  h = 1/nCells;
  nGhostCells=4;
  nDim=1;
  //Adding 4 ghost cells, 2 on each side
  TCC = nCells+nGhostCells;
  grid = (double*) new double[(TCC)];
  firstCell = grid[nGhostCells/2];
  //
  // 2nd to last cell, one shifted over;
  lastCell = grid[nCells-(nGhostCells/2)-1]; 
  initializeGrid(initialCond); 
  
  //Periodic boundary conditions
  
  grid[0] = U[N];
  grid[N+1] = grid[1];
}

Geometry::~Geometry() {
  delete[] grid;
}

void Geometry::initializeGrid(string initialCond) {
  int i = 0; 
  currentCell = firstCell; 
  if ( initialCond.compare("Square Wave") == 0) {
    int lCut = nCells/4;
    int rCut = nCells*3/4;
    for(currentCell = firstCell, i=0; currentCell != (lastCell+1); currentCell++, i++){
      *currentCell = 0;
      if ( (i-1)/M >= 0.25 && (i-1)/M <= 0.75 )
        grid[i] = grid[N+1-i] = 1;
    }
  }
  else if ( initialCond.compare("Semicircle") == 0 ) {
    double tmp;
    for(int i =1; i <= N/2;  ++i) {
      x = (i-1)/M;
      tmp = (x-0.5)*(x-0.5);
      grid[i] = grid[N+1-i] = sqrt(0.25 - tmp);
    }
  }
  else if ( initialCond.compare("Gaussian Pulse") == 0 ) {
     double tmp;
     for(int i = 1; i <= N/2; i++) { 
      x = (i-1)/M;
      tmp = -256 * (x-0.5)*(x-0.5);
      grid[i] = grid[N-i+1] = exp(tmp);
     }
  }
  else {
   for(int i=0; i < N; i++) {
     x = i/M;
     grid[i+1]=cos(2*PI*x);
   }
  } 
}
double* Geometry::dispU() {
//  double* disp_u = new double[(int) M];
  double *full_u = this->getU();
  return full_u+1;
/*  for ( int j=1; j<=M; j++) 
    disp_u[j-1] = full_u[j];
  return disp_u;*/
}
