#include <cmath>
#include <string>

class Geometry {
  public:
 //Adding 4 ghost cells, 2 on each side
    Geometry(int nCells, std::string initialCond);
    ~Geometry();
 
    //getter methods

    double getM() const;
    double getDh() const;
    double* getGrid();
    //returns a pointer that points to a copy of the grid without the extra 2 cells
    double* dispU();
//    double* getTB();
    //setter methods

    private:
    int nCells;
    // Ghosts Ooooooooooo
    int nGhostCells;
    // T.C.C = Total Cell Count
    int TCC;
    // Only one D for now.... muhahahah
    int nDim;
    double deltaH; 
    double* grid;
    double* firstCell;
    double* currentCell;
    double* lastCell;
   
   //Types include the Square Wave, Semicircle, and the Gaussian Pulse
    initializeGrid(string type); 
};   

inline double Geometry::getM() const { return M; }
inline double Geometry::getDh() const { return h; }
inline double* Geometry::getU() { return U; }
//inline double* Geometry::getT() { return T; }

//inline double* Geometry::getTB() { return TB; }
