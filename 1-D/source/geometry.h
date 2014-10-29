#include <cmath>
#include <string>

class Geometry {
  public:
    Geometry(double numberOfCells, double h, std::string initialCond);
    ~Geometry();
  
    //getter methods
    //Returns number of cells (1-D)
    double getM() const;
    //Returns delta h
    double getDh() const;
    //returns array U
    double* getU();
//    double* getTB();
    //setter methods

    private:
    double M;
    double h; 
    double* U;
};   

inline double Geometry::getM() const { return M; }
inline double Geometry::getDh() const { return h; }
inline double* Geometry::getU() { return U; }
//inline double* Geometry::getT() { return T; }

//inline double* Geometry::getTB() { return TB; }
