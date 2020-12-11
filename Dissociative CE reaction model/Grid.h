#include <string>
using namespace std;

class Grid{
    int n;
public:
    double g;
	double gB;
    double* x;
    double* conc;
	double* concA;
	double* concB;
	double* concY;
	double* concZ;
    Grid(int nn);
    ~Grid();
    void grid(double dx,double gamma);
    void init_c(double concA, double concB, double concY,double concZ);
    void grad();
	void gradB();
	void updateAll();
	void saveA(std::string filename);
	void saveB(std::string filename);
	void saveY(std::string filename);
	void saveZ(std::string filename);

};