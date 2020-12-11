#include <string>
using namespace std;

class Grid{
    int n;
public:
    double g;
    double *x;
    double *conc;
	double* concA;
	double* concB;
	double* concY;
    Grid(int nn);
    ~Grid();
    void grid(double dx,double gamma);
    void init_c(double concA, double concB, double concY);
    void grad();
	void updateAll();
	void saveA(std::string filename);
	void saveB(std::string filename);
	void saveY(std::string filename);
	

};