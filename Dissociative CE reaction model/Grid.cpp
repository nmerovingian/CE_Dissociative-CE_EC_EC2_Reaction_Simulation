#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include<assert.h>
#include "Grid.h"
using namespace std;

Grid::Grid(int nn){
    n = nn;
	std::cout << n << "\n";
    x = new double[n];
    conc = new double[4*n];
	concA = new double[n];
	concB = new double[n];
	concY = new double[n];
	concZ = new double[n];
	g = 0.0;
	gB = 0.0;
}

Grid::~Grid(){
    delete[] x;
    delete[] conc;
}

void Grid::grid(double dx, double gamma){
    x[0] = 0.0;
    for(int i = 1; i < n; i++){
        x[i] = x[i-1] + dx;
		dx *= (1 + gamma);
    }
}

/*void Grid::init_c(double sc){
    for(int i = 0; i < n; i++){
        conc[i] = sc;
    }
}*/

void Grid::init_c(double A, double B, double Y,double Z) {
	for (int i = 0; i < 4 * n; i = i + 4) {
		conc[i] = A;
		conc[i + 1] = B;
		conc[i + 2] = Y;
		conc[i + 3] = Z;
	}
	for (int i = 0; i < n; i++) {
		concA[i] = A;
		concB[i] = B;
		concY[i] = Y;
		concZ[i] = Z;
	}
}

void Grid::grad(){
    g = -(conc[5] - conc[1]) / (x[1] - x[0]);
	//std::cout << conc[5] << "\t" << conc[1] << "\t" << x[1] << "\t" << x[0]<< "\n";
}
void Grid::gradB() {
	gB = (conc[5] - conc[1]) / (x[1] - x[0]);
}
void Grid::updateAll() {
	for (int j = 0, i = 0; j <n*4; j = j + 4, i++) {
		concA[i] = conc[j];
		concB[i] = conc[j + 1];
		concY[i] = conc[j + 2];
		concZ[i] = conc[j + 3];

	}
	/*for (int i = 3; i < n * 3 - 3; i++) {
		assert(conc[i] <= 1.000000001 && conc[i] >= -0.001 && "Concentration of species out of bound");
		if (i >= 6 && i < n - 3) {
			assert(fabs(conc[i] - conc[i - 3]) < 0.1 && "Large variation with previous concentrtaion");
			assert(fabs(conc[i] - conc[i + 3]) < 0.1 && "Large variation from next concentration");
		}
	}*/
}

void Grid::saveA(string filename)
{
    ofstream file(filename.c_str(), ios::out);
    for(int i = 1; i < n; i++)
    {
        file << x[i] << "," << concA[i] << endl;
    }
    file.close();
}
void Grid::saveB(string filename)
{
	ofstream file(filename.c_str(), ios::out);
	for (int i = 1; i < n; i++)
	{
		file << x[i] << "," << concB[i] << endl;
	}
	file.close();
}
void Grid::saveY(string filename)
{
	ofstream file(filename.c_str(), ios::out);
	for (int i = 1; i < n; i++)
	{
		file << x[i] << "," << concY[i] << endl;
	}
	file.close();
}

void Grid::saveZ(string filename)
{
	ofstream file(filename.c_str(), ios::out);
	for (int i = 1; i < n; i++)
	{
		file << x[i] << "," << concZ[i] << endl;
	}
	file.close();
}
