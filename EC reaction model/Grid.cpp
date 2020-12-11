#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "Grid.h"
using namespace std;

Grid::Grid(int nn){
    n = nn;
    x = new double[n];
    conc = new double[3*n];
	concA = new double[n];
	concB = new double[n];
	concY = new double[n];
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

void Grid::init_c(double A, double B, double Y) {
	for (int i = 0; i < 3 * n; i = i + 3) {
		conc[i] = A;
		conc[i + 1] = B;
		conc[i + 2] = Y;
	}
	for (int i = 0; i < n; i++) {
		concA[i] = A;
		concB[i] = B;
		concY[i] = Y;
	}
}

void Grid::grad(){
    g = (conc[3] - conc[0]) / (x[1] - x[0]);
}
void Grid::updateAll() {
	for (int j = 0, i = 0; j <n*3; j = j + 3, i++) {
		concA[i] = conc[j];
		concB[i] = conc[j + 1];
		concY[i] = conc[j + 2];
	} 
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
