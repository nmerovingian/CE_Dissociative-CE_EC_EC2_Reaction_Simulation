#include<iostream>
#include "helper.h"

double convergence(double* conc, double* space, int spaceSteps) {
	double massFinal = 0.0;
	double temp = 0.0;
	for (int i = 0; i < spaceSteps - 1; i++) {
		temp = (conc[i] + conc[i + 1]) / 2.0;
		massFinal += (temp * (space[i + 1] - space[i]));
	}
	double massInitial = 0.0;
	for (int i = 0; i < spaceSteps - 1; i++) {
		massInitial += (space[i + 1] - space[i]) * 1.0;
	}
	double massDifference = massInitial - massFinal;
	return massDifference; 
}
double conservation(double CA, double CB, double CY,double CZ,double* concA, double* concB, double* concY,double* concZ, double* grid, int n) {
	double massInitial = 0.0;
	double massA = 0.0;
	double massB = 0.0;
	double massY = 0.0;
	double massZ = 0.0;
	double tempA = 0.0;
	double tempB = 0.0;
	double tempY = 0.0;
	double tempZ = 0.0;
	//Print1DArray(concA, n);
	//Print1DArray(concZ, n);
	for (int i = 0; i < n-1; i++) {
		double h = grid[i + 1] - grid[i];
		massInitial += h * (2*CA + CB + CY + CZ);
		tempA = (concA[i] + concA[i + 1]) / 2.0;
		massA += 2.0*tempA * h;
		tempB = (concB[i] + concB[i + 1]) / 2.0;
		massB += tempB * h;
		tempY = (concY[i] + concY[i + 1]) / 2.0;
		massY += tempY * h;
		tempZ = (concZ[i] + concZ[i + 1]) / 2.0;
		massZ += tempZ * h;
	}
	std::cout << "Initial mass=" << massInitial << "\n";

	double massFinal = massA + massB +  massY+massZ;
	std::cout << "finalMassA=" << massA << "\tfinalMassB=" << massB << "\tfinalMassY=" << massY << "\tFinal MassZ=" << massZ  << "\tmassFinal = " << massFinal << "\n";
	return (massFinal - massInitial) / massInitial * 100;
	
}
