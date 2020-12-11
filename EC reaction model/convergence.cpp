#include<iostream>

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
double conservation(double CA, double CB, double CY,double* concA, double* concB, double* concY, double* grid, int n) {
	double massInitial = 0.0;
	double massA = 0.0;
	double massB = 0.0;
	double massY = 0.0;
	double tempA = 0.0;
	double tempB = 0.0;
	double tempY = 0.0;
	for (int i = 0; i < n-1; i++) {
		double h = grid[i + 1] - grid[i];
		massInitial += h * (CA + CB + CY);
		tempA = (concA[i] + concA[i + 1]) / 2.0;
		massA += tempA * h;
		tempB = (concB[i] + concB[i + 1]) / 2.0;
		massB += tempB * h;
		tempY = (concY[i] + concY[i + 1]) / 2.0;
		massY += tempY * h;
	}
	std::cout << "Initial mass=" << massInitial << "\n";

	double massFinal = massA + massB +  massY;
	std::cout << "finalMassA=" << massA << "\t finalMassB=" << massB << "\t finalMassY=" << massY << "\t" << "massFinal = " << massFinal << "\n";
	return (massFinal - massInitial) / massInitial * 100;
	
}
