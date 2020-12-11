#include <iostream>
#include <cmath>
#include <fstream>
#include <windows.h>
#include <thread>
#include <vector>
#include <thread>
#include <iterator>
#include <mutex>
#include <assert.h>
#include "Grid.h"
#include "Trid.h"
#include "Coeff.h"
#include "convergence.h"
#include "ThomasXL.h"
#include "helper.h"

using namespace std;
std::mutex mlock;
void spawnThread(int n1, double variable1[],int n2, double variable2[], int n3, double variable3[]);
void foo(double variable1, double variable2, double variable3);
double DA = 1e-9;
double E0f = 0.01;
double dElectrode = 0.005641895835478;
void foo(double variable1, double variable2,double variable3) {
	//variable 1 is scan rate, Variable 3 is K0. Now change Varible 2 
	//Input parameters
	double concA = 1.0;
	double concY = 0.0;
	double cAstar = 1e-3;
	double cAstarSI = cAstar * 1000;


	//double theta_i = 20.0;
	double theta_i = 40.0;
	double theta_v = -40.0;
	double sigma = variable1;


	double deltaX = 1e-6;
	double deltaTheta = 1e-2;
	double K0 = 1e-1;
	//K0 /= cAstar;
	//double dimKeq = 10; //unit is M^(-1)
	//double Kf = 0.0;
	//double dimKf = Kf * DA / (cAstar * dElectrode * dElectrode);
	double dimKeq = variable2;
	double dimKf = variable3;
	double dimKb =  dimKf / dimKeq;


	double Kf = dimKf * dElectrode * dElectrode / DA;


	double Kb = dimKb * dElectrode * dElectrode / DA;
	std::cout << "Kf is " << Kf << "Kb is " << Kb << "\n";

	double concB = dimKeq;


	//	std::cout << "dimKf is " << dimKf << " dimKb is " << dimKb << " Kf is " << Kf << " Kb is " << Kb << "\n";

	double alpha = 0.5;
	double gamma = 0.05 ;
	//if (sigma > 38943) {
	//	gamma = 0.02;
	//}

	double dB = 1.0; 
	double dY = 1.0;
	double epsilon = 1e-2;
	int numOfIteration = 2 ;
	int maxNumOfIteration = 8;

	std::vector<double> CV_theta = {};
	std::vector<double> CV_flux = {};
	double calculated_alpha = 0.0;
	
	bool outputConcProfile = false;
	bool Print = false; // if true, will not print concentration profiles. 
	bool printA = true;  //if print forward half peak max profile.
	bool printB = true;	 //if print forward  peak max profile.
	bool printC = true;	 //if print theta_v peak max profile.	
	bool printD = true;  //if print backward half peak max profile
	bool printE = true;  // if print backward peak max profile
	if (!Print) {
		printA = false;
		printB = false;
		printC = false;
		printD = false;
		printE = false;
	}
	//the point where needs to be printed
	double pointA = 0.0;
	double pointB = -50;
	double pointC = theta_v;
	double pointD = 0.0;
	double pointE = 30;

	//Calcute the maximums of time and position 
	double deltaT = deltaTheta / sigma;
	double maxT = 2.0 * fabs(theta_v - theta_i) / sigma;
	double maxX = 6.0 * sqrt(maxT); 
	Coeff coeff(maxX, K0,Kf, Kb, deltaT, alpha,cAstar, gamma, dB, dY);
	mlock.lock();
	coeff.calc_n(deltaX);
	mlock.unlock();
	int m = (int)(maxT / deltaT);
	coeff.ini_jacob();
	coeff.ini_fx();
	coeff.ini_dx();
	Grid grid(coeff.n); 
	grid.grid(deltaX, gamma);
	grid.init_c(concA, concB, concY);
	coeff.get_XX(grid.x);
	coeff.update(grid.conc, theta_i, deltaX, concA, concB, concY);
	coeff.Allcalc_abc(deltaT, theta_i, deltaX);
	coeff.calc_jacob(grid.conc, theta_i);
	coeff.calc_fx(grid.conc, theta_i);
	ThomasXL Trid(coeff.n * 3,7);
	Trid.solve(coeff.J, coeff.fx, coeff.dx);
	coeff.xupdate(grid.conc);
	for (int i = 0; i < 10; i++) {

		coeff.calc_jacob(grid.conc, theta_i);
		coeff.calc_fx(grid.conc, theta_i);
		Trid.solve(coeff.J, coeff.fx, coeff.dx);
		coeff.xupdate(grid.conc);
	}
	//double logVar2 = log10(variable2);

	std::string BKineticsLocation = "Concentration B C var1=" + to_string(variable1) + "var2=" + to_string(variable2) + "var3=" + to_string(variable3);
	//std::string CVLocation = "Variable =" + to_sci(variable2) + " Sigma="; //+"K1=" + to_sci(K1) + "Sigma=" + to_sci(sigma) + "K0=" + to_sci(K0) + "dB=" + to_sci(dB) + "dY=" + to_sci(dY);
	std::string CVLocation = "var1=" + to_string(variable1) + "var2=" + to_string(dimKf) + "var3=" + to_string(dimKeq); //+"K1=" + to_sci(K1) + "Sigma=" + to_sci(sigma) + "K0=" + to_sci(K0) + "dB=" + to_sci(dB) + "dY=" + to_sci(dY);

	std::string ConcALocation = "concA=" + to_sci(concA) + "var1=" + to_string(variable1) + "var2=" + to_string(dimKf) + "var3=" + to_string(dimKeq);
	std::string ConcBLocation = "concB=" + to_sci(concB) + "var1=" + to_string(variable1) + "var2=" + to_string(dimKf) + "var3=" + to_string(dimKeq);
	std::string ConcYLocation = "concY=" + to_sci(concY) + "var1=" + to_string(variable1) + "var2=" + to_string(dimKf) + "var3=" + to_string(dimKeq);


	//Calculate the first flux and write it in output file
	grid.grad();
	grid.gradB();

	//std::cout << "first flux is " << grid.g << "\n";
	std::ofstream myFile;
	myFile.open(genAddress(CVLocation));
	myFile << theta_i << "," << grid.g << "\n"; //Dimensionless Form 
	//myFile << (theta_i / (96385 / (8.314 * 298)) + E0f) << "," << (grid.g * 3.1415926 * dElectrode * 96485 * DA * cAstar * 1000.0) << "\n"; //Dimensional Form
	CV_theta.push_back(theta_i);
	CV_flux.push_back(grid.g);
	mlock.lock();
	cout << "Time step is " << m << "\n";
	mlock.unlock();
	std::ofstream bFile;
	bFile.open(genAddress(BKineticsLocation));
	bFile << theta_i << "," << concB <<"," <<concY << "\n";


	double fluxTotal = 0.0;
	int startTime = static_cast<int>(GetTickCount64());
	//Simulate through the given time
	double Theta = theta_i;
	double max_flux = 0.0;
	double min_flux = 0.0;
	double forwardPeak = 0.0;
	double backwardPeak = 0.0;
	std::thread::id this_id = std::this_thread::get_id();
	for (int i = 0; i < m; i++) {
		if (i < m / 2) {
			Theta -= deltaTheta;
		}
		else {
			Theta += deltaTheta;
		}
		if (i % 10000 == 0) {
			std::cout << "Thread id " << this_id<< "Time Step" << i << '\n';
		}
		if (i == int(m * .001)) {
			estimateRunTime(startTime);
		}
		/*coeff.Allcalc_abc(deltaT, Theta, deltaX);
		coeff.update(grid.conc, Theta, deltaX);
		Trid.solve(coeff.A, coeff.d, grid.conc);
		grid.grad();*/
		coeff.update(grid.conc, Theta, deltaX, concA,concB, concY); //get new d from grid.conc;
		coeff.Allcalc_abc(deltaT, Theta, deltaX);
		for (int ii = 0; ii < numOfIteration; ii++) {

			coeff.calc_jacob(grid.conc, Theta);
			coeff.calc_fx(grid.conc, Theta);
			Trid.solve(coeff.J, coeff.fx, coeff.dx);
			coeff.xupdate(grid.conc);
		}
		//std::cout << "max dx is " << "\n";
		//std::cout << coeff.max_dx() << "\n";
		if (numOfIteration < maxNumOfIteration && coeff.max_dx() > 1e-11) {
			
			if (numOfIteration < maxNumOfIteration) {
				numOfIteration++;
				std::cout << "/n/nNot converging!/n/n" << Kf << "\t" << "avg_dx()is" << coeff.avg_dx() << "\n" << "max_dX is" << coeff.max_dx() << "\n";
				std::cout << "Num of Iteration Now is: " << numOfIteration << "\n";
			}
			
		}
		grid.updateAll();
		grid.grad();
		grid.gradB();
		//std::cout << "\n\n\n\n\n\n" << grid.g << "\n\n\n";
		//std::cout << "grid" << grid.g << "\t";
		fluxTotal += grid.g * deltaT;
 		myFile << Theta << "," << grid.g << "\n";  //Dimensionless Form
		//myFile << (theta_i / (96385 / (8.314 * 298)) + E0f) << "," << (grid.g * 3.1415926 * dElectrode * 96485 * DA * cAstar * 1000.0) << "\n"; //Dimensional Form
		CV_theta.push_back(Theta);
		CV_flux.push_back(grid.g);
		bFile << Theta << "," << grid.concB[0] <<"," <<grid.concY[0] << "\n";
		// find the maximum flux
		if (grid.g > max_flux) {
			max_flux = grid.g;
			forwardPeak = Theta;

		}
		//find the minimum flux
		if (grid.g < min_flux) {
			min_flux = grid.g;
			backwardPeak = Theta;
		}
		if (printA && (fabs(Theta - pointA) < epsilon && i < m / 2)) {
			std::cout << "Saving PointA at " << Theta << "\n";
			//Point A Forward scan Half Max
			std::string  str = "Point=A,Theta=" + to_string(Theta);
			grid.saveA(genAddress( str + ConcALocation));
			grid.saveB(genAddress( str + ConcBLocation));
			grid.saveY(genAddress(str +  ConcYLocation));
			printA = false;
		}
		if (printB && fabs(Theta - pointB) < epsilon && i < m / 2) {
			std::cout << "Saving PointB at " << Theta << "\n";
			std::string  str = "Point=B,Theta=" + to_string(Theta);
			grid.saveA(genAddress(str + ConcALocation));
			grid.saveB(genAddress(str + ConcBLocation));
			grid.saveY(genAddress(str + ConcYLocation));
			printB = false;
		}
		if (printC && fabs(Theta - pointC) < epsilon) {
			std::string  str = "Point=C,Theta=" + to_string(Theta);
			std::cout << "Saving PointC at " << Theta << "\n";
			grid.saveA(genAddress(str + ConcALocation));
			grid.saveB(genAddress(str + ConcBLocation));
			grid.saveY(genAddress(str + ConcYLocation));
			printC = false;
		}
		if (printD && (fabs(Theta - pointD) < epsilon && i > m / 2)) {
			std::string  str = "Point=D,Theta=" + to_string(Theta);
			std::cout << "Saving PointD at " << Theta << "\n";
			grid.saveA(genAddress(str + ConcALocation));
			grid.saveB(genAddress(str + ConcBLocation));
			grid.saveY(genAddress(str + ConcYLocation));
			printD = false;
		}
		if (printE && (fabs(Theta - pointE) < epsilon && i > m / 2)) {
			std::string  str = "Point=E,Theta=" + to_string(Theta);
			std::cout << "Saving PointE at " << Theta << "\n";
			grid.saveA(genAddress(str + ConcALocation));
			grid.saveB(genAddress(str + ConcBLocation));
			grid.saveY(genAddress(str + ConcYLocation));
			printE = false;
		}
		//Debug 
		//std::cout.precision(3);
		//std::cout << Theta << " " << grid.conc[0] << " " << grid.conc[1] << " " << grid.conc[2] << " " << grid.conc[3] << " " << grid.conc[4] << " " << grid.conc[5] << " " << grid.conc[0] + grid.conc[1] * 2 + grid.conc[2] * 3 << " " << grid.conc[3] + 2 * grid.conc[4] + grid.conc[5] * 3 << "\n";  //<< grid.conc[6] <<" " << grid.conc[7] << " " << grid.conc[8] << " " << grid.conc[9] << " " << grid.conc[10] << " " << grid.conc[11] << " " << grid.conc[12] << " " << grid.conc[13] << " " << grid.conc[14] << " " << grid.conc[15] << " " << grid.conc[16] << " " << grid.conc[17] << " " << grid.conc[18] << " " << grid.conc[19] << " " << grid.conc[20] << " Gradient: " << -grid.g <<"\n";
		//assert(grid.conc[0] < 1.000000001 && "Error of mass conservation");
		//assert(grid.conc[3] < 1.000000001 && "Error of mass conservation");
		//assert(grid.conc[6] < 1.000000001 && "Error of mass conservation");

		//std::cout << grid.g << " "  <<grid.gB * (2.0) << "  " <<grid.conc[0] << " " <<  grid.conc[1] << "\n";

	}
	//std::cout.precision(9);
	double predIrrFlux = 0.496 * sqrt(alpha) * sqrt(sigma);
	std::cout<<"Kf= " << Kf << "\t" << "Kb= " << Kb << "\n";
	std::cout << "dB= " << dB << "\t" << "dY= " << dY << "\n";
	double diffFromPredIrrFlux = (max_flux - predIrrFlux) / predIrrFlux * 100;
	std::cout << "Difference from Predicted Irr Flux is " << diffFromPredIrrFlux << "% \n";
	double predRevFlux = 0.446 * sqrt(sigma);
	double diffFromPredRevFlux = (max_flux - predRevFlux) / predRevFlux * 100;
	std::cout << "Difference from Predicted Rev Flux is " << diffFromPredRevFlux << "% \n";
	std::cout << "Min Flux is " << min_flux << "\n" << "Forward Scan Peak at " << forwardPeak << "\n";
	std::cout << "Max Flux is " << max_flux << "\n" << "Backward Scan Peak at " << backwardPeak << "\n";
	std::cout << "Peak Separation is " << fabs(forwardPeak - backwardPeak) << " \n";
	std::cout << "flux_integrate is: " << fluxTotal << "\n";
	double massDifference = convergence(grid.concA, grid.x, coeff.n);
	std::cout << "mass difference is " << massDifference << "\n";
	std::cout << "percentage" << (1.0 + fluxTotal / massDifference) * 100.0 << "%" << "\n";
	double massConservation = conservation(concA, concB, concY, grid.concA, grid.concB, grid.concY, grid.x, coeff.n);
	std::cout << "mass Conservation is: " << massConservation << "%\n";

	

	calculated_alpha = cal_alpha(CV_theta, CV_flux);
	std::cout << "calculated alpha is: " << calculated_alpha << "\t" << "alpha is: " << alpha << "\n";
	std::cout << "\n\n\n";

	myFile.close();
	bFile.close();


	if (outputConcProfile) {
		std::string  str = "Point=F,Theta=" + to_string(Theta);
		grid.saveA(genAddress(str + ConcALocation));
		grid.saveB(genAddress(str + ConcBLocation));
		grid.saveY(genAddress(str + ConcYLocation));
	}

	// Write log file in here 
	mlock.lock();
	std::ofstream logFile("C:/Users/nmero/OneDrive - Nexus365/Log/log.txt", std::ios::app); //append to logfile
	logFile << "\n";
	logFile << "Thread ID=" << this_id << "\n";
	logFile << "Sigma = " << to_sci(sigma) << "\t" << "Kf = " << to_sci(Kf) << "\t" << "Kb = " << to_sci(Kb) << "\n";
	logFile << "dB= " << dB << "\t" << "dY= " << dY << "\n";
	logFile << "MaxX = " << maxX << "\t" <<"MaxT = " << maxT << "\n";
	logFile << "Min Flux is " << min_flux << "\t" << "Forward Scan Peak at " << forwardPeak << "\n";
	logFile << "Max Flux is " << max_flux << "\t" << "Backward Scan Peak at " << backwardPeak << "\n";
	logFile << "Difference from Predicted Irr Flux is " << diffFromPredIrrFlux << "% \n";
	logFile << "Difference from Predicted Rev Flux is " << diffFromPredRevFlux << "% \n";
	logFile << "Peak Separation is " << fabs(forwardPeak-backwardPeak) << " \n";
	logFile << "flux_integrate is: " << fluxTotal << "\n";
	logFile << "mass Conservation is: " << massConservation << "%\n";
	logFile << "calculated alpha is: " << calculated_alpha << "\t" << "alpha is: " << alpha << "\n";
	logFile << "deltaX is: " << deltaX << "\n";
	logFile << "deltaTheta is: " << deltaTheta << "\n";
	logFile << "gamma is: " << gamma << "\n";
	logFile << "\n";
	logFile.close();
	std::ofstream alphaFile("C:/Users/nmero/OneDrive - Nexus365/Log/alpha.txt", std::ios::app);
	alphaFile << sigma << "\t" << Kf << "\t" << calculated_alpha << "\n";
	alphaFile.close();
	mlock.unlock();
}
void spawnThread(int n1, double variable1[],int n2,double variable2[], int n3, double variable3[] ) {
	std::cout <<"Maximum Concurrency is: " << (int)std::thread::hardware_concurrency <<"\n";
	assert((int)(std::thread::hardware_concurrency) >= static_cast<int>(n1 * n2 * n3) && "Threads exceeding Maximum Threads of Hardware!");
	std::vector<std::thread> threads(n1*n2*n3);
	//spawn n1*n2 threads
	int index = 0;

	for (int i = 0; i < n1; i++) {
		for (int j = 0; j < n2; j++) {
			for (int k = 0; k < n3; k++) {
				threads[index] = std::thread(foo, variable1[i], variable2[j], variable3[k]);
				index++;
			}
		}
	}
	for (auto& th : threads) {
		th.join();
	}
}
int main() {


	double dimScanRate = 1.0; //dimensional scanRate

	double scanRate = (dElectrode) * (dElectrode) / (DA) * (96485 / (8.314 * 298)) * dimScanRate; 
	std::cout << "Scan Rate is " << scanRate << "\n";                                                            
	//double sigma[] = {scanRate*0.001,scanRate*0.002,scanRate*0.003,scanRate*0.004,scanRate*0.005,scanRate*0.006,scanRate*0.007,scanRate*0.008,scanRate*0.009,scanRate*0.01,scanRate*0.02, scanRate*0.03,scanRate*0.04,scanRate*0.05,scanRate*0.06,scanRate*0.07,scanRate*0.08,scanRate*0.09,scanRate*0.10,scanRate*0.20,scanRate*0.30,scanRate*0.40,scanRate*0.50,scanRate*0.60,scanRate*0.70,scanRate*0.80,scanRate*0.90,scanRate*1.0,scanRate*2.0,scanRate*3.0,scanRate*4.0,scanRate*5.0,scanRate*6.0,scanRate*7.0,scanRate*8.0,scanRate*9.0,scanRate*10,scanRate*20,scanRate*30,scanRate*40,scanRate*50,scanRate*60,scanRate*70,scanRate*80,scanRate*90,scanRate*100};
	double sigma[] = { scanRate*1e-3,scanRate*1e-2,scanRate*1e-1,scanRate*1e0,scanRate*1e1,scanRate*1e2};
	//double K1[] = {0,1e0,pow(10,0.5),1e1,pow(10,1.5),1e2,pow(10,2.5),1e3,pow(10,3.5),1e4,pow(10,4.5),1e5,pow(10,5.5),1e6,pow(10,6.5),1e7,pow(10,7.5),1e8,pow(10,8.5),1e9,pow(10,9.5),1e10};
	double variables3[] = { 1e4 };  //dimensional Kf
	//double K1[] = { 0.0,1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e 9,1e10,1e11,1e12,1e13,1e14,1e15,1e16,1e17,1e18,1e19,1e20,1e21,1e22,1e23,1e24,1e25,1e26,1e27,1e28,1e29,1e30 };
	//double K1[]{0,1e0,1e1,1e3,1e5,1e6,1e7,1e8,1e9,1e10 }; 
	double variables2[] = {1e-5}; //dimensional keq
	//double Vary[] = { 1e-1,1e-2,1e-3,1e-4 }; //Vary CAsyar
	//double Vary[] = { 0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.15,0.20}; // now vary expanding factor 
	//double Vary[] = { 5e-3,1e-3,5e-4,1e-4,5e-5,1e-5,5e-6,1e-6 }; //now vary dX
	//double Vary[] = { 1e-1,8e-2,5e-2,3e-2,1e-2,8e-3,5e-3,3e-3,1e-3,8e-4,5e-4,3e-4,1e-4 }; // vary DTheta

	/*double Vary[120];
	int count = 0;
	for (double i = 4.05; i < 10.01; i = i + 0.05) {
		Vary[count] = pow(10, i);
		count++;
	}*/


	//double dB[]{ 0.1,1,10 };
	//double dY[]{ 0.1,1,10 };
	writeStartTimeToLog();
	spawnThread(std::size(sigma), sigma, std::size(variables2),variables2, size(variables3), variables3);
	//spawnThread(std::size(sigma), sigma, std::size(dB), dB, std::size(dY), dY);
	writeEndTimeToLog();
	endAlphaToLog();
	
	/*for (double i = 2; i > -4.1;i=i-0.5) {
		sigma[0] = scanRate*pow(10, i);
		writeStartTimeToLog();
		startAlphaToLog();
		spawnThread(std::size(sigma), sigma, std::size(K1), K1, size(K0), K0);
		//spawnThread(std::size(sigma), sigma, std::size(dB), dB, std::size(dY), dY);
		writeEndTimeToLog();
		endAlphaToLog();
	}*/






	std::cout << "\a";
	return 0;
}