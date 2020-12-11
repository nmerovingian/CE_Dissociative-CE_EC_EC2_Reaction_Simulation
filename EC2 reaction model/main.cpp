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
double dElectrode = 0.005641895835478;
double DA = 1e-9;
void spawnThread(int n1, double variable1[],int n2, double variable2[], int n3, double variable3[]);
void foo(double variable1, double variable2, double variable3);
void foo(double variable1, double variable2,double variable3) {
	//variable 1 is scan rate, Variable 3 is K0. Now change Varible 2 
	//Input parameters
	double concA = 1.0;
	double concB = 0.0;
	double concY = 0.0;

	double cAstar = 1e-3;

	//double theta_i = 20.0;
	double theta_i = 40.0;
	double theta_v = -40.0;
	double sigma = variable1;


	double deltaX = 1e-6;
	double deltaTheta = 1e-2;
	double K0 = 1e-1;

	double dimKf = variable3;
	double dimKeq = variable2;
	double dimKb = dimKf / dimKeq;
	double K1 = dimKf*dElectrode*dElectrode*cAstar / DA;
	double Kb = dimKb * dElectrode * dElectrode / DA;
	double alpha = 0.5;
	double gamma = 0.05;
	//if (sigma > 38943) {
	//	gamma = 0.02;
	//}

	double dB = 1.0;
	double dY = 1.0;
	double epsilon = 1e-3;
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
	double pointB = -20.0;
	double pointC = theta_v;
	double pointD = 0.0;
	double pointE = 25;

	//Calcute the maximums of time and position 
	double deltaT = deltaTheta / sigma;
	double maxT = 2.0 * fabs(theta_v - theta_i) / sigma;
	double maxX = 6.0 * sqrt(maxT); 
	Coeff coeff(maxX, K0, K1,Kb, deltaT, alpha, gamma, dB, dY);
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
	for (int i = 0; i < 30; i++) {

		coeff.calc_jacob(grid.conc, theta_i);
		coeff.calc_fx(grid.conc, theta_i);
		Trid.solve(coeff.J, coeff.fx, coeff.dx);
		coeff.xupdate(grid.conc);
	}
	//double logVar2 = log10(variable2);
	std::string BKineticsLocation  = "Kinetics concB =" + to_sci(concB) + "K1 = " + to_sci(K1) + "Sigma=" + to_sci(sigma) +"K0=" + to_sci(K0) + "dB=" +to_sci(dB) + "dY=" + to_sci(dY);
	//std::string CVLocation = "Variable =" + to_sci(variable2) + " Sigma="; //+"K1=" + to_sci(K1) + "Sigma=" + to_sci(sigma) + "K0=" + to_sci(K0) + "dB=" + to_sci(dB) + "dY=" + to_sci(dY);
	std::string CVLocation = "var1=" + to_string(variable1) + "var2=" + to_string(dimKf) + "var3=" + to_string(dimKeq); //+"K1=" + to_sci(K1) + "Sigma=" + to_sci(sigma) + "K0=" + to_sci(K0) + "dB=" + to_sci(dB) + "dY=" + to_sci(dY);


	std::string ConcALocation = "concA =" + to_string(variable1) + "var2=" + to_string(dimKf) + "var3=" + to_string(dimKeq);
	std::string ConcBLocation = "concB =" + to_string(variable1) + "var2=" + to_string(dimKf) + "var3=" + to_string(dimKeq);
	std::string ConcYLocation = "concY =" + to_string(variable1) + "var2=" + to_string(dimKf) + "var3=" + to_string(dimKeq);


	//Calculate the first flux and write it in output file
	grid.grad();
	//std::cout << "first flux is " << grid.g << "\n";
	std::ofstream myFile;
	myFile.open(genAddress(CVLocation));
	myFile << theta_i << "," << -grid.g << "\n";
	CV_theta.push_back(theta_i);
	CV_flux.push_back(-grid.g);
	mlock.lock();
	cout << "Time step is " << m << "\n";
	mlock.unlock();
	std::ofstream bFile;
	bFile.open(genAddress(BKineticsLocation));
	bFile << 0 << "," << concB << "\n";


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
		if (numOfIteration < maxNumOfIteration && coeff.max_dx() > 1e-12) {
			
			if (numOfIteration < maxNumOfIteration) {
				numOfIteration++;
				//std::cout << "/n/nNot converging!/n/n" << K1 << "\t" << "avg_dx()is" << coeff.avg_dx() << "\n" << "max_dX is" << coeff.max_dx() << "\n";
				//std::cout << "Num of Iteration Now is: " << numOfIteration << "\n";
			}
			
		}
		grid.updateAll();
		grid.grad();
		//std::cout << "grid" << grid.g << "\t";
		fluxTotal += -grid.g * deltaT;
		myFile << Theta << "," << -grid.g << "\n";
		CV_theta.push_back(Theta);
		CV_flux.push_back(-grid.g);
		bFile << double (i)/double(m)*maxT << "," << grid.concB[1] << "\n";
		// find the maximum flux
		if (-grid.g > max_flux) {
			max_flux = -grid.g;
			backwardPeak = Theta;

		}
		//find the minimum flux
		if (-grid.g < min_flux) {
			min_flux = -grid.g;
			forwardPeak = Theta;
		}
		if (printA && (fabs(Theta - pointA) < epsilon && i < m / 2)) {
			std::string str = "Point A Forward Scan Half Max Theta=";
			str += to_string(Theta);

			grid.saveA(genAddress("concA" + str));
			grid.saveB(genAddress("concB" + str));
			grid.saveY(genAddress("ConcY" + str));
			printA = false;
		}
		if (printB && fabs((Theta - (pointB)) < epsilon && i < m / 2)) {
			std::string str = "Point B Forward Scan Peak Theta=";
			str += to_string(Theta);

			grid.saveA(genAddress("concA" + str));
			grid.saveB(genAddress("concB" + str));
			grid.saveY(genAddress("concY" + str));
			printB = false;
		}
		if (printC && fabs(Theta - pointC) < epsilon) {
			std::string str = "Point C Theta_v Theta=";
			str += to_string(Theta);

			grid.saveA(genAddress("concA" + str));
			grid.saveB(genAddress("concB" + str));
			grid.saveY(genAddress("concY" + str));
			printC = false;
		}
		if (printD && (fabs(Theta - pointD) < epsilon && i > m / 2)) {
			std::string str = "Point D Backward Scan Half Peak Theta=";
			str += to_string(Theta);

			grid.saveA(genAddress("concA" + str));
			grid.saveB(genAddress("concB" + str));
			grid.saveY(genAddress("concY" + str));
			printD = false;
		}
		if (printE && (fabs(Theta - pointE) < epsilon && i > m / 2)) {
			std::string str = "Point E Backward Scan Peak Theta=";
			str += to_string(Theta);

			grid.saveA(genAddress("concA" + str));
			grid.saveB(genAddress("concB" + str));
			grid.saveY(genAddress("concY" + str));
			printE = false;
		}
	}
	//std::cout.precision(9);
	double predIrrFlux = -0.496 * sqrt(alpha) * sqrt(sigma);
	std::cout<<"K0= " << K0 << "\t" << "K1= " << K1 << "\n";
	std::cout << "dB= " << dB << "\t" << "dY= " << dY << "\n";
	double diffFromPredIrrFlux = (min_flux - predIrrFlux) / predIrrFlux * 100;
	std::cout << "Difference from Predicted Irr Flux is " << diffFromPredIrrFlux << "% \n";
	double predRevFlux = -0.446 * sqrt(sigma);
	double diffFromPredRevFlux = (min_flux - predRevFlux) / predRevFlux * 100;
	std::cout << "Difference from Predicted Rev Flux is " << diffFromPredRevFlux << "% \n";
	std::cout << "Min Flux is " << min_flux << "\n" << "Forward Scan Peak at " << forwardPeak << "\n";
	std::cout << "Max Flux is " << max_flux << "\n" << "Backward Scan Peak at " << backwardPeak << "\n";
	std::cout << "Peak Separation is " << fabs(forwardPeak) + fabs(backwardPeak) << " \n";
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
		grid.saveA(genAddress(ConcALocation));
		grid.saveB(genAddress(ConcBLocation));
		grid.saveY(genAddress(ConcYLocation));
	}

	// Write log file in here 
	mlock.lock();
	std::ofstream logFile("Your local address/log.txt", std::ios::app); //append to logfile
	logFile << "\n";
	logFile << "Thread ID=" << this_id << "\n";
	logFile << "Sigma = " << to_sci(sigma) << "\t" << "K0 = " << to_sci(K0) << "\t" << "K1 = " << to_sci(K1) << "\n";
	logFile << "dB= " << dB << "\t" << "dY= " << dY << "\n";
	logFile << "MaxX = " << maxX << "\t" <<"MaxT = " << maxT << "\n";
	logFile << "Min Flux is " << min_flux << "\t" << "Forward Scan Peak at " << forwardPeak << "\n";
	logFile << "Max Flux is " << max_flux << "\t" << "Backward Scan Peak at " << backwardPeak << "\n";
	logFile << "Difference from Predicted Irr Flux is " << diffFromPredIrrFlux << "% \n";
	logFile << "Difference from Predicted Rev Flux is " << diffFromPredRevFlux << "% \n";
	logFile << "Peak Separation is " << fabs(forwardPeak) + fabs(backwardPeak) << " \n";
	logFile << "flux_integrate is: " << fluxTotal << "\n";
	logFile << "mass Conservation is: " << massConservation << "%\n";
	logFile << "calculated alpha is: " << calculated_alpha << "\t" << "alpha is: " << alpha << "\n";
	logFile << "deltaX is: " << deltaX << "\n";
	logFile << "deltaTheta is: " << deltaTheta << "\n";
	logFile << "gamma is: " << gamma << "\n";
	logFile << "\n";
	logFile.close();
	std::ofstream alphaFile("Your local address/alpha.txt", std::ios::app);
	alphaFile << sigma << "\t" << K1 << "\t" << calculated_alpha << "\n";
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


	double dimScanRate = 1.0; 

	double scanRate = (dElectrode) * (dElectrode) / (DA) * (96485 / (8.314 * 298)) * dimScanRate;
	std::cout << "Scan Rate is " << scanRate << "\n";
	double sigma[] = {scanRate,scanRate*0.1,scanRate*1e-2,scanRate*1e-3,scanRate*1e1,scanRate*1e2};

	double variables3[] = { 1e9 }; 
	double variables2[] = { 1e9 }; 
	writeStartTimeToLog();
	spawnThread(std::size(sigma), sigma, std::size(variables2),variables2, size(variables3), variables3);

	writeEndTimeToLog();
	endAlphaToLog();
	

	std::cout << "\a";
	return 0;
}