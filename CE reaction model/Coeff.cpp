#include <iostream>
#include <cmath>
#include<boost/multiprecision/cpp_dec_float.hpp>
#include "Coeff.h"
#include "ThomasXL.h"
#include "helper.h"

using namespace std;
void Print2dArray(double** array, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			std::cout << array[i][j] << "\t";
		} 
		std::cout << "\n\n\n";
	}
}
void Print1dArray(double* array, int n) {
	for (int i = 0; i < n; i++) {
		std::cout << array[i] << "\t";
	}
	std::cout << "\n";
}
Coeff::Coeff(double maxX, double K0_input,double Kf_input, double Kb_input,double deltaT, double alpha_input, double cAstar_input, double gamma_input, double dB_input, double dY_input) {
    n = 0;
    xi = 0.0;
	dB = dB_input;
	dY = dY_input;
    xm = maxX;
	K0 = K0_input;
	Kf = Kf_input * deltaT;
	Kb = Kb_input * deltaT;
	alpha = alpha_input;
	cAstar = cAstar_input;
	gamma = gamma_input;


}

Coeff::~Coeff(){
    delete[] aA;
    delete[] bA;
    delete[] cA;
	delete[] aB;
	delete[] bB;
	delete[] cB;
	delete[] aY;
	delete[] bY;
	delete[] cY;
    delete[] d;
    delete[] XX;
	delete[] J;
	delete[] fx;
}
//Calculate the maximum spaceSteps
void Coeff::calc_n(double dx){
    while(xi < xm){
        xi += dx;
		dx *= (1 + gamma);
        n++;
    }
    n += 1;
	std::cout << "n is " << n << "\n";
	/*A = new double* [3*n];
	for (int i = 0; i < 3*n; i++) {
		A[i] = new double[3*n];
	}
	for (int i = 0; i < 3*n; i++) {
		for (int j = 0; j < 3*n; j++) {
			A[i][j] = 0.0;
		}
	}
	std::cout << "Initialize A in calc_nY_ThomasXL_expanding sucess" << " \n";*/
	aA = new double[n];
	bA = new double[n];
	cA = new double[n];
	aB = new double[n];
	bB = new double[n];
	cB = new double[n];
	aY = new double[n];
	bY = new double[n];
	cY = new double[n];
    d = new double[3*n];


	XX = new double[n];
}

void Coeff::get_XX(double *xx){
    for(int i = 0; i < n; i++){
        XX[i] = xx[i];
    }
	//Print1dArray(XX, n);
	std::cout << "Get X success!" << "\n";
}//get space steps

/*void Coeff::calc_abc(double deltaT, double Theta, double deltaX){
	f_theta = exp(-alpha * Theta);
  

	A[0][0] = 1.0 + deltaX * f_theta * K0 * (1.0 + exp(Theta)); //b
	A[0][1] = -1.0; //c
	// std::cout << "x[1] is: " << x[1];
	for (int i = 1; i < n - 1; i++) {
		deltaX_m = x[i] - x[i - 1];
		deltaX_p = x[i + 1] - x[i];
		A[i][i - 1] = -(2.0 * deltaT) / (deltaX_m * (deltaX_m + deltaX_p)); //a
		A[i][i + 1] = -(2.0 * deltaT) / (deltaX_p * (deltaX_m + deltaX_p)); //c
		A[i][i] = 1.0 - A[i][i - 1] - A[i][i + 1]; //b

	}
	A[n - 1][n - 2] = 0.0;
	A[n - 1][n - 1] = 1.0;
}*/
void Coeff::Acalc_abc(double deltaT, double Theta, double deltaX) {
	aA[0] = 0.0;
	bA[0] = 0.0;
	cA[0] = 0.0;
	for (int i = 1; i < n-1; i++) {
		deltaX_m = XX[i] - XX[i - 1];
		deltaX_p = XX[i + 1] - XX[i];
		aA[i] = -(2.0 * deltaT) / (deltaX_m * (deltaX_m + deltaX_p));
		cA[i] = -(2.0 * deltaT) / (deltaX_p * (deltaX_m + deltaX_p));
		bA[i] = 1.0 - aA[i] - cA[i];
	}
	aA[n - 1] = 0.0;
	bA[n - 1] = 0.0;
	cA[n - 1] = 0.0;
	//Print1dArray(aA, n);
	//std::cout << "calculate Aabc success!" << "\n";

}
void Coeff::Bcalc_abc(double deltaT, double Theta, double deltaX) {
	aB[0] = 0.0;
	bB[0] = 0.0;
	cB[0] = 0.0;
	for (int i = 1; i < n - 1; i++) {
		deltaX_m = XX[i] - XX[i - 1];
		deltaX_p = XX[i + 1] - XX[i];
		aB[i] = dB*(-(2.0 * deltaT) / (deltaX_m * (deltaX_m + deltaX_p)));
		cB[i] = dB*(-(2.0 * deltaT) / (deltaX_p * (deltaX_m + deltaX_p)));
		bB[i] = 1.0 - aB[i] - cB[i];
	}
	aB[n - 1] = 0.0;
	bB[n - 1] = 0.0;
	cB[n - 1] = 0.0;
	//Print1dArray(aB, n);
	//std::cout << "calculate Babc success!" << "\n";

}
void Coeff::Ycalc_abc(double deltaT, double Theta, double deltaX) {
	aY[0] = 0.0;
	bY[0] = 0.0;
	cY[0] = 0.0;
	for (int i = 1; i < n - 1; i++) {
		deltaX_m = XX[i] - XX[i - 1];
		deltaX_p = XX[i + 1] - XX[i];
		aY[i] = dY*(-(2.0 * deltaT) / (deltaX_m * (deltaX_m + deltaX_p)));
		cY[i] = dY*(-(2.0 * deltaT) / (deltaX_p * (deltaX_m + deltaX_p)));
		bY[i] = 1.0 - aY[i] - cY[i];
	}
	aY[n - 1] = 0.0;
	bY[n - 1] = 0.0;
	cY[n - 1] = 0.0;
	//Print1dArray(aY, n);
	//std::cout << "calculate Yabc success!" << "\n";

}
void Coeff::Allcalc_abc(double deltaT, double Theta, double deltaX) {
	Acalc_abc( deltaT, Theta, deltaX);
	Bcalc_abc( deltaT,  Theta,  deltaX);
	Ycalc_abc( deltaT,  Theta,  deltaX);
}

void Coeff::ini_jacob() {
	J = new double* [3*n];
	for (int i = 0; i < 3*n; i++) {
		J[i] = new double[3*n];
	}
	for (int i = 0; i < 3*n; i++) {
		for (int j = 0; j < 3*n; j++) {
			J[i][j] = 0.0;
		}
	}
	//Print2dArray(J, 3 * n);
	//std::cout << "Initialize J success" << " \n";
}
void Coeff::ini_fx() {
	fx = new double [3*n];
	for (int i = 0; i < 3 * n; i=i+3) {
		fx[i] = 0.0;
		fx[i + 1] = 0.0;
		fx[i + 2] = 0.0;
	}
	//std::cout << "Initialize fx Success!" << "\n";
}
void Coeff::ini_dx() {
	dx = new double[3 * n];
	for (int i = 0; i <3*n; i++) {
		dx[i] = 0.0;
	}
}
void Coeff::calc_fx(double* x, double Theta) {
	//have to solve a quadratic equation to solve the equation.
	//boost::multiprecision::cpp_dec_float_100 termA = 1.0;
	//boost::multiprecision::cpp_dec_float_100 termB = exp(-2.0 * Theta) * (0.5) * (1.0 / cAstar);
	//boost::multiprecision::cpp_dec_float_100 termC = -exp(-2.0 * Theta) * (0.5) * (1.0 / cAstar);
	//boost::multiprecision::cpp_dec_float_100 root_high_precision = (-termB + sqrt(termB * termB - (boost::multiprecision::cpp_dec_float_100)4.0 * termA * termC)) / ((boost::multiprecision::cpp_dec_float_100)2.0 * termA);
	//double root = root_high_precision.convert_to<double>();


	//std::cout << "\n\n\n concentration of A at surface of electrode  " << root << "\n";
	//double flux = -0.5 * (x[3] - x[1])/dB;
	//x[0] = root;
	double Kred = K0*exp(-alpha*Theta); //reduction
	double Kox = K0 * exp((1.0 - alpha) * Theta);
	double h = XX[1] - XX[0];
	//double flux = x[0] * x[0] * h * Kox - x[1] * h * Kred;
	// double fluxB = x[0] * x[0] *h* -0.5 * 1.0 / dB * Kox + x[1] * (0.5 * 1.0 / dB * h * Kred);
	//std::cout << "Kox" << Kox << " Kred" << Kred <<"fluxA is " << flux << " fluxB is" << fluxB << "\n";
	fx[0] = x[3]-x[0];
	fx[1] = x[1] + (1.0/dB)*Kred*h*x[1]-(1.0/dB)*Kox*h*x[2]-x[4];
	fx[2] = (-1.0 / dY) * Kred * h * x[1] + (1.0 / dY * Kox * h + 1.0) * x[2] - x[5];


	for (int j = 3, i = 1; j < 3 * n - 3; j = j + 3, i++) {
		fx[j] = aA[i] * x[3 * i - 3] + bA[i] * x[3 * i] + cA[i] * x[3 * i + 3] + Kf*x[3*i]  - Kb *x[3*i+1]-  d[3*i];
		fx[j + 1] = aB[i] * x[3 * i - 2] + bB[i] * x[3 * i + 1] + cB[i] * x[3 * i + 4] - Kf*x[3*i] + Kb*x[3*i+1] - d[3 * i + 1]; // 
		fx[j + 2] = aY[i] * x[3 * i - 1] + bY[i] * x[3 * i + 2] + cY[i] * x[3 * i + 5] - d[3 * i + 2]; //
	}
	fx[3 * n - 3] = x[3 * n - 3] - d[3*n-3];
	fx[3 * n - 2] = x[3 * n - 2] - d[3 * n - 2];
	fx[3 * n - 1] = x[3 * n - 1]-d[3*n-1];

	negative_fx();
	//std::cout << "fX is: " << '\n';
	//Print1dArray(fx, 3 * n);
}
void Coeff::negative_fx() {
	for (int i = 0; i < 3 * n; i++) {
		fx[i] = -fx[i];
	}
}
void Coeff::calc_jacob(double* x,double Theta) {

	//Initialize The First Three Rows of Jacobian
	double Kred = K0 * exp(-alpha * Theta); //reduction
	double Kox = K0 * exp((1.0 - alpha) * Theta);
	double h = XX[1] - XX[0];
	//Initialize The First Three Rows of Jacobian
	J[0][0] = -1.0;
	//J[0][1] = - Kred*h ;
	J[0][3] = 1.0;
	//J[1][0] = (1.0+(1.0/dB)*Kox*h);
	//J[1][1] = -(1.0 / dB) * Kred * h;
	J[1][1] = (1.0 + (1.0 / dB) * Kred * h);
	J[1][2] = -(1.0 / dB) * Kox * h;
	J[1][4] = -1.0;
	J[2][1] = (-1.0 / dY) * Kred * h;
	J[2][2] = (1.0 / dY * Kox * h + 1.0);
	J[2][5] = -1.0;
	for (int row = 3, i = 1; row < 3 * n-3; row = row + 3, i++) {
		//Initialzie Species A;
		J[row][row - 3] = aA[i];
		J[row][row] = bA[i] +Kf; 
		J[row][row + 1] = -Kb;
		J[row][row + 3] = cA[i];
		J[row + 1][row - 2] = aB[i];
		J[row + 1][row] = -Kf;
		J[row + 1][row + 1] = bB[i] + Kb;
		J[row + 1][row + 4] = cB[i];
		J[row + 2][row - 1] = aY[i];
		J[row + 2][row + 2] = bY[i];
		J[row + 2][row + 5] = cY[i];
	}
	J[3 * n - 3][3 * n - 3] = 1.0;
	J[3 * n - 2][3 * n - 2] = 1.0;
	J[3 * n - 1][3 * n - 1] = 1.0;
	//std::cout << "Jacobian is " << "\n";
	//Print2dArray(J, 3 * n);
}
//Update the d array
/*void Coeff::update(double *conc, double Theta, double deltaX){
	f_theta = exp(-alpha * Theta);
    for(int i = 0; i < n; i++){
        d[i] = conc[i];
    }
	d[0] = deltaX * f_theta * K0 *exp(Theta);
    d[n - 1] = 1.0;
}*/

void Coeff::update(double* x,double Theta, double deltaX, double concA,double concB,double concY) {
	for (int i = 0; i < n*3-3; i++) {
		d[i] = x[i];
	}
	d[0] = 1.0;
	d[1] = 0.0;
	d[2] = 0.0;
	d[3 * n - 3] = concA;
	d[3 * n - 2] = concB;
	d[3 * n - 1] = concY;

	//std::cout << "deltaArray is " << "\n";
	//Print1dArray(d, 3 * n);
}
void Coeff::xupdate(double* x) {
	//std::cout << "dX is" << "\n";
	//Print1dArray(dx, 3*n);
	for (int i = 0; i < 3 * n; i++) {
		x[i] += (dx[i]);
		//if (x[i] < -0.01 || x[i] >1.01) {
		//	std::string str = to_sci(x[i]) + " " + to_sci(dx[i]) + "\n";
		//	colorP(str,244);
		//}
		//assert(x[i] >= -0.01 && x[i] <= 1.01 && "Wrong Concentration");
		/*if (x[i] < 0.0) {
			x[i] = 0.0;
		}else if(x[i]>1.0) {
			x[i] = 1.0;
		}*/

	}
	//std::cout << "X array is " << "\n";
	//Print1dArray(x, 3*n);
	/*for (int i = 0; i < 3 * n; i++) {
		if (x[i] < 0.0) {
			x[i] = 0.0;
		}  
	}*/
}

double Coeff::avg_dx() {
	double avg = 0.0;
	for (int i = 0; i < 3 * n; i++) {
		avg += dx[i];
	}
	return avg / (double(n)*3);
}
double Coeff::max_dx(){
	double max = 0.0;
	for (int i = 0; i < 3 * n; i++) {
		if (fabs(dx[i]) > max) {
			max = fabs(dx[i]);
		}
	}
	return fabs(max);
}

