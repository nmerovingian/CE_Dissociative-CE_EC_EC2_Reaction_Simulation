#include<iostream>

class ThomasXL {
public:
	int n;
	int nn; //bandidth
	double** L;
	double** U;
	double* Y;

	ThomasXL(int n, int numOfDiag);
	void LUDecomposition(double** A);
	void ForwardSubstitution(double* b);
	void BackwardSubstitution(double* x);
	void solve(double** A, double* b, double* x);


};