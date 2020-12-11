#include <iostream>
#include "Trid.h"

Trid::Trid(int nn){
	n = nn;
	cc = new double[n-1];
    dd = new double[n];
}

Trid::~Trid(){
	delete[] cc;
	delete[] dd;
}

void Trid::trid(double *a, double *b, double *c, double *d, double *x){
	cc[0] = c[0]/b[0];
	dd[0] = d[0]/b[0];
	for(int i = 1; i < n-1; i++){
		cc[i] = c[i]/(b[i]-a[i]*cc[i-1]);
		dd[i] = (d[i]-a[i]*dd[i-1])/(b[i]-a[i]*cc[i-1]);
	} 
	dd[n-1] = (d[n-1]-a[n-1]*dd[n-2])/(b[n-1]-a[n-1]*cc[n-2]);

	x[n-1] = dd[n-1];
	for(int i = n-2; i >= 0; i--){
		x[i] = dd[i] - cc[i]*x[i+1];
    }	
}