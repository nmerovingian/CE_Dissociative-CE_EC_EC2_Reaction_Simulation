#include <vector>
void Print2DArray(double** array, int n);
void Print1DArray(double* array, int n);
void estimateRunTime(int startTime);
std::string genAddress(std::string filename);
std::string to_sci(double num);
std::string to_sci(int num);
bool dcomp(double a, double b, double epsilon);
//bool dcomp(double a, double b);
void writeStartTimeToLog();
void writeEndTimeToLog();
void endAlphaToLog();
void startAlphaToLog();
double cal_alpha( const std::vector<double> &CV_theta, const std::vector<double> &CV_flux );
double cal_slope(const std::vector<double>& x, const std::vector<double> &y);
void colorP(std::string str, int color);

