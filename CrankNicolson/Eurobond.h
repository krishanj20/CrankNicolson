#pragma once
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <chrono> 
using namespace std::chrono;
using namespace std;

void sorSolve_EURO(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& rhs,
	std::vector<double>& x, int iterMax, double tol, double omega, int& sor)
{
	// assumes vectors a,b,c,d,rhs and x are same size (doesn't check)
	int n = a.size() - 1;
	// sor loop
	for (sor = 0; sor < iterMax; sor++)
	{
		double error = 0.;
		// implement sor in here
		{
			double y = (rhs[0] - c[0] * x[1]) / b[0];
			x[0] = x[0] + omega * (y - x[0]);
		}
		for (int j = 1; j < n; j++)
		{
			double y = (rhs[j] - a[j] * x[j - 1] - c[j] * x[j + 1]) / b[j];
			x[j] = x[j] + omega * (y - x[j]);
		}
		{
			double y = (rhs[n] - a[n] * x[n - 1]) / b[n];
			x[n] = x[n] + omega * (y - x[n]);
		}
		// calculate residual norm ||r|| as sum of absolute values
		error += std::fabs(rhs[0] - b[0] * x[0] - c[0] * x[1]);
		for (int j = 1; j < n; j++)
			error += std::fabs(rhs[j] - a[j] * x[j - 1] - b[j] * x[j] - c[j] * x[j + 1]);
		error += std::fabs(rhs[n] - a[n] * x[n - 1] - b[n] * x[n]);
		// make an exit condition when solution found
		if (error < tol)
			break;
	}
}

/* Solution code for the Crank Nicolson Finite Difference
search for COURSEWORK EDIT for parts that needed to be altered for the coursework
 */

//This contains linear interpolation at the end
double crank_nicolson_E_LINEAR(double S0, double X, double F, double T, double r, double sigma,
	double R, double kappa, double mu, double C, double alpha, double beta, int iMax, int jMax, int S_max, double tol, double omega, int iterMax)
{
	// declare and initialise local variables (ds,dt)
	double dS = S_max / jMax;
	double dt = T / iMax;
	// create storage for the stock price and option price (old and new)
	vector<double> S(jMax + 1), vOld(jMax + 1), vNew(jMax + 1);
	// setup and initialise the stock price
	for (int j = 0; j <= jMax; j++)
	{
		S[j] = j * dS;
	}
	// setup and initialise the final conditions on the option price
	for (int j = 0; j <= jMax; j++)
	{
		vOld[j] = max(F, R * S[j]);
		vNew[j] = max(F, R * S[j]);
	}
	// start looping through time levels
	for (int i = iMax - 1; i >= 0; i--)
	{
		// declare vectors for matrix equations
		vector<double> a(jMax + 1), b(jMax + 1), c(jMax + 1), d(jMax + 1);
		// set up matrix equations a[j]=
		double theta = (1 + mu) * X * exp(mu * i * dt);
		a[0] = 0;
		b[0] = (-1 / dt) - (r / 2) - (kappa * theta / dS);
		c[0] = (kappa * theta / dS);
		d[0] = (-C * exp(-alpha * i * dt)) + (vOld[0] * (-(1 / dt) + (r / 2)));

		for (int j = 1; j <= jMax - 1; j++)
		{
			a[j] = (pow(sigma, 2) * pow(j * dS, 2 * beta) / (4 * pow(dS, 2))) - (kappa * (theta - j * dS) / (4 * dS));
			b[j] = (-1 / dt) - ((pow(sigma, 2.) * pow(j * dS, 2. * beta)) / (2. * pow(dS, 2))) - (r / 2.);
			c[j] = ((pow(sigma, 2.) * pow(j * dS, 2. * beta)) / (4. * pow(dS, 2.))) + ((kappa * (theta - j * dS)) / (4. * dS));
			d[j] = (-vOld[j] / dt) - ((pow(sigma, 2.) * pow(j * dS, 2. * beta) / (4. * pow(dS, 2.))) * (vOld[j + 1] - 2. * vOld[j] + vOld[j - 1])) - (((kappa * (theta - j * dS)) / (4. * dS)) * (vOld[j + 1] - vOld[j - 1])) + ((r / 2.) * vOld[j]) - (C * exp(-alpha * dt * i));
		}
		double A = R * exp((kappa + r) * (i * dt - T));
		double B = X * R * (1 - exp((kappa + r) * (i * dt - T))) + (C / (alpha + r)) * (exp(-alpha * i * dt) - exp(-alpha * T));
		B = X * R * exp((kappa + r) * (dt * i - T)) + C * exp(-alpha * i * dt) / (alpha + r) + R * exp(r * (i * dt - T)) - C * exp(-(alpha + r) * T + r * i * dt);
		B = -X * A + C * exp(-alpha * i * dt) / (alpha + r) + X * R * exp(r * (i * dt - T)) - C * exp(-(alpha + r) * T + r * i * dt) / (alpha + r);
		a[jMax] = 0;
		b[jMax] = 1;
		c[jMax] = 0;
		d[jMax] = dS * jMax * A + B;
		int sor;
		// solve matrix equations with SOR
		sorSolve_EURO(a, b, c, d, vNew, iterMax, tol, omega, sor);
		if (sor == iterMax)
			cout << "\n NOT CONVERGED \n";

		// set old=new
		vOld = vNew;
	}
	// finish looping through time levels

	// output the estimated option price
	double sum;
	{
		int jStar = S0 / dS;
		sum = 0.;
		if (jStar > 0 && jStar < jMax) {
			sum += ((S0 - S[jStar]) * (S0 - S[jStar + 1]) / (2 * dS * dS)) * vNew[jStar - 1];
			sum -= ((S0 - S[jStar - 1]) * (S0 - S[jStar + 1]) / (dS * dS)) * vNew[jStar];
			sum += ((S0 - S[jStar - 1]) * (S0 - S[jStar]) / (2 * dS * dS)) * vNew[jStar + 1];
		}
		else {
			sum += (S0 - S[jStar]) / dS * vNew[jStar + 1];
			sum += (S[jStar + 1] - S0) / dS * vNew[jStar];
		}
	}
	return sum;
}


/* Template code for the Crank Nicolson Finite Difference
 */
//This contains quadratic interpolation at the end
double crank_nicolson_E_QUAD(double S0, double X, double F, double T, double r, double sigma,
	double R, double kappa, double mu, double C, double alpha, double beta, int iMax, int jMax, int S_max, double tol, double omega, int iterMax)
{
	// declare and initialise local variables (ds,dt)
	double dS = S_max / jMax;
	double dt = T / iMax;
	// create storage for the stock price and option price (old and new)
	vector<double> S(jMax + 1), vOld(jMax + 1), vNew(jMax + 1);
	// setup and initialise the stock price
	for (int j = 0; j <= jMax; j++)
	{
		S[j] = j * dS;
	}
	// setup and initialise the final conditions on the option price
	for (int j = 0; j <= jMax; j++)
	{
		vOld[j] = max(F, R * S[j]);
		vNew[j] = max(F, R * S[j]);
	}
	// start looping through time levels
	for (int i = iMax - 1; i >= 0; i--)
	{
		// declare vectors for matrix equations
		vector<double> a(jMax + 1), b(jMax + 1), c(jMax + 1), d(jMax + 1);
		// set up matrix equations a[j]=
		double theta = (1 + mu) * X * exp(mu * i * dt);
		a[0] = 0;
		b[0] = (-1 / dt) - (r / 2) - (kappa * theta / dS);
		c[0] = (kappa * theta / dS);
		d[0] = (-C * exp(-alpha * i * dt)) + (vOld[0] * (-(1 / dt) + (r / 2)));

		for (int j = 1; j <= jMax - 1; j++)
		{
			//
			a[j] = (pow(sigma, 2) * pow(j * dS, 2 * beta) / (4 * pow(dS, 2))) - (kappa * (theta - j * dS) / (4 * dS));
			b[j] = (-1 / dt) - ((pow(sigma, 2.) * pow(j * dS, 2. * beta)) / (2. * pow(dS, 2))) - (r / 2.);
			c[j] = ((pow(sigma, 2.) * pow(j * dS, 2. * beta)) / (4. * pow(dS, 2.))) + ((kappa * (theta - j * dS)) / (4. * dS));
			d[j] = (-vOld[j] / dt) - ((pow(sigma, 2.) * pow(j * dS, 2. * beta) / (4. * pow(dS, 2.))) * (vOld[j + 1] - 2. * vOld[j] + vOld[j - 1])) - (((kappa * (theta - j * dS)) / (4. * dS)) * (vOld[j + 1] - vOld[j - 1])) + ((r / 2.) * vOld[j]) - (C * exp(-alpha * dt * i));
			//
		}
		double A = R * exp((kappa + r) * (i * dt - T));
		double B = -X * A + C * exp(-alpha * i * dt) / (alpha + r) + X * R * exp(r * (i * dt - T)) - C * exp(-(alpha + r) * T + r * i * dt) / (alpha + r);
		a[jMax] = 0;
		b[jMax] = 1;
		c[jMax] = 0;
		d[jMax] = jMax * dS * A + B;
		int sor;
		// solve matrix equations with SOR
		sorSolve_EURO(a, b, c, d, vNew, iterMax, tol, omega, sor);
		//vNew = thomasSolve(a, b, c, d);

		if (sor == iterMax)
			return -1;

		// set old=new
		vOld = vNew;
	}
	// finish looping through time levels

	// output the estimated option price
	double optionValue;
	{
		int jStar = S0 / dS;
		double sum = 0.;
		sum += (S0 - S[jStar]) / (dS)* vNew[jStar + 1];
		sum += (S[jStar + 1] - S0) / (dS)* vNew[jStar];
		optionValue = sum;
	}
	return optionValue;
}

//This function creates a txt file (sigmabeta.txt) Which explores the effect of changinging sigma and beta on option value
void Getsigmabeta() {
	//Creates three csv files with increasing sigma at a given beta for a fixed s0
	//	double S0 = X;
	double T = 3., F = 56., R = 1., r = 0.0038, kappa = 0.083333333333,
		mu = 0.0073, X = 56.47, C = 0.106, alpha = 0.01, beta = 1., sigma = 3.73, tol = 1.e-7, omega = 1., S_max = 10 * X;
	//
	int iterMax = 10000, iMax = 100, jMax = 200;
	double S0 = X;
	jMax = 40;
	iMax = 26;
	double sMax = 20 * X;
	beta = 0.425;
	sigma = 3.73;

	//Explore effect of Smax
	//Look at given imax and jmax, then increase Smax
	std::ofstream sigmabeta("./sigmabeta.txt");
	for (int i = 0; i < 11; i++) {
		sigmabeta << i * 0.5 << " , " <<
			crank_nicolson_E_LINEAR(X, X, F, T, r, sigma = i*0.5, R, kappa, mu, C, alpha, 0.2, iMax, jMax, sMax, tol, omega, iterMax) << " , " <<
			crank_nicolson_E_LINEAR(X, X, F, T, r, sigma = i * 0.5, R, kappa, mu, C, alpha, 0.8, iMax, jMax, sMax, tol, omega, iterMax) << " , " <<
			crank_nicolson_E_LINEAR(X, X, F, T, r, sigma = i * 0.5, R, kappa, mu, C, alpha, 1.2, iMax, jMax, sMax, tol, omega, iterMax) <<
			"\n";
	}
	cout << "DONE V FUNC S_MAX" << endl;
}

//This attempts to help find the most efficient value for imax, jmax and Smax
void GetEfficientResult() {
	// declare and initialise Black Scholes parameters - Currently looking at a solution we can get a definite answer for
	double T = 3., F = 56., R = 1., r = 0.0038, kappa = 0.083333333333,
		mu = 0.0073, X = 56.47, C = 0.106, alpha = 0.01, beta = 1., sigma = 3.73, tol = 1.e-7, omega = 1., S_max = 10 * X;
	//
	int iterMax = 10000, iMax = 100, jMax = 200;
	double S0 = X;
	beta = 0.425;
	sigma = 3.73;
	jMax = 40*10;
	iMax = 26*10;
	double sMax = 200 * X;
	cout <<
		"imax  = " << iMax << "   " <<
		"jmax  = " << jMax << "   " <<
		"Smax  = " << sMax << "   ";
	auto start = high_resolution_clock::now();
	double V1 = crank_nicolson_E_QUAD(S0, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, iMax, jMax, sMax, tol, omega, iterMax);
	auto stop = high_resolution_clock::now();
	auto duration1 = duration_cast<microseconds>(stop - start);
	cout << "OPTION VALUE = " << V1 << "  ";
	cout << "DURATION (microseconds): " << duration1.count() << endl;
}

//This code creates csv files to allow us to explore the effect of imax, jmax and smax on time
void GetTimeData() {
	// declare and initialise Black Scholes parameters - Currently looking at a solution we can get a definite answer for
	double T = 3., F = 56., R = 1., r = 0.0038, kappa = 0.08333333,
		mu = 0.0073, X = 56.47, C = 0.106, alpha = 0.01, beta = 1., sigma = 3.73, tol = 1.e-7, omega = 1., S_max = 10 * X;
	//
	int iterMax = 10000, iMax = 100, jMax = 200;
	double S0 = X;
	beta = 0.425;
	sigma = 3.73;
	//Get data for increasing smax with time
	//Look at given imax and jmax, then increase Smax
	std::ofstream time_sMax("./time_sMax.txt");
	for (int i = 6; i < 20; i++) {
		time_sMax << i * X << " , ";
			auto start = high_resolution_clock::now();
			crank_nicolson_E_LINEAR(X, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, iMax, jMax, i * X, tol, omega, iterMax);
			auto stop = high_resolution_clock::now();
			auto duration = duration_cast<milliseconds>(stop - start);
			time_sMax << duration.count() << "\n";
	}
	//Get data for increasing imax with time
		//Look at given imax and jmax, then increase Smax
	std::ofstream time_iMax("./time_iMax.txt");
	for (int i = 0; i < 200; i++) {
		time_iMax << i << " , ";
		auto start = high_resolution_clock::now();
		crank_nicolson_E_LINEAR(S0, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, i, jMax, S_max, tol, omega, iterMax);
		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<milliseconds>(stop - start);
		time_iMax << duration.count() << "\n";
	}
	cout << "DONE iMax as function of time" << endl;
	//Get data for increasing jmax with timw
	std::ofstream time_jMax("./time_jMax.txt");
	for (int i = 1; i < 200; i++) {
		time_jMax << i << " , ";
		auto start = high_resolution_clock::now();
		crank_nicolson_E_LINEAR(S0, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, iMax, i, S_max, tol, omega, iterMax);
		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<milliseconds>(stop - start);
		time_jMax << duration.count() << "\n";
	}
}


//This populates the rest of the csv files on european bond data
void getEurobondData() {
	// declare and initialise Black Scholes parameters - Currently looking at a solution we can get a definite answer for
	double T = 3., F = 56., R = 1., r = 0.0038, kappa = 0.08333333,
		mu = 0.0073, X = 56.47, C = 0.106, alpha = 0.01, beta = 1., sigma = 3.73, tol = 1.e-7, omega = 1., S_max = 10 * X;
	//
	int iterMax = 10000, iMax = 100, jMax = 25;
	beta = 0.425;
	sigma = 3.73;
	double S0 = X;
	double V1 = 0, V2 = 0;

	//Accuracy code
	for (int i = 1; i < 400; i++) {

		auto start = high_resolution_clock::now();
		V1 = crank_nicolson_E_QUAD(S0, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, 10 * 26, 10 * 40, 10 * 20 * X, tol, omega, iterMax);
		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(stop - start);
		if (abs(V2 - V1) < 1e-5) {
			cout << "S_max = " << i << endl;
			cout << "OPTION VALUE = " << V1 << endl;
			cout << "DURATION (microseconds):" << duration.count() << endl;
			break;
		}
		else {
			V2 = V1;
		}
	}
	//



	// declare and initialise grid paramaters 
	//int iMax = 100, jMax = 25;
	// declare and initialise local variables (ds,dt)
	//double S_max = 10 * X;
	//int iterMax = 5000000;
	//double tol = 1.e-7, omega = 1.;

	//Checking value against theory
	std::ofstream analytical("./analytical.txt");
	for (int s = 1; s <= 300; s++) {
		analytical << s << " , " << crank_nicolson_E_LINEAR(s, X, F, T, r, 3.73, R, 0, mu, C, alpha, 1., iMax, jMax, S_max, tol, omega, iterMax) << "\n";
	}

	//Create graph of varying  S and optionvalue
	int length = 50;
	double S_maxi = 4 * X;
	std::ofstream V_function_S("./V_function_S.txt");
	std::ofstream V_function_beta("./V_function_beta.txt");

	for (int j = 1; j <= length - 1; j++) {
		//Plottting V(S,t) as a function of S (t=0) for two cases
		beta = 1.;
		sigma = 0.369;
		//Puts value of S, and value of stock for different parameters into csv file
		V_function_S << j * S_maxi / length << " , " << crank_nicolson_E_LINEAR(j * S_maxi / length, X, F, T, r, 0.369, R, kappa, mu, C, alpha, 1., iMax, jMax,
			S_max, tol, omega, iterMax) << " , " << crank_nicolson_E_LINEAR(j * S_maxi / length, X, F, T, r, 3.73, R, kappa, mu, C, alpha, 0.425, iMax, 35,
				S_max, tol, omega, iterMax) << "\n";


		//"You may wish to explore different values of beta and sigma?" Do research before implementing this. But preliminary (BETA =1)
		V_function_beta << j * S_maxi / length << " , " <<
			crank_nicolson_E_LINEAR(j * S_maxi / length, X, F, T, r, 0.369, R, kappa, mu, C, alpha, 3., iMax, jMax, S_max, tol, omega, iterMax) << " , " <<
			crank_nicolson_E_LINEAR(j * S_maxi / length, X, F, T, r, 0.369, R, kappa, mu, C, alpha, 2, iMax, jMax, S_max, tol, omega, iterMax) << " , " <<
			crank_nicolson_E_LINEAR(j * S_maxi / length, X, F, T, r, 0.369, R, kappa, mu, C, alpha, 1, iMax, jMax, S_max, tol, omega, iterMax) << " , " <<
			crank_nicolson_E_LINEAR(j * S_maxi / length, X, F, T, r, 0.369, R, kappa, mu, C, alpha, 0.5, iMax, jMax, S_max, tol, omega, iterMax) << " , " <<
			//crank_nicolson2(j * S_maxi / length, X, F, T, r, 0.369, R, kappa, mu, C, alpha, 0, iMax, jMax, S_max, tol, omega, iterMax) <<
			"\n";

	}
	cout << "DONE V FUNC S" << endl;

	//Assume now,
//	double S0 = X;
	beta = 0.425;
	sigma = 3.73;

	//Explore effect of Smax
	//Look at given imax and jmax, then increase Smax
	std::ofstream V_function_sMax("./V_function_sMax.txt");
	for (int i = 6; i < 20; i++) {
		V_function_sMax << i * X << " , " <<
			crank_nicolson_E_LINEAR(X, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, 100, 50 * i, i * X, tol, omega, iterMax) << " , " <<
			crank_nicolson_E_LINEAR(X, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, 100, 2 * i, i * X, tol, omega, iterMax) << " , " <<
			crank_nicolson_E_LINEAR(X, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, 100, 2 * i, i * X, tol, omega, iterMax) <<
			"\n";
	}
	cout << "DONE V FUNC S_MAX" << endl;
	//Explore effect of imax
	//Look at given imax and jmax, then increase Smax
	std::ofstream V_function_iMax("./V_function_iMax.txt");
	for (int i = 0; i < 200; i++) {
		V_function_iMax << i << " , " <<
			crank_nicolson_E_LINEAR(S0, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, i, 15, S_max, tol, omega, iterMax) << " , " <<
			crank_nicolson_E_LINEAR(S0, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, i, 25, S_max, tol, omega, iterMax) << " , " <<
			crank_nicolson_E_LINEAR(S0, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, i, 50, S_max, tol, omega, iterMax) << " , " <<
			"\n";
	}
	cout << "DONE V FUNC iMax" << endl;
	//Explore effect of jmax
	//Look at given imax and jmax, then increase Smax
	std::ofstream V_function_jMax("./V_function_jMax.txt");
	for (int i = 1; i < 200; i++) {
		V_function_jMax << i << " , " <<
			crank_nicolson_E_LINEAR(S0, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, 20, i, S_max, tol, omega, iterMax) <<
			"\n";
	}
	cout << "DONE V FUNC jMax" << endl;

}