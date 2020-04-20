#pragma once
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <chrono> 
using namespace std::chrono;
using namespace std;


void sorSolve_AM(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& rhs,
	std::vector<double>& x, int iterMax, double tol, double omega, int& sor, double dS, double cp, double t0, int i, double dt)
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
			error += std::fabs(rhs[0] - b[0] * x[0] - c[0] * x[1]);
			if (i * dt < t0) {
				x[0] = max(0., min(cp, x[0] + omega * (y - x[0])));
			}
			else {
				x[0] = max(x[0] + omega * (y - x[0]), 0.);
			}
		}
		for (int j = 1; j < n; j++)
		{
			double y = (rhs[j] - a[j] * x[j - 1] - c[j] * x[j + 1]) / b[j];
			error += std::fabs(rhs[j] - a[j] * x[j - 1] - b[j] * x[j] - c[j] * x[j + 1]);
			if (i * dt < t0) {
				x[j] = max(min(x[j] + omega * (y - x[j]),cp), j * dS);
			}
			else {
				x[j] = max(x[j] + omega * (y - x[j]), j * dS);
			}
		}
		{
			double y = (rhs[n] - a[n] * x[n - 1]) / b[n];
			error += std::fabs(rhs[n] - a[n] * x[n - 1] - b[n] * x[n]);
			if (i * dt < t0) {
				x[n] = max(min(x[n] + omega * (y - x[n]),cp), n * dS);
			}
			else {
				x[n] = max(x[n] + omega * (y - x[n]), n * dS);
			}
		}
		// make an exit condition when solution found
		if (error < tol)
			break;
	}
}

/* Solution code for the Crank Nicolson Finite Difference
search for COURSEWORK EDIT for parts that needed to be altered for the coursework
 */

double crank_nicolson_AM_LINEAR(double S0, double X, double F, double T, double r, double sigma,
	double R, double kappa, double mu, double C, double alpha, double beta, int iMax, int jMax, int S_max, double tol, double omega, int iterMax, double cp, double t0)
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
		sorSolve_AM(a, b, c, d, vNew, iterMax, tol, omega, sor, dS, cp, t0, i, dt);
		if (sor == iterMax)
			//cout << "NOT SOLVED" << endl;
			//break;

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

//Get data for ameribond
void getAmeribondData() {
	// declare and initialise Black Scholes parameters - Currently looking at a solution we can get a definite answer for
	double T = 3., F = 56., R = 1., r = 0.0038, kappa = 0.0833333333,
		mu = 0.0073, X = 56.47, C = 0.106, alpha = 0.01, beta = 0.425, sigma = 3.73, tol = 1.e-8, omega = 1., S_max = 10 * X;
	//
	int iterMax = 10000, iMax = 40, jMax = 25;
	beta = 0.425;
	sigma = 3.73;
	double S0 = X;
	double t0 = 1.2448, cp = 67;

	//Checking value against theory
	std::ofstream analytical("./Ameribond1.txt");
	for (int s = 1; s <= 100; s++) {
		analytical << s << " , " << crank_nicolson_AM_LINEAR(s, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax, cp, t0) << "\n";
	}
	cout << "AMERICAN BOND PART 1 DONE" << endl;

	//Checking how the value changes with different values of r
	std::ofstream r_file("./Ameribond_r.txt");
	for (int s = 1; s <= 66; s++) {
		r_file << s << " , " <<
			crank_nicolson_AM_LINEAR(s, X, F, T, 0.0019, sigma, R, kappa, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax,cp, t0) << " , " <<
			crank_nicolson_AM_LINEAR(s, X, F, T, 0.0038, sigma, R, kappa, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax,cp, t0) << " , " <<
			crank_nicolson_AM_LINEAR(s, X, F, T, 0.0057, sigma, R, kappa, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax, cp, t0) <<
			"\n";

	}
	cout << "AMERICAN BOND PART 2 DONE" << endl;
		//auto start = high_resolution_clock::now();
		//cout << crank_nicolson_AM_LINEAR(14, X, F, T, 0.0057, sigma, R, kappa, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax) << endl;
		//auto stop = high_resolution_clock::now();
		//auto duration = duration_cast<milliseconds>(stop - start);
		//cout << duration.count() << endl;

}