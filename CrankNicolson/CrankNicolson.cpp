#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <chrono> 
#include "Eurobond.h"
#include "A.h"
using namespace std::chrono;
using namespace std;

void sorSolve_AM(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& rhs,
	std::vector<double>& x, int iterMax, double tol, double omega, int& sor, double dS)
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
			x[0] = max(x[0] + omega * (y - x[0]), 0.);
		}
		for (int j = 1; j < n; j++)
		{
			double y = (rhs[j] - a[j] * x[j - 1] - c[j] * x[j + 1]) / b[j];
			x[j] = max(x[j] + omega * (y - x[j]), j * dS);
		}
		{
			double y = (rhs[n] - a[n] * x[n - 1]) / b[n];
			x[n] = max(x[n] + omega * (y - x[n]), n * dS);
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

double crank_nicolson_AM_LINEAR(double S0, double X, double F, double T, double r, double sigma,
	double R, double kappa, double mu, double C, double alpha, double beta, int iMax, int jMax, int S_max, double tol, double omega, int iterMax, double t0,double cp)
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
		for (sor = 0; sor < iterMax; sor++)
		{
			// SOR equations in here
			{
				double y = (d[0] - c[0] * vNew[1]) / b[0];
				//if (i * dt <= t0 ){
				//	vNew[0] = max(0., min(cp, vNew[0] + omega * (y - vNew[0])));
				//}
				//else {
				vNew[0] = max(vNew[0] + omega * (y - vNew[0]), 0.);
				//}
			}
			for (int j = 1; j < jMax; j++)
			{
				double y = (d[j] - a[j] * vNew[j - 1] - c[j] * vNew[j + 1]) / b[j];
				if (i * dt <= t0) {
					vNew[j] = max(R*j*dS, min(cp, vNew[j] + omega * (y - vNew[j]) ));
				}
				else {
					vNew[j] = max(vNew[j] + omega * (y - vNew[j]), R*j*dS);
				}
			}
			{
				double y = (d[jMax] - a[jMax] * vNew[jMax - 1]) / b[jMax];
				if (i * dt <= t0) {
					vNew[jMax] = max(R * jMax * dS, min(cp, vNew[jMax] + omega * (y - vNew[jMax])));
				}
				else {
					vNew[jMax] = max(vNew[jMax] + omega * (y - vNew[jMax]), R * jMax * dS);
				}
			}
			// calculate residual
			double error = 0.;
			error += fabs(d[0] - b[0] * vNew[0] - c[0] * vNew[1]);
			for (int j = 1; j < jMax; j++)
				error += fabs(d[j] - a[j] * vNew[j - 1] - b[j] * vNew[j] - c[j] * vNew[j + 1]);
			error += fabs(d[jMax] - a[jMax] * vNew[jMax - 1] - b[jMax] * vNew[jMax]);
			// check for convergence and exit loop if converged
			if (error < tol)
				break;
		}
		if (sor == iterMax)
			//cout << "\n NOT CONVERGED \n";
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
void main() {
	// declare and initialise Black Scholes parameters - Currently looking at a solution we can get a definite answer for
	double T = 3., F = 56., R = 1., r = 0.0038, kappa = 0.08333333,
		mu = 0.0073, X = 56.47, C = 0.106, alpha = 0.01, beta = 0.425, sigma = 3.73,tol = 1.e-1, omega = 1., S_max = 10 * X;
	//
	int iterMax = 10000, iMax = 40, jMax = 25;
	beta = 0.425;
	sigma = 3.73;
	double S0 = X;
	double t0 = 1.2448, cp = 67;
	/*
	//Checking value against theory
	std::ofstream analytical("./Ameribond1.txt");
	for (int s = 1; s <= 100; s++) {
		analytical << s << " , " << crank_nicolson_AM_LINEAR(s, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax, t0, cp) << "\n";
	}
	cout << "AMERICAN BOND PART 1 DONE" << endl;

	//Checking value against theory
	std::ofstream r_file("./Ameribond_r.txt");
	for (int s = 1; s <= 66; s++) {
		r_file << s << " , " <<
			crank_nicolson_AM_LINEAR(s, X, F, T,0.0019, sigma, R, kappa, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax, t0, cp) << " , " <<
			crank_nicolson_AM_LINEAR(s, X, F, T,0.0038, sigma, R, kappa, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax, t0, cp) << " , " <<
			crank_nicolson_AM_LINEAR(s, X, F, T,0.0057, sigma, R, kappa, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax, t0, cp) <<
			"\n";

	}
	cout << "AMERICAN BOND PART 2 DONE" << endl;
	*/
	auto start = high_resolution_clock::now();
	cout << crank_nicolson_AM_LINEAR(56.47, X, F, T, 0.0038, sigma, R, kappa, mu, C, alpha, beta, iMax = 320, jMax = 50, S_max = 56.47*6.7, tol, omega, iterMax, t0, cp) << endl;
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << duration.count() << endl;
}
