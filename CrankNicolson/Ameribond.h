#pragma once
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <chrono> 
#include <iomanip>
using namespace std::chrono;
using namespace std;

std::vector<double> thomasSolve(const std::vector<double>& a, const std::vector<double>& b_, const std::vector<double>& c, std::vector<double>& d)
{
	int n = a.size();
	std::vector<double> b(n), temp(n);
	// initial first value of b
	b[0] = b_[0];
	for (int j = 1; j < n; j++)
	{
		b[j] = b_[j] - c[j - 1] * a[j] / b[j - 1];
		d[j] = d[j] - d[j - 1] * a[j] / b[j - 1];
	}
	// calculate solution
	temp[n - 1] = d[n - 1] / b[n - 1];
	for (int j = n - 2; j >= 0; j--)
		temp[j] = (d[j] - c[j] * temp[j + 1]) / b[j];
	return temp;
}

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
			if (i * dt < t0) {
				x[0] = max(0., min(cp, x[0] + omega * (y - x[0])));
			}
			else {
				x[0] = max(x[0] + omega * (y - x[0]), 0.);
			}
			error += (y - x[0]) * (y - x[0]);
		}
		for (int j = 1; j < n; j++)
		{
			double y = (rhs[j] - a[j] * x[j - 1] - c[j] * x[j + 1]) / b[j];
			if (i * dt < t0) {
				x[j] = max(min(x[j] + omega * (y - x[j]),cp), j * dS);
			}
			else {
				x[j] = max(x[j] + omega * (y - x[j]), j * dS);
			}
			error += (y - x[j]) * (y - x[j]);
		}
		{
			double y = (rhs[n] - a[n] * x[n - 1]) / b[n];
			if (i * dt < t0) {
				x[n] = max(min(x[n] + omega * (y - x[n]),cp), n * dS);
			}
			else {
				x[n] = max(x[n] + omega * (y - x[n]), n * dS);
			}
			error += (y - x[n]) * (y - x[n]);
		}
		// make an exit condition when solution found
		if (error < tol) 
			break;
	}
	if (sor >= iterMax)
	{
		std::cout << " Error NOT converging within required iterations\n";
	}
}

/* Solution code for the Crank Nicolson Finite Difference
search for COURSEWORK EDIT for parts that needed to be altered for the coursework
 */
double crank_nicolson_AM_LINEAR(double S0, double X, double F, double T, double r, double sigma,
	double R, double kappa, double mu, double C, double alpha, double beta, int iMax, int jMax, int S_max, double tol, double omega, int iterMax, double cp, double t0)
{
	// declare and initialise local variables (ds,dt)
	cp = 67;
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
		}
		double A = R * exp((kappa + r) * (i * dt - T));
		double B = -X * A + C * exp(-alpha * i * dt) / (alpha + r) + X * R * exp(r * (i * dt - T)) - C * exp(-(alpha + r) * T + r * i * dt) / (alpha + r);
		a[jMax] = 0;
		b[jMax] = 1;
		c[jMax] = 0;
		d[jMax] = jMax * dS * A + B;
		// solve matrix equations with SOR
		int sor;
		for (sor = 0; sor < iterMax; sor++)
		{
			double error = 0.;
			// implement sor in here
			{
				double y = (d[0] - c[0] * vNew[1]) / b[0];
				y = vNew[0] + omega * (y - vNew[0]);
				if (i * dt < t0)
				{
					y = max(0., min(cp, y));
				}
				else
				{
					y = std::max(y, R * S[0]);
				}
				error += (y - vNew[0]) * (y - vNew[0]);
				vNew[0] = y;
			}
			for (int j = 1; j < jMax; j++)
			{
				double y = (d[j] - a[j] * vNew[j - 1] - c[j] * vNew[j + 1]) / b[j];
				y = vNew[j] + omega * (y - vNew[j]);
				if (i * dt < t0)
				{
					y = max(min(y, cp), j * dS);
				}
				else
				{
					y = std::max(y, R * j * dS);
				}
				error += (y - vNew[j]) * (y - vNew[j]);
				vNew[j] = y;
			}
			{
				double y = (d[jMax] - a[jMax] * vNew[jMax - 1]) / b[jMax];
				y = vNew[jMax] + omega * (y - vNew[jMax]);
				if (i * dt < t0)
				{
					y = max(min(y, cp), jMax * dS);
				}
				else
				{
					y = std::max(y, R * jMax * dS);
				}
				error += (y - vNew[jMax]) * (y - vNew[jMax]);
				vNew[jMax] = y;
			}
			// make an exit condition when solution found
			if (error < tol)
				break;
		}
		if (sor >= iterMax)
		{
			std::cout << " Error NOT converging within required iterations\n";
			std::cout.flush();
			throw;
		}

		if (sor == iterMax)
			return -1;

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

/* This code seems to run much faster when in a separate solution without header files, maybe just copy and paste this
 */
double crank_nicolson_AM_FAST(double S0, double X, double F, double T, double r, double sigma,
	double R, double kappa, double mu, double C, double alpha, double beta, int iMax, int jMax, int S_max, double tol, double omega, int iterMax, int& sorCount, double t0)
{
	// declare and initialise local variables (ds,dt)
	double cp = 67.;

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
		//if (i * dt < t0) { dt = t0 / iMax; }
		//if (i * dt >= t0) { dt = (T - t0) / iMax; }
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
		}
		double A = R * exp((kappa + r) * (i * dt - T));
		double B = -X * A + C * exp(-alpha * i * dt) / (alpha + r) + X * R * exp(r * (i * dt - T)) - C * exp(-(alpha + r) * T + r * i * dt) / (alpha + r);
		a[jMax] = 0;
		b[jMax] = 1;
		c[jMax] = 0;
		d[jMax] = jMax * dS * A + B;
		// solve matrix equations with SOR
		int sor;
		for (sor = 0; sor < iterMax; sor++)
		{
			double error = 0.;
			// implement sor in here
			{
				double y = (d[0] - c[0] * vNew[1]) / b[0];
				y = vNew[0] + omega * (y - vNew[0]);
				if (i * dt < t0)
				{
					y = max(0., min(cp, y));
				}
				else
				{
					y = std::max(y, R * S[0]);
				}
				error += (y - vNew[0]) * (y - vNew[0]);
				vNew[0] = y;
			}
			for (int j = 1; j < jMax; j++)
			{
				double y = (d[j] - a[j] * vNew[j - 1] - c[j] * vNew[j + 1]) / b[j];
				y = vNew[j] + omega * (y - vNew[j]);
				if (i * dt < t0)
				{
					y = max(min(y, cp), j * dS);
				}
				else
				{
					y = std::max(y, R * j * dS);
				}
				error += (y - vNew[j]) * (y - vNew[j]);
				vNew[j] = y;
			}
			{
				double y = (d[jMax] - a[jMax] * vNew[jMax - 1]) / b[jMax];
				y = vNew[jMax] + omega * (y - vNew[jMax]);
				if (i * dt < t0)
				{
					y = max(min(y, cp), jMax * dS);
				}
				else
				{
					y = std::max(y, R * jMax * dS);
				}
				error += (y - vNew[jMax]) * (y - vNew[jMax]);
				vNew[jMax] = y;
			}
			// make an exit condition when solution found
			if (error < tol)
				break;
		}
		if (sor >= iterMax)
		{
			std::cout << " Error NOT converging within required iterations\n";
			std::cout.flush();
			throw;
		}

		if (sorCount == iterMax)
			return -1;

		// set old=new
		vOld = vNew;
	}
	// finish looping through time levels

	// output the estimated option price
	double optionValue;

	int jStar = S0 / dS;
	double sum = 0.;
	sum += (S0 - S[jStar]) / (dS)* vNew[jStar + 1];
	sum += (S[jStar + 1] - S0) / (dS)* vNew[jStar];
	optionValue = sum;
	//optionValue = lagrangeInterpolation(vNew, S, S0, vNew.size());

	return optionValue;
}

/* Template code for the Crank Nicolson Finite Difference
 */
double crank_nicolson_penalty(double S0, double X, double F, double T, double r, double sigma,
	double R, double kappa, double mu, double C, double alpha, double beta, int iMax, int jMax, int S_max, double tol, double omega, int iterMax, int& sorCount, double t0)
{
	// declare and initialise local variables (ds,dt)
	double cp = 67.;
	//dS calculated as before
	double dS = S_max / jMax;
	//What is f used for?
	double f = (T - t0) / T;
	//STILL T/iMax
	double dt = (T - t0) / (iMax * f);
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
	for (int i = iMax; i >= 0; i--)
	{
		//If you are before t-, change increments to
		if (i * dt < t0)
		{
			dt = t0 / (iMax * (1 - f));
		}

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
		}
		double A = R * exp((kappa + r) * (i * dt - T));
		double B = -X * A + C * exp(-alpha * i * dt) / (alpha + r) + X * R * exp(r * (i * dt - T)) - C * exp(-(alpha + r) * T + r * i * dt) / (alpha + r);
		a[jMax] = 0;
		b[jMax] = 1;
		c[jMax] = 0;
		d[jMax] = jMax * dS * A + B;
		double penalty = 1.e8;
		int q;
		for (q = 0; q < 100000; q++)
		{
			vector<double> bHat(b), dHat(d);
			for (int j = 1; j < jMax; j++)
			{
				if (i * dt < t0)
				{
					if (vNew[j] > max(R * S[j], cp))
					{
						bHat[j] = b[j] - penalty;
						dHat[j] = d[j] - penalty * max(R * S[j], cp);

					}
				}
				else
				{
					// turn on penalty if V < RS
					if (vNew[j] < R * S[j])
					{
						bHat[j] = b[j] - penalty;
						dHat[j] = d[j] - penalty * R * S[j];
					}
				}
			}
			// solve matrix equations with SOR
			vector<double> y = thomasSolve(a, bHat, c, dHat);
			// calculate difference from last time
			double error = 0.;
			for (int j = 0; j <= jMax; j++)
				error += fabs(vNew[j] - y[j]);
			vNew = y;
			if (error < 1.e-8)
				break;
		}
		if (q == 100000)
		{
			std::cout << " Error NOT converging within required iterations\n";
			std::cout.flush();
			throw;
		}

		// set old=new
		vOld = vNew;
	}
	// finish looping through time levels

	// output the estimated option price
	//Why 4?
	return lagrangeInterpolation(vNew, S, S0, 4);
}


//ALlows you to extract information in csv file for efficiency of penalty method
void GetPenaltyEfficiency()
{
	// declare and initialise Black Scholes parameters - Currently looking at a solution we can get a definite answer for
	double T = 3., F = 56., R = 1., r = 0.0038, kappa = 0.0833333333,
		mu = 0.0073, X = 56.47, C = 0.106, alpha = 0.01, beta = 0.425, sigma = 3.73, tol = 1.e-7, omega = 1., S_max = 10 * X;
	double t0 = 1.2448;
	//

	int iterMax = 100000;
	//Create graph of varying S0 and beta and bond
	int length = 300;
	int sorCount;
	double S0 = X;
	std::ofstream outFile("ameribond_eff.txt");
	double oldResult = 0, oldDiff = 0;
	double S = X;
	int iMax = 100;
	int jMax = 100;
	S_max = 200 * X;
	//cout << crank_nicolson1(X, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax, sorCount, t0) << endl;

	for (int n = 100; n <= 10000; n *= 2)
	{
		//Set aprameters for iteration
		iMax = n;
		jMax = n;
		S_max = n / 20 * X;

		auto t1 = std::chrono::high_resolution_clock::now();
		double result = crank_nicolson_penalty(X, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax, sorCount, t0);
		double diff = result - oldResult;
		auto t2 = std::chrono::high_resolution_clock::now();
		auto time_taken = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

		outFile << n << "," << setprecision(10) << result << "," << time_taken << "," << setprecision(3) << oldDiff / diff << "\n";
		cout << "RESULT : " << result << "  TIME: " << time_taken << "," << setprecision(3) << oldDiff / diff << "\n";

		oldDiff = diff;
		oldResult = result;
	}

}

void getAmeribondEfficiency() {
	// declare and initialise Black Scholes parameters - Currently looking at a solution we can get a definite answer for
	double T = 3., F = 56., R = 1., r = 0.0038, kappa = 0.0833333333,
		mu = 0.0073, X = 56.47, C = 0.106, alpha = 0.01, beta = 0.425, sigma = 3.73, tol = 1.e-7, omega = 1., S_max = 20 * X;
	//
	int iterMax = 100000, iMax = 40, jMax = 25;
	beta = 0.425;
	sigma = 3.73;
	double S0 = X;
	double t0 = 1.2448, cp = 67;

	std::ofstream outFile5("american_varying_smax.txt");
	tol = 1.e-7;
	for (int i = 1; i <= 1; i++)
	{
		double jMax = 300;
		double S = X;
		int sorCount;
		auto t1 = std::chrono::high_resolution_clock::now();
		double result = crank_nicolson_AM_FAST(S, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, iMax = 1000, jMax = 300 * i, S_max = S * cp, tol = 1e-8, omega, iterMax, sorCount, t0 = 0);
		auto t2 = std::chrono::high_resolution_clock::now();
		auto time_taken =
			std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1)
			.count();
		cout << result << "," << time_taken << "\n";
		t1 = std::chrono::high_resolution_clock::now();
		result = crank_nicolson_AM_LINEAR(S, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, iMax = 1000, jMax = 300 * i, S_max = S * cp, tol = 1e-8, omega, iterMax, sorCount, t0 = 0);
		t2 = std::chrono::high_resolution_clock::now();
		time_taken =
			std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1)
			.count();
		cout << result << "," << time_taken << "\n";
	}
	outFile5.close();

}


//Get data for ameribond. Gets information for first part of American Bond graphs
void getAmeribondData() {
	// declare and initialise Black Scholes parameters - Currently looking at a solution we can get a definite answer for
	double T = 3., F = 56., R = 1., r = 0.0038, kappa = 0.0833333333,
		mu = 0.0073, X = 56.47, C = 0.106, alpha = 0.01, beta = 0.425, sigma = 3.73, tol = 1.e-7, omega = 1., S_max = 15 * X;
	//
	int iterMax = 10000, iMax = 700, jMax = 700;
	beta = 0.425;
	sigma = 3.73;
	double S0 = X;
	double t0 = 1.2448, cp = 67;
	int sor;
	//Checking value against theory

	std::ofstream analytical("./Ameribond1.txt");
	for (int s = 1; s <= 100; s++) {
		analytical << s << " , " << crank_nicolson_penalty(s, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax,sor, t0) << "\n";
	}
	cout << "AMERICAN BOND PART 1 DONE" << endl;

	//Checking how the value changes with different values of r
	std::ofstream r_file("./Ameribond_r.txt");
	for (int s = 1; s <= 67; s++) {
		r_file << s << " , " <<
			crank_nicolson_penalty(s, X, F, T, 0.0019, sigma, R, kappa, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax,sor, t0) << " , " <<
			crank_nicolson_penalty(s, X, F, T, 0.0038, sigma, R, kappa, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax,sor, t0) << " , " <<
			crank_nicolson_penalty(s, X, F, T, 0.0057, sigma, R, kappa, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax,sor, t0) <<
			"\n";

	}
	
	cout << "AMERICAN BOND PART 2 DONE" << endl;

}