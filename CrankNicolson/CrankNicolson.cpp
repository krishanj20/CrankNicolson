#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;

void sorSolve(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& rhs,
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
double crank_nicolson1(double S0, double X, double F, double T, double r, double sigma,
	double R, double kappa, double mu, double C, double alpha, double beta,int iMax, int jMax, double S_max, double tol, double omega, int iterMax)
{
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
		//COURSEWORK EDIT - Edited to equation 1
		vOld[j] = max(F, R*S[j]);
		vNew[j] = max(F, R*S[j]);
	}
	// start looping through time levels
	for (int i = iMax - 1; i >= 0; i--)
	{
		// declare vectors for matrix equations
		vector<double> a(jMax + 1), b(jMax + 1), c(jMax + 1), d(jMax + 1);
		// set up matrix equations a[j]=
		//COURSEWORK EDIT
		//This is boundary condition at S=0: edited
		double theta = (1 + mu) * X * exp(mu * i * dt);
		a[0] = 0.;
		b[0] = -1 / dt - kappa * theta / dS - r/2;
		c[0] = kappa * theta / dS;
		d[0] = -C * exp(-alpha * i * dt) - vOld[0] / dt + r*vOld[0]/2;


		for (int j = 1; j <= jMax - 1; j++)
		{
			//COURSEWORK EDIT
			//This is from geenral PDE, edited
			a[j] = sigma * sigma * pow(j * dS, 2 * beta) / (4 * dS * dS) - kappa * (theta - j * dS) / 4 * dS;
			b[j] = -1 / dt - sigma * sigma * pow(j * dS, 2 * beta) / 2 * dS * dS - r / 2;
			c[j] = sigma * sigma * pow(j * dS, 2 * beta) / (4 * dS * dS) + kappa * (theta - j * dS) / 4 * dS;
			d[j] = -a[j] * vOld[j - 1] - (b[j] + 2. / dt) * vOld[j] - c[j] * vOld[j + 1] - C * exp(-alpha * i * dt);
		}

		//COURSEWORK EDIT
		//THis is from calculating the general solution
		double A = R * exp((kappa + r) * (i * dt - T));
		//double B = X * A + C * exp(-alpha * i * dt) / (alpha + r) + X * R * exp(r * (i * dt - T)) - C * exp(-(alpha + r) * T + r * i * dt) / (alpha + r);
		double B = X * R * (1 - exp((kappa + r) * (i * dt - T))) + C * (exp(-alpha * i * dt) - exp(-alpha * T)) / (alpha + r);
		a[jMax] = 0.;
		b[jMax] = 1.;
		c[jMax] = 0.;
		d[jMax] = A * jMax * dS + B;

		// solve matrix equations with SOR
		int sor;
		for (sor = 0; sor < iterMax; sor++)
		{
			// SOR equations in here
			{
				double y = (d[0] - c[0] * vNew[1]) / b[0];
				vNew[0] = vNew[0] + omega * (y - vNew[0]);
			}
			for (int j = 1; j < jMax; j++)
			{
				double y = (d[j] - a[j] * vNew[j - 1] - c[j] * vNew[j + 1]) / b[j];
				vNew[j] = vNew[j] + omega * (y - vNew[j]);
			}
			{
				double y = (d[jMax] - a[jMax] * vNew[jMax - 1]) / b[jMax];
				vNew[jMax] = vNew[jMax] + omega * (y - vNew[jMax]);
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
			cout << "\n NOT CONVERGED \n";
		// set old=new 
		vOld = vNew;
	}
	// finish looping through time levels

	// output the estimated option price
	double optionValue;
	{
		int jStar = S0 / dS;
		double sum = 0.;
		sum += (S0 - S[jStar]) / dS * vNew[jStar + 1];
		sum += (S[jStar + 1] - S0) / dS * vNew[jStar];
		optionValue = sum;
	}
	return optionValue;

}

double crank_nicolson2(double S0, double X, double F, double T, double r, double sigma,
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
		sorSolve(a, b, c, d, vNew, iterMax, tol, omega, sor);
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
		sum += (S0 - S[jStar]) / dS * vNew[jStar + 1];
		sum += (S[jStar + 1] - S0) / dS * vNew[jStar];
		optionValue = sum;
	}
	return optionValue;
}


int main() {
	// declare and initialise Black Scholes parameters - Currently looking at a solution we can get a definite answer for
	double T = 3., F = 56., R = 1., r = 0.0038, kappa = 0.08333333,
		mu = 0.0073, X = 56.47, C = 0.106, alpha = 0.01, beta = 1., sigma = 3.73;
	// declare and initialise grid paramaters 
	int iMax = 100, jMax = 25;
	// declare and initialise local variables (ds,dt)
	double S_max = 10 * X;
	int iterMax = 5000000;
	double tol = 1.e-7, omega = 1.;

	//Checking value against theory
	std::ofstream analytical ("./analytical.txt");
	for (int s = 1; s <= 300; s++) {
		analytical << s << " , " << crank_nicolson2(s, X, F, T, r, 3.73, R, 0, mu, C, alpha, 1., iMax, jMax, S_max, tol, omega, iterMax) << "\n";
	}
	
	//Create graph of varying  S and optionvalue
	int length = 50;
	double S_maxi = 4 * X;
	std::ofstream V_function_S("./V_function_S.txt");
	std::ofstream V_function_sigma("./V_function_sigma.txt");
	
	for (int j = 1; j <= length -1; j++) {
		//Plottting V(S,t) as a function of S (t=0) for two cases
		beta = 1.;
		sigma = 0.369;
		//Puts value of S, and value of stock for different parameters into csv file
		V_function_S << j * S_maxi / length << " , " << crank_nicolson2(j * S_maxi / length, X, F, T, r, 0.369, R, kappa, mu, C , alpha, 1., iMax, jMax,
			S_max, tol, omega, iterMax) << " , " << crank_nicolson2(j * S_maxi / length, X, F, T, r, 3.73, R, kappa, mu, C, alpha, 0.425, iMax, 35,
				S_max, tol, omega, iterMax) << "\n";
		
		//"You may wish to explore different values of beta and sigma?" Do research before implementing this. But preliminary (BETA =1)
		//V_function_sigma << j * S_maxi / length << " , " <<
		//	CrankNicolson(j * S_maxi / length, X, F, T, r, 1, R, kappa, mu, C, alpha, 1., iMax, jMax,S_max, tol, omega, iterMax) << " , " <<
		//	CrankNicolson(j * S_maxi / length, X, F, T, r, 3., R, kappa, mu, C, alpha, 1., iMax, jMax, S_max, tol, omega, iterMax) << " , " <<
		//	CrankNicolson(j * S_maxi / length, X, F, T, r, 5., R, kappa, mu, C, alpha, 1., iMax, jMax, S_max, tol, omega, iterMax) << " , " <<
		//	CrankNicolson(j * S_maxi / length, X, F, T, r, 7., R, kappa, mu, C, alpha, 1., iMax, jMax, S_max, tol, omega, iterMax) << " , " <<
		//	CrankNicolson(j * S_maxi / length, X, F, T, r, 9., R, kappa, mu, C, alpha, 1., iMax, jMax, S_max, tol, omega, iterMax) << " , " <<
		//	"\n";
			
	}
	cout << "DONE V FUNC S" << endl;
	
	//Assume now,
	double S0 = 56.47;
	beta = 0.425;
	sigma = 3.73;

	//Explore effect of Smax
	//Look at given imax and jmax, then increase Smax
	std::ofstream V_function_sMax("./V_function_sMax.txt");
	for (int i = 2; i < 30; i++) {
		V_function_sMax << i * X << " , " << crank_nicolson2(X, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, 100, 25,
			i*X, tol, omega, iterMax) << "\n";
	}
	cout << "DONE V FUNC S_MAX" << endl;
	//Explore effect of imax
	//Look at given imax and jmax, then increase Smax
	std::ofstream V_function_iMax("./V_function_iMax.txt");
	for (int i = 0; i < 200; i++) {
		V_function_iMax << i << " , " <<
			crank_nicolson2(S0, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, i, 15,S_max, tol, omega, iterMax) << " , " <<
			crank_nicolson2(S0, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, i, 25, S_max, tol, omega, iterMax) << " , " <<
			crank_nicolson2(S0, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, i, 50, S_max, tol, omega, iterMax) << " , " <<
			"\n";
	}
	cout << "DONE V FUNC iMax" << endl;
	//Explore effect of jmax
	//Look at given imax and jmax, then increase Smax
	std::ofstream V_function_jMax("./V_function_jMax.txt");
	for (int i = 1; i < 200; i++) {
		V_function_jMax << i<< " , " <<
			crank_nicolson2(S0, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, 20, i, S_max, tol, omega, iterMax) << " , " <<
			crank_nicolson2(S0, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, 35, i, S_max, tol, omega, iterMax) << " , " <<
			crank_nicolson2(S0, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, 50, i, S_max, tol, omega, iterMax) << " , " <<
			"\n";
	}
	cout << "DONE V FUNC jMax" << endl;
}
