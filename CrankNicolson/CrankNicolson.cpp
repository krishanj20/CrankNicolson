#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;

/* Solution code for the Crank Nicolson Finite Difference
search for COURSEWORK EDIT for parts that needed to be altered for the coursework
 */
double CrankNicolson(double S0, double X, double F, double T, double r, double sigma,
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
		double theta = (1 + mu) * X * exp(-mu * i * dt);
		a[0] = 0.;
		b[0] = (-1/dt - kappa* theta / dS - r/2);
		c[0] = kappa* theta / dS;
		d[0] = -(1/dt - r/2)*vOld[1] - C*exp(-alpha*i*dt);

		for (int j = 1; j <= jMax - 1; j++)
		{
			//COURSEWORK EDIT
			//This is from geenral PDE, edited
			a[j] = (1/(4*dS*dS)) * sigma * sigma* pow(j*dS, 2*beta) - (1/(4*dS)) * kappa * ( theta - j * dS);
			b[j] = -1/dt - (1/(dS*dS)) *sigma * sigma * pow(j*dS,2*beta) - r/2;
			c[j] = (1 / (4 * dS * dS)) * sigma * sigma * pow(j * dS, 2 * beta) + (1 / (4 * dS)) * kappa * (theta - j * dS);
			d[j] = -a[j] * vOld[j - 1] - (b[j] + 2. / dt) * vOld[j] - c[j] * vOld[j + 1];
		}

		//COURSEWORK EDIT
		//THis is from calculating the general solution
		double A = R * exp((kappa + r) * (i*dt - T));
		double B = X * A - C / (alpha + r) * exp(-alpha * i * dt) + C / (alpha + r) * exp(-(alpha * r) * T) - X * R * exp(-r * T);
		a[jMax] = 0.;
		b[jMax] = 1.;
		c[jMax] = 0.;
		d[jMax] = A*jMax*dS + B;

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

int main() {
	// declare and initialise Black Scholes parameters - Currently looking at a solution we can get a definite answer for
	double T = 3., F = 56., R = 1., r = 0.0038, kappa = 0.0833333333,
		mu = 0.0073, X = 56.47, C = 0.106, alpha = 0.01, beta = 1., sigma = 3.73;
	// declare and initialise grid paramaters 
	int iMax = 10, jMax = 10;
	// declare and initialise local variables (ds,dt)
	double S_max = 5 * X;
	int iterMax = 10000;
	double tol = 1.e-8, omega = 1.;

	//Create graph of varying  S and optionvalue
	int length = 50;

	std::ofstream outFile1("./Varying_S_beta_1.txt");
	std::ofstream outFile2("./Varying_S_beta_425.txt");
	for (int j = 0; j <= length -1; j++) {
		//Case 1
		beta = 1.;
		sigma = 0.369;
		outFile1 << j * S_max / length << " , " << CrankNicolson(j * S_max / length, X, F, T, r, sigma, R, kappa, mu, C , alpha, beta, iMax, jMax,
			S_max, tol, omega, iterMax) << "\n";

		//Case 2
		beta = 0.425;
		sigma = 3.73;
		outFile2 << j * S_max / length << " , " << CrankNicolson(j * S_max / length, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, iMax, jMax,
			S_max, tol, omega, iterMax) << "\n";
	}
}
