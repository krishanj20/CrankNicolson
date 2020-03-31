/*
 *  A simple example of the sor method for you to fill in in class.
 *  The ideas here are easily extendable to Finite Difference methods
 */
#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

int main()
{
	// declare vector to store matrix A
	vector<double> a = { 0.,-1,-1,-1 };
	vector<double> b = { 1.,2,2,1 };
	vector<double> c = { 0.,-1,-1,0. };
	// declare vector for d and solution x 
	vector<double> d = { 1.,0.25,0.5,0. };
	vector<double> x;
	// first try
	// reset initial guess for x
	x = { 0.,0.,0.,0. };

	// MATRIX SOLVE
	// sor variables
	int iterMax = 100;
	double error, tol = 1.e-8, omega = 1;
	// sor loop
	for (int sor = 0; sor < iterMax; sor++)
	{
		// implement sor in here
		// x[0]=...

		// output guess at (sor+1)th iteration
		cout << sor + 1 << " x = { ";
		for (auto xi : x)cout << xi << " ";
		cout << "} \n";
		// make an exit condition when solution found
		// ...

	}
}