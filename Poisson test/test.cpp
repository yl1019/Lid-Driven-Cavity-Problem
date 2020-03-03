#include"PoissonSolver.h"
#include<fstream>
#include<iomanip>
#include<iostream>

using namespace std;

int main()
{
	double Lx = 1;
	double Ly = 1;
	int Nx = 10;
	int Ny = 10;
	Nx -= 2;
	Ny -= 2;
	double dx = Lx/(Nx+1);
	double dy = Ly/(Ny+1);
	int size = Nx*Ny;

	double *top = new double[Nx]();
	double *bot = new double[Nx]();
	double *left = new double[Ny]();
	double *right = new double[Ny]();
	
	double *x = new double[size]();
	double *f = new double[size]();
	
	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			x[i*Ny+ j] = 5*i+j;
			// f = 2x + 3y
			f[i*Ny + j] = 2*i*dx + 3*j*dy;
		}
	}

	PoissonSolver test(Nx, Ny, dx, dy);
	test.SetBoundary(top, left, bot, right);
	test.Solve(x, f);

	
	// output
	ofstream out("Output.txt", ios::out|ios::trunc);

	out.precision(5);
	
	for (int i = 0; i < size; i++)
	{
		out << setw(15) << x[i];
	}
	out.close();

	delete[] top;
	delete[] bot;
	delete[] left;
	delete[] right;
	delete[] x;
	delete[] f;
	
	return 0;
}
