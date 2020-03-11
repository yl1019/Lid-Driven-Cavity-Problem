#include"PoissonSolver.h"
#include<fstream>
#include<iomanip>
#include<iostream>
#include<cmath>

using namespace std;

int main()
{
	double Lx = 1;
	double Ly = 1;
	int Nx = 100;
	int Ny = 50;
	Nx -= 2;
	Ny -= 2;
	double dx = Lx/(Nx+1);
	double dy = Ly/(Ny+1);
	int size = Nx*Ny;

	double *top = new double[Nx]();
	double *bot = new double[Nx]();
	double *left = new double[Ny]();
	double *right = new double[Ny]();

	for (int i = 0; i < Nx; i++)
	{
		top[i] = sin(i*dx);
		bot[i] = 0;
	}
	for (int j = 0; j < Ny; j++)
	{
		left[j]= 2;
		right[j] = 5;
	}

	double *x = new double[size]();
	double *f = new double[size]();
	
	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			x[i*Ny+ j] = 0.0;
			
			f[i*Ny + j] = 10*sin(i*dx)*cos(j*dy);
		}
	}

	PoissonSolver test(Nx, Ny, dx, dy);
	test.SetBoundary(top, left, bot, right);
	test.Solve_Chol(x, f);
	//test.Solve_Conj(x, f);

	
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
