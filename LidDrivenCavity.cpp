#include "LidDrivenCavity.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>

#define F77NAME(x) x##_
extern "C"
{
	void F77NAME(dgemv) (const char& trans, const int& m,
			const int& n, const double& alpha,
			const double* A, const int& lda,
			const double* x, const int& incx,
			const double& beta, double* y,
			const int& incy);
	void F77NAME(daxpy) (const int &n, const double &alpha, 
			const double *x, const int &incx,
			double *y, const int &incy);	
	void F77NAME(dscal) (const int &n, const double &alpha, 
		        double *x, const int &incx);
	void F77NAME(dgesv) (const int &n, const int &nrhs, const double *A,
			const int &lda, int *ipiv, double *B, 
			const int &ldb, int &info);
}

LidDrivenCavity::LidDrivenCavity()
{
}

LidDrivenCavity::~LidDrivenCavity()
{
	delete[] v;
	delete[] s;
	delete[] v_top;
	delete[] v_bot;
	delete[] v_left;
	delete[] v_right;
	delete[] s_top;
	delete[] s_bot;
	delete[] s_left;
	delete[] s_right;

}

void LidDrivenCavity::SetDomainSize(double xlen, double ylen)
{
	Lx = xlen;
	Ly = ylen;
}
void LidDrivenCavity::SetGridSize(int nx, int ny)
{
	/// Nx and Ny is the number of the interior points in x and y direction
	Nx = nx;
	Ny = ny;
}
void LidDrivenCavity::SetTimeStep(double deltat)
{
	dt = deltat;
}
void LidDrivenCavity::SetFinalTime(double finalt)
{
	T = finalt;
}
void LidDrivenCavity::SetReynoldsNumber(double re)
{
	Re = re;
}
void LidDrivenCavity::GridSpace()
{
	dx = Lx / (Nx + 1);
	dy = Ly / (Ny + 1);
}

/**
 * @brief Construct constant matrices for linear systems
 * @param A matrix for calculate v at time t, i.e. v = As + b
 * @param B matrix for calculate v at time t+dt, with only subdiagonal values
 * @param C matrix for calculate v at time t+dt, with other off-diagonal values
 */
void LidDrivenCavity::LinearMatrices()
{
	static int size = Nx * Ny;        ///< size of matrices
	A = new double[size * size]();      ///< only consider inner points
	B = new double[size * size]();      ///< only consider inner points
	C = new double[size * size]();      ///< only consider inner points
	double A_DiagVal = 2/(dx*dx) + 2/(dy*dy); ///< diagonal entries
	double A_SubDiagVal = -1/(dy*dy);     ///< entries above and below the diagonal 
	double A_OffDiagVal = -1/(dx*dx);     ///< entries in the off-diagonal blocks 
	double B_SubDiagVal = 1/(2*dy);
	double C_OffDiagVal = 1/(2*dx);
	/// Diagonal entries
	for (int i = 0; i < size; i++)
	{
		A[i*size + i] = A_DiagVal;
	}
	/// Sub-diagonal entries
	for (int i = 0; i < size-2; i++)
	{
		if ( (i+1) % Ny  != 0)
		{
			A[i*size + i+1] = A_SubDiagVal;
			A[(i+1)*size + i] = A_SubDiagVal;
			B[i*size + i+1] = B_SubDiagVal;
			B[(i+1)*size + i] = -B_SubDiagVal;
		}
	}
	/// Other off-diagonal entries
	for (int i = 0; i < size - Ny; i++)
	{
		A[i*size + (i+Ny)] = A_OffDiagVal;
		A[(i+Ny)*size + i] = A_OffDiagVal;
		C[i*size + (i+Ny)] = C_OffDiagVal;
		C[(i+Ny)*size + i] = -C_OffDiagVal;
	}
}

/**
 * @brief construct the boundary vector b in linear system y = Ax + b
 * @param b the vector b
 * @param matrix determine the matrix of the linear system (A, B or C)
 * @param x determine the solution vector (v or s)
 */
void LidDrivenCavity::BoundaryVector(double *b, char matrix, char x)
{	
	static int size = Nx * Ny;           ///< size of vector b 
	static double A_SubDiagVal = -1/(dy*dy);     ///< entries above and below the diagonal 
	static double A_OffDiagVal = -1/(dx*dx);     ///< entries in the off-diagonal blocks 
	static double B_SubDiagVal = 1/(2*dy);  
	static double C_OffDiagVal = 1/(2*dx);  
	switch(matrix)
	{
		case 'A':
		{	
			if (x == 's')
			{
				for (int i = 0; i < Ny; i++)
				{
					b[i] += s_left[i] * A_OffDiagVal;
					b[i+(Nx-1)*Ny] += s_right[i] * A_OffDiagVal;
				}
				for (int i = 0; i < Nx; i++)
				{
					b[i*Ny] += s_bot[i] * A_SubDiagVal;
					b[i*Ny + Ny-1] += s_top[i] * A_SubDiagVal;
				}
			}
			else if (x == 'v')
			{
				for (int i = 0; i < Ny; i++)
				{
					b[i] += v_left[i] * A_OffDiagVal;
					b[i+(Nx-1)*Ny] += v_right[i] * A_OffDiagVal;
				}
				for (int i = 0; i < Nx; i++)
				{
					b[i*Ny] += v_bot[i] * A_SubDiagVal;
					b[i*Ny + Ny-1] += v_top[i] * A_SubDiagVal;
				}
			} break;
		}
		case 'B':
		{
			if (x == 's')
			{
				for (int i = 0; i < Nx; i++)
				{
					b[i*Ny] += s_bot[i] * B_SubDiagVal;
					b[i*Ny + Ny-1] += s_top[i] * (-B_SubDiagVal);
				}
			}
			else if (x == 'v')
			{
				for (int i = 0; i < Nx; i++)
				{
					b[i*Ny] += v_bot[i] * B_SubDiagVal;
					b[i*Ny + Ny-1] += v_top[i] * (-B_SubDiagVal);
				}
			} break;
		}
		case 'C':
		{	
			if (x == 's')
			{
				for (int i = 0; i < Ny; i++)
				{
					b[i] += s_left[i] * C_OffDiagVal;
					b[i+(Nx-1)*Ny] += s_right[i] * C_OffDiagVal;
				}
			}
			else if (x == 'v')
			{
				for (int i = 0; i < Ny; i++)
				{
					b[i] += v_left[i] * C_OffDiagVal;
					b[i+(Nx-1)*Ny] += v_right[i] * C_OffDiagVal;
				}
			} break;
		}
		default: break;
	}	
}

void LidDrivenCavity::Initialise()
{
	/// Allocating interior memory, initialise to zero
	v = new double[Nx*Ny]();
	s = new double[Nx*Ny]();
	/// Allocate the boundary memory, initialise to zero
	v_top = new double[Nx]();
	v_bot = new double[Nx]();
	v_left = new double[Ny]();
	v_right = new double[Ny]();
	s_top = new double[Nx]();
	s_bot = new double[Nx]();
	s_left = new double[Ny]();
	s_right = new double[Ny]();
}

void LidDrivenCavity::Integrate()
{
}

/**
 * @brief Calculate Vorticity boundary conditions at time t
 *        Message passing and receive here
 */
void LidDrivenCavity::VorticityBCs()
{
	static double U = 1.0;
	/// Top boundary
	for (int i = 0; i < Nx; i++)
	{
		v_top[i] = (s_top[i] - s[i*Ny+(Ny-1)]) * 2
		      / (dy*dy) - 2*U/dy;
	}
	/// Bottom boundary
	for (int i = 0; i < Nx; i++)
	{
		v_bot[i] = (s_bot[i] - s[i*Ny]) * 2 / (dy*dy);
	}

	/// Left boundary
	for (int j = 0; j < Ny; j++)
	{
		 v_left[j] =(s_left[j] - s[j]) * 2 / (dx*dx); 
	}
	/// Right boundary
	for (int j = 0; j < Ny; j++)
	{
		 v_right[j] = (s_right[j] - s[(Nx-1)*Ny+j]) * 2
			 / (dx*dx);
	}
}

/**
 * @brief Calculate interior vorticity at time t, i.e. calculate v = As + b;
 */
void LidDrivenCavity::VorticityInterior()
{
	/// Construct RHS vector b using boundary conditions
	static int size = Nx * Ny;           ///< size of vector b 
	double *b =  new double[size]();    ///< allocate memory and initialise to zero
	BoundaryVector(b, 'A', 's');

	/// Calculate v = As + b;
	F77NAME(dgemv) ('T', size, size, 1.0, A, size, s, 1, 1.0, b, 1);
	memcpy(v, b, size * sizeof(double));
	delete[] b;	
}

/**
 * @brief Calculate interior vorticity at time t+dt, i.e. calculate v = v + f, where b includes three parts
 */
void LidDrivenCavity::VorticityUpdate()
{
	static int size =  Nx * Ny;   ///< size of vector f
	/// three parts of vector f
	double *b1 = new double[size]();   ///< viscosity term
	double *b2 = new double[size]();   ///< advection term
	double *b3 = new double[size]();   ///< also advection term
	
	/// First calculate b1-----------------------------------
	BoundaryVector(b1, 'A', 'v');
	F77NAME(dgemv) ('T', size, size, 1.0, A, size, v, 1, 1.0, b1, 1 );
	/// Calculate b2-----------------------------------------
	double *temp1 = new double[size](); ///< temp1 = Cs + b, b is the boundary term
	double *temp2 = new double[size](); ///< temp2 = Bv + b
	// calculate temp1
	BoundaryVector(temp1, 'C', 's');
	F77NAME(dgemv) ('T', size, size, 1.0, C, size, s, 1, 1.0, temp1, 1 );
	// calculate temp2
	BoundaryVector(temp2, 'B', 'v');
	F77NAME(dgemv) ('T', size, size, 1.0, B, size, v, 1, 1.0, temp2, 1 );
	// calculate b2 by multiply each vector element i.e. b2 = temp1 .* temp2
	for (int i = 0; i < size; i++)
	{
		b2[i] = temp1[i] * temp2[i];
	}
	// release memory for temp1 and temp2
	delete[] temp1;
	delete[] temp2;
	
	/// Calculate b3-----------------------------------------
	double *temp3 = new double[size](); ///< temp3 = Cv + b, b is the boundary term
	double *temp4 = new double[size](); ///< temp4 = Bs + b
	// calculate temp3
	BoundaryVector(temp3, 'C', 'v');
	F77NAME(dgemv) ('T', size, size, 1.0, C, size, v, 1, 1.0, temp3, 1 );
	// calculate temp4
	BoundaryVector(temp4, 'B', 's');
	F77NAME(dgemv) ('T', size, size, 1.0, B, size, s, 1, 1.0, temp4, 1 );
	// calculate b3 by multiply each vector element i.e. b3 = temp3 .* temp4
	for (int i = 0; i < size; i++)
	{
		b3[i] = temp3[i] * temp4[i];
	}
	// release memory for temp3 and temp4
	delete[] temp3;
	delete[] temp4;

	/// Update interior vorticity-------------------------------
	F77NAME(daxpy) (size, -dt/Re, b1, 1, v, 1);
	F77NAME(daxpy) (size, dt, b2, 1, v, 1);
	F77NAME(daxpy) (size, -dt, b3, 1, v, 1);

	/// release memory
	delete[] b1;
	delete[] b2;
	delete[] b3;
}

/**
 * @brief Possion solver to update stream function at time t+dt
 */
void LidDrivenCavity::PossionSolver()
{
	static int size = Nx * Ny;
	double *b =  new double[size]();    ///< allocate memory and initialise to zero
	BoundaryVector(b, 'A', 's');

	/// Solving As = v - b
	// Calculate  RHS =  v - b
	F77NAME(daxpy) (size, -1.0, v, 1, b, 1);
	F77NAME(dscal) (size, -1.0, b, 1);
	// Solving As = RHS (b)
	int *ipiv = new int[size];   ///< vector for pivots
	int info = 0;
	// A is symmetric so doesn't need to be transpose
	F77NAME(dgesv) (size, 1, A, size, ipiv, b, size, info);
	memcpy (s, b, size*sizeof(double));
	delete[] ipiv;
	delete[] b;	
}

/**
 * @brief Output the whole domain vorticity and stream function as matrix form
 */
void LidDrivenCavity::Output()
{
	/// here only output interior values of each process
	ofstream vOut("Vorticity.txt", ios::out | ios::trunc);
	ofstream sOut("StreamFunction.txt", ios::out | ios::trunc);
	vOut.precision(5);
	sOut.precision(5);

	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			vOut << setw(15) << v[i*Ny+j];
			sOut << setw(15) << s[i*Ny+j];
		}
		vOut << endl;
		sOut << endl;
	}
	vOut.close();
	sOut.close();
}
