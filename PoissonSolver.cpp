#include<cmath>
#include<cstring>
#include<iostream>
using namespace std;

#include "PoissonSolver.h"
#include "cblas.h"
/**
 * @brief Constructor
 * @param Nx	number of grid points in x direction
 * @param Ny	number of grid points in y direction
 * @param dx	element size in x direction
 * @param dy	element size in y direction
 */
PoissonSolver::PoissonSolver(const int &Nx, const int &Ny, const double &dx, const double &dy)
{
	/// accept the input parameter
	this->Nx = Nx;
	this->Ny = Ny;
	this->dx = dx;
	this->dy = dy;
	/// build matrix A of the system in banded and packed storage
	lda = Ny + 1;
	size = Nx * Ny;
	A = new double[size * lda];
  	memset(A, 0.0, size*lda*sizeof(double));	
	A_DiagVal = 2/(dx*dx) + 2/(dy*dy);	///< diagonal entries
	A_SubDiagVal = -1/(dy*dy);	///< entries above and below the diagonal 
	A_OffDiagVal = -1/(dx*dx);	///< entries in the off-diagonal block	
	for (int j = 0; j < size; j++)
	{
		A[lda-1 + j*lda] = A_DiagVal;
		if (j % Ny != 0)   A[lda-2 + j*lda] = A_SubDiagVal;
		if (j >= Ny)       A[j * lda] = A_OffDiagVal;
	}	
}

PoissonSolver::~PoissonSolver()
{
	delete[] b;
	delete[] A;
}

/**
 * @brief Accept boundary conditions from 4 direction and construct boundary vector b
 * @param top	top boundary vector
 * @param left	left boundary vector
 * @param bot	bottom boundary vector
 * @param right	right boundary vector
 */
void PoissonSolver::SetBoundary(const double *top, const double *left, const double *bot, const double *right)
{	
	b = new double[size];
	memset(b, 0.0, size*sizeof(double));
	for (int i = 0; i < Ny; i++)
	{
		b[i] += left[i] * A_OffDiagVal;
		b[i+(Nx-1)*Ny] += right[i] * A_OffDiagVal;
	}
	for (int i = 0; i < Nx; i++)
	{
		b[i*Ny] += bot[i] * A_SubDiagVal;
		b[i*Ny + Ny-1] += top[i] * A_SubDiagVal;
	}
}


/**
 * @brief Solve Poisson equation Ax = f - b using Conjugate Gradient Method, store solution in x
 * @param x	initial guess of the solution vector
 * @param f	Source vector 
 */
void PoissonSolver::Solve(double *x, const double *f)
{	
	double *r = new double[size];
	double *p = new double[size];
	double *t = new double[size]; ///< for temporary use
	double alpha;
	double beta;
	double eps;
	double tol = 0.00001;	///< tolerance for termination
	cblas_dcopy(size, f, 1, r, 1);	///< r0 = f
	cblas_daxpy(size, -1.0, b, 1, r, 1);	///< r0 = r0 - b
	cblas_dsbmv(CblasColMajor, CblasUpper, size, Ny,-1.0, A, lda, x, 1, 1.0, r, 1);	///< r0 = b - Ax0
	cblas_dcopy(size, r, 1, p, 1);	///< p0 = r0
	int k = 0;
	do
	{
		cblas_dsbmv(CblasColMajor, CblasUpper, size, Ny, 1.0, A, lda, p, 1, 0.0, t, 1);	///< t = A*p_k
		alpha = cblas_ddot(size, t, 1, p, 1);	///< alpha = p_k^T * A * p_k
		beta = cblas_ddot(size, r, 1, r, 1);	///< beta = r_k^T * r_k 
		alpha = beta / alpha;	///< alpha_k = beta / alpha
		cblas_daxpy(size, alpha, p, 1, x, 1); 	///< x_{k+1} = x_k + alpha_k * p_k
		cblas_daxpy(size, -alpha, t, 1, r, 1); 	///< r_{k+1} = r_k - alpha_k * t
		/// Take norm2 as the criteria of iteration
		eps = cblas_dnrm2(size, r, 1);
		if (eps < tol*tol) break;
		beta = cblas_ddot(size, r, 1, r, 1) / beta;	///< beta_k = r_{k+1}^T * r_{k+1} / dbeta
		cblas_dcopy(size, r, 1, t, 1);	///< t = r_{k+1}, avoid overwritting r 
		cblas_daxpy(size, beta, p, 1, t, 1);	///< t = t + beta_k * p_k
		cblas_dcopy(size, t, 1, p, 1);	///< p_{k+1} = t;
		
		k++;
	}while(k < 500);	///< maximum iteration time
	
	delete[] r;
	delete[] p;
	delete[] t;
}


