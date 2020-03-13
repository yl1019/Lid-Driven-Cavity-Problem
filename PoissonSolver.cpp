#include<cmath>
#include<cstring>
#include<iostream>
using namespace std;

#include "PoissonSolver.h"
#include "cblas.h"
<<<<<<< HEAD
#include "mpi.h"

#define F77NAME(x) x##_
extern "C"
{
	void F77NAME(dpbtrf) (const char &UPLO, const int &N, const int &KD, const double *AB,
		       const int &LDAB,	int &INFO);
	void F77NAME(dpbtrs) (const char &UPLO, const int &N, const int &KD, const int &NRHS,
		       	const double *AB, const int &LDAB, double *B, const int &LDB, int &INFO);
}

=======
>>>>>>> parent of 940f1f8... (1) Parallel program works fine, but really slow.
/**
 * @brief Constructor
 * @param Nx	number of grid points in x direction
 * @param Ny	number of grid points in y direction
 * @param dx	element size in x direction
 * @param dy	element size in y direction
 */
PoissonSolver::PoissonSolver(MPI_Comm mygrid, int *neighbor, int rank, const int &Nx,
	       	const int &Ny, const double &dx, const double &dy)
{
	/// accept the input parameter
	this->mygrid =mygrid;
	memcpy(this->neighbor, neighbor, 4*sizeof(double));
	this->rank = rank;
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
<<<<<<< HEAD
	}
	info = 1;
	x_top = new double[Nx]();
	x_bot = new double[Nx]();
	x_left = new double[Ny]();
	x_right = new double[Ny]();
	flag = 0;
	FLAG = 0;
=======
	}	
>>>>>>> parent of 940f1f8... (1) Parallel program works fine, but really slow.
}

PoissonSolver::~PoissonSolver()
{
	delete[] b;
	delete[] A;
	delete[] x_left;
	delete[] x_right;
	delete[] x_top;
	delete[] x_bot;
}

/**
<<<<<<< HEAD
 * @brief Cholesky Factorization of matrix A, to avoid doing this multiple times
 */
void PoissonSolver::CholeskyFactor()
{
	F77NAME(dpbtrf) ('U', size , Ny, A, lda, info);
}

/**
 * @brief Send and receive stream function to neighbor domains
 */
void PoissonSolver::SendAndRecv(double *x)
{
	/** Left to right */
	MPI_Sendrecv(x + (Nx-1)*Ny, Ny, MPI_DOUBLE, neighbor[3], 3, x_left, Ny,
				MPI_DOUBLE, neighbor[1], 3, mygrid, MPI_STATUS_IGNORE);

	/** Right to left */
	MPI_Sendrecv(x, Ny, MPI_DOUBLE, neighbor[1], 1, x_right, Ny,
			MPI_DOUBLE, neighbor[3], 1, mygrid, MPI_STATUS_IGNORE);

	/** Top to bottom */
	// extract temporary bottom boundary vector
	double temp_bot[Nx];
	cblas_dcopy(Nx, x, Ny, temp_bot, 1);
	MPI_Sendrecv(temp_bot, Nx, MPI_DOUBLE, neighbor[2], 2, x_top, Nx,
			MPI_DOUBLE, neighbor[0], 2, mygrid, MPI_STATUS_IGNORE);

	/** Bottom to top */
	// extract temporary top boundary vector
	double temp_top[Nx];
	cblas_dcopy(Nx, x + Ny-1, Ny, temp_top, 1);
	MPI_Sendrecv(temp_top, Nx, MPI_DOUBLE, neighbor[0], 0, x_bot, Nx,
			MPI_DOUBLE, neighbor[2], 0, mygrid, MPI_STATUS_IGNORE);

}

/**
 * @brief Setter for accepting boundary conditions from 4 direction 
=======
 * @brief Accept boundary conditions from 4 direction and construct boundary vector b
>>>>>>> parent of 940f1f8... (1) Parallel program works fine, but really slow.
 * @param top	top boundary vector
 * @param left	left boundary vector
 * @param bot	bottom boundary vector
 * @param right	right boundary vector
 */
void PoissonSolver::SetBoundary(const double *top, const double *left, const double *bot, const double *right)
{	
	cblas_dcopy(Nx, top, 1, x_top, 1);
	cblas_dcopy(Nx, bot, 1, x_bot, 1);
	cblas_dcopy(Ny, left, 1, x_left, 1);
	cblas_dcopy(Ny, right, 1, x_right, 1);
}

/**
 * @brief Build boundary vector b
 */   
void PoissonSolver::BoundaryVector()
{
	b = new double[size];
	memset(b, 0.0, size*sizeof(double));
	for (int i = 0; i < Ny; i++)
	{
		b[i] += x_left[i] * A_OffDiagVal;
		b[i+(Nx-1)*Ny] += x_right[i] * A_OffDiagVal;
	}
	for (int i = 0; i < Nx; i++)
	{
		b[i*Ny] += x_bot[i] * A_SubDiagVal;
		b[i*Ny + Ny-1] += x_top[i] * A_SubDiagVal;
	}
}


/**
 * @brief Solve Poisson equation Ax = f - b using Conjugate Gradient Method, store solution in x
 * @param x	solution vector
 * @param f	Source vector 
 */
void PoissonSolver::Solve(double *x, const double *f)
{
	if (cblas_dnrm2(size, f, 1) == 0 && cblas_dnrm2(size, f, 1) == 0)
	{
		memset(x, 0.0, size*sizeof(double));
		return;
	}	
	double *r = new double[size];
	double *p = new double[size];
	double *t = new double[size]; ///< for temporary use
	double alpha;
	double beta;
	double tol = 0.0001;	///< tolerance for termination
	int k = 0;

	do
	{
		BoundaryVector();
		cblas_dcopy(size, f, 1, r, 1);	///< r0 = f
		cblas_daxpy(size, -1.0, b, 1, r, 1);	///< r0 = r0 - b
		cblas_dsbmv(CblasColMajor, CblasUpper, size, Ny,-1.0, A, lda, x, 1, 1.0, r, 1);	///< r0 = b - Ax0
		cblas_dcopy(size, r, 1, p, 1);	///< p0 = r0

		cblas_dsbmv(CblasColMajor, CblasUpper, size, Ny, 1.0, A, lda, p, 1, 0.0, t, 1);	///< t = A*p_k
		alpha = cblas_ddot(size, t, 1, p, 1);	///< alpha = p_k^T * A * p_k
		beta = cblas_ddot(size, r, 1, r, 1);	///< beta = r_k^T * r_k 
		alpha = beta / alpha;	///< alpha_k = beta / alpha
		cblas_daxpy(size, alpha, p, 1, x, 1); 	///< x_{k+1} = x_k + alpha_k * p_k
		cblas_daxpy(size, -alpha, t, 1, r, 1); 	///< r_{k+1} = r_k - alpha_k * t
		/// Take norm2 as the criteria of iteration
		eps = cblas_dnrm2(size, r, 1);
		if (eps < tol*tol)
		{
			flag = 1;
			MPI_Reduce(&flag, &FLAG, 1, MPI_INT, MPI_PROD, 0, mygrid);
			if (rank == 0) flag = FLAG;
			MPI_Bcast(&flag, 1, MPI_INT, 0, mygrid);	
		}
		SendAndRecv(x);
		if (flag == 1) break;
		beta = cblas_ddot(size, r, 1, r, 1) / beta;	///< beta_k = r_{k+1}^T * r_{k+1} / dbeta
		cblas_dcopy(size, r, 1, t, 1);	///< t = r_{k+1}, avoid overwritting r 
		cblas_daxpy(size, beta, p, 1, t, 1);	///< t = t + beta_k * p_k
		cblas_dcopy(size, t, 1, p, 1);	///< p_{k+1} = t;
		
		k++;
	}while(1);	///< maximum iteration time
	
	delete[] r;
	delete[] p;
	delete[] t;
}

<<<<<<< HEAD
void PoissonSolver::Solve_Chol(double *x, const double *f)
{
	BoundaryVector();
	/// Only do the factorization once
	if (info != 0)	CholeskyFactor();
	/// Set RHS vector
	cblas_daxpy(size, -1.0, f, 1, b, 1);
	cblas_dscal(size, -1.0, b, 1);
	/// Solve Ax = RHS
	F77NAME(dpbtrs) ('U', size, Ny, 1, A, lda, b, size, info);
	cblas_dcopy(size, b, 1, x, 1);
}

=======
>>>>>>> parent of 940f1f8... (1) Parallel program works fine, but really slow.

