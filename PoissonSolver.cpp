#include<cmath>
#include<cstring>
using namespace std;

#include "LidDrivenCavity.h"
#include "PoissonSolver.h"
#include "cblas.h"



PoissonSolver::PoissonSolver()
{}
PoissonSolver::~PoissonSolver()
{}
/// Set the linear system
void PoissonSolver::SetSize(int &size)
{
	this->size = size;
}

void PoissonSolver::SetBoundary (double *b)
{
	if (b == nullptr) this->b = new double[size];
	memcpy(this->b, b, size*sizeof(double));
}

void PoissonSolver::SetSource(double *f)
{
	if (f == nullptr) this->f = new double[size];
	memcpy(this->f, f, size*sizeof(double));
}

/**
 * @brief Passing the established matrix A to the solver
 * @param A the symmetric banded P.D. matrix
 * @param bwidth banded width of A
 */
void PoissonSolver::SetMatrix(double *A, int &bwidth)
{
	/// No need to construct A twice
	if (A == nullptr)
	{
		this->A = new double[size*size];
		memcpy(this->A, A, size*size*sizeof(double));
		this->bwidth = bwidth;
	}
}

/// Initial guess from the last iteration
void PoissonSolver::InitialGuess(double *x)
{
	if (x == nullptr) this->x = new double[size];
	memcpy(this->x, x, size*sizeof(double));
}

/// Solve the linear system
void PoissonSolver::Solve()
{
	/// if no initial guess inputed, set to be zero
	if (x == nullptr) memset(x, 0, size*sizeof(double));
	int lda = bwidth + 1;
	double *r = new double[size];
	double *p = new double[size];
	double *t = new double[size]; /// for temporary use
	double alpha;
	double beta;
	double eps;
	double tol = 0.001;

	cblas_dcopy(size, b, 1, r, 1); ///< r0 = b
	cblas_dsbmv(CblasColMajor, CblasUpper, size, bwidth ,-1.0, A, lda, x, 1, 1.0, r, 1);         ///< r0 = b - Ax0
	cblas_dcopy(size, r, 1, p, 1); ///< p0 = r0
	int k = 0;
	do
	{
		cblas_dsbmv(CblasColMajor, CblasUpper, size, bwidth,1.0, A, lda, p, 1, 0.0, t, 1);   ///< t = A*p_k
		alpha = cblas_ddot(size, t, 1, p, 1); ///< alpha = p_k^T * A * p_k
		beta = cblas_ddot(size, r, 1, r, 1); ///< beta = r_k^T * r_k 
		alpha = beta / alpha;              ///< alpha_k = beta / alpha
		cblas_daxpy(size, alpha, p, 1, x, 1); ///x_{k+1} = x_k + alpha_k * p_k
		cblas_daxpy(size, -alpha, t, 1, r, 1); ///r_{k+1} = r_k - alpha_k * t
		/// take norm2 as the criteria of iteration
		eps = cblas_dnrm2(size, r, 1);
		/// avoid calling sqrt
		if (eps < tol*tol) break;
		beta = cblas_ddot(size, r, 1, r, 1) / beta; ///beta_k = r_{k+1}^T * r_{k+1} / beta
		cblas_dcopy(size, r, 1, t, 1);  /// t = r_{k+1}, avoid overwritting r 
		cblas_daxpy(size, beta, p, 1, t, 1); /// t = t + beta_k * p_k
		cblas_dcopy(size, t, 1, p, 1); /// p_{k+1} = t;
		
		k++;
	}while(k < 500);// maximum number of iteration
	
	delete[] r;
	delete[] p;
	delete[] t;
}
