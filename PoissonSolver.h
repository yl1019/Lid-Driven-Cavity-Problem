#ifndef POISSON_SOLVER
#define POISSON_SOLVER

/**
 * @class PoissonSolver
 * @brief Solving the poisson problem Ax + b = f using Conjugate Gradient algorithm / Cholesky factorization
 */
class PoissonSolver
{
public:
	PoissonSolver(const int &Nx, const int &Ny, const double &dx, const double &dy);
	~PoissonSolver();
	void SetBoundary(const double *top, const double *left, const double *bot, const double *right);
	void Solve_Conj(double *x, const double *f);
	void CholeskyFactor();
	void Solve_Chol(double *x, const double *f);
private:
	/// Size of the linear system
	int Nx;
	int Ny;
	double dx;
	double dy;
	int size, info;     

	double *A = nullptr;	///< store matrix A
	double A_DiagVal;
	double A_SubDiagVal;
	double A_OffDiagVal;	
	int lda;	///< leading dimension of A	
	double *b = nullptr;	///< store boundary vector
};

#endif
