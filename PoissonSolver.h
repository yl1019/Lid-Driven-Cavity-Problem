#ifndef POISSON_SOLVER
#define POISSON_SOLVER

/**
 * @class PoissonSolver
 * @brief Solving the poisson problem Ax + b = f using Conjugate Gradient algorithm
 */
class PoissonSolver
{
public:
	PoissonSolver(const int &Nx, const int &Ny, const double &dx, const double &dy);
	~PoissonSolver();
	void SetBoundary(const double *top, const double *left, const double *bot, const double *right);
	void Solve(double *x, const double *f);

private:
	/// Size of the linear system
	int Nx;
	int Ny;
	double dx;
	double dy;
	int size;     

	double *A = nullptr;	///< store matrix A
	double A_DiagVal;
	double A_SubDiagVal;
	double A_OffDiagVal;	
	int lda;	///< leading dimension of A	
	double *b = nullptr;	///< store boundary vector
};

#endif
