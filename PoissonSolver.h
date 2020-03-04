#ifndef POISSON_SOLVER
#define POISSON_SOLVER

/**
 * @class PoissonSolver
 * @brief Solving the poisson problem Ax + b = f using Conjugate Gradient algotithm
 */
class PoissonSolver
{
public:
	PoissonSolver(const int &Nx, const int &Ny, const double &dx, const double &dy);
	~PoissonSolver();
	void SetBoundary(double *top, double *left, double *bot, double *right);
	void Solve(double *x, double *f);

private:
	/// size of the linear system problem
	int Nx;
	int Ny;
	double dx;
	double dy;
	int size;     

	double *A = nullptr;        ///< store matrix A
	double A_DiagVal;
	double A_SubDiagVal;
	double A_OffDiagVal;	
	int lda;                    ///< leading dimension of A	
	double *b = nullptr;        ///< store boundary vector
};

#endif
