#ifndef POISSON_SOLVER
#define POISSON_SOLVER

#include "mpi.h"
/**
 * @class PoissonSolver
 * @brief Solving the poisson problem Ax + b = f using Conjugate Gradient algorithm
 */
class PoissonSolver
{
public:
	PoissonSolver(MPI_Comm mygrid, int *neighbor, int rank, const int &Nx,
		       	const int &Ny, const double &dx, const double &dy);
	~PoissonSolver();
<<<<<<< HEAD
	void BoundaryVector();
	void SendAndRecv(double *x);
	void SetBoundary(const double *top, const double *left,
		       	const double *bot, const double *right);
	void Solve_Conj(double *x, const double *f);
	void CholeskyFactor();
	void Solve_Chol(double *x, const double *f);
=======
	void SetBoundary(const double *top, const double *left, const double *bot, const double *right);
	void Solve(double *x, const double *f);

>>>>>>> parent of 940f1f8... (1) Parallel program works fine, but really slow.
private:
	/** MPI configuration */
	MPI_Comm mygrid;
	int neighbor[4];
	int rank;
	

	/** Size of the linear system */
	int Nx;
	int Ny;
	double dx;
	double dy;
	int size;     

	int FLAG;	///< FLAG for root process
	int flag;	///< flag for stop iteration, 1 for stop and 0 for keep going
	double eps;
	double *A = nullptr;	///< store matrix A
	double A_DiagVal;
	double A_SubDiagVal;
	double A_OffDiagVal;	
	int lda;	///< leading dimension of A	
	double *b = nullptr;	///< store boundary vector


	double *x_left = nullptr;
	double *x_right = nullptr;
	double *x_top = nullptr;
	double *x_bot = nullptr;
};

#endif
