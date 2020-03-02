/**
 * @class PoissonSolver
 * @brief Solving the poisson problem Ax + b = f using Conjugate Gradient algotithm
 */

class PoissonSolver
{
public:
	PoissonSolver();
	~PoissonSolver();
	
	void SetSize(int &size);
	void SetBoundary(double *b);
	void SetSourse(double *f);
	void SetMatrix(double *A, int &bwidth);
	void InitialGuess(double *x);
	void Solve();

private:
	int size = 0;         ///< size of the linear system
	double *x = nullptr;        ///< store solution vector
	double *b = nullptr;        ///< store boundary vector
	double *f = nullptr;        ///< store sourse vector
	double *A = nullptr;        ///< store matrix A	
	int bwidth;             ///< store banded width of A
}

