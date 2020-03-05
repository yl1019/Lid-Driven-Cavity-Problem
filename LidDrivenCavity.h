#ifndef LID_DRIVEN_CAVITY
#define LID_DRIVEN_CAVITY

#include <string>
#include "PoissonSolver.h"
using namespace std;

/**
 * @class LidDrivenCavity
 * @brief Contains information in 2d flow domain
 */
class LidDrivenCavity
{
public:
    /// Default Constructor
    LidDrivenCavity();
    /// Clean up and deallocate memory
    ~LidDrivenCavity();

    /// Set Domain size of each process
    void SetDomainSize(double xlen, double ylen);
    /// Set the grid of each subspace
    void SetGridSize(int nx, int ny);
    /// Set time step with respect to the restricion
    bool SetTimeStep(double deltat);
    /// Set terminal time
    void SetFinalTime(double finalt);
    /// Set Reynolds number for the problem
    void SetReynoldsNumber(double Re);
    /// Calculate grid space
    void GridSpace();
    /// Construct constant matrices for computing
    void LinearMatrices();	
    /// Generate vector b of linear system using boundary conditions
    void BoundaryVector(double *b, char matrix, char x);

    /// Initialise the flow quantities
    void Initialise();
    /// Gather all information to rank 0 for output
    void Integrate();

    /// Calculate vorticity in ghost cells
    void VorticityBCs();
    /// Calculate interior vorticity at time t
    void VorticityInterior();
    /// Update interior vorticity
    void VorticityUpdate();
    /// Solver poisson equation using the class PossionSolver
    void SolvePoisson();
    /// Output the solutions in different files
    void Output();
    
private: 
    double *v;	///< vorticity vector
    double *s;  ///< stream function vector
    /// Vorticity and streamfunction boundaries
    double *v_top, *v_left, *v_bot, *v_right;
    double *s_top, *s_left, *s_bot, *s_right;
    /// User input parameters for each partition
    double dt;	///< time step size
    double T;	///< terminal time
    int    Nx;  ///< number of interior grid points in x-direction
    int    Ny; 	///< number of interior grid points in y-direction
    double Lx;  ///< length of subdomain in x-direction
    double Ly;  ///< length of subdomain in y-direction
    double Re;  ///< Reynolds number
    double dx;  ///< grid spacing in x-direction
    double dy;  ///< grid spacing in y-direction

    double *A;	///< store a symmetric banded matrix of Laplacian equation
    double *B;  ///< store a banded matrix with only subdiagonal values
    double *C;  ///< store a banded matrix with other off-diagonal values
    /// leading dimension of these matrices (for packed storage)
    int lda, ldb, ldc;
    int size;	///< size of the linear matrices

    PoissonSolver *ps;	///< poisson solver object
};
    
#endif
