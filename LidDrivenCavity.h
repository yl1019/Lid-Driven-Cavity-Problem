#ifndef LID_DRIVEN_CAVITY
#define LID_DRIVEN_CAVITY

#include <string>
#include <mpi.h>

#include "PoissonSolver.h"

using namespace std;

/**
 * @class LidDrivenCavity
 * @brief Contains information in 2d flow domain
 */
class LidDrivenCavity
{
public:
    /** Constructor and Destructor */
    LidDrivenCavity(MPI_Comm mygrid, int rank, int *coords, int *neighbor, int nx,
	       	int ny, double deltat, double finalt, double re, double dx, double dy, bool &dt_flag);
    ~LidDrivenCavity();

    /** Setters */
    void SetGridSize(int nx, int ny);
    bool SetTimeStep(double deltat);
    void SetFinalTime(double finalt);
    void SetReynoldsNumber(double Re);
   
    /** Methods for solver */
    void LinearMatrices();	
    void BoundaryVector(double *b, char matrix, char x);
    void Initialise();
    void Integrate();

    void VorticityBCs();
    void VorticityInterior();
    void VorticityUpdate();
    void SolvePoisson();
    void Output();
    void Solve();

private: 
    /** MPI configuration */
    MPI_Comm mygrid;
    int rank;
    int coords[2];
    int neighbor[4];

    /** Vorticity and streamfunction interior and boundary vectors */
    double *v = nullptr;
    double *s = nullptr;
    double *v_top = nullptr;
    double *v_left = nullptr;
    double *v_bot = nullptr;
    double *v_right = nullptr;
    double *s_top = nullptr;
    double *s_left = nullptr;
    double *s_bot = nullptr;
    double *s_right = nullptr;

    /** User input parameters for each partition */
    double dt;	///< time step size
    double T;	///< terminal time
    int    Nx;  ///< number of interior grid points in x-direction
    int    Ny; 	///< number of interior grid points in y-direction
    double Lx;  ///< length of subdomain in x-direction
    double Ly;  ///< length of subdomain in y-direction
    double Re;  ///< Reynolds number
    double dx;  ///< grid spacing in x-direction
    double dy;  ///< grid spacing in y-direction

    /** Matrices of linear system */
    double *A = nullptr;	///< store a symmetric banded matrix of Laplacian equation
    double *B = nullptr;	///< store a banded matrix with only subdiagonal values
    double *C = nullptr;	///< store a banded matrix with other off-diagonal values
    int lda, ldb, ldc;
    int size;

    PoissonSolver *ps = nullptr;	///< poisson solver object
};
    
#endif
