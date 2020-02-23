#pragma once

#include <string>
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
    void SetTimeStep(double deltat);
    /// Set terminal time
    void SetFinalTime(double finalt);
    /// Set Reynolds number for the problem
    void SetReynoldsNumber(double Re);
    /// Calculate grid space
    void GridSpace();
    /// Construct constant matrices for computing
    void LinearMatrices();	
    /// extract interior vector from whole matrix for solving linear system 
    void Extract();


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
    /// Possion solver, to be separated later
    void PossionSolver();

    
private:
    double *v = nullptr;   ///< vorticity matrix stored in row major
    double *s = nullptr;   ///<	stream function matrix stored in row major
    double *v_int = nullptr; ///< interior vorticity matrix
    double *s_int = nullptr; ///< interior stream function matrix

    double dt;             ///< time step size
    double T;              ///< terminal time
    int    Nx;             ///< number of grid points in x-direction
    int    Ny; 	           ///< number of grid points in y-direction
    double Lx;		   ///< length of subdomain in x-direction
    double Ly;             ///< length of subdomain in y-direction
    double Re;             ///< Reynolds number
    double dx;             ///< grid spacing in x-direction
    double dy;             ///< grid spacing in y-direction

    double *A = nullptr;   ///< store constant linear matrix to compute v
};
    double *B = nullptr;   ///< store constant matrix with only subdiagonal values
    double *C = nullptr;   ///< store constant matrix with other off-diagonal values
    

