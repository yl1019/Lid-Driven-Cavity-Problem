#include <iostream>
#include <mpi.h>

using namespace std;
#include "MpiConfig.h"

/**
 * @brief Check if the number of workers is compatible with number of partitions
 * @param np	number of processores
 * @param Px	number of partitions in x-direction
 * @param Py	number of partitions in x-direction
 * @returns	true if compatible, false if compatible
 */
bool CheckWorkers(const int &np, const int &Px, const int &Py)
{
	if (np == Px * Py) return true;
	else	return false;
}

/**
 * @brief Find the neighborhood in four directions for current rank
 * @param mygrid	Cartesian communicator
 * @param neighbor	array to store ranks of neighbor. sequence: top, left, bot, right
 */
void FindNeighbor(MPI_Comm mygrid, int *neighbor)
{
	MPI_Cart_shift(mygrid, 0, 1, neighbor + 2, neighbor);	///< bot and top
	MPI_Cart_shift(mygrid, 1, 1, neighbor + 1, neighbor + 3);	///< left and right
}

/**
 * @brief Distribute work as evenly as possible to each worker
 * @param Nx	global number of grids in x direction
 * @param Ny	global number of grids in y direction
 * @param Px	number of partitions in x-direction
 * @param Py	number of partitions in x-direction
 * @param dx	element size in x direction
 * @param dy	element size in y direction
 * @param coords	coordinate of current rank
 * @param start	starting point in global coordinate, used to output
 * @param nx	local number of grids in x direction
 * @param ny	local number of grids in y direction
 */
void DistributeWork(const int &rank, const int &Nx, const int &Ny, const int &Px, const int &Py,
	       	const double &dx, const double &dy, int *coords, double *start, int &nx, int &ny)
{
	/** ranks with smaller coordinates have 1 larger size in grids */
	nx = (Nx-2) / Px;
	ny = (Ny-2) / Py;
	int remainder_x = (Nx-2) % Px;
	int remainder_y = (Ny-2) % Py;

	// Note here start[i, j] and coords[j, i] !!!
	// x direction
	if (coords[1] < remainder_x)
	{	
		nx++;
		start[0] = (coords[1]*nx + 1)*dx;
	}
	else
	{
		start[0] = (Nx-1)*dx - (Px - coords[1])*dx*nx;
	}
	// y direction
	if (coords[0] < remainder_y)
	{	
		ny++;
		start[1] = (coords[0]*ny + 1)*dy;
	}
	else
	{
		start[1] = (Ny-1)*dy - (Py - coords[0])*dy*ny;
	}
	
	// anouncement before working
	cout << "This is rank " << rank << ", my work load is nx = " <<  nx <<  ", ny = " <<  ny
		<< "; start from global coords(" << start[0] << ", " << start[1] <<
	       	") with Cartesian coords ["<< coords[0] << ", " << coords[1] << "]" << endl;
}


