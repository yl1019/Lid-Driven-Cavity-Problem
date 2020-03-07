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
 * @param neighbor	array to store ranks of neighbor. sequece: top, left, bot, right
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
 * @param coords	coordinate of current rank
 * @param nx	local number of grids in x direction
 * @param ny	local number of grids in y direction
 */
void DistributeWork(const int &Nx, const int &Ny, const int &Px, const int &Py, int *coords, 
		int &nx, int &ny)
{
	/** ranks with smaller coordinates have 1 larger size in grids */
	nx = Nx / Px;
	ny = Ny / Py;
	int remainder_x = Nx % Px;
	int remainder_y = Ny % Py;
	if (coords[0] < remainder_x)	nx++;
	if (coords[1] < remainder_y)	ny++;

}
