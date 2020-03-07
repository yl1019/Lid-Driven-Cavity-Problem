#ifndef MPI_CONFIG
#define MPI_CONFIG

#include <mpi.h>

bool CheckWorkers(const int &np, const int &Px, const int &Py);
void FindNeighbor(MPI_Comm mygrid, int *neighbor);
void DistributeWork(const int &Nx, const int &Ny, const int &Px, const int &Py, int *coords, 
		int &nx, int &ny);

#endif
