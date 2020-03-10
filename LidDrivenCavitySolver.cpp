#include <iostream>
#include <exception>
#include <iomanip>
#include <mpi.h>
using namespace std;

#include "LidDrivenCavity.h"
#include "ProgramOptions.h"
#include "MpiConfig.h"
int main(int argc, char **argv)
{
	/** Initialize MPI and retrieve rank & size */
	MPI_Init(&argc, &argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

    	/** Parameters from command input */
    	double Lx, Ly;
    	int Nx, Ny;
    	int Px, Py;
   	double dt;
    	double T;
   	double Re;
    	po::variables_map vm;
    	bool status;	///< to decide the input status, 1 for successful input
    	status = OptionStatus(argc, argv, vm);
    	/// If help is called or error occurs, terminate the program
    	if (!status)
	{
		MPI_Finalize();
		return 0;
	}
    	/// Read all parameters
    	ReadVals(vm, Lx, Ly, Nx, Ny, Px, Py, dt, T, Re);

	/** Check if the number of processors are compatible with the size */ 
    	double dx = Lx/(Nx - 1);
	double dy = Ly/(Ny - 1);
	if (!CheckWorkers(size, Px, Py))
	{
		if (rank == 0)
		{
			cout << "Number of processores doe not match input Px and Py." << endl;
			cout << "Please ensure np = Px * Py." << endl;
		}
		MPI_Finalize();
		return 0;
	}

	/** Work before operating */
	MPI_Comm mygrid;
	const int dims = 2;
	int sizes[dims] = {Py, Px};
	int periods[dims] = {0, 0};
	int reorder = 0;
	/// Create a new communicator based on Cartesian Topology
	MPI_Cart_create(MPI_COMM_WORLD, dims, sizes, periods, reorder, &mygrid);
	int coords[dims];
	/// Give each rank a coordinate
	MPI_Cart_coords(mygrid, rank, dims, coords);
	int neighbor[4] = {0};
	/// Obtain neighborhood information for each rank
	FindNeighbor(mygrid, neighbor);
	
//	cout << "rank " << rank << " 's neighbor is " << neighbor[0] << "  "
//	    <<   	neighbor[1] << "  " << neighbor[2] << " " << neighbor[3] << endl;

	int nx, ny;	///< number of grids for each rank
	double start[2] = {0.0, 0.0};
	/// Distribute work to each process
	DistributeWork(rank, Nx, Ny, Px, Py, dx, dy, coords, start, nx, ny);
		
//	MPI_Finalize();
//	return 0;

	/** Solving the main problem */
	bool dt_flag;
   	/// Create a new instance of the LidDrivenCavity class for each rank
   	LidDrivenCavity* solver = new LidDrivenCavity(mygrid, rank, coords, start, neighbor,
		       	nx, ny, dt, T, Re, dx, dy, dt_flag);
	if (!dt_flag)
	{	
		MPI_Finalize();
		return 0;
	}
	/// Run the solver
        solver->Solve();
	/// Output the result
	solver->OUTPUT(Px, Py, Lx, Ly);
	//solver->Output(Lx, Ly, Px, Py, Nx-2, Ny-2);


	/// Exit
	MPI_Finalize();
	return 0;
}
