#include <iostream>
#include <exception>
using namespace std;

#include "LidDrivenCavity.h"

int main(int argc, char **argv)
{
    /// Parameters from command input
    double Lx = 1, Ly = 1;
    int Nx = 21, Ny = 21;
    int Px = 1, Py = 1;
    double dt = 0.01;
    double T = 5.0;
    double Re = 100;

    /// Parameters for each partition
    double xlen = Lx/Px;
    double ylen = Ly/Py;
    int nx = (Nx-2)/Px + 2;
    int ny = (Ny-2)/Py + 2;
    
   //try
	 
   	 /// Create a new instance of the LidDrivenCavity class
   	 LidDrivenCavity* solver = new LidDrivenCavity();

   	 solver->SetDomainSize(xlen, ylen);
	 solver->SetGridSize(nx, ny);
	 solver->SetTimeStep(dt);
	 solver->SetFinalTime(T);
	 solver->SetReynoldsNumber(Re);
	 solver->GridSpace();
	 solver->LinearMatrices();
	 solver->Extract();

   	 solver->Initialise();
   	 // solver->Integrate();
	 double t = 0.0;
	 while (t < T)
	 {
		 solver->VorticityBCs();
		 solver->VorticityInterior();
		 solver->VorticityUpdate();
		 solver->PossionSolver();
		 t = t + dt;
		 cout << "t = " << t << endl;
	 }
	return 0;
}
