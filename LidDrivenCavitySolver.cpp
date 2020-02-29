#include <iostream>
#include <exception>
#include <iomanip>
using namespace std;

#include "LidDrivenCavity.h"

int main(int argc, char **argv)
{
    /// Parameters from command input
    double Lx = 1, Ly = 1;
    int Nx = 10, Ny = 10;
    int Px = 1, Py = 1;
    double dt = 0.1;
    double T = 10.0;
    double Re = 100;

    /// Parameters for each partition
    double xlen = Lx/Px;
    double ylen = Ly/Py;
    int nx = (Nx-2)/Px;
    int ny = (Ny-2)/Py;
    
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
		 // print the matrix out to debug
	       	 cout.precision(3);

		 cout << "vorticity: " << endl;
		 for (int i = 0; i < nx; i++)
		 {
		 	 for (int j = 0; j < ny; j++)
		 	 {
				 cout << setw(10) << solver->v[i*ny+j];
			 }
			 cout << endl;
		 }
		 cout << "stream function: " << endl;
		 for (int i = 0; i < nx; i++)
		 {
		 	 for (int j = 0; j < ny; j++)
		 	 {
				 cout << setw(10) << solver->s[i*ny+j];
			 }
			 cout << endl;
		 }
	 }
	 solver->Output();

	 return 0;
}
