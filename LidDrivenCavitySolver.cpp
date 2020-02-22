#include <iostream>
#include <exception>
using namespace std;

#include "LidDrivenCavity.h"

int main(int argc, char **argv)
{
    /// Parameters from command input
    double Lx = 1, Ly = 1;
    int Nx = 161; Ny = 161;
    int Px = 1; Py = 1;
    double dt = 0.01;
    double T = 5.0;
    double Re = 100;
    try
    {
    

   	 /// Create a new instance of the LidDrivenCavity class
   	 LidDrivenCavity* solver = new LidDrivenCavity();

   	 /// ...
   	 solver->Initialise();

    	/// Run the solver
   	 solver->Integrate();
    }
    catch
    {
	    
    }
	return 0;
}
