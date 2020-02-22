#include "LidDrivenCavity.h"
#include "CLAlib.h"

LidDrivenCavity::LidDrivenCavity()
{
}

LidDrivenCavity::~LidDrivenCavity()
{
	delete[] v;
	delete[] s;
}

void LidDrivenCavity::SetDomainSize(double xlen, double ylen)
{
	Lx = xlen;
	Ly = ylen;
}
void LidDrivenCavity::SetGridSize(int nx, int ny)
{
	Nx = nx;
	Ny = ny;
}
void LidDrivenCavity::SetTimeStep(double deltat)
{
	dt = deltat;
}
void LidDrivenCavity::SetFinalTime(double finalt)
{
	T = finalt;
}
void LidDrivenCavity::SetReynoldsNumber(double re)
{
	Re = re;
}
void LidDrivenCavity::GridSpace()
{
	dx = Lx / (Nx - 1);
	dy = Ly / (Ny - 1);
}

/**
 * @brief Construct constant matrices for linear systems
 * @param A matrix for calculate v at time t, i.e. v = As + b
 */
void LidDrivenCavity::LinearMatrices(double *A, double *B, double *C)
{
	int size = (Nx-2)*(Ny-2);        ///< size of matrices
	A = new double[size * size];      ///< only consider inner points
	B = new double[size * size];      ///< only consider inner points
	C = new double[size * size];      ///< only consider inner points
	double A_DiagVal = 2/(dx*dx) + 2/(dy*dy); ///< diagonal entries
	double A_SubDiagVal = -1/(dy*dy);     ///< entries above and below the diagonal 
	double A_OffDiagVal = -1/(dx*dx);     ///< entries in the off-diagonal blocks 
	double B_SubDiagVal = 1/(2*dy);
	double C_OffDiagVal = 1/(2*dx);
	/// Diagonal entries
	for (int i = 0; i < size-1; i++)
	{
		A[i*size + i] = A_DiagVal;
	}
	/// Sub-diagonal entries
	for (int i = 0; i < size-2; i++)
	{
		if ((i+1) % (Ny-2) != 0)
		{
			A[i*size + i+1] = A_SubDiagVal;
			A[(i+1)*size + i] = A_SubDiagVal;
			B[i*size + i+1] = B_SubDiagVal;
			B[(i+1)*size + i] = -B_SubDiagVal;
		}
	}
	/// Other off-diagonal entries
	for (int i = 0; i < size-1 - (Ny-2); i++)
	{
		A[i*size + (i+Ny-2)] = A_OffDiagVal;
		A[(i+Ny-2)*size + i] = A_OffDiagVal;
		C[i*size + (i+Ny-2)] = C_OffDiagVal;
		C[(i+Ny-2)*size + i] = -C_OffDiagVal;
	}
}

void LidDrivenCavity::Initialise()
{
	/// Allocating memory including ghost cells, initialise to zero
	v = new double[Nx*Ny]();
	s = new double[Nx*Ny]();
}

void LidDrivenCavity::Integrate()
{
}

/**
 * @brief Calculate Vorticity boundary conditions at time t
 * @param U Lid velocity
 */
void LidDrivenCavity::VorticityBCs(int U)
{
	/// Top boundary
	for (int i = 1; i < Nx-1; i++)
	{
		v[i*Ny+(Ny-1)] = (s[i*Ny+(Ny-1)] - s[i*Ny+(Ny-2)]) * 2
		      / (dy*dy) - 2*U/dy;
	}
	/// Bottom boundary
	for (int i = 1; i < Nx-1; i++)
	{
		v[i*Ny] = (s[i*Ny] - s[i*Ny+1]) * 2 / (dy*dy);
	}

	/// Left boundary
	for (int j = 1; j < Ny-1; j++)
	{
		 v[j] =(s[j] - s[Ny+j]) * 2 / (dx*dx); 
	}
	/// Right boundary
	for (int j = 1; j < Ny-1; j++)
	{
		 v[(Nx-1)*Ny+j] = (s[(Nx-1)*Ny+j] - s[(Nx-2)*Ny+j]) * 2
			 / (dx*dx);
	}
}

/**
 * @brief Calculate interior vorticity at time t, i.e. calculate v = As + b;
 */
void LidDrivenCavity::VorticityInterior();
{
	/// Construct RHS vector b using boundary conditions
	static int size = (Nx-2)*(Ny-2);           ///< size of vector b 
	double *b =  new double[size]();    ///< allocate memory and initialise to zero
	static double SubDiagVal = -1/(dy*dy);     ///< entries above and below the diagonal 
	static double OffDiagVal = -1/(dx*dx);     ///< entries in the off-diagonal blocks 
	for (int i = 0; i < Ny-2; i++)
	{
		b[i] += s[i+1] * OffDiagVal;
		b[i+(Nx-1)*(Ny-2)] += s[(Nx-1)*Ny + i+1] * OffDiagVal;
		b[i*(Ny-2)] = s[Ny*(i+1)] * SubDiagVal;
		b[Ny-3 + i*(Ny-2)] = s[Ny*(i+1) + Ny-1] * SubDiagVal;
	}

	/// Calculate v = As + b;
	F77NAME(dgemv) ('T', size, size, 1.0, A, size, &s[1], 1, 1.0, b, 1 );
	memcpy(&v[1], b, size * sizeof(double));
	delete[] b;	
}

/**
 * @brief Calculate interior vorticity at time t+dt, i.e. calculate v = v + b, where b includes three parts
 */
void LidDrivenCavity::VoriticityUpdate();
{
	static int size = (Nx-2)*(Ny-2);   ///< size of vector b
	/// three parts of vector b
	double *b1 = new double[size]();   ///< viscosity term
	double *b2 = new double[size]();   ///< advection term
	double *b3 = new double[size]();   ///< also advection term
	double A_SubDiagVal = -1/(dy*dy);     ///< entries above and below the diagonal 
	double A_OffDiagVal = -1/(dx*dx);     ///< entries in the off-diagonal blocks 

	static double B_SubDiagVal = 1/(2*dy);  
	static double C_OffDiagVal = 1/(2*dx);  
	
	/// First calculate b1-----------------------------------
	for (int i = 0; i < Ny-2; i++)
	{
		b1[i] += -v[i+1] * A_OffDiagVal;
		b1[i+(Nx-1)*(Ny-2)] += -v[(Nx-1)*Ny + i+1] * A_OffDiagVal;
		b1[i*(Ny-2)] = -v[Ny*(i+1)] * A_SubDiagVal;
		b1[Ny-3 + i*(Ny-2)] = -v[Ny*(i+1) + Ny-1] * A_SubDiagVal;
	}
	F77NAME(dgemv) ('T', size, size, -dt/Re, A, size, &v[1], 1, 1.0, b1, 1 );
	/// Calculate b2-----------------------------------------
	double *temp1 = new double[size](); ///< temp1 = Cs + f, f is boundary term
	double *temp2 = new double[size](); ///< temp2 = Bv + f
	// boundary term for temp1
	for (int i = 0; i < Ny-2; i++)
	{
		temp1[i] += -s[i+1] * C_OffDiagVal;
		temp1[i+(Nx-1)*(Ny-2)] += s[(Nx-1)*Ny + i+1] * C_OffDiagVal;
	}
	// boundary term for temp2
	for (int i = 0; i < Ny-2; i++)
	{
		temp2[i*(Ny-2)] = -v[Ny*(i+1)] * B_SubDiagVal;
		temp2[Ny-3 + i*(Ny-2)] = v[Ny*(i+1) + Ny-1] * B_SubDiagVal;
	}
	// calculate temp1 and temp2	
	F77NAME(dgemv) ('T', size, size, 1.0, C, size, &s[1], 1, 1.0, temp1, 1 );

	F77NAME(dgemv) ('T', size, size, 1.0, B, size, &v[1], 1, 1.0, temp2, 1 );
	// calculate b2 by multiply each vector element
	for (int i = 0; i < size; i++)
	{
		b2[i] = temp1[i] * temp2[i];
	}
	// release memory for temp1 and temp2
	delete[] temp1;
	delete[] temp2;

	/// Calculate b3--------------------------------------------
	double *temp1 = new double[size](); ///< temp1 = Bs + f, f is boundary term
	double *temp2 = new double[size](); ///< temp2 = Cv + f
	// boundary term for temp1
	for (int i = 0; i < Ny-2; i++)
	{
		temp1[i*(Ny-2)] = -s[Ny*(i+1)] * B_SubDiagVal;
		temp1[Ny-3 + i*(Ny-2)] = s[Ny*(i+1) + Ny-1] * B_SubDiagVal;
	}
	// boundary term for temp2
	for (int i = 0; i < Ny-2; i++)
	{
		temp2[i] += -v[i+1] * C_OffDiagVal;
		temp2[i+(Nx-1)*(Ny-2)] += v[(Nx-1)*Ny + i+1] * C_OffDiagVal;
	}
	// calculate temp1 and temp2	
	F77NAME(dgemv) ('T', size, size, 1.0, B, size, &s[1], 1, 1.0, temp1, 1 );

	F77NAME(dgemv) ('T', size, size, 1.0, C, size, &v[1], 1, 1.0, temp2, 1 );
	// calculate b2 by multiply each vector element
	for (int i = 0; i < size; i++)
	{
		b3[i] = temp1[i] * temp2[i];
	}
	// release memory for temp1 and temp2
	delete[] temp1;
	delete[] temp2;

	/// Update interior vorticity-------------------------------
	F77NAME(daxpy) (size, 1.0, b1, 1, s_int, 1);
	F77NAME(daxpy) (size, 1.0, b2, 1, s_int, 1);
	F77NAME(daxpy) (size, -1.0, b3, 1, s_int, 1);

	/// release memory
	delete[] b1;
	delete[] b2;
	delete[] b3;
}
