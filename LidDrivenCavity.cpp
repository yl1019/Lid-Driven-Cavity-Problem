#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <mpi.h>

#include "LidDrivenCavity.h"
#include "cblas.h"
#include "PoissonSolver.h"

/* Constructor and destructor
======================================================================== */
LidDrivenCavity::LidDrivenCavity(MPI_Comm mygrid, int rank, int *coords, double *start, int *neighbor,
	       	int nx,	int ny, double deltat, double finalt, double re, double dx, double dy,
	       	bool &dt_flag)
{
	this->mygrid = mygrid;
	this->rank = rank;
	memcpy(this->coords, coords, 2*sizeof(int));
	memcpy(this->neighbor, neighbor, 4*sizeof(int));
	memcpy(this->start, start, 2*sizeof(double));
	Nx = nx;
	Ny = ny;
	size = Nx * Ny;
	lda = Ny + 1;              
	ldb = 3;
	ldc = 2*Ny + 1;
	T = finalt;
	Re = re;
	this->dx = dx;
	this->dy = dy;
	dt_flag = SetTimeStep(deltat);	///< flag to determine the stability restriction on dt
}
LidDrivenCavity::~LidDrivenCavity()
{
	delete[] v;
	delete[] s;
	delete[] v_top;
	delete[] v_bot;
	delete[] v_left;
	delete[] v_right;
	delete[] s_top;
	delete[] s_bot;
	delete[] s_left;
	delete[] s_right;
	delete[] A;
	delete[] B;
	delete[] C;
	delete[] ps;

}

/* Setters
======================================================================== */
/**
 * @brief Set the number of interior grids in x and y direction
 * @param nx	number of interior grids in x direction
 * @param nx	number of interior grids in x direction
 */
void LidDrivenCavity::SetGridSize(int nx, int ny)
{
	Nx = nx;
	Ny = ny;
}
/**
 * @brief Set the time step and check if it satisfies the restriction
 * @returns return 0 if the time step is too large; return 1 if it's OK
 */
bool LidDrivenCavity::SetTimeStep(double deltat)
{
   try
   {
	dt = deltat;
	/// check the restriction on time step dt < Re*dx*dy/4
	if (dt >= Re*dx*dy/4) 
	{
		throw std::logic_error("dt input is too large.");	
	}
   }
   catch (const std::logic_error &e)
   {
	if (rank == 0) std::cout << "An error occured: " << e.what() << std::endl;
	return 0;
   }
   return 1;
}
/**
 * @brief Set final time
 * @param finalt	final time
 */
void LidDrivenCavity::SetFinalTime(double finalt)
{
	T = finalt;
}
/**
 * @brief Set Reynolds number
 * @param re	Reynolds number
 */
void LidDrivenCavity::SetReynoldsNumber(double re)
{
	Re = re;
}

/* Methods for solver
======================================================================== */
/**
 * @brief Construct constant matrices for linear systems
 * @param A	matrix for calculate v at time t(also for v at t+dt), i.e. v = As + b
 * @param B	matrix for calculate v at time t+dt, with only subdiagonal values
 * @param C	matrix for calculate v at time t+dt, with other off-diagonal values
 */
void LidDrivenCavity::LinearMatrices()
{
	A = new double[size * lda];	///< store as symmetric banded matrix with bwidth = Ny
	B = new double[size * ldb];    	///< store as banded matrix with bwidth = 1
	C = new double[size * ldc];	///< store as banded matrix with bwidth = Ny  
	double A_DiagVal = 2/(dx*dx) + 2/(dy*dy);	///< diagonal entries
	double A_SubDiagVal = -1/(dy*dy);	///< entries above and below the diagonal 
	double A_OffDiagVal = -1/(dx*dx);	///< entries in the off-diagonal blocks 
	double B_SubDiagVal = 1/(2*dy);
	double C_OffDiagVal = 1/(2*dx);
	/// Matrix A
	for (int j = 0; j < size; j++)
	{
		A[lda-1 + j*lda] = A_DiagVal;
		if (j % Ny != 0)   A[lda-2 + j*lda] = A_SubDiagVal;
		if (j >= Ny)       A[j * lda] = A_OffDiagVal;
		
	}
	/// Matrix B
	for (int j = 0; j < size; j++)
	{
		if (j % Ny != 0)     B[j*ldb] = B_SubDiagVal;
		if (j % Ny != Ny-1)  B[2 + j*ldb] = -B_SubDiagVal;
	}
	/// Matrix C
	for (int j = 0; j < size; j++)
	{
		if (j >= Ny)      C[j*ldc] = C_OffDiagVal;
		if (j < size-Ny)  C[ldc-1 + j*ldc] = -C_OffDiagVal;
	}
}

/**
 * @brief Construct the boundary vector b in linear system y = Ax + b
 * @param b	the vector b
 * @param matrix	determine the matrix of the linear system (A, B or C)
 * @param x	determine the solution vector (v or s)
 */
void LidDrivenCavity::BoundaryVector(double *b, char matrix, char x)
{	
	static double A_SubDiagVal = -1/(dy*dy); 
	static double A_OffDiagVal = -1/(dx*dx);
	static double B_SubDiagVal = 1/(2*dy);  
	static double C_OffDiagVal = 1/(2*dx);  
	switch(matrix)
	{
		case 'A':
		{	
			if (x == 's')
			{
				for (int i = 0; i < Ny; i++)
				{
					b[i] += s_left[i] * A_OffDiagVal;
					b[i+(Nx-1)*Ny] += s_right[i] * A_OffDiagVal;
				}
				for (int i = 0; i < Nx; i++)
				{
					b[i*Ny] += s_bot[i] * A_SubDiagVal;
					b[i*Ny + Ny-1] += s_top[i] * A_SubDiagVal;
				}
			}
			else if (x == 'v')
			{
				for (int i = 0; i < Ny; i++)
				{
					b[i] += v_left[i] * A_OffDiagVal;
					b[i+(Nx-1)*Ny] += v_right[i] * A_OffDiagVal;
				}
				for (int i = 0; i < Nx; i++)
				{
					b[i*Ny] += v_bot[i] * A_SubDiagVal;
					b[i*Ny + Ny-1] += v_top[i] * A_SubDiagVal;
				}
			} break;
		}
		case 'B':
		{
			if (x == 's')
			{
				for (int i = 0; i < Nx; i++)
				{
					b[i*Ny] = s_bot[i] * (-B_SubDiagVal);
					b[i*Ny + Ny-1] = s_top[i] * B_SubDiagVal;
				}
			}
			else if (x == 'v')
			{
				for (int i = 0; i < Nx; i++)
				{
					b[i*Ny] = v_bot[i] * (-B_SubDiagVal);
					b[i*Ny + Ny-1] = v_top[i] * B_SubDiagVal;
				}
			} break;
		}
		case 'C':
		{	
			if (x == 's')
			{
				for (int i = 0; i < Ny; i++)
				{
					b[i] = s_left[i] * (-C_OffDiagVal);
					b[i+(Nx-1)*Ny] = s_right[i] * C_OffDiagVal;
				}
			}
			else if (x == 'v')
			{
				for (int i = 0; i < Ny; i++)
				{
					b[i] = v_left[i] * (-C_OffDiagVal);
					b[i+(Nx-1)*Ny] = v_right[i] * C_OffDiagVal;
				}
			} break;
		}
		default: break;
	}	
}

/**
 * @brief Initialize all interior and boundary vectors
 */
void LidDrivenCavity::Initialise()
{
	/// Allocating interior memory and initialize to zero
	v = new double[size];
	s = new double[size];	
	memset(v, 0.0, size*sizeof(double));	
	memset(s, 0.0, size*sizeof(double));	
	/// Allocate the boundary memory and initialize to zero
	v_top = new double[Nx];
	v_bot = new double[Nx];
	v_left = new double[Ny];
	v_right = new double[Ny];
	s_top = new double[Nx];
	s_bot = new double[Nx];
	s_left = new double[Ny];
	s_right = new double[Ny];
	memset(v_top, 0.0, Nx*sizeof(double));	
	memset(v_bot, 0.0, Nx*sizeof(double));	
	memset(v_left, 0.0, Ny*sizeof(double));	
	memset(v_right, 0.0, Ny*sizeof(double));	
	memset(s_top, 0.0, Nx*sizeof(double));	
	memset(s_bot, 0.0, Nx*sizeof(double));	
	memset(s_left, 0.0, Ny*sizeof(double));	
	memset(s_right, 0.0, Ny*sizeof(double));	
}



/**
 * @brief Calculate Vorticity boundary conditions at time t
 *        Message passing and receive here
 */
void LidDrivenCavity::VorticityBCs()
{
	/** Top boundary */
	if (neighbor[0] == -2)
	{
		for (int i = 0; i < Nx; i++)
		{
			v_top[i] = (s_top[i] - s[i*Ny+(Ny-1)]) * 2
		     	 / (dy*dy) - 2*U/dy;
		}
	}
	/** Bottom boundary */
	if (neighbor[2] == -2)
	{
		for (int i = 0; i < Nx; i++)
		{
			v_bot[i] = (s_bot[i] - s[i*Ny]) * 2 / (dy*dy);
		}
	}

	/** Left boundary */
	if (neighbor[1] == -2)
	{
		for (int j = 0; j < Ny; j++)
		{	
			 v_left[j] =(s_left[j] - s[j]) * 2 / (dx*dx); 
		}
	}

	/** Right boundary */
	if (neighbor[3] == -2)
	{
		for (int j = 0; j < Ny; j++)
		{
			 v_right[j] = (s_right[j] - s[(Nx-1)*Ny+j]) * 2
			 / (dx*dx);
		}
	}
}

/**
 * @brief Calculate interior vorticity at time t, i.e. calculate v = As + b;
 */
void LidDrivenCavity::VorticityInterior()
{
	/// Construct RHS vector b using boundary conditions
	double *b =  new double[size];
	BoundaryVector(b, 'A', 's');
	cblas_dsbmv (CblasColMajor, CblasUpper, size, Ny, 1.0, A, lda, s, 1, 1.0, b, 1);	///< b = As + b
	cblas_dcopy (size, b, 1, v, 1);	///< v = b
	delete[] b;	
}

/**
 * @brief Calculate interior vorticity at time t+dt, i.e. calculate v = v + f, where b includes three parts
 */
void LidDrivenCavity::VorticityUpdate()
{
	/// Three parts of vector f
	double *b1 = new double[size];	///< viscosity term
	double *b2 = new double[size];	///< advection term
	double *b3 = new double[size];	///< also advection term
	memset(b1, 0.0, size*sizeof(double));	
	memset(b2, 0.0, size*sizeof(double));	
	memset(b3, 0.0, size*sizeof(double));	

	/// First calculate b1 = Av + b
	BoundaryVector(b1, 'A', 'v');
	cblas_dsbmv (CblasColMajor, CblasUpper, size, Ny, 1.0, A, lda, v, 1, 1.0, b1, 1); 

	/// Calculate b2
	double *temp1 = new double[size];	///< temp1 = Cs + b, b is the boundary term
	double *temp2 = new double[size];	///< temp2 = Bv + b
	memset(temp1, 0.0, size*sizeof(double));	
	memset(temp2, 0.0, size*sizeof(double));	
	// calculate temp1
	BoundaryVector(temp1, 'C', 's');
	cblas_dgbmv (CblasColMajor, CblasNoTrans, size, size, Ny, Ny, 1.0, C, ldc, s, 1, 1.0, temp1, 1);
	// calculate temp2
	BoundaryVector(temp2, 'B', 'v');
	cblas_dgbmv (CblasColMajor, CblasNoTrans, size, size, 1, 1, 1.0, B, ldb, v, 1, 1.0, temp2, 1);
	// calculate b2 by multiply each vector element, i.e. b2 = temp1 .* temp2
	for (int i = 0; i < size; i++)
	{
		b2[i] = temp1[i] * temp2[i];
	}
	// release memory for temp1 and temp2
	delete[] temp1;
	delete[] temp2;
	
	/// Calculate b3
	double *temp3 = new double[size];	///< temp3 = Cv + b, b is the boundary term
	double *temp4 = new double[size];	///< temp4 = Bs + b
	memset(temp3, 0.0, size*sizeof(double));	
	memset(temp4, 0.0, size*sizeof(double));	
	// calculate temp3
	BoundaryVector(temp3, 'C', 'v');
	cblas_dgbmv (CblasColMajor, CblasNoTrans, size, size, Ny, Ny, 1.0, C, ldc, v, 1, 1.0, temp3, 1);
	// calculate temp4
	BoundaryVector(temp4, 'B', 's');
	cblas_dgbmv (CblasColMajor, CblasNoTrans, size, size, 1, 1, 1.0, B, ldb, s, 1, 1.0, temp4, 1);
	// calculate b3 by multiply each vector element, i.e. b3 = temp3 .* temp4
	for (int i = 0; i < size; i++)
	{
		b3[i] = temp3[i] * temp4[i];
	}
	// release memory for temp3 and temp4
	delete[] temp3;
	delete[] temp4;

	/// Update interior vorticity, i.e. v = v -(dt/Re)*b1 + dt*b2 - dt*b3
	cblas_daxpy (size, -dt/Re, b1, 1, v, 1);
	cblas_daxpy (size, dt, b2, 1, v, 1);
	cblas_daxpy (size, -dt, b3, 1, v, 1);

	/// release memory
	delete[] b1;
	delete[] b2;
	delete[] b3;
}


/**
 * @brief Call PossionSolver to update stream function at time t+dt
 */
void LidDrivenCavity::SolvePoisson()
{
	/// Create a new PoissonSolver instance and solve Poisson problem
	if (ps == nullptr)    ps = new PoissonSolver(Nx, Ny, dx, dy);
	ps->SetBoundary(s_top, s_left, s_bot, s_right);
	ps->Solve_Chol(s, v);
	//ps->Solve_Conj(s, v);
}

/**
 * @brief Send and receive stream function to neighbor domains
 */
void LidDrivenCavity::SendAndRecv_s()
{
	/** Left to right */
	MPI_Sendrecv(s + (Nx-1)*Ny, Ny, MPI_DOUBLE, neighbor[3], 3, s_left, Ny,
			MPI_DOUBLE, neighbor[1], 3, mygrid, MPI_STATUS_IGNORE);

	/** Right to left */
	MPI_Sendrecv(s, Ny, MPI_DOUBLE, neighbor[1], 1, s_right, Ny,
			MPI_DOUBLE, neighbor[3], 1, mygrid, MPI_STATUS_IGNORE);

	/** Top to bottom */
	// extract temporary bottom boundary vector
	double temp_bot[Nx];
	cblas_dcopy(Nx, s, Ny, temp_bot, 1);

	MPI_Sendrecv(temp_bot, Nx, MPI_DOUBLE, neighbor[2], 2, s_top, Nx,
			MPI_DOUBLE, neighbor[0], 2, mygrid, MPI_STATUS_IGNORE);

	/** Bottom to top */
	// extract temporary top boundary vector
	double temp_top[Nx];
	cblas_dcopy(Nx, s + Ny-1, Ny, temp_top, 1);
	MPI_Sendrecv(temp_top, Nx, MPI_DOUBLE, neighbor[0], 0, s_bot, Nx,
			MPI_DOUBLE, neighbor[2], 0, mygrid, MPI_STATUS_IGNORE);

}

/**
 * @brief Send and receive vorticity to neighbor domains
 */
void LidDrivenCavity::SendAndRecv_v()
{
	/** Left to right */
	MPI_Sendrecv(v + (Nx-1)*Ny, Ny, MPI_DOUBLE, neighbor[3], 3, v_left, Ny,
			MPI_DOUBLE, neighbor[1], 3, mygrid, MPI_STATUS_IGNORE);

	/** Right to left */
	MPI_Sendrecv(v, Ny, MPI_DOUBLE, neighbor[1], 1, v_right, Ny,
			MPI_DOUBLE, neighbor[3], 1, mygrid, MPI_STATUS_IGNORE);

	/** Top to bottom */
	// extract temporary bottom boundary vector
	double temp_bot[Nx];
	cblas_dcopy(Nx, v, Ny, temp_bot, 1);

	MPI_Sendrecv(temp_bot, Nx, MPI_DOUBLE, neighbor[2], 2, v_top, Nx,
			MPI_DOUBLE, neighbor[0], 2, mygrid, MPI_STATUS_IGNORE);

	/** Bottom to top */
	// extract temporary top boundary vector
	double temp_top[Nx];
	cblas_dcopy(Nx, v + Ny-1, Ny, temp_top, 1);
	MPI_Sendrecv(temp_top, Nx, MPI_DOUBLE, neighbor[0], 0, v_bot, Nx,
			MPI_DOUBLE, neighbor[2], 0, mygrid, MPI_STATUS_IGNORE);

}


/* Solver 
======================================================================== */
/**
 * @brief Solve the whole problem
 */
void LidDrivenCavity::Solve()
{
	LinearMatrices();
	Initialise();
	double t = 0.0;
	while (t < T)
	{
		if (rank == 0) cout << "t = " << t << endl;
		VorticityBCs();
		VorticityInterior();

		SendAndRecv_v();
		VorticityUpdate();
		
		for (int i = 0; i < 5; i++)
		{
			SolvePoisson();
			SendAndRecv_s();
		}
		t += dt;
	}
}

/* Output 
======================================================================== */
/**
 * @brief Output the result in each process
 * @param Px	number of processes in x direction
 * @param Py	number of processes in y direction
 * @param Lx	length of domain in x direction	
 * @param Ly	length of domain in y direction	
 */
void LidDrivenCavity::Output(int Px, int Py, double Lx, double Ly)
{
	for (int k = 0; k < Px *Py; k++)
	{
		if (k == rank)
		{	
			/** Write header and corner values */
			if (k == 0)
			{
				ofstream vOut("Output.txt", ios::out | ios::trunc);
				double x, y;
				vOut << setw(12) << "x" << setw(12) << "y" << setw(20) << "vorticity"
				<< setw(20) << "stream function" << endl;

				vOut << setw(12) << 0.0 << setw(12) << 0.0 << setw(20) << 0.0 
					<< setw(20) << 0.0 << endl;
				vOut << setw(12) << 0.0 << setw(12) << Ly << setw(20) << 0.0 
					<< setw(20) << 0.0 << endl;
				vOut << setw(12) << Lx << setw(12) << 0.0 << setw(20) << 0.0 
					<< setw(20) << 0.0 << endl;
				vOut << setw(12) << Lx << setw(12) << Ly << setw(20) << 0.0 
					<< setw(20) << 0.0 << endl;
				vOut.close();
			}

			ofstream vOut("Output.txt", ios::out | ios::app);
			double x, y;
			/** Write interior values */
			vOut.precision(6);
			for (int i = 0; i < Nx; i++)
			{
				for (int j = 0; j < Ny; j++)
				{	
					x = start[0] + i*dx;
					y = start[1] + j*dy;
					vOut << setw(12) << x << setw(12) << y << setw(20) << v[i*Ny + j] 
						<< setw(20) << s[i*Ny + j] << endl;
				}
			}
			/** Write boundary values */
			// top 
			if (coords[0] == Py - 1)
			{
				y = Ly;
				for (int i = 0; i < Nx; i++)
				{
					x = start[0] + i*dx;
					vOut << setw(12) << x << setw(12) << y << setw(20) << v_top[i] 
					<< setw(20) << 0.0 << endl;
				}		
			}
			// left
			if (coords[1] == 0)
			{
				x = 0;
				for (int j = 0; j < Ny; j++)
				{
					y = start[1] + j*dy;
					vOut << setw(12) << x << setw(12) << y << setw(20) << v_left[j] 
					<< setw(20) << 0.0 << endl;
				}		
			}
			// bottom
			if (coords[0] == 0)
			{
				y = 0;
				for (int i = 0; i < Nx; i++)
				{
					x = start[0] + i*dx;
					vOut << setw(12) << x << setw(12) << y << setw(20) << v_bot[i] 
					<< setw(20) << 0.0 << endl;
				}		
			}
			// left
			if (coords[1] == Px - 1)
			{
				x = Lx;
				for (int j = 0; j < Ny; j++)
				{
					y = start[1] + j*dy;
					vOut << setw(12) << x << setw(12) << y << setw(20) << v_right[j] 
					<< setw(20) << 0.0 << endl;
				}		
			}
			vOut.close();
		}
		MPI_Barrier(mygrid);
	}
}
