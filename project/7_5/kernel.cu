#include <iostream>
#include <vector>
#include <time.h> 
#include <cuda.h>//need cuda.
#include "cuda_runtime.h"//need cuda.
#include "device_launch_parameters.h"//need cuda.
#include <device_functions.h>//need cuda. 
#include "classes\simec.h" //attached classes for simulation. 

using std::vector;
# define PI  3.14159265358979323846 /* pi */

# define cv    //choose "cv" for voltametry or "ca" for amperometry.
# define inputxt
typedef double real; //define all varibles for double precision.

//----------------------------------------------------------------------------------------------
//prototype for sending to help function.
cudaError_t addWithCuda(real *u, real *z, real *r, real *j, unsigned int size, unsigned int sizeY, unsigned int blocksize, unsigned int iter, int shared_mem, real dT, real gammaT, real Zn, real Rn, real Rd, real dV, real Vf, real Vi, real K_bv, real a_bv);
//allocate various memories on the gpu device.
__constant__ real gamma_d, dV_d, Vf_d, Vi_d, K_d, a_d, Rd_d; //send constants to device : time expansion,delta V,max V, min V,kinetic const,alfa coef, and electrode size (respectively).
__device__ unsigned int  n_d = 0 /* time step solution */, double V_d = -20/* initial potential (only voltametry)*/, int S_d = -1/* votametry direction scan */, double d_index = 0/* atomic varible for current*/, dt_d/* delta time */;//send varibles to device: atomic addition, 
//include file of the kernels
#include "kernel_gpu.h" // gpu 4 kernels. 
//----------------------------------------------------------------------------------------------


int main()
{
	//Parameters voltammetry-------------------------------
#ifdef cv
		real SR = 100000.0, Vi = -10, Vf =10, V = -10, K_bv = 1e10, a_bv = 0.5;  int S = -1; //voltammetry parameters: scan rate, min V, max V, initial V, kinetic konst, alpha coef, dierection of scan. 
		real dV = 0.05; int niter = 2.0*(Vf - Vi) / dV; /* number of iterations*/;// dV small for accuracy, very performance costely but must be decreased for low scan rates.  
		real T = 0/* time */, dt = dV / SR, dt_init = dt, Da = 1,/*(diffusion coef), is always 1 if dimensionless is applied*/  gammaT = 1.0;/*no expansion time for voltammetry,thus 1.0*/;
		real dR = 1e-4, gammaR = 1.05, dZ = 1e-4, gammaZ = 1.05; //mesh parameters: minimum distances and expansion factors.
		real Rmax = 500, Zmax = 500;/* Max 2D cell sizes in dimensinless units */ 
#ifdef inputxt//reading parameters from input.txt for a test version.
		cout << "Taking these data from input file for a simple test-run:"<<'\n';
		std::fstream input("input.txt", std::ios_base::in);
		input >> dR >> gammaR >> dZ >> gammaZ >> K_bv >> SR >> niter;
		if (niter==0) niter = 2.0*(Vf - Vi) / dV;
		printf("%s%f\t%s%f\t%s%f\t%s%f\t%s%f\t%s%f\t%s%d\n", "dR=", dR, "gammaR=", gammaR, "dZ =", dZ, "gammaZ =", gammaZ, "Kinetic const=", K_bv, "scan rate=", SR, "n. of time steps=",niter);
#endif 

#endif 	
   //Parameters amperometry-------------------------------
#ifdef ca
		real dV = 0.00;//no change in potential in amperometry 
		real T = 0, dt_init = 1e-5/* intial small time step */, dt = dt_init, Da = 1;
		real dR = 0.125 / 64, gammaR = 1.2, dZ = 0.125 / 64, gammaZ = 1.2, gammaT = 1.05 /* expending time step */;
		real Rmax = 40, Zmax = 40; int niter = 500;
		real SR = 20.0, Vi = 20, Vf = 20, V = 20, K_bv = 1e8, a_bv = 0.5;  int S = -1;/* not important for amperometry unless potential is not high enough*/
#endif
		
	//------------MESH(using the Grid class)--------------------
	//using a mesh class to create a grid on cpu.
	//R vector 
	Grid Rgrid(dR, gammaR); //create grid for z direction.
	Rgrid.new_node(1); //add new node of electrode edge at point 1.0.
	Rgrid.open_node(Rmax);//add another node for end, open node means without high density of points around the node.
	int Rd = Rgrid.get_node(1) - 1, Rn = Rgrid.get_node(2) - 1;//get the index of Rd (electrode)and of the cell size, Rn;
	//Z vector
	Grid Zgrid(dZ, gammaZ);
	Zgrid.new_node(1.0);//unnecessary point, but makes Z grid and R grid equals.
	Zgrid.open_node(Zmax);
	int Zn = Zgrid.get_node(2) - 1;
	if (Rn < Zn) { std::cout << "sorry, only works for Rn>Zn"; return 0; }//warning: because of the GPU algorithm made only for the case where Rn>Zn.
	//---------------------------------------------------------
	
	specie A(Rgrid, Zgrid, 1.0);
	PreTomas phys(Rgrid, Zgrid);
	
	//----------file and clock------------------------------------
	const char* PATH1("CPU-GPU.txt");  //saving output comparison cpu to gpu.
	ofstream myfile;  myfile.open(PATH1);
	cout << '\n' << "matrix size" << '\t' << Rn << '\t' << Zn << '\n';
	//-----------------------------------------------------------

	
	
	//*****************************************************************************************************************
	////---------------Simulation on the cpu--------------------
	cout << "press enter to start CPU simulation"<<'\n';
	cin.get();
	//-----------World-----------------------------------------
	int Domain = A.world.make_domain(0, 0, Rn, Zn);//assign domain space to 2D grid
	int Electrode = A.world.make_edge(0, 0, Rd, 0);//assign  electrode boundary to 2D grid
	int BulkZn = A.world.make_edge(0, Zn, Rn - 1, Zn);//assign  cell edge z to 2D grid
	int BulkRn = A.world.make_edge(Rn, 0, Rn, Zn);//assign  cell edge r to 2D grid
	int Symmetry = A.world.make_edge(0, 0, 0, Zn);//assign  axes symmetry edge to 2D grid
	int Sheath = A.world.make_edge(Rd + 1, 0, Rn, 0); //assign supporting sheath next to the electrode to 2D grid
	//--------------------------------------------------------------
	//---cpu variables----------------------------------------------
	real *h_u/*concentration data */, *h_j /* current data */, *h_z /* grid z */, *h_r /* grid r */;
	real *j_cpu = (real *)malloc(sizeof(real)*niter * 2) /* allocate for current respons solution on cpu */;
	//-------------------------------------------------------------
	clock_t begin = clock();//start clock for the cpu calculation.
		for (int iter = 0; iter < niter; iter++)
		{
			//------------------define Z elements-----------------
			for (int n = 0; n <= Rn - 1; n++)
			{
				for (int m = 0; m <= Zn; m++)
				{	   
					//domain equations Z directon, assign coefficients: alpha, beta, gamma, delta.
					if (A.world.is(n, m, Domain) || A.world.is(n, m, Symmetry) || A.world.is(n, m, BulkRn))
					{
							*A.set_az(n, m) = phys.Z_a(n, m); //alpha 
						*A.set_bz(n, m) = phys.Z_b(n, m, dt, Da);     //beta
						 *A.set_gz(n, m) = phys.Z_g(n, m); //gamma
						*A.set_dz(n, m) = phys.Z_d(&A, n, m, dt, Da); //delta
					}
					//boundary conditions Z directon
					if (A.world.is(n, m, Electrode)) phys.butlervolmer_boundary1(&A, 'Z', n, m, V, a_bv, K_bv);//alpha/beta, gamma,and delta on electrode.
					if (A.world.is(n, m, Sheath)) phys.insulation(&A, 'Z', n, m);//alpha/beta, gamma,and delta on supporting sheath.
					if (A.world.is(n, m, BulkZn)) phys.bulk(&A, 'Z', n, m);//alpha/beta, gamma,and delta on cell edges.
				}
			}
	     //----------solve Z elements using thomas algorithm (half time step)---------
			for (int n = 1; n <= Rn - 1; n++)
			{
				tomas * tomA = new tomas(&A, 'Z', n); //create object of class thomas
				(*tomA).modify_g(); (*tomA).modify_d(); (*tomA).solve_c(); //solve thomas with modify  gamma -> modify delta -> backward sloving.
				delete tomA;
			}
         //---------------------------------------------------------------------------
		
		//----------Flux-time-voltage- output---------------------------
			printf("\r %d %s",  int(100*(iter+1)/niter), "%");
		    j_cpu[iter] = 2 * PI*A.current(); j_cpu[niter + iter] = T; //save to file (2pi is to get the flux from all the disc).
		//--------------------------------------------------------------

		//------------------define R elements-----------------
		for (int n = 0; n <= Rn; n++)
		 {
			for (int m = 0; m < Zn; m++)
			{	//domain equations Z directon, assign coefficients: alpha, beta, gamma, delta.
				if (A.world.is(n, m, Domain) || A.world.is(n, m, Electrode) || A.world.is(n, m, Sheath))
				{
					 *A.set_ar(n, m) = phys.R_a(n, m); //alpha 
					*A.set_br(n, m) = phys.R_b(n, m, dt, Da);//beta
					 *A.set_gr(n, m) = phys.R_g(n, m);//gamma
					*A.set_dr(n, m) = phys.R_d(&A, n, m, dt, Da);//delta
				}
					if (A.world.is(n, m, BulkRn)) phys.bulk(&A, 'R', n, m);//alpha/beta, gamma, delta for the edges of the cell
					if (A.world.is(n, m, Symmetry)) phys.insulation(&A, 'R', n, m);//alpha/beta, gamma, delta for axial symmetry at r=0
			}
		}
		
		//----------solve R direction-----------------------------------------------
		for (int m = 1; m < Zn; m++)
		{
			tomas * tomA = new tomas(&A, 'R', m);
			(*tomA).modify_g(); (*tomA).modify_d(); (*tomA).solve_c();
			delete tomA;
		}
		//--------------------------------------------------------------------------

	
		if (V >= Vf - dV) { S = -1; }//change direction of scan if gets to max value
		if (V <= Vi + dV) { S = 1; } //change direction of scan if gets to min value
	
		V = V + S*dV; //scan to the next potentnial step only relevant for voltammetry
		dt *= gammaT;// expand the dt only relvant for amperometery
 		T = T + dt; // the next time step
	}

    //---------------------------------------------------------------------------	
	clock_t end = clock(); //finish clock for performance estimation
	//---------------------------------------------------------------------------
	//*****************************************************************************************************************

	
	
	//-------prepare again concentration and grid matrix for gpu calculations------------------------------------------ 	
	//re-assigning NX and Ny as the number of points on the grid for GPU
	    int NX = Rn + 1, BLK_size =  (Rn + 1);
		int NY = Zn + 1;
		// allocate memories on host
		h_u = (real *)malloc(sizeof(real)*NX*BLK_size); //concentration size. (NX>NY is necessary).
		h_z = (real *)malloc(sizeof(real)*NY); //later for sending to the device the grid points on z direction .
		h_r = (real *)malloc(sizeof(real)*NX); //later for sending to the device the grid points on r direction.
		h_j = (real *)malloc(sizeof(real)*niter*2);//for saving the flux solution.
		int shared_mem_size = 3 * sizeof(real) * NX; //define alfa  gamma delta for each block later will allocate the shared memory,  beta normalized to 1 so no need to save.
		vector<double>Z = Zgrid.get_vector();
		vector<double>R = Rgrid.get_vector();
		for (int i = 0; i < NX ; i++)
		{
			for (int j = 0; j < BLK_size; j++)
				h_u[i + j*NX] = 1;
			
			h_r[i] = R[i];
			if (i <= Zn)
				h_z[i] = Z[i];
		}
		//------------------------------------------------------------------------------------------------------------ 


  //**********************************send stuff to helper function to use cuda *************************************	
		cudaError_t cudaStatus = addWithCuda(h_u, h_z, h_r, h_j, NX,NY, BLK_size, niter, shared_mem_size, dt_init, gammaT, Zn, Rn, Rd, dV,  Vf,  Vi,  K_bv,  a_bv);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addWithCuda failed!");
        return 1;
    }
	cout << '\n' << "Finished cuda." ;
	cout << "\n" << "\n" << "CPU CV:" << (end - begin) << "(ms)" << "\n" << "GPU CV:" << h_j[niter] << "(ms)" << "\n" << (end - begin) / h_j[niter] << " ratio" << "\n" << '\t' << "\n";
	cout << '\n' << "press enter to print and save: t, V,flux, and deviation from cpu(cpu-gpu)"<<'\n' ;
	cin.get();
  //*****************************************************************************************************************	

	//-----------------print voltammetry/amperometry to file---------------------------------------
	T = 0; dt = dt_init; real error=0.0;
	for (int i = 1; i < niter; i++)
	{	
		printf("%f %f %f %f %f \n", T, h_j[niter + i], h_j[i], j_cpu[i], j_cpu[i] - h_j[i]);
		myfile << T << '\t' << h_j[niter + i] << '\t' << h_j[i] << '\t' << j_cpu[i] << '\n';
		dt *= gammaT;T = T + dt;
		if (error < j_cpu[i] - h_j[i])error = j_cpu[i] - h_j[i];
	}
	//---------------------------------------------------------------------------------------------

	// cudaDeviceReset must be called before exiting in order for profiling and
    // tracing tools such as Nsight and Visual Profiler to show complete traces.
    cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceReset failed!");
        return 1;
    }
	cout << "\n" << "CPU CV:" << (end - begin) << "(ms)" << "\n" << "GPU CV:" << h_j[niter] << "(ms)" << "\n" << " ratio" << '\t' << (end - begin) / h_j[niter] << "\n" << "deviation" << '\t' << error << "\n";
	cin.get();
	myfile.close();
    return 0;
 }
//****************************************************************************************************************************************

 
 
 
 
 
 
 
 
 //************************************ Helper function for using CUDA*********************************************************************
	 cudaError_t addWithCuda(real *u, real *z, real *r, real *j, unsigned int size, unsigned int sizeY, unsigned int blocksize, unsigned int iter, int shared_mem, real dt_init, real gammaT, real Zn, real Rn, real Rd, real dV, real Vf, real Vi, real K_bv, real a_bv)
{
	// Choose which GPU to run on, change this on a multi-GPU system.
	cudaError_t cudaStatus;
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA- capable GPU installed?");
		goto Error;
	}
	
	
	// Allocate GPU buffers for all vectors vectors (two input, one output)   
	real *dev_u, *dev_r, *dev_z, *dev_j, *dev_u1;
	 
	cudaStatus = cudaMalloc((void**)&dev_u, sizeof(real)*size*blocksize);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	cudaStatus = cudaMalloc((void**)&dev_u1, sizeof(real)*size*blocksize);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	cudaStatus = cudaMalloc((void**)&dev_r, sizeof(real)*(Rn + 1));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	cudaStatus = cudaMalloc((void**)&dev_z, sizeof(real)*(Zn + 1));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	cudaStatus = cudaMalloc((void**)&dev_j, sizeof(real)*iter*2);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}



	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_u, u, sizeof(real)*size*blocksize, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}
	cudaStatus = cudaMemcpy(dev_u1, u, sizeof(real)*size*blocksize, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}
	cudaStatus = cudaMemcpy(dev_r, r, sizeof(real)*(Rn + 1), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}
	cudaStatus = cudaMemcpy(dev_z, z, sizeof(real)*(Zn + 1), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}
	cudaStatus = cudaMemcpy(dev_j, j, iter * sizeof(real), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}

	
	// Copy global varibles to device.
	cudaMemcpyToSymbol(Vf_d, &Vf, sizeof(real));
	cudaMemcpyToSymbol(Vi_d, &Vi, sizeof(real));
	cudaMemcpyToSymbol(a_d, &a_bv, sizeof(real));
	cudaMemcpyToSymbol(K_d, &K_bv, sizeof(real));
	cudaMemcpyToSymbol(dV_d, &dV, sizeof(real));
	cudaMemcpyToSymbol(Rd_d, &Rd, sizeof(real));
	cudaMemcpyToSymbol(dt_d, &dt_init, sizeof(real));
	cudaMemcpyToSymbol(gamma_d, &gammaT, sizeof(real));
	unsigned int delta = 0;
	int current0 = 0;
	unsigned int n_0= 0;
	int S_h = -1;
	cudaMemcpyToSymbol(V_d, &Vi, sizeof(real), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(S_d, &S_h, sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(n_d, &n_0, sizeof(unsigned int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_index, &current0, sizeof(double), 0, cudaMemcpyHostToDevice);
	
	
	cout << '\n'<< "Now with GPU, press enter to start" << '\n';
	cin.get();
	float milli;
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	
	cudaEventRecord(start, 0);

    // Launch the kernels on the GPU.
	for (int N = 0; N < iter; N++)
	{
		addKernel1 << <blocksize, size, shared_mem >> >(size,sizeY, iter, dev_u, dev_u1, dev_r, dev_z);
		addKernel2 << <blocksize, size, shared_mem >> >(size, sizeY, iter, dev_u, dev_u1, dev_r, dev_z);
		addKernel3 << <blocksize, size, shared_mem >> >(size, iter, dev_j);
		printf("\r %d %s", int(100 * (N + 1) / iter), "%");
	}
	cudaEventRecord(stop, 0);
	
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&milli, start, stop);


	
	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
		goto Error;
	}

	// Copy output vector from GPU buffer to host memory.
	cudaStatus = cudaMemcpy(j, dev_j, sizeof(real)*iter*2, cudaMemcpyDeviceToHost);
	j[iter] = milli;
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}
	
Error:
	cudaFree(dev_u);
	cudaFree(dev_u1);
	cudaFree(dev_z);
	cudaFree(dev_r);
	cudaFree(dev_j);
	return cudaStatus;
  }
  //****************************************************************************************************************************************