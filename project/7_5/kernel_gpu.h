
//-----------------------------------------------------------------------------------------------------------------------------------
//-------Routine for making double precision of atomic varible - only for the summation integration of flux (name: index_d)----------
__device__ double atomicAdd(double* address, double val)
{
	unsigned long long int* address_as_ull =
		(unsigned long long int*)address;
	unsigned long long int old = *address_as_ull, assumed;
	do {
		assumed = old;
		old = atomicCAS(address_as_ull, assumed,
			__double_as_longlong(val +
			__longlong_as_double(assumed)));
	} while (assumed != old);
	return __longlong_as_double(old);
}
//----------------------------------------------------------------------------------------------------------------------------------

//-----------------------------------Karnel for calculating Z direction-------------------------------------------------------------
__global__ void addKernel1(int NX0, int NY0, int niter, real *u, real *u1, real * const __restrict__ X, real * const __restrict__ Z)
{
	extern __shared__  real coef[];//shared memory for each block that solved with PCR algo.
	const unsigned int NX = NX0;  const unsigned int NY = NY0;//number of elements
	const unsigned	int Rn = NX - 1; const unsigned  int Zn = NY - 1; //end indexes
	
	//will be used to save coefficients in the blocks.
	real* a = &coef[0], *c = &coef[NX],* d = &coef[2 * NX];
	//temporal coefficient for each thread.
	real at = 0.0, bt = 1.0, ct = 0, dt = 0, bti;	
	  
	int n = threadIdx.x; //threads in a block are assigned to one column in the concentration array.
	int m = blockIdx.x; //the block index is the index of the row in the concentration array.

	unsigned int nNX = n*NX;//line counter in array for itterating over a pseudo 2D array.

	// defines the difference above and below  and the difference to the right and left of a certain point in the grid X,Z. 
	//--------------------------------------------------------------------------------------------------------------------
	real dZ_ms, dZ_ps, dZ_msps, dZ_psms; 
	real dR_ms, dR_ps, dR_msps, dR_psms; 

	if (m > 0 && m < Rn)//everywhere beside the left and right edges
	{
		dR_ms = X[m] - X[m - 1]; dR_ps = X[m + 1] - X[m];
		dR_msps = dR_ms*dR_ms + dR_ps*dR_ms; dR_psms = dR_ps*dR_ps + dR_ps*dR_ms; //also reduce arithmetic operations later.
		if (n > 0 && n < Zn)//everywhere beside the edges
		{
			dZ_ms = Z[n] - Z[n - 1]; dZ_ps = Z[n + 1] - Z[n];
			dZ_msps = dZ_ms*dZ_ms + dZ_ps*dZ_ms; dZ_psms = dZ_ps*dZ_ps + dZ_ps*dZ_ms;//also reduce arithmetic operations later.
		}
		else
		{
				dZ_ps = Z[1] - Z[0];
		}			
	   //--------------------------------------------------------------------------------------------------------------------
	

		//--------------------------------------setting coefficients--------------------------------------------------------
		if (n > 0 && n < Zn)// medium (rate determining step of performance since it includes many arithmetic operations, espacially setting the delta coefficient.
			{
				at = 2.0 / (dZ_msps);
				bt = -2.0 / (dZ_msps)-2 / (dZ_psms)-2 / (dt_d);
				ct = 2.0 / (dZ_psms);
				dt = (-2.0 / (dR_msps)+(1.0 / ((dR_ps + dR_ms)* X[m])))*u[m - 1 + nNX] + (2.0 / (dR_msps)+2 / (dR_psms)-2 / (dt_d))*u[m + nNX] + (-2.0 / (dR_psms)-(1.0 / ((dR_ps + dR_ms)* X[m])))*u[m + 1 + nNX];
			}
	

			if (n == Zn) // boundary condition of the edge (bulk concentration)
			{
				at = 0;
				bt = 1;
				ct = 0;
				dt = 1;
			}
		
			if (n == 0 && m <= Rd_d)// boundary condition of electrode, kinetics using butler volmer model.
			{
				at= 0;
				bt = 1 + K_d*(dZ_ps)*(exp(a_d*V_d))*(1 + exp(-V_d));
				ct = -1;
				dt = K_d*(dZ_ps)*(exp(a_d*V_d))*(exp(-V_d));				
			}
	
			if (n == 0 && m > Rd_d)// insulation/symmetry boundray condition on the supporting sheath and on the axial symmetry.
			{
				at = 0;
				bt =-1;
				ct = 1;
				dt = 0;
			}
			__syncthreads();//waiting the branch to finish in the block 
			
			//normailzation to get beta equals 1.
			bti = 1 / bt;
			at = at*bti;
			ct = ct*bti;
			dt = dt*bti;
	
			a[n] = at;
			c[n] = ct;
			d[n] = dt;

			//------starting PCR algorithm for Z direction (block syncs are must)-------
			for (int nt = 1; nt < NY; nt = 2 * nt) {
				__syncthreads();
				bt = 1.0;

				if (n - nt >= 0) {
					dt = dt - at*d[n - nt];
					bt = bt - at*c[n - nt];
					at = -at*a[n - nt];
			
				}
				__syncthreads();
				if (n + nt < NY) {
					dt = dt - ct*d[n + nt];
					bt = bt - ct*a[n + nt];
					ct = -ct*c[n + nt];
				}
				__syncthreads();
				bti = 1.0 / bt;
				at = at*bti;
				ct = ct*bti;
				dt = dt*bti;

				a[n] = at;
				c[n] = ct;
				d[n] = dt;
			}
			//--------------------------------------------------------------
		    	__syncthreads(); 
			
			    u1[m + nNX] = dt;// copy concentration solution from final delta to u1 (not to old u since u might be still active in other blocks! see setting delta coefficient above.
		}
				d_index = 0;	//initializing flux calculation, this will be now calculate in the next kernel to ensure that all concentration written to u1.
}
//----------------------------------------------------------------------------------------------------------------------



//-----------------------------------Karnel for calculating R direction and flux summation------------------------------------------
__global__ void addKernel2(int NX0, int NY0, int niter, real *u, real *u1, real * const __restrict__ X, real * const __restrict__ Z)
{	
	extern __shared__  real coef[];	
	real *P = &coef[0];// here we use another way to save coefficeint: a0 b0 c0 d0 a1 b1 c1 d1 .....an bn cn dn. instead of a0 a1 a2...an b0 b1 b2....bn. this saves a bit of performance but not necessary.  
	//temporal coefficient for each thread.
	real at = 0.0, bt = 1.0, ct = 0, dt = 0, bti;
	
	const unsigned int NX = NX0;  const unsigned int NY = NY0;
	const unsigned	int Rn = NX - 1; const unsigned  int Zn = NY - 1;
	 
	int n = threadIdx.x;
	int m = blockIdx.x;

	//finding the flux at each point on the electrode (m<Rd_d) and make integration by intergral summation (tarpez method) over all the electrode.
	if (n == 1 && m < Rd_d){

		atomicAdd(&d_index, (X[m + 1] - X[m])*((u1[m + NX] - u1[m])*X[m] + (u1[m + 1 + NX] - u1[m + 1])*X[m + 1]) / ((Z[1] - Z[0]) * 2));
	}

	// defines the difference above and below  and the difference to the right and left of a certain point in the grid X,Z. 
	//--------------------------------------------------------------------------------------------------------------------
	real dZ_ms, dZ_ps, dZ_msps, dZ_psms;
	real dR_ms, dR_ps, dR_msps, dR_psms;

	if (m > 0 && m <  Zn )
	{
		if (n > 0 && n < Rn)
		{
			dR_ms = X[n] - X[n - 1]; dR_ps = X[n + 1] - X[n];
			dR_msps = dR_ms*dR_ms + dR_ps*dR_ms; dR_psms = dR_ps*dR_ps + dR_ps*dR_ms;
		}
		dZ_ms = Z[m] - Z[m - 1]; dZ_ps = Z[m + 1] - Z[m];
		dZ_msps = dZ_ms*dZ_ms + dZ_ps*dZ_ms; dZ_psms = dZ_ps*dZ_ps + dZ_ps*dZ_ms;

	}	
	//--------------------------------------------------------------------------------------------------------------------

		if (m > 0 && m <  Zn )
		{
			//--------------------------------------setting coefficients-------------------------------------------------
			if (n > 0 && n < Rn)//medium
			{
				at = (1.0 / (dR_ps + dR_ms))*(2 / dR_ms - 1 / X[n]);
				bt = -2.0*(1 / dR_msps+1 / dR_psms+1 / dt_d);
				ct = (1.0 / (dR_ps + dR_ms))*((2 / dR_ps) + (1 / X[n]));
				dt = (-2.0 / dZ_msps)*u1[n + (m - 1)*NX] + (2.0 / dZ_msps+2 / dZ_psms-2 / dt_d)*u1[n + m*NX] + (-2.0 / dZ_psms)*u1[n + (m + 1)*NX];

			}
			else
			{
				if (n == Rn) //bulk boundary conditions.
				{
					at = 0;
					bt = 1;
					ct = 0;
					dt = 1;
				}
				if (n == 0)//axial symmetry boundary conditions.
				{
					at = 0;
					bt =-1;
					ct = 1;
					dt = 0;
				}
			}
			__syncthreads();

			bti = 1 / bt;
			at = at*bti;
			ct = ct*bti;
			dt = dt*bti;
		
			P[3*n] = at;
			P[3 * n + 1] = ct;
			P[3 * n + 2] = dt;

			//------starting PCR algorithm for R direction (block syncs are must)-------
			for (int nt = 1; nt < NX; nt = 2 * nt) {
				int k = 3 * (n - nt);
				int j = 3 * (n + nt);
				__syncthreads();

				bt = 1;

				if (n - nt >= 0) {
					dt = dt - at*P[k + 2];
					bt = bt - at*P[k + 1];
					at = -at*P[k];
				}
				__syncthreads();
				if (n + nt < NX) {
					dt = dt - ct*P[j + 2];
					bt = bt - ct*P[j];
					ct = -ct*P[j + 1];
				}
				__syncthreads();
				bti = 1 / bt;
				at = at*bti;
				ct = ct*bti;
				dt = dt*bti;
		//		__syncthreads();
				P[3*n] = at;
				P[3 * n + 1] = ct;
				P[3 * n + 2] = dt;

			}
			//-----------------------------------------------------------------------
			__syncthreads();
			u[n + m*NX] = dt;//save solution  to the first concentration array.
		}	
}
//----------------------------------------------------------------------------------------------------------------------------


//-----------------------------------Karnel for update parameters-------------------------------------------------------------
__global__ void addKernel3(int NX, int niter, real *J)
{
		int n = threadIdx.x;
    	int m = blockIdx.x;
	if (m == 0 && n == 0)
	{
		J[n_d] = 2 * PI*d_index; //save result to global.
		if (n_d > 0)
			J[n_d + niter] = V_d; //save potential to global.
		else
			J[n_d + niter] = 0;

		dt_d = dt_d*gamma_d;//amperometry expand time step (in voltammetry gamma_d=1).
		n_d++;
		if (V_d >= Vf_d - dV_d) { S_d = -1; } //voltammetry change directions.
		if (V_d <= Vi_d + dV_d) { S_d = +1; }
		V_d = V_d + S_d*dV_d;
	}
}
//----------------------------------------------------------------------------------------------------------------------------

