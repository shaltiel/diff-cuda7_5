#include <vector>

using namespace std;

class specie 
{	
	friend class Grid;
	friend class World;
	//friend class PreTomas;
public:
	vector<double> ar, br, gr, dr, az, bz, gz, dz, C, C0;
	int Zsize, Rsize; 
	double init_c=1;
	Grid Rg, Zg; 
	World world;
	
	//PreTomas phys;
public:

	specie(Grid, Grid, double);

	~specie(){};

	void set_grid(Grid, Grid);
	double current();
	double current0();
	double current( int);
	double current_density(int, int);
	double current_range(int N_int, int N_end, char U);
	double current_range_Z(int N_int, int N_end, char U);
	int where_am_I(int n, int m) { return (n + m*Rsize); }
	double* set_ar(int n, int m) { return &ar[n + m*Rsize];}
	double* set_br(int n, int m) { return &br[n + m*Rsize]; }
	double* set_gr(int n, int m) { return &gr[n + m*Rsize];}
	double* set_dr(int n, int m) { return &dr[n + m*Rsize];}
	double* set_C(int n, int m) {  return &C[n + m*Rsize];}
	double* set_C0(int n, int m) { return &C0[n + m*Rsize]; }
	double* set_az(int n, int m) { return &az[n + m*Rsize]; }
	double* set_bz(int n, int m) { return &bz[n + m*Rsize]; }
	double* set_gz(int n, int m) { return &gz[n + m*Rsize]; }
	double* set_dz(int n, int m) { return &dz[n + m*Rsize]; }

	double get_ar(int n, int m)  { return ar[n + m*Rsize]; }
	double get_br(int n, int m)  { return br[n + m*Rsize]; }
	double get_gr(int n, int m)  { return gr[n + m*Rsize]; }
	double get_dr(int n, int m)  { return dr[n + m*Rsize]; }
	double get_az(int n, int m)  { return az[n + m*Rsize]; }
	double get_bz(int n, int m)  { return bz[n + m*Rsize]; }
	double get_gz(int n, int m)  { return gz[n + m*Rsize]; }
	double get_dz(int n, int m)  { return dz[n + m*Rsize]; }
	double get_C(int n, int m)  { return C[n + m*Rsize]; }
	double get_C0(int n, int m)  { return C0[n + m*Rsize]; }
	double get_sizear(){  return ar.size(); } 
};


specie::specie(Grid GrdR, Grid GrdZ, double initial_c)
{
	init_c = initial_c;
	Rsize = GrdR.V.size();
	Rg = GrdR;
	
	world.setworld(GrdR, GrdZ);
	
	Zsize = GrdZ.V.size();
	Zg = GrdZ;

	//PreTomas ph_temp(Rg,Zg, this);
	//phys = ph_temp;

	ar.resize(Rsize*Zsize, 0);
	br.resize(Rsize*Zsize, 0);
	gr.resize(Rsize*Zsize, 0);
	dr.resize(Rsize*Zsize, 0);
	
	az.resize(Rsize*Zsize, 0);
	bz.resize(Rsize*Zsize, 0);
	gz.resize(Rsize*Zsize, 0);
	dz.resize(Rsize*Zsize, 0);
	
	C.resize(Rsize*Zsize, initial_c);
	C0.resize(Rsize*Zsize, initial_c);
}


//function implemetation 
void specie::set_grid(Grid GrdR, Grid GrdZ)
{
	 Rsize = GrdR.V.size();
	 
	 Zsize = GrdZ.V.size();	

	 ar.resize(Rsize*Zsize,0);
	 br.resize(Rsize*Zsize,0);
	 gr.resize(Rsize*Zsize,0);
	 dr.resize(Rsize*Zsize,0);

	 az.resize(Rsize*Zsize, 0);
	 bz.resize(Rsize*Zsize, 0);
	 gz.resize(Rsize*Zsize, 0);
	 dz.resize(Rsize*Zsize, 0);
	 C.resize(Rsize*Zsize, init_c);
	 C0.resize(Rsize*Zsize, init_c);
}
//all current functions are actually flux and not current
double specie::current()
{
	double J = 0.0;
	//----flux integration-----
	for (int n = 0; n < Rg.get_node(1) - 1; n++)
	{

		J = J + (Rg.V[n + 1] - Rg.V[n])*(((*this).get_C(n, 1) - (*this).get_C(n, 0))*Rg.V[n] + ((*this).get_C(n + 1, 1) - (*this).get_C(n + 1, 0))*Rg.V[n + 1]) / ((Zg.V[1] - Zg.V[0]) * 2);
		
	}
	return J;
}
double specie::current0()
{
	double J = 0.0;
	//----flux integration-----
	for (int n = 0; n < Rg.get_node(1) - 1; n++)
	{

		J = J + (Rg.V[n + 1] - Rg.V[n])*(((*this).get_C0(n, 1) - (*this).get_C0(n, 0))*Rg.V[n] + ((*this).get_C0(n + 1, 1) - (*this).get_C0(n + 1, 0))*Rg.V[n + 1]) / ((Zg.V[1] - Zg.V[0]) * 2);

	}
	return J;
}

	double specie::current(int N)
	{
		double J = 0.0;
		//----flux integration by layer---------
		for (int n = 0; n < Rg.get_node(1) - 1; n++)
		{
			J = J + (Rg.V[n + 1] - Rg.V[n])*(((*this).get_C(n, N + 1) - (*this).get_C(n, N))*Rg.V[n] + ((*this).get_C(n + 1, N + 1) - (*this).get_C(n + 1, N))*Rg.V[n + 1]) / ((Zg.V[N+1] - Zg.V[N]) * 2);

		}
		return J;
	}
	double specie::current_range(int N_int, int N_end,char U)
	{
		double J = 0.0;
		if (U == 'R')//radial coordinates
		{
		
			//----flux integration by layer---------
			for (int n = N_int; n < N_end; n++)
			{
				J = J + (Rg.V[n + 1] - Rg.V[n])*(((*this).get_C(n, 1) - (*this).get_C(n, 0))*Rg.V[n] + ((*this).get_C(n + 1, 1) - (*this).get_C(n + 1, 0))*Rg.V[n + 1]) / ((Zg.V[1] - Zg.V[0]) * 2);

			}
		}
		if (U == 'C') //caterzian cube
		{
			
			//----flux integration by layer---------
			for (int n = N_int; n < N_end; n++)
			{
				J = J + (Rg.V[n + 1] - Rg.V[n])*(((*this).get_C(n, 1) - (*this).get_C(n, 0))+ ((*this).get_C(n + 1, 1) - (*this).get_C(n + 1, 0))) / ((Zg.V[1] - Zg.V[0]) * 2);

			}
		}
		return J;
	}
	double specie::current_range_Z(int N_int, int N_end,char U)
	{
		double J = 0.0;
		
		if (U == 'R') //caterzian cube
		{
			for (int n = N_int; n < N_end; n++)
			{
				J = J + (Zg.V[n + 1] - Zg.V[n])*(((*this).get_C(1, n) - (*this).get_C(0, n))*Zg.V[n] + ((*this).get_C(1, n + 1) - (*this).get_C(0, n + 1))*Zg.V[n + 1]) / ((Rg.V[1] - Rg.V[0]) * 2);

			}
		}
		if (U == 'C') //caterzian cube
		{
			for (int n = N_int; n < N_end; n++)
			{
				J = J + (Zg.V[n + 1] - Zg.V[n])*(((*this).get_C(1, n) - (*this).get_C(0, n)) + ((*this).get_C(1, n + 1) - (*this).get_C(0, n + 1))) / ((Rg.V[1] - Rg.V[0]) * 2);

			}
		}
		return J;
	}

	double specie::current_density(int N, int n) //layer number (N) and node number (n)
	{
		double J = 0.0;
		
	
			J = ((*this).get_C(n, N + 1) - (*this).get_C(n, N))/ ((Zg.V[N+1] - Zg.V[N]) );

		return J;
	}
	//for (int m = 0; m < 1; m++)//for calculations of current in different layers
	//{
	//	J_s[m] = 0;
	//	for (int n = node1 + 1; n < Re_cell; n++)
	//	{
	//		J_s[n] = (R[n + 1] - R[n])*((C[cA + n - (m + 1)*raw] - C[cA + n - m*raw])*R[n] + (C[cA + n + 1 - (m + 1)*raw] - C[cA + n + 1 - m*raw])*R[n + 1]) / ((Z[m + 1] - Z[m]) * 2);
	//	}
	//}

