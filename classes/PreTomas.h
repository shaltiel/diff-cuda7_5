#include <vector>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <fstream>

using namespace std;

class PreTomas
{
	friend class Grid;
	friend class specie;
private:
	Grid R; Grid Z;
	vector<double> *pRm; vector<double>* pRp;
	vector<double> *pZm; vector<double>* pZp ;
	
public:
	PreTomas();

	PreTomas(Grid, Grid);
	

	 
	//double cell_ap(int n , int m)       { return (1.0 / ((*pRp)[n] + (*pRm)[n]))*((2 / (*pRm)[n]) - (1 / R.V[n]));}
	
	//physics-operation on specie -domain for a cylinder coordinates
	double R_a(int n , int m)                       { return (1.0 / (R.dp[n] + R.dm[n]))*((2 / R.dm[n]) - (1 / R.V[n]));}
	double R_b(int n , int m, double dT, double D)  { return (2.0 / (R.dp[n] + R.dm[n]))*((-1 / R.dp[n]) + (-1 / R.dm[n])) - 2 /(D*dT); }
	double R_g(int n , int m)                       { return (1.0 / (R.dp[n] + R.dm[n]))*((2 / R.dp[n]) + (1 / R.V[n])); }
	double R_d(specie *, int, int, double, double);

	double Z_a(int n, int m)                       { return  2.0 / ((Z.dp[m] + Z.dm[m])*Z.dm[m]); }
	double Z_b(int n, int m, double dT, double D)  { return  (2.0 / (Z.dp[m] + Z.dm[m]))*((-1 / Z.dp[m]) + (-1 / Z.dm[m])) - 2 / (D*dT); }
	double Z_g(int n, int m)                       { return  2.0 / ((Z.dp[m] + Z.dm[m])*Z.dp[m]); }
	double Z_d(specie *, int, int, double, double);
	
	//physics-operation on specie -domain for a caterzian coordinates
	double X_a(int n, int m)                       { return  2.0 / ((R.dp[n] + R.dm[n])*R.dm[n]); }
	double X_b(int n, int m, double dT, double D)  { return (2.0 / (R.dp[n] + R.dm[n]))*((-1 / R.dp[n]) + (-1 / R.dm[n])) - 2 / (D*dT); }
	double X_g(int n, int m)                      { return  2.0 / ((R.dp[n] + R.dm[n])*R.dp[n]); }
	double X_d(specie *, int, int, double, double);

	double Y_a(int n, int m)                       { return  2.0 / ((Z.dp[m] + Z.dm[m])*Z.dm[m]); }
	double Y_b(int n, int m, double dT, double D)  { return  (2.0 / (Z.dp[m] + Z.dm[m]))*((-1 / Z.dp[m]) + (-1 / Z.dm[m])) - 2 / (D*dT); }
	double Y_g(int n, int m)                       { return  2.0 / ((Z.dp[m] + Z.dm[m])*Z.dp[m]); }
	double Y_d(specie *, int, int, double, double);

	//1D
	double plane1D_ar(int n, int m)                       { return (-2.0 / (R.dp[n] * R.dm[n] + R.dm[n] * R.dm[n])); }
	double plane1D_br(int n, int m, double dT, double D)  { return (2.0 / (R.dp[n] * R.dm[n] + R.dm[n] * R.dm[n])) + (2.0 / (R.dp[n] * R.dm[n] + R.dp[n] * R.dp[n])) + 1 / (D*dT); }
	double plane1D_gr(int n, int m)                       { return (-2.0 / (R.dp[n] * R.dm[n] + R.dp[n] * R.dp[n])); }
	double plane1D_dr(specie * S, int n, int m, double dT, double D)   { return (*S).get_C(n, m) / (D*dT); }
	
	double sphere1D_ar(int n, int m)                       { return (-2.0 / (R.dp[n] * R.dm[n] + R.dm[n] * R.dm[n])) + (1 / R.V[n])*(1 / (R.dp[n] +R.dm[n])); }
	double sphere1D_br(int n, int m, double dT, double D)  { return (2.0 / (R.dp[n] * R.dm[n] + R.dm[n] * R.dm[n])) + (2.0 / (R.dp[n] * R.dm[n] + R.dp[n] * R.dp[n])) + 1 / (D*dT); }
	double sphere1D_gr(int n, int m)                       { return (-2.0 / (R.dp[n] * R.dm[n] + R.dp[n] * R.dp[n])) - (1 / R.V[n])*(1 / (R.dp[n] + R.dm[n])); }
	double sphere1D_dr(specie * S, int n, int m, double dT, double D)   { return (*S).get_C(n, m) / (D*dT); }
	
	double plane1D_az(int n, int m)                       { return (-2.0 / (Z.dp[m] * Z.dm[m] + Z.dm[m] * Z.dm[m])); }
	double plane1D_bz(int n, int m, double dT, double D)  { return (2.0 / (Z.dp[m] * Z.dm[m] + Z.dm[m] * Z.dm[m])) + (2.0 / (Z.dp[m] * Z.dm[m] + Z.dp[m] * Z.dp[m])) + 1 / (D*dT); }
	double plane1D_gz(int n, int m)                       { return (-2.0 / (Z.dp[m] * Z.dm[m] + Z.dp[m] * Z.dp[m])); }
	double plane1D_dz(specie * S, int n, int m, double dT, double D)   { return (*S).get_C(n, m) / (D*dT); }

	double sphere1D_az(int n, int m)                       { return (-2.0 / (Z.dp[m] * Z.dm[m] + Z.dm[m] * Z.dm[m])) + (1 / Z.V[m])*(1 / (Z.dp[m] + Z.dm[m])); }
	double sphere1D_bz(int n, int m, double dT, double D)  { return (2.0 / (Z.dp[m] * Z.dm[m] + Z.dm[m] * Z.dm[m])) + (2.0 / (Z.dp[m] * Z.dm[m] + Z.dp[m] * Z.dp[m])) + 1 / (D*dT); }
	double sphere1D_gz(int n, int m)                       { return (-2.0 / (Z.dp[m] * Z.dm[m] + Z.dp[m] * Z.dp[m])) - (1 / Z.V[m])*(1 / (Z.dp[m] + Z.dm[m])); }
	double sphere1D_dz(specie * S, int n, int m, double dT, double D)   { return (*S).get_C(n, m) / (D*dT); }
	//
	
	
	
	//physics-operation on specie - common 

	void insulation(specie * ,char, int , int );
	void bulk(specie *, char, int, int);
	void butlervolmer_boundary(specie *, specie *, char, int, int, double V, double a_bv, double K_bv, double Db);
	void butlervolmer_boundary1(specie *, char, int, int, double V, double a_bv, double K_bv);
	void langmuir_boundary(specie *, char, int, int, double , double , double, double);
	void nernst(specie *, char, int, int, double);
	void continous1D(specie *, int, int, double, double );
	~PreTomas(){};

};

    PreTomas::PreTomas(){} 

	PreTomas::PreTomas(Grid a, Grid b)
	{
		R = a; 
		Z = b;
		pRm = R.VectorMinus();
		pRp = R.VectorPlus();
		pZm = Z.VectorMinus();
		pZp = Z.VectorPlus();
		
		R.VectorMin();
		R.VectorPlu();
		Z.VectorMin();
		Z.VectorPlu();
		//pS = S_temp;
	}
	
	void PreTomas::continous1D(specie * S, int n, int m,double D, double dT)
	{
		*(S->set_ar(n, m)) = (-2.0 / (Z.dm[m] * R.dm[n] + Z.dm[m] * Z.dm[m]));	
		*(S->set_br(n, m)) = (2.0 / (Z.dm[m] * R.dm[n] + R.dm[n] * R.dm[n])) + (2.0 / (Z.dm[m] * R.dm[n] + Z.dm[m] * Z.dm[m])) + 1 / (D*dT) + 1 / (R.V[n] * R.dm[n]);
		*(S->set_gr(n, m)) = -(2.0 / (Z.dm[m] * R.dm[n] + R.dm[n] * R.dm[n])) - 1 / (R.V[n] * R.dm[n]);
		*(S->set_dr(n, m)) = (*S).get_C(n, m) / (D*dT);
		
	}
	
	void PreTomas::insulation(specie * S,char d, int n, int m)
	{
		if (d == 'R')
		{
			if (n == 0)
			{
				*(S->set_ar(n, m)) = 0;
				*(S->set_br(n, m)) = -1;
				*(S->set_gr(n, m)) = 1;
				*(S->set_dr(n, m)) = 0;
			}
			else
			{
				*(S->set_ar(n, m)) = 1;
				*(S->set_br(n, m)) = -1;
				*(S->set_gr(n, m)) = 0;
				*(S->set_dr(n, m)) = 0;
			}
		}
		if (d == 'Z')
		{
			if (m == 0)//depends if its in the begining or in the end, for the tomas alg to fit to the number of elements
			{
				*(S->set_az(n, m)) = 0;
				*(S->set_bz(n, m)) = -1;
				*(S->set_gz(n, m)) = 1;
				*(S->set_dz(n, m)) = 0;
			}
			else
			{
			*(S->set_az(n, m)) = 1;
			*(S->set_bz(n, m)) = -1;
			*(S->set_gz(n, m)) = 0;
			*(S->set_dz(n, m)) = 0; 
			}
		}
	}

	void PreTomas::bulk(specie * S, char d, int n, int m)
	{
		if (d == 'R')
		{
			*(S->set_ar(n, m)) = 0;
			*(S->set_br(n, m)) = 1;
			*(S->set_gr(n, m)) = 0;
			*(S->set_dr(n, m)) = (*S).init_c;
		}
		else
		if (d == 'Z')
		{
			*(S->set_az(n, m)) = 0;
			*(S->set_bz(n, m)) = 1;
			*(S->set_gz(n, m)) = 0;
			*(S->set_dz(n, m)) = (*S).init_c;
		}
		
	}
	
	void PreTomas::nernst(specie * S, char d, int n, int m, double etha)
	{
		if (d == 'R')
		{
			*(S->set_ar(n, m)) = 0;
			*(S->set_br(n, m)) = 1;
			*(S->set_gr(n, m)) = 0;
			*(S->set_dr(n, m)) = 1/(1+exp(etha));
		}
		else
		if (d == 'Z')
		{
			*(S->set_az(n, m)) = 0;
			*(S->set_bz(n, m)) = 1;
			*(S->set_gz(n, m)) = 0;
			*(S->set_dz(n, m)) = 1 / (1 + exp(etha));
		}
	}
	
	
	void PreTomas::butlervolmer_boundary(specie * s1, specie* s2, char d, int n, int m, double V, double a_bv, double K_bv, double Db)
	{
		if (d == 'Z')
		{
			*(s1->set_gz(n, m)) = (-(Z.dp[m]) * K_bv*exp(-a_bv*V)); 
			*(s1->set_bz(n, m)) = (1 + (Z.dp[m]) * K_bv*exp(V - a_bv*V));
			*(s1->set_az(n, m)) = (-1);
			*(s1->set_dz(n, m)) = 0;

			*(s2->set_az(n, m)) = (-(Z.dp[m] / Db) * K_bv*exp(V - a_bv*V));
			*(s2->set_bz(n, m)) = (1 + (Z.dp[m] / Db) * K_bv*exp(-a_bv*V));
			*(s2->set_gz(n, m)) = (-1);
			*(s2->set_dz(n, m)) = 0;
		}
		
		if (d == 'R')//need to be checked if its true for R too.
		{
			*(s1->set_gr(n, m)) = (-(R.dp[m]) * K_bv*exp(-a_bv*V)); 
			*(s1->set_br(n, m)) = (1 + (R.dp[m]) * K_bv*exp(V - a_bv*V));
			*(s1->set_ar(n, m)) = (-1);
			*(s1->set_dr(n, m)) = 0;

			*(s2->set_ar(n, m)) = (-(R.dp[m] / Db) * K_bv*exp(V - a_bv*V));
			*(s2->set_br(n, m)) = (1 + (R.dp[m] / Db) * K_bv*exp(-a_bv*V));
			*(s2->set_gr(n, m)) = (-1);
			*(s2->set_dr(n, m)) = 0;
		}
	}
	void PreTomas::butlervolmer_boundary1(specie * s1, char d, int n, int m, double V, double a_bv, double K_bv)
	{
		
		if (d == 'Z')
		{
			*(s1->set_az(n, m)) = 0; 
			*(s1->set_bz(n, m)) = 1 + K_bv*(Z.dp[m]) *(exp(a_bv*V))*(1 + exp(-V));
			*(s1->set_gz(n, m)) =-1 ;
			*(s1->set_dz(n, m)) = K_bv*(Z.dp[m]) *(exp(a_bv*V))*(exp(-V));
		}
	
	}
	void PreTomas::langmuir_boundary(specie * S, char d, int n, int m, double K0, double K1,double Db, double theta)
	{
		if (d == 'R')
		{
			*(S->set_ar(n, m)) = 0;
			*(S->set_br(n, m)) = (K0 / Db)*R.dp[m] * (1 - theta);
			*(S->set_gr(n, m)) = -1;
			*(S->set_dr(n, m)) = (K1 / Db)*R.dp[m] * (theta);
		}
		else
		if (d == 'Z')
		{
			*(S->set_az(n, m)) = 0;
			*(S->set_bz(n, m)) = 1 + (K0 / Db)*Z.dp[m] * (1 - theta);
			*(S->set_gz(n, m)) =-1;
			*(S->set_dz(n, m)) = (K1 / Db)*Z.dp[m] * (theta);
		}
		


	}

	double PreTomas::R_d(specie * S, int n, int m, double dT, double D)  
	
	{ 
		double Temp;
			if (m > 0){
			Temp = -(2.0 / ((Z.dp[m] + Z.dm[m])*Z.dm[m]))*(*S).get_C(n, m - 1);
			Temp += (-(2.0 / (Z.dp[m] + Z.dm[m]))*((-1 / Z.dp[m]) + (-1 / Z.dm[m])) - 2 / (D*dT))*(*S).get_C(n, m);
			Temp += -(2.0 / ((Z.dp[m] + Z.dm[m])*Z.dp[m]))*(*S).get_C(n, m + 1);
			}
			else
		{
				Temp = (-(4.0 / (Z.dp[m] + Z.dp[m]))*(-1 / Z.dp[m]) - 2 / (D*dT))*(*S).get_C(n, m);
				Temp += -((4.0 / (Z.dp[m] + Z.dp[m]))*(1 / Z.dp[m]))*(*S).get_C(n , m+1);
		}
		return   Temp;
	}
	double PreTomas::Z_d(specie * S, int n, int m, double dT, double D)	
	{
		double Temp;
		if (n > 0){
			Temp = ((-1 / (R.dp[n] + R.dm[n]))*((2 / R.dm[n]) - (1 / R.V[n])))* (*S).get_C(n - 1, m);
			Temp += ((2.0 / (R.dp[n] + R.dm[n]))*((1 / R.dp[n]) + (1 / R.dm[n])) - 2 / (D*dT))*(*S).get_C(n, m);
			Temp += ((-1 / (R.dp[n] + R.dm[n]))*((2 / R.dp[n]) + (1 / R.V[n])))*(*S).get_C(n + 1, m);
		}
		else
		{
			Temp = (-(4.0 / (R.dp[n] + R.dp[n]))*(-1 / R.dp[n]) - 2 / (D*dT))*(*S).get_C(n, m);
			Temp += -((4.0 / (R.dp[n] + R.dp[n]))*(1 / R.dp[n]))*(*S).get_C(n+1, m);
		}
		return   Temp;
	}
	
	double PreTomas::Y_d(specie * S, int n, int m, double dT, double D)

	{
		double Temp = -(2.0 / ((Z.dp[m] + Z.dm[m])*Z.dm[m]))*(*S).get_C(n, m - 1);
		Temp += (-(2.0 / (Z.dp[m] + Z.dm[m]))*((-1 / Z.dp[m]) + (-1 / Z.dm[m])) - 2 / (D*dT))*(*S).get_C(n, m);
		Temp += -(2.0 / ((Z.dp[m] + Z.dm[m])*Z.dp[m]))*(*S).get_C(n, m + 1);
		return   Temp;
	}
	double PreTomas::X_d(specie * S, int n, int m, double dT, double D)
	{
		double Temp = -(2.0 / ((R.dp[n] + R.dm[n])*R.dm[n]))*(*S).get_C(n-1, m);
		Temp += (-(2.0 / (R.dp[n] + R.dm[n]))*((-1 / R.dp[n]) + (-1 / R.dm[n])) - 2 / (D*dT))*(*S).get_C(n, m);
		Temp += -(2.0 / ((R.dp[n] + R.dm[n])*R.dp[n]))*(*S).get_C(n+1, m);
		return   Temp;
	}



	
	/*PreTomas::~PreTomas()
	{
	}*/

