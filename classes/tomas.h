#include <vector>

class tomas
{
	friend class specie;
private:
	
	vector<double> a ;
	vector<double> b ; 
	vector<double> g ; 
	vector<double> d ; 
	vector<double> c; int size; int I; specie *S; char direction; 

public:
	tomas(specie*, char, int);
	tomas(vector<double>*, vector<double>*, vector<double>*, vector<double>*);
	tomas(vector<double>*, vector<double>*, vector<double>*, vector<double>*, int, int);
	void modify_g(); 
	void modify_d();
	void solve_c();
	void reduce_size(int red) { size = size - red; }
	vector<double> get_c(){	return c;}
	vector<double>* ads_c(){ return &c; }
	~tomas(){};
};

tomas::tomas(specie *s, char direct , int nm)
{
	direction = direct;
	I = nm;
	S = s;

	if (direction == 'R')
	{
		size = (*S).Rsize;
	(a).resize(size);
	(b).resize(size);
	(g).resize(size);
	(d).resize(size);
	}
	if (direction == 'Z')
	{
		size = (*S).Zsize;
	(a).resize(size);
	(b).resize(size);
	(g).resize(size);
	(d).resize(size);
	}
}
tomas::tomas(vector<double> *A, vector<double> *B, vector<double> *G, vector<double> *D)
{
	int size = (*A).size();
	(c).resize(size);
	(*G)[0] /= (*B)[0];

	for (int i = 1; i < size; i++)
	{
		(*G)[i] /= ((*B)[i] - (*G)[i - 1] * (*A)[i]);
	}
	//Modify delta coefficients
	(*D)[0] /= (*B)[0];
	for (int i = 1; i < size; i++)
	{
		(*D)[i] = ((*D)[i] - (*D) [i - 1] * (*A)[i]) / ((*B)[i] - (*G)[i - 1] * (*A)[i]);
	}
	//Back Substitution
	c[size - 1] = (*D)[size - 1];
	for (int i = size - 2; i >= 0; i--)		
	{
		c[i] = (*D)[i] - (*G)[i] * (c)[i + 1];
	}
	////calculating new thet a
}
tomas::tomas(vector<double> *A, vector<double> *B, vector<double> *G, vector<double> *D,int start, int end)
{
	int size = end - start;
	(c).resize(size);
	(*G)[start] /= (*B)[start];

	for (int i = 1; i < size; i++)
	{
		(*G)[start + i] /= ((*B)[start + i] - (*G)[start + i - 1] * (*A)[start+i]);
	}
	//Modify delta coefficients
	(*D)[start + 0] /= (*B)[start + 0];
	for (int i = 1; i < size; i++)
	{
		(*D)[start + i] = ((*D)[start + i] - (*D)[start + i - 1] * (*A)[start + i]) / ((*B)[start + i] - (*G)[start + i - 1] * (*A)[start + i]);
	}
	//Back Substitution
	c[size - 1] = (*D)[start + end - 1];
	for (int i = size - 2; i >= 0; i--)
	{
		c[i] = (*D)[start + i] - (*G)[start + i] * (c)[i + 1];
	}
	////calculating new thet a
}
void tomas::modify_g()
{

	if (direction == 'R')
	{
		g[0] =     *(*S).set_gr(0, I)    /    *(*S).set_br(0, I);
		for (int i = 1; i < size; i++)
		{
			g[i] = *(*S).set_gr(i, I)    / (    *(*S).set_br(i, I) - g[i-1]    *      *(*S).set_ar(i, I)    );
		}
	}


	if (direction == 'Z')
	{
		g[0] = *(*S).set_gz(I, 0)        /      *(*S).set_bz(I, 0);
		for (int i = 1; i < size; i++)
		{
			g[i] = *(*S).set_gz(I, i) /  (   *(*S).set_bz(I, i) - g[i-1]     *     *(*S).set_az(I, i)  );
		}
	}

}
void tomas::modify_d()
{

	if (direction == 'R')
	{
		d[0] = *(*S).set_dr(0, I) / *(*S).set_br(0, I);
		for (int i = 1; i < size; i++)
		{
			d[i] = (*(*S).set_dr(i, I) - d[i - 1] * *(*S).set_ar(i, I)) / (*(*S).set_br(i, I) - g[i - 1] * *(*S).set_ar(i, I));
		}
	}
	if (direction == 'Z')
	{
		d[0] =  *(*S).set_dz(I, 0) / *(*S).set_bz(I, 0);
		for (int i = 1; i < size; i++)
		{
			d[i] = (*(*S).set_dz(I, i) - d[i - 1] * *(*S).set_az(I, i)) / (*(*S).set_bz(I, i) - g[i - 1] * *(*S).set_az(I, i));
		}
	}
}
void tomas::solve_c()
	//Back Substitution
{	

	if (direction == 'Z')
	{
		*(*S).set_C(I, size - 1) = d[size - 1];
		for (int i = size - 2; i >= 0; i--)
		{
			*(*S).set_C(I, i) = d[i] - g[i] * *(*S).set_C(I, i + 1);
		}
	}
	
	if (direction == 'R')
	{
		*(*S).set_C(size-1, I) = d[size - 1];
		for (int i = size - 2; i >= 0; i--)
		{
			*(*S).set_C(i, I) = d[i] - g[i] * *(*S).set_C(i + 1,I);
		}
	}
}
