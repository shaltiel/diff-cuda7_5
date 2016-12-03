
#include <vector>

using namespace std;

class Grid {
	

public:
	vector<double> V; double h; double gamma; vector<int> node;
	vector<double> *dRp = new vector<double>(); vector<double> *dRm = new vector<double>();
	vector<double> dp; vector<double> dm;

	Grid(double, double);
	 Grid();
	~Grid(){};
	void new_node(double);
	void open_node(double);
	vector<double> get_vector(){ return V; }
	int get_node(int n){ return node[n]; }
	void VectorPlusMinus();
	vector<double>* VectorMinus();
	vector<double>* VectorPlus();
	void VectorMin();
	void VectorPlu();
};

//Constructor 
Grid::Grid(){}
Grid::Grid(double a, double b)
{
	h = a;
	gamma = b;
	node.push_back(0);
	V.push_back(0.0);
}

//function implementation 
void Grid::new_node(double max){
	double h_temp = h;
	//inital is r0=
	if (max > V.back())
	{
		V.push_back(max); //inital is r0=1
		int	mirrormesh = V.size() - 2;
		while (V[mirrormesh] < V[mirrormesh + 1]) //writing the distances
		{
			mirrormesh++;
			V.insert(V.begin() + mirrormesh, V[mirrormesh - 1] + h);
			V.insert(V.begin() + mirrormesh + 1, V[mirrormesh + 1] - h);
			h *= gamma; //growing h
		}

		V.erase(V.begin() + mirrormesh);
	//	V.erase(V.begin() + mirrormesh);
		node.push_back(V.size());
		h = h_temp;
	}
	if (max < V.back())
	{
		V.push_back(max); //inital is r0=1
		int	mirrormesh = V.size() - 2;
		while (V[mirrormesh] > V[mirrormesh + 1]) //writing the distances
		{
			mirrormesh++;
			V.insert(V.begin() + mirrormesh, V[mirrormesh - 1] - h);
			V.insert(V.begin() + mirrormesh + 1, V[mirrormesh + 1] + h);
			h *= gamma; //growing h
		}

		V.erase(V.begin() + mirrormesh);
		V.erase(V.begin() + mirrormesh);
		node.push_back(V.size());
		h = h_temp;
	}
}



void Grid::open_node(double max){
	double h_temp = h;
	//inital is r0=1
	V.push_back(max); //inital is r0=1
	int	mirrormesh = V.size() - 2;
	while (V[mirrormesh] < V[mirrormesh + 1]) //writing the distances
	{
		mirrormesh++;
		V.insert(V.begin() + mirrormesh, V[mirrormesh - 1] + h);
		h *= gamma; //growing h
	}

	//R.erase(R.begin() + mirrormesh);
	V.erase(V.begin() + mirrormesh);
	node.push_back(V.size());
	h = h_temp;
}


void Grid::VectorPlusMinus()
{
	dRp->resize(V.size(),0.0);
	dRm->resize(V.size(),0.0);
	for (size_t n = 0; n < V.size(); n++)
	{
		if (n < (V.size() - 1)){ (*dRp)[n] = V[n + 1] - V[n]; }
		if (n > 0)          { (*dRm)[n] = V[n] - V[n - 1]; }
	}
}


vector<double>* Grid::VectorPlus()
{
	dRp->resize(V.size(), 0.0);	
	for (std::size_t n = 0; n < V.size(); n++)
	{
		if (n < (V.size() - 1)){ (*dRp)[n] = V[n + 1] - V[n]; }
	}
	return (dRp);
}

vector<double>* Grid::VectorMinus()
{
	dRm->resize(V.size(), 0.0);
	for (std::size_t n = 0; n < V.size(); n++)
	{
		if (n > 0)          { (*dRm)[n] = V[n] - V[n - 1]; }
	}
	return (dRm);
}

void Grid::VectorPlu()
{
	dp.resize(V.size(), 0.0);
	for (std::size_t n = 0; n < V.size(); n++)
	{
		if (n < (V.size() - 1)){ dp[n] = V[n + 1] - V[n]; }
	}
}

void Grid::VectorMin()
{
	dm.resize(V.size(), 0.0);
	for (std::size_t n = 0; n < V.size(); n++)
	{
		if (n > 0)          { dm[n] = V[n] - V[n - 1]; }
	}
	
}