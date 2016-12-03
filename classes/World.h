#include <vector>

using namespace std;

 int edge_n = 10, domain_n = 1, point_n = -1;
class World{
private:
	int Rsize, Zsize;
	vector<int> w;  
	friend class Grid;
public:
	int make_edge(int,int,int,int);
	int make_domain(int, int, int, int);
	int make_point(int, int);
	bool is(int, int,int);
	void setworld(Grid, Grid);
	World();
	World(Grid, Grid);
	~World(){};
};

World::World(Grid GrR, Grid GrZ){
	Rsize = GrR.V.size();
	Zsize = GrZ.V.size();
	w.resize(Rsize*Zsize, 0);	
}
World::World(){}
 int World::make_edge(int n1, int m1,int n2,int m2)
{
	for (int n = n1; n <= n2; n++)
	{
	 for (int m = m1; m <= m2; m++)
		{
		 if (w[n + m*Rsize]>=0)
			w[n + m*Rsize] = edge_n;
		 else std::cout << "warning: you can't write edge over defined point" << '\n';
		}
	}	
	return edge_n++;
}
 int World::make_domain(int n1, int m1, int n2, int m2)
 {
	 for (int n = n1; n <= n2; n++)
	 {
		 for (int m = m1; m <= m2; m++)
		 {
			 if (w[n + m*Rsize]==0)
			 w[n + m*Rsize] = domain_n;
			 else std::cout << "warning: you can't write domain over defined point/edge" << '\n';
		 }
	 }
	 return domain_n++;
 }
 int World::make_point(int n1, int m1)
 {
			 w[n1 + m1*Rsize] = point_n;
	 return point_n--;
 }
bool World::is(int n1, int m1 ,int prop )
 {
	return  (w[n1 + m1*Rsize] == prop);
 }

void World::setworld(Grid GrR, Grid GrZ){
	Rsize = GrR.V.size();
	Zsize = GrZ.V.size();
	w.resize(Rsize*Zsize, 0);
}