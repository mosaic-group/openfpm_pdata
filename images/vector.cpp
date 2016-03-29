#include "Vector/vector_dist.hpp"
#include "Decomposition/CartDecomposition.hpp"


template<typename T> class Particle
{
public:

	typedef boost::fusion::vector<T,T[3],T[3][3]> type;

	type data;

	static const unsigned int s = 0;
	static const unsigned int v = 1;
	static const unsigned int t = 2;
	static const unsigned int max_prop = 3;
};

int main(int argc, char* argv[])
{
	//
	// ### WIKI 2 ###
	//
	// Here we Initialize the library, than we create a uniform random generator between 0 and 1 to to generate particles
	// randomly in the domain, we create a Box that define our domain, boundary conditions, and ghost
	//
	init_global_v_cluster(&argc,&argv);
	Vcluster & v_cl = *global_v_cluster;
	
	typedef Point<2,float> s;

	// set the seed
	// create the random generator engine
	std::default_random_engine eg(v_cl.getProcessUnitID()*100);
	std::uniform_real_distribution<float> ud(0.0f, 1.0f);

	Box<2,float> domain({0.0,0.0},{1.0,1.0});
    size_t bc[2]={PERIODIC,PERIODIC};
	Ghost<2,float> g(0.01);
	
	vector_dist<2,float, Particle<float>, CartDecomposition<2,float> > vd(4096,domain,bc,g);

	auto it = vd.getIterator();

	while (it.isNext())
	{
		auto key = it.get();

		vd.template getPos<s::x>(key)[0] = ud(eg);
		vd.template getPos<s::x>(key)[1] = ud(eg);

		++it;
	}

	vd.write("Vector/vector_before_map");

	vd.map();

	vd.write("Vector/vector_after_map");
	
	vd.ghost_get<0,1,2>();

	delete_global_v_cluster();
}
