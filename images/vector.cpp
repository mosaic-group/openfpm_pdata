#include "Vector/vector_dist.hpp"
#include "Decomposition/CartDecomposition.hpp"


template<typename T> class Particle
{
public:

#ifdef SE_CLASS3

	typedef boost::fusion::vector<T,T[3],T[3][3],SE3_ADD_PROP(3)> type;

#else

	typedef boost::fusion::vector<T,T[3],T[3][3]> type;

#endif

	type data;

	static const unsigned int s = 0;
	static const unsigned int v = 1;
	static const unsigned int t = 2;

#ifdef SE_CLASS3

	static const unsigned int max_prop = SE3_MAX_PROP(3);
	static const unsigned int max_prop_real = 3;

#else

	static const unsigned int max_prop = 3;
	static const unsigned int max_prop_real = 3;

#endif

	static inline bool noPointers()
	{
		return true;
	}

};

int main(int argc, char* argv[])
{
	//
	// ### WIKI 2 ###
	//
	// Here we Initialize the library, than we create a uniform random generator between 0 and 1 to to generate particles
	// randomly in the domain, we create a Box that define our domain, boundary conditions, and ghost
	//
	openfpm_init(&argc,&argv);
	Vcluster<> & v_cl = create_vcluster();

	// set the seed
	// create the random generator engine
	std::default_random_engine eg(v_cl.getProcessUnitID()*100);
	std::uniform_real_distribution<float> ud(0.0f, 1.0f);

	Box<2,float> domain({0.0,0.0},{1.0,1.0});
    size_t bc[2]={PERIODIC,PERIODIC};
	Ghost<2,float> g(0.01);
	
	vector_dist<2,float, Particle<float> > vd(4096,domain,bc,g);

	auto it = vd.getIterator();

	while (it.isNext())
	{
		auto key = it.get();

		vd.getPos(key)[0] = ud(eg);
		vd.getPos(key)[1] = ud(eg);

		vd.template getProp<1>(key)[0] = sin(10.0*vd.getPos(key)[0]);
		vd.template getProp<1>(key)[1] = sin(10.0*vd.getPos(key)[1]);

		++it;
	}

	vd.write("Vector/vector_before_map",CSV_WRITER);

	vd.map();

	vd.write("Vector/vector_after_map",CSV_WRITER);
	
	vd.ghost_get<0>();

	vd.write("Vector/vector_ghost_fill",CSV_WRITER);

	vd.getDecomposition().write("Vector/vect_decomposition");

	// move the particles

	for (size_t i = 0 ; i < 100 ; i++)
	{
		auto it = vd.getDomainIterator();

		while (it.isNext())
		{
			auto key = it.get();

			vd.getPos(key)[0] += 0.005;
			vd.getPos(key)[1] += 0.005;
                        
            vd.template getProp<1>(key)[0] = 0.005;
            vd.template getProp<1>(key)[1] = 0.005;

			++it;
		}
		vd.write_frame("Vector/vector_move_before_map",i,CSV_WRITER);
		vd.map();
		vd.write_frame("Vector/vector_move",i,CSV_WRITER);
	}

	openfpm_finalize();
}
