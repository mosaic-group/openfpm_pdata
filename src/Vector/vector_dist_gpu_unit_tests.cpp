
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "VCluster/VCluster.hpp"
#include <Vector/vector_dist.hpp>

BOOST_AUTO_TEST_SUITE( vector_dist_gpu_test )

void print_test(std::string test, size_t sz)
{
	if (create_vcluster().getProcessUnitID() == 0)
		std::cout << test << " " << sz << "\n";
}

BOOST_AUTO_TEST_CASE( vector_dist_gpu_test)
{
	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	// set the ghost based on the radius cut off (make just a little bit smaller than the spacing)
	Ghost<3,float> g(0.01);

	// Boundary conditions
	size_t bc[3]={NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

	vector_dist_gpu<3,float,aggregate<float,float[3]>> vd(1000,domain,bc,g);

	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		vd.getPos(p)[0] = (float)rand() / RAND_MAX;
		vd.getPos(p)[1] = (float)rand() / RAND_MAX;
		vd.getPos(p)[2] = (float)rand() / RAND_MAX;

		++it;
	}

	vd.map();
}

BOOST_AUTO_TEST_SUITE_END()
