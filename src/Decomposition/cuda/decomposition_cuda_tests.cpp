#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "VCluster/VCluster.hpp"
#include "Decomposition/CartDecomposition.hpp"

#define SUB_UNIT_FACTOR 1024

BOOST_AUTO_TEST_SUITE( decomposition_to_gpu_test )

BOOST_AUTO_TEST_CASE( decomposition_to_gpu_test_use )
{
	// Vcluster
	Vcluster & vcl = create_vcluster();

	//! [Create CartDecomposition]
	CartDecomposition<3, float> dec(vcl);

	// Physical domain
	Box<3, float> box( { 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
	size_t div[3];

	// Get the number of processor and calculate the number of sub-domain
	// for each processor (SUB_UNIT_FACTOR=64)
	size_t n_proc = vcl.getProcessingUnits();
	size_t n_sub = n_proc * SUB_UNIT_FACTOR;

	// Set the number of sub-domains on each dimension (in a scalable way)
	for (int i = 0; i < 3; i++)
	{	div[i] = openfpm::math::round_big_2(pow(n_sub,1.0/3));}

	// Define ghost
	Ghost<3, float> g(0.01);

	// Boundary conditions
	size_t bc[] = { PERIODIC, PERIODIC, PERIODIC };

	// Decompose
	dec.setParameters(div,box,bc,g);
	dec.decompose();

	dec.toKernel()
}

BOOST_AUTO_TEST_SUITE_END()
