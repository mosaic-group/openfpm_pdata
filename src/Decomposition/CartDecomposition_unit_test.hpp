#ifndef CARTDECOMPOSITION_UNIT_TEST_HPP
#define CARTDECOMPOSITION_UNIT_TEST_HPP

#include "CartDecomposition.hpp"
#include "mathutil.hpp"

BOOST_AUTO_TEST_SUITE( CartDecomposition_test )

#define SUB_UNIT_FACTOR 64

BOOST_AUTO_TEST_CASE( CartDecomposition_test_use)
{
	// Vcluster
	Vcluster & vcl = *global_v_cluster;

	// Initialize the global VCluster
	init_global_v_cluster(&boost::unit_test::framework::master_test_suite().argc,&boost::unit_test::framework::master_test_suite().argv);

	CartDecomposition<3,float> dec(vcl);

	// Physical domain
	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});
	size_t div[3];

	// Get the number of processor and calculate the number of sub-domain
	// for decomposition
	size_t n_proc = vcl.getProcessingUnits();
	size_t n_sub = n_proc * SUB_UNIT_FACTOR;

	// Calculate the number of sub-domain on each dimension
	for (int i = 0 ; i < 3 ; i++)
	{div[i] = round_big_2(pow(n_sub,1.0/3));}

	// Decompose
	dec.setParameters(div,box);
}

BOOST_AUTO_TEST_SUITE_END()


#endif
