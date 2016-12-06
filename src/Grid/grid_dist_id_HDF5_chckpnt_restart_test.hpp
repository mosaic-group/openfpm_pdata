/*
 * grid_dist_id_HDF5_chckpnt_restart_test.hpp
 *
 *  Created on: Nov 9, 2016
 *      Author: Yaroslav Zaluzhnyi
 */

#ifndef SRC_GRID_GRID_DIST_ID_HDF5_CHCKPNT_RESTART_TEST_HPP_
#define SRC_GRID_GRID_DIST_ID_HDF5_CHCKPNT_RESTART_TEST_HPP_

#include "Grid/grid_dist_id.hpp"

BOOST_AUTO_TEST_SUITE( gd_hdf5_chckpnt_rstrt_test )

BOOST_AUTO_TEST_CASE( grid_dist_id_hdf5_save_test )
{

	// Input data
	size_t k = 1000;

	size_t ghost_part = 0.02;

	// Domain
	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	Vcluster & v_cl = create_vcluster();

	// Skip this test on big scale
	if (v_cl.getProcessingUnits() >= 32)
		return;

	if (v_cl.getProcessUnitID() == 0)
			std::cout << "Saving Distributed 2D Grid..." << std::endl;

	// grid size
	size_t sz[2];
	sz[0] = k;
	sz[1] = k;

	// Ghost
	Ghost<2,float> g(ghost_part);

	// Distributed grid with id decomposition
	grid_dist_id<2, float, scalar<float>, CartDecomposition<2,float>> g_dist(sz,domain,g);

	// get the decomposition
	auto & dec = g_dist.getDecomposition();

	// check the consistency of the decomposition
	bool val = dec.check_consistency();
	BOOST_REQUIRE_EQUAL(val,true);

	timer t;
	t.start();
	// Save the grid
    g_dist.save("grid_dist_id.h5");
	t.stop();

	std::cout << "Saving time: " << t.getwct() << std::endl;
}

BOOST_AUTO_TEST_CASE( grid_dist_id_hdf5_load_test )
{

	// Input data
	size_t k = 1000;

	size_t ghost_part = 0.02;

	// Domain
	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	Vcluster & v_cl = create_vcluster();

	// Skip this test on big scale
	if (v_cl.getProcessingUnits() >= 32)
		return;

	if (v_cl.getProcessUnitID() == 0)
			std::cout << "Loading Distributed 2D Grid..." << std::endl;

	// grid size
	size_t sz[2];
	sz[0] = k;
	sz[1] = k;

	// Ghost
	Ghost<2,float> g(ghost_part);

	// Distributed grid with id decomposition
	grid_dist_id<2, float, scalar<float>, CartDecomposition<2,float>> g_dist(sz,domain,g);

	timer t;
	t.start();
	// Save the grid
	g_dist.load("grid_dist_id.h5");
	t.stop();

	std::cout << "Loading time: " << t.getwct() << std::endl;

	auto it = g_dist.getOldDomainIterator();

	size_t count = 0;

	while (it.isNext())
	{
		//key
		grid_dist_key_dx<2> key = it.get();

		//g_dist.get(key);

		++it;
		count++;
	}

	openfpm::vector<size_t> count_total;
	v_cl.allGather(count,count_total);
	v_cl.execute();

	size_t sum = 0;

	for (size_t i = 0; i < count_total.size(); i++)
		sum += count_total.get(i);

	BOOST_REQUIRE_EQUAL(sum, (size_t)k*k);
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* SRC_GRID_GRID_DIST_ID_HDF5_CHCKPNT_RESTART_TEST_HPP_ */
