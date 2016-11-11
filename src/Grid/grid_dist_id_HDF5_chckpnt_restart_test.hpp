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
	size_t k = 100;

	size_t ghost_part = 0.01;

	/////////////////
	size_t bc[3] = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};

	// Domain
	Box<3,float> domain({-0.3,-0.3,-0.3},{1.0,1.0,1.0});

	Vcluster & v_cl = create_vcluster();

	// Skip this test on big scale
	if (v_cl.getProcessingUnits() >= 32)
		return;

	if (v_cl.getProcessUnitID() == 0)
			std::cout << "Testing 3D grid HDF5 save/load" << std::endl;

	// grid size
	size_t sz[3];
	sz[0] = k;
	sz[1] = k;
	sz[2] = k;

	// Ghost
	Ghost<3,float> g(ghost_part);

	// Distributed grid with id decomposition
	grid_dist_id<3, float, scalar<float>, CartDecomposition<3,float>> g_dist(sz,domain,g);

	// get the decomposition
	auto & dec = g_dist.getDecomposition();

	// check the consistency of the decomposition
	bool val = dec.check_consistency();
	BOOST_REQUIRE_EQUAL(val,true);

	// for each local volume
	// Get the number of local grid needed
	size_t n_grid = dec.getNSubDomain();

	size_t vol = 0;

	// vector of boxes
	openfpm::vector<Box<3,size_t>> vb;

	// Allocate the grids
	for (size_t i = 0 ; i < n_grid ; i++)
	{
		// Get the local hyper-cube
		SpaceBox<3,float> sub = dec.getSubDomain(i);
		sub -= domain.getP1();

		Box<3,size_t> g_box = g_dist.getCellDecomposer().convertDomainSpaceIntoGridUnits(sub,bc);

		vb.add(g_box);

		vol += g_box.getVolumeKey();
	}


	// Save the vector
    g_dist.save("grid_dist_id.h5");
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* SRC_GRID_GRID_DIST_ID_HDF5_CHCKPNT_RESTART_TEST_HPP_ */
