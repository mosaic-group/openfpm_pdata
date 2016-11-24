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
	size_t k = 10;

	size_t ghost_part = 0.01;

	// Domain
	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	Vcluster & v_cl = create_vcluster();

	// Skip this test on big scale
	if (v_cl.getProcessingUnits() >= 32)
		return;

	if (v_cl.getProcessUnitID() == 0)
			std::cout << "Testing 2D grid HDF5 save" << std::endl;

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

	// Save the vector
    g_dist.save("grid_dist_id.h5");
}

BOOST_AUTO_TEST_CASE( grid_dist_id_hdf5_load_test )
{

	// Input data
	size_t k = 100;

	size_t ghost_part = 0.01;

	// Domain
	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	Vcluster & v_cl = create_vcluster();

	// Skip this test on big scale
	if (v_cl.getProcessingUnits() >= 32)
		return;

	if (v_cl.getProcessUnitID() == 0)
			std::cout << "Testing 2D grid HDF5 save/load" << std::endl;

	// grid size
	size_t sz[2];
	sz[0] = k;
	sz[1] = k;

	// Ghost
	Ghost<2,float> g(ghost_part);

	// Distributed grid with id decomposition
	grid_dist_id<2, float, scalar<float>, CartDecomposition<2,float>> g_dist(sz,domain,g);

	g_dist.load("grid_dist_id.h5");
	/*
	auto NN = vd.getCellList(0.5);

	auto it_v = vd.getDomainIterator();

	while (it_v.isNext())
	{
		//key
		vect_dist_key_dx key = it_v.get();

		size_t count = 0;

		// Get the position of the particles
		Point<dim,float> p = vd.getPos(key);

		// Get the neighborhood of the particle
		auto cell_it = NN.template getNNIterator<NO_CHECK>(NN.getCell(p));

		while(cell_it.isNext())
		{
			//Next particle in a cell
			++cell_it;
			count++;
		}

		std::cout << "Count: " << count << std::endl;

		//Next particle in cell list
		++it_v;
	}
*/
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* SRC_GRID_GRID_DIST_ID_HDF5_CHCKPNT_RESTART_TEST_HPP_ */
