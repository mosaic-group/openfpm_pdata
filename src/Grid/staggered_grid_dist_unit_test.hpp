/*
 * staggered_grid_unit_test.hpp
 *
 *  Created on: Aug 20, 2015
 *      Author: i-bird
 */

#ifndef SRC_GRID_STAGGERED_GRID_DIST_UNIT_TEST_HPP_
#define SRC_GRID_STAGGERED_GRID_DIST_UNIT_TEST_HPP_

#include "staggered_dist_grid.hpp"
#include "Point_test.hpp"

BOOST_AUTO_TEST_SUITE( staggered_grid_dist_id_test )


BOOST_AUTO_TEST_CASE( staggered_grid_dist_unit_test)
{
	typedef Point2D_test<float> p;

	Vcluster & v_cl = *global_v_cluster;

	// Domain
	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	size_t k = 1024;

/*	for (size_t k = 1024 ; k >= 2 ; k--)
	{*/
		BOOST_TEST_CHECKPOINT( "Testing grid k=" << k );

		// grid size
		size_t sz[2];
		sz[0] = k;
		sz[1] = k;

		// Ghost
		Ghost<2,float> g(0.01);

		staggered_grid_dist<2,float,Point2D_test<float>,CartDecomposition<2,float>> sg(sz,domain,g);

		sg.setDefaultStagPosition();

		auto it = sg.getDomainIterator();

		while (it.isNext())
		{
			auto key = it.get();

			grid_key_dx<2> keyg = sg.getGKey(key);

			sg.template get<p::s>(key) = 1;

			sg.template get<p::v>(key)[0] = 0;
			sg.template get<p::v>(key)[1] = 1;

			sg.template get<p::t>(key)[0][0] = 0;
			sg.template get<p::t>(key)[0][1] = 1;
			sg.template get<p::t>(key)[1][0] = 2;
			sg.template get<p::t>(key)[1][1] = 3;

			++it;
		}

		sg.write("stag_test.vtk");
/*	}*/
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* SRC_GRID_STAGGERED_GRID_DIST_UNIT_TEST_HPP_ */
