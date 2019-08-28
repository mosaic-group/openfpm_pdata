/*
 * amr_base_gpu_unit_tests.cu
 *
 *  Created on: Aug 28, 2019
 *      Author: i-bird
 */

/*
 * amr_base_unit_test.cpp
 *
 *  Created on: Oct 5, 2017
 *      Author: i-bird
 */
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Amr/grid_dist_amr.hpp"
#include "Point_test.hpp"
#include "Grid/tests/grid_dist_id_util_tests.hpp"

BOOST_AUTO_TEST_SUITE( amr_grid_dist_id_test )


BOOST_AUTO_TEST_CASE( grid_dist_id_amr_gpu )
{
	// Domain
	Box<3,float> domain3({0.0,0.0,0.0},{1.0,1.0,1.0});


	Ghost<3,float> g(0.05);
	sgrid_dist_amr_gpu<3,float,aggregate<float>> amr_g(domain3,g);

	size_t g_sz[3] = {4,4,4};

	size_t n_lvl = 10;

//	amr_g.initLevels(n_lvl,g_sz);


/*	for (size_t i = 0 ; i < amr_g.getNLvl() ; i++)
	{
		// Fill the AMR with something

		size_t count = 0;

		auto it = amr_g.getGridIterator(i);

		while (it.isNext())
		{
			auto key = it.get_dist();
			auto akey = amr_g.getAMRKey(i,key);

			amr_g.template insert<0>(akey) = 3.0;

			count++;

			++it;
		}
	}

	// Iterate across all the levels initialized
	auto it = amr_g.getDomainIterator();

	size_t count = 0;

	while (it.isNext())
	{
		count++;

		++it;
	}

	Vcluster<> & v_cl = create_vcluster();

	v_cl.sum(count);
	v_cl.execute();

	BOOST_REQUIRE_EQUAL(count,correct_result);

	auto itc = amr_g.getDomainIteratorCells();

	size_t count_c = 0;

	while (itc.isNext())
	{
		count_c++;

		++itc;
	}

	v_cl.sum(count_c);
	v_cl.execute();

	auto it_level = amr_g.getDomainIteratorCells(3);

	while (it_level.isNext())
	{
		auto key = it_level.get();

		amr_g.template get<0>(3,key);

		++it_level;
	}

	BOOST_REQUIRE_EQUAL(count_c,correct_result_cell);*/
}

BOOST_AUTO_TEST_SUITE_END()
