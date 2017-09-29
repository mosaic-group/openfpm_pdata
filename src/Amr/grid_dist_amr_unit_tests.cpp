/*
 * grid_dist_amr_dist_unit_tests.cpp
 *
 *  Created on: Sep 21, 2017
 *      Author: i-bird
 */


#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include "grid_dist_amr.hpp"

BOOST_AUTO_TEST_SUITE( grid_dist_amr_test )

/*! \brief Coarsest levels of the grid
 *
 * \param domain Simulation domain
 * \param coars_g coarsest grid resolution
 * \param n_lvl number of levels
 *
 */
void Test3D_amr_create_levels(Box<3,float> & domain, size_t coars_g, size_t n_lvl)
{
	Ghost<3,float> g(0.05);

	size_t g_sz[3] = {coars_g,coars_g,coars_g};

	size_t tot = coars_g*coars_g*coars_g;
	size_t correct_result = 0;
	size_t fact = 1;

	for (size_t i = 0 ; i <  n_lvl ; i++)
	{
		correct_result += tot*fact;
		fact *= 8;
	}


	grid_dist_amr<3,float,aggregate<float>> amr_g(domain,g);

	amr_g.initLevels(n_lvl,g_sz);

	// Iterate across all the levels initialized
	auto it = amr_g.getDomainIterator();

	size_t count = 0;

	while (it.isNext())
	{
		auto key = it.get();

		count++;

		++it;
	}

	Vcluster & v_cl = create_vcluster();

	v_cl.sum(count);
	v_cl.execute();

	BOOST_REQUIRE_EQUAL(count,correct_result);
}

BOOST_AUTO_TEST_CASE( grid_dist_amr_test )
{
	// Domain
	Box<3,float> domain3({0.0,0.0,0.0},{1.0,1.0,1.0});

	long int k = 16*16*16*create_vcluster().getProcessingUnits();
	k = std::pow(k, 1/3.);

	Test3D_amr_create_levels(domain3,k,4);
}

BOOST_AUTO_TEST_SUITE_END()
