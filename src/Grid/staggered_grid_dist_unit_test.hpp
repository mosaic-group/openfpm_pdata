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
	// Domain
	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	size_t k = 1024;

	// grid size
	size_t sz[2] = {k,k};

	// Ghost
	Ghost<2,float> g(0.0);

	staggered_grid_dist<2,float,Point2D_test<float>,CartDecomposition<2,float>> sg(sz,domain,g);
	sg.setDefaultStagPosition();

	// We check that the staggered position is correct
	const openfpm::vector<comb<2>> (& cmbs)[6] = sg.getStagPositions();


	BOOST_REQUIRE_EQUAL(cmbs[0].size(),1ul);
	BOOST_REQUIRE_EQUAL(cmbs[1].size(),1ul);
	BOOST_REQUIRE_EQUAL(cmbs[2].size(),1ul);
	BOOST_REQUIRE_EQUAL(cmbs[3].size(),1ul);
	BOOST_REQUIRE_EQUAL(cmbs[4].size(),2ul);
	BOOST_REQUIRE_EQUAL(cmbs[5].size(),4ul);

	BOOST_REQUIRE(cmbs[0].get(0) == comb<2>({0,0}));
	BOOST_REQUIRE(cmbs[1].get(0) == comb<2>({0,0}));
	BOOST_REQUIRE(cmbs[2].get(0) == comb<2>({0,0}));
	BOOST_REQUIRE(cmbs[3].get(0) == comb<2>({0,0}));
	BOOST_REQUIRE(cmbs[4].get(0) == comb<2>({0,-1}));
	BOOST_REQUIRE(cmbs[4].get(1) == comb<2>({-1,0}));
	BOOST_REQUIRE(cmbs[5].get(0) == comb<2>({0,0}));
	BOOST_REQUIRE(cmbs[5].get(1) == comb<2>({-1,-1}));
	BOOST_REQUIRE(cmbs[5].get(2) == comb<2>({-1,-1}));
	BOOST_REQUIRE(cmbs[5].get(3) == comb<2>({0,0}));
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* SRC_GRID_STAGGERED_GRID_DIST_UNIT_TEST_HPP_ */
