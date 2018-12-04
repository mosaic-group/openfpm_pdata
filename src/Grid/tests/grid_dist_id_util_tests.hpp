/*
 * grid_dist_id_util_tests.hpp
 *
 *  Created on: Oct 15, 2017
 *      Author: i-bird
 */

#ifndef SRC_GRID_TESTS_GRID_DIST_ID_UTIL_TESTS_HPP_
#define SRC_GRID_TESTS_GRID_DIST_ID_UTIL_TESTS_HPP_


static void print_test(std::string test, size_t sz)
{
	if (create_vcluster().getProcessUnitID() == 0)
		std::cout << test << " " << sz << "\n";
}


static void Test2D_core(grid_dist_id<2, float, scalar<float>> & g_dist, const size_t (& sz)[2], size_t k)
{
	// check the consistency of the decomposition
	bool val = g_dist.getDecomposition().check_consistency();
	BOOST_REQUIRE_EQUAL(val,true);

	// Grid sm
	grid_sm<2,void> info(sz);

	// get the domain iterator
	size_t count = 0;

	auto dom = g_dist.getDomainIterator();

	while (dom.isNext())
	{
		auto key = dom.get();
		auto key_g = g_dist.getGKey(key);

		g_dist.template get<0>(key) = info.LinId(key_g);

		// Count the point
		count++;

		++dom;
	}

	//! [Create and access a distributed grid]

	// Get the virtual cluster machine
	Vcluster & vcl = g_dist.getVC();

	// reduce
	vcl.sum(count);
	vcl.execute();

	// Check
	BOOST_REQUIRE_EQUAL(count,(size_t)k*k);

	auto dom2 = g_dist.getDomainIterator();

	grid_key_dx<2> start = dom2.getStart();
	grid_key_dx<2> stop = dom2.getStop();

	BOOST_REQUIRE_EQUAL((long int)stop.get(0),(long int)g_dist.size(0)-1);
	BOOST_REQUIRE_EQUAL((long int)stop.get(1),(long int)g_dist.size(1)-1);

	BOOST_REQUIRE_EQUAL(start.get(0),0);
	BOOST_REQUIRE_EQUAL(start.get(1),0);

	bool match = true;

	// check that the grid store the correct information
	while (dom2.isNext())
	{
		auto key = dom2.get();
		auto key_g = g_dist.getGKey(key);

		match &= (g_dist.template get<0>(key) == info.LinId(key_g))?true:false;

		++dom2;
	}

	BOOST_REQUIRE_EQUAL(match,true);

	g_dist.template ghost_get<0>();

	// check that the communication is correctly completed

	auto domg = g_dist.getDomainGhostIterator();

	// check that the grid with the ghost past store the correct information
	while (domg.isNext())
	{
		auto key = domg.get();
		auto key_g = g_dist.getGKey(key);

		// In this case the boundary condition are non periodic
		if (g_dist.isInside(key_g))
		{
			match &= (g_dist.template get<0>(key) == info.LinId(key_g))?true:false;
		}

		++domg;
	}

	BOOST_REQUIRE_EQUAL(match,true);
}

#endif /* SRC_GRID_TESTS_GRID_DIST_ID_UTIL_TESTS_HPP_ */
