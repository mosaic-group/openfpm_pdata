#ifndef GRID_DIST_UNIT_TEST_HPP
#define GRID_DIST_UNIT_TEST_HPP

#include "grid_dist_id.hpp"
#include "data_type/scalar.hpp"

BOOST_AUTO_TEST_SUITE( grid_dist_id_test )

template<typename iterator> void jacobi_iteration(iterator g_it, grid_dist_id<2, scalar<float>, CartDecomposition<2,size_t>> & g_dist)
{
	// scalar
	typedef scalar<size_t> S;

	// iterator

	while(g_it.isNext())
	{
		// Jacobi update

		auto pos = g_it.get();

		g_dist.template get<S::ele>(pos) = (g_dist.template get<S::ele>(pos.move(0,1)) +
	                             g_dist.template get<S::ele>(pos.move(0,-1)) +
	                             g_dist.template get<S::ele>(pos.move(1,1)) +
	                             g_dist.template get<S::ele>(pos.move(1,-1)) / 4.0);

		++g_it;
	}
}

BOOST_AUTO_TEST_CASE( grid_dist_id_iterator_test_use)
{
	// Initialize the global VCluster
	init_global_v_cluster(&boost::unit_test::framework::master_test_suite().argc,&boost::unit_test::framework::master_test_suite().argv);

	// grid size
	size_t sz[2] = {1024,1024};

	// Distributed grid with id decomposition

	grid_dist_id<2, scalar<float>, CartDecomposition<2,size_t>> g_dist(sz);

	// get the domain iterator

	size_t count = 0;

	auto dom = g_dist.getDomainIterator();

	while (dom.isNext())
	{
		auto key = dom.get();

		g_dist.template get<0>(key) = count;

		// Count the point
		count++;

		++dom;
	}

	size_t count_check = 0;
	dom = g_dist.getDomainIterator();

	while (dom.isNext())
	{
		auto key = dom.get();

		BOOST_REQUIRE_EQUAL(g_dist.template get<0>(key),count_check);

		count_check++;
		++dom;
	}

	// Get the virtual cluster machine
	Vcluster & vcl = g_dist.getVC();

	// reduce
	vcl.reduce(count);
	vcl.execute();

	// Check
	BOOST_REQUIRE_EQUAL(count,1024*1024);

/*	auto g_it = g_dist.getIteratorBulk();

	auto g_it_halo = g_dist.getHalo();

	// Let try to solve the poisson equation d2(u) = f with f = 1 and computation
	// comunication overlap (100 Jacobi iteration)

	for (int i = 0 ; i < 100 ; i++)
	{
		g_dist.ghost_get();

		// Compute the bulk

		jacobi_iteration(g_it);

		g_dist.ghost_sync();

		// Compute the halo

		jacobi_iteration(g_it_halo);
	}*/
}

BOOST_AUTO_TEST_CASE( grid_dist_id_poisson_test_use)
{
	// grid size
/*	size_t sz[2] = {1024,1024};

	// Distributed grid with id decomposition

	grid_dist_id<2, scalar<float>, CartDecomposition<2,size_t>> g_dist(sz);

	// Create the grid on memory

	g_dist.Create();*/

/*	auto g_it = g_dist.getIteratorBulk();

	auto g_it_halo = g_dist.getHalo();

	// Let try to solve the poisson equation d2(u) = f with f = 1 and computation
	// comunication overlap (100 Jacobi iteration)

	for (int i = 0 ; i < 100 ; i++)
	{
		g_dist.ghost_get();

		// Compute the bulk

		jacobi_iteration(g_it);

		g_dist.ghost_sync();

		// Compute the halo

		jacobi_iteration(g_it_halo);
	}*/
}

BOOST_AUTO_TEST_SUITE_END()

#endif
