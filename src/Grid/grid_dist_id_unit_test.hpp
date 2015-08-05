#ifndef GRID_DIST_UNIT_TEST_HPP
#define GRID_DIST_UNIT_TEST_HPP

#include "grid_dist_id.hpp"
#include "data_type/scalar.hpp"

BOOST_AUTO_TEST_SUITE( grid_dist_id_test )

template<typename iterator> void jacobi_iteration(iterator g_it, grid_dist_id<2, float, scalar<float>, CartDecomposition<2,float>> & g_dist)
{
	// scalar
	typedef scalar<float> S;

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

BOOST_AUTO_TEST_CASE( grid_dist_id_domain_grid_unit_converter_test)
{
	// Domain
	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	// Initialize the global VCluster
	init_global_v_cluster(&boost::unit_test::framework::master_test_suite().argc,&boost::unit_test::framework::master_test_suite().argv);

	Vcluster & v_cl = *global_v_cluster;

	// Test several grid dimensions

	for (size_t k = 1024 ; k >= 2 ; k--)
	{
		std::cout << "Testing: " << k << "\n";

		// grid size
		size_t sz[2];
		sz[0] = k;
		sz[1] = k;

		// Ghost
		Ghost<2,float> g(0.01);

		// Distributed grid with id decomposition
		grid_dist_id<2, float, scalar<float>, CartDecomposition<2,float>> g_dist(sz,domain,g);

		// get the decomposition
		auto & dec = g_dist.getDecomposition();

		// check the consistency of the decomposition
		bool val = dec.check_consistency();
		BOOST_REQUIRE_EQUAL(val,true);

		// for each local volume
		// Get the number of local grid needed
		size_t n_grid = dec.getNLocalHyperCube();

		size_t vol = 0;

		// Allocate the grids
		for (size_t i = 0 ; i < n_grid ; i++)
		{
			// Get the local hyper-cube
			SpaceBox<2,float> sub = dec.getLocalHyperCube(i);

			Box<2,size_t> g_box = g_dist.getCellDecomposer().convertDomainSpaceIntoGridUnits(sub);

			vol += g_box.getVolumeKey();
		}

		v_cl.reduce(vol);
		v_cl.execute();

		BOOST_REQUIRE_EQUAL(vol,sz[0]*sz[1]);
	}
}

BOOST_AUTO_TEST_CASE( grid_dist_id_iterator_test_use)
{
	// Domain
	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	// Initialize the global VCluster
	init_global_v_cluster(&boost::unit_test::framework::master_test_suite().argc,&boost::unit_test::framework::master_test_suite().argv);

	for (long int k = 1024 ; k >= 2 ; k-= (k >= 66)?33:1 )
	{
		std::cout << "Testing: " << k << "\n";

		// grid size
		size_t sz[2];
		sz[0] = k;
		sz[1] = k;

		// Ghost
		Ghost<2,float> g(0.01);

		// Distributed grid with id decomposition
		grid_dist_id<2, float, scalar<float>, CartDecomposition<2,float>> g_dist(sz,domain,g);

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

		// Get the virtual cluster machine
		Vcluster & vcl = g_dist.getVC();

		// reduce
		vcl.reduce(count);
		vcl.execute();

		// Check
		BOOST_REQUIRE_EQUAL(count,k*k);

		auto dom2 = g_dist.getDomainIterator();

		// check that the grid store the correct information
		while (dom2.isNext())
		{
			auto key = dom2.get();
			auto key_g = g_dist.getGKey(key);

			BOOST_REQUIRE_EQUAL(g_dist.template get<0>(key),info.LinId(key_g));

			++dom2;
		}

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
				if (g_dist.template get<0>(key) != info.LinId(key_g))
				{
					int debug = 0;
					debug++;
				}

				BOOST_REQUIRE_EQUAL(g_dist.template get<0>(key),info.LinId(key_g));
			}

			++domg;
		}


	}

//	g_dist.write("");

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
