/*
 * grid_dist_id_unit_test_unb_ghost.hpp
 *
 *  Created on: Jul 11, 2016
 *      Author: i-bird
 */

#ifndef SRC_GRID_GRID_DIST_ID_UNIT_TEST_UNB_GHOST_HPP_
#define SRC_GRID_GRID_DIST_ID_UNIT_TEST_UNB_GHOST_HPP_

void Test3D_unb_ghost(const Box<3,float> & domain, long int k)
{
	long int big_step = k / 30;
	big_step = (big_step == 0)?1:big_step;
	long int small_step = 21;

	if (create_vcluster().getProcessingUnits() > 48)
		return;

	print_test( "Testing 3D grid unbound ghost k<=",k);

	// 3D test
	for ( ; k >= 2 ; k-= (k > 2*big_step)?big_step:small_step )
	{
		BOOST_TEST_CHECKPOINT( "Testing 3D grid k=" << k );

		// grid size
		size_t sz[3];
		sz[0] = k;
		sz[1] = k;
		sz[2] = k;

		// Ghost
		Ghost<3,float> g(0.49);

		// Distributed grid with id decomposition
		grid_dist_id<3, float, aggregate<float>, CartDecomposition<3,float>> g_dist(sz,domain,g);

		g_dist.getDecomposition().write("no_bound_decomposition");

		// check the consistency of the decomposition
		bool val = g_dist.getDecomposition().check_consistency();
		BOOST_REQUIRE_EQUAL(val,true);

		// Grid sm
		grid_sm<3,void> info(sz);

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
		vcl.sum(count);
		vcl.execute();

		// Check
		BOOST_REQUIRE_EQUAL(count,(size_t)k*k*k);

		bool match = true;

		auto dom2 = g_dist.getDomainIterator();

		// check that the grid store the correct information
		while (dom2.isNext())
		{
			auto key = dom2.get();
			auto key_g = g_dist.getGKey(key);

			match &= (g_dist.template get<0>(key) == info.LinId(key_g))?true:false;

			++dom2;
		}

		BOOST_REQUIRE_EQUAL(match,true);

		//! [Synchronize the ghost and check the information]

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

		//! [Synchronize the ghost and check the information]
	}
}


// Test grid periodic
void Test3D_unb_ghost_periodic(const Box<3,float> & domain, long int k)
{
	Vcluster & v_cl = create_vcluster();

	if ( v_cl.getProcessingUnits() > 24 )
		return;

	long int big_step = k / 30;
	big_step = (big_step == 0)?1:big_step;
	long int small_step = 21;

	print_test( "Testing grid periodic unbound ghost k<=",k);

	// 3D test
	for ( ; k >= 2 ; k-= (k > 2*big_step)?big_step:small_step )
	{
		BOOST_TEST_CHECKPOINT( "Testing grid unbound ghost periodic k<=" << k );

		// grid size
		size_t sz[3];
		sz[0] = k;
		sz[1] = k;
		sz[2] = k;

		// Ghost
		Ghost<3,float> g(0.491);

		// periodicity
		periodicity<3> pr = {{PERIODIC,PERIODIC,PERIODIC}};

		// Distributed grid with id decomposition
		grid_dist_id<3, float, aggregate<long int>, CartDecomposition<3,float>> g_dist(sz,domain,g,pr);

		// check the consistency of the decomposition
		bool val = g_dist.getDecomposition().check_consistency();
		BOOST_REQUIRE_EQUAL(val,true);

		// Grid sm
		grid_sm<3,void> info(sz);

		size_t count = 0;

		// Set to zero the full grid

		auto dom1 = g_dist.getDomainGhostIterator();

		while (dom1.isNext())
		{
			auto key = dom1.get();

			g_dist.template get<0>(key) = -1;

			++dom1;
		}

		auto dom = g_dist.getDomainIterator();

		while (dom.isNext())
		{
			auto key = dom.get();
			auto key_g = g_dist.getGKey(key);

			g_dist.template get<0>(key) = info.LinId(key_g);

			// Count the points
			count++;

			++dom;
		}

		// Get the virtual cluster machine
		Vcluster & vcl = g_dist.getVC();

		// reduce
		vcl.sum(count);
		vcl.execute();

		// Check
		BOOST_REQUIRE_EQUAL(count,(size_t)k*k*k);

		// sync the ghosts
		g_dist.ghost_get<0>();

		bool match = true;

		// Domain + Ghost iterator
		auto dom_gi = g_dist.getDomainGhostIterator();

		size_t out_cnt = 0;

		while (dom_gi.isNext())
		{
			bool out_p = false;

			auto key = dom_gi.get();
			auto key_g = g_dist.getGKey(key);

			// Return the boxes containing the grids
			auto & gb = dom_gi.getGBoxes();

			// transform the key to be periodic
			for (size_t i = 0 ; i < 3 ; i++)
			{
				if (key_g.get(i) < 0)
				{key_g.set_d(i,key_g.get(i) + k);out_p = true;}
				else if (key_g.get(i) >= k)
				{key_g.set_d(i,key_g.get(i) - k);out_p = true;}
			}

			if (g_dist.template get<0>(key) != -1 && out_p == true)
				out_cnt++;

			// The last points can be invalid because of rounding off problems
			bool can_invalid = false;
			if (key.getKey().get(0) == 0 || key.getKey().get(1) == 0 || key.getKey().get(2) == 0)
				can_invalid = true;
			else if (key.getKey().get(0) == gb.get(key.getSub()).GDbox.getHigh(0) ||
					 key.getKey().get(1) == gb.get(key.getSub()).GDbox.getHigh(1) ||
					 key.getKey().get(2) == gb.get(key.getSub()).GDbox.getHigh(2))
				can_invalid = true;

			if (can_invalid == true)
			{
				if ( g_dist.template get<0>(key) != -1 && info.LinId(key_g) != g_dist.template get<0>(key) )
					match &= false;
			}
			else
			{
				if (info.LinId(key_g) != g_dist.template get<0>(key) )
					match &= false;
			}

			++dom_gi;
		}

		BOOST_REQUIRE_EQUAL(match, true);
		if (k > 83)
		{
			vcl.sum(out_cnt);
			vcl.execute();
			BOOST_REQUIRE(out_cnt != 0ul);
		}
	}
}


#endif /* SRC_GRID_GRID_DIST_ID_UNIT_TEST_UNB_GHOST_HPP_ */
