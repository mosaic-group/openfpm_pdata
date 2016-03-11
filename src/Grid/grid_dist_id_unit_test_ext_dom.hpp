/*
 * grid_dist_id_unit_test_ext_dom.hpp
 *
 *  Created on: Feb 24, 2016
 *      Author: i-bird
 */

#ifndef SRC_GRID_GRID_DIST_ID_UNIT_TEST_EXT_DOM_HPP_
#define SRC_GRID_GRID_DIST_ID_UNIT_TEST_EXT_DOM_HPP_

#include "Decomposition/CartDecomposition_ext.hpp"

// Test duplicated topology

void Test3D_extended_grid(const Box<3,float> & domain, long int k)
{
	long int big_step = k / 30;
	big_step = (big_step == 0)?1:big_step;
	long int small_step = 21;

	Vcluster & v_cl = *global_v_cluster;

	if ( v_cl.getProcessingUnits() > 32 )
		return;

	print_test( "Testing 3D extended grid k<=",k);

	// 3D test
	for ( ; k >= 2 ; k-= (k > 2*big_step)?big_step:small_step )
	{
		BOOST_TEST_CHECKPOINT( "Testing 3D extended grid k=" << k );

		// grid size
		size_t sz[3];
		sz[0] = k;
		sz[1] = k;
		sz[2] = k;

		// factor
		float factor = pow(global_v_cluster->getProcessingUnits()/2.0f,1.0f/3.0f);

		// Ghost
		Ghost<3,float> g(0.01 / factor);

		//! [Construct an extended grid]

		// Distributed grid with id decomposition
		grid_dist_id<3, float, aggregate<size_t[3],size_t>, CartDecomposition<3,float>> g_dist1(sz,domain,g);

		// Extend the grid by 2 points
		Box<3,size_t> ext({2,2,2},{2,2,2});

		// another grid perfectly overlapping the previous, extended by 2 points
		grid_dist_id<3, float, aggregate<size_t[3],size_t>, CartDecomposition_ext<3,float>> g_dist2(g_dist1,ext);

		// Given an iterator on grid 1
		auto dom_g1 = g_dist1.getDomainIterator();
		// And a sub-iterator on grid 2 overlapping grid 1
		auto dom_g2 = g_dist2.getSubDomainIterator({0,0,0},{k-1,k-1,k-1});

		// the 2 iterator must match

		bool check = true;

		while (dom_g2.isNext())
		{
			auto key1 = dom_g1.get();
			auto key2 = dom_g2.get();

			grid_key_dx<3> g1_k = g_dist1.getGKey(key1);
			grid_key_dx<3> g2_k = g_dist2.getGKey(key2);

			check &= (g1_k == g2_k)?true:false;

			++dom_g1;
			++dom_g2;
		}

		BOOST_REQUIRE_EQUAL(check,true);

		//! [Construct an extended grid]

		bool ret = g_dist2.getDecomposition().check_consistency();
		BOOST_REQUIRE_EQUAL(ret,true);

		// Get domain iterator

		grid_sm<3,void> info = g_dist2.getGridInfo();

		size_t cnt = 0;
		auto dom_g3 = g_dist2.getDomainIterator();

		check = false;

		while (dom_g3.isNext())
		{
			auto key1 = dom_g3.get();

			auto keyg = g_dist2.getGKey(key1);

			g_dist2.template get<0>(key1)[0] = keyg.get(0);
			g_dist2.template get<0>(key1)[1] = keyg.get(1);
			g_dist2.template get<0>(key1)[2] = keyg.get(2);

			g_dist2.template get<1>(key1) = info.LinId(keyg);

			++dom_g3;
		}

		g_dist2.ghost_get<0,1>();

		auto dom_g4 = g_dist2.getSubDomainIterator({-1,-1,-1},{(long int) sz[0]+2-2, (long int) sz[1]+2-2, (long int) sz[2]+2-2});

		check = true;

		while (dom_g4.isNext())
		{
			auto key1 = dom_g4.get();

			key1 = key1.move(0,1);
			key1 = key1.move(1,1);
			key1 = key1.move(2,1);

			auto key2 = g_dist2.getGKey(key1);

			check &= g_dist2.template get<0>(key1)[0] == key2.get(0);
			check &= g_dist2.template get<0>(key1)[1] == key2.get(1);
			check &= g_dist2.template get<0>(key1)[2] == key2.get(2);

			if (check == false)
			{
				int debug = 0;
				debug++;
			}

			auto key3 = dom_g4.get();

			key3 = key3.move(0,-1);
			key3 = key3.move(1,-1);
			key3 = key3.move(2,-1);

			auto key4 = g_dist2.getGKey(key3);

			check &= g_dist2.template get<0>(key3)[0] == key4.get(0);
			check &= g_dist2.template get<0>(key3)[1] == key4.get(1);
			check &= g_dist2.template get<0>(key3)[2] == key4.get(2);

			if (check == false)
			{
				int debug = 0;
				debug++;
			}

			++dom_g4;
		}

		std::cout << "k=" << k << "\n";
		BOOST_REQUIRE_EQUAL(check,true);
	}
}


#endif /* SRC_GRID_GRID_DIST_ID_UNIT_TEST_EXT_DOM_HPP_ */
