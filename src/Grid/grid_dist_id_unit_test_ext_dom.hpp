/*
 * grid_dist_id_unit_test_ext_dom.hpp
 *
 *  Created on: Feb 24, 2016
 *      Author: i-bird
 */

#ifndef SRC_GRID_GRID_DIST_ID_UNIT_TEST_EXT_DOM_HPP_
#define SRC_GRID_GRID_DIST_ID_UNIT_TEST_EXT_DOM_HPP_


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
		grid_dist_id<3, float, Point_test<float>, CartDecomposition<3,float>> g_dist1(sz,domain,g);

		// Extend the grid by 2 points
		Box<3,size_t> ext({2,2,2},{2,2,2});

		// another grid perfectly overlapping the previous, extended by 2 points
		grid_dist_id<3, float, Point_test<float>, CartDecomposition<3,float>> g_dist2(g_dist1,ext);

		//! [Construct an extended grid]

		bool ret = g_dist2.getDecomposition().check_consistency();
		BOOST_REQUIRE_EQUAL(ret,true);

		// Given an iterator on grid 1
		auto dom_g1 = g_dist1.getDomainIterator();
		// And a sub-iterator on grid 2 overlapping grid 1
		auto dom_g2 = g_dist2.getSubDomainIterator({2,2,2},{k+2-1,k+2-1,k+2-1});

		grid_key_dx<3> kb({2l,2l,2l});

		// the 2 iterator must match

		bool check = true;

		while (dom_g2.isNext())
		{
			auto key1 = dom_g1.get();
			auto key2 = dom_g2.get();

			grid_key_dx<3> g1_k = g_dist1.getGKey(key1);
			grid_key_dx<3> g2_k = g_dist2.getGKey(key2);

			g2_k = g2_k - kb;

			check &= (g1_k == g2_k)?true:false;

			std::cout << "KEY: " << g1_k.to_string() << "   " << g2_k.to_string() << "\n";

			if (check == false)
			{
				std::cout << "ERROR: " << g1_k.to_string() << "   " << g2_k.to_string() << "\n";
				break;
			}

			++dom_g1;
			++dom_g2;
		}

		BOOST_REQUIRE_EQUAL(check,true);
	}
}


#endif /* SRC_GRID_GRID_DIST_ID_UNIT_TEST_EXT_DOM_HPP_ */
