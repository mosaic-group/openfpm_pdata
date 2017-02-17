/*
 * grid_dist_id_iterators_unit_tests.hpp
 *
 *  Created on: Jan 4, 2017
 *      Author: i-bird
 */

#ifndef SRC_GRID_ITERATORS_GRID_DIST_ID_ITERATORS_UNIT_TESTS_HPP_
#define SRC_GRID_ITERATORS_GRID_DIST_ID_ITERATORS_UNIT_TESTS_HPP_

#include "grid_dist_id_iterator_dec_skin.hpp"

BOOST_AUTO_TEST_SUITE( grid_dist_id_iterators_test )

void print_test(std::string test, size_t sz)
{
	if (create_vcluster().getProcessUnitID() == 0)
		std::cout << test << " " << sz << "\n";
}

void Test2D_sub(const Box<2,float> & domain, long int k)
{
	long int big_step = k / 30;
	big_step = (big_step == 0)?1:big_step;
	long int small_step = 21;

	// this test is only performed when the number of processor is <= 32
	if (create_vcluster().getProcessingUnits() > 32)
		return;

	print_test( "Testing 2D grid sub iterator k<=",k);

	// 2D test
	for ( ; k >= 2 ; k-= (k > 2*big_step)?big_step:small_step )
	{
		BOOST_TEST_CHECKPOINT( "Testing 2D grid k=" << k );

		// grid size
		size_t sz[2];
		sz[0] = k;
		sz[1] = k;

		float factor = pow(create_vcluster().getProcessingUnits()/2.0f,1.0f/2.0f);

		// Ghost
		Ghost<2,float> g(0.01 / factor);

		// Distributed grid with id decomposition
		grid_dist_id<2, float, scalar<float>> g_dist(sz,domain,g);

		// check the consistency of the decomposition
		bool val = g_dist.getDecomposition().check_consistency();
		BOOST_REQUIRE_EQUAL(val,true);

		size_t count;

		// Grid sm
		grid_sm<2,void> info(sz);

		{
		//! [Usage of a sub_grid iterator]

		grid_key_dx<2> one(1,1);
		grid_key_dx<2> one_end(k-2,k-2);

		bool check = true;
		count = 0;

		// get the sub-domain iterator
		auto dom = g_dist.getSubDomainIterator(one,one_end);

		while (dom.isNext())
		{
			auto key = dom.get();
			auto key_g = g_dist.getGKey(key);

			// key_g should never be 1 or k-1
			check &= (key_g.get(0) == 0 || key_g.get(0) == k-1)?false:true;
			check &= (key_g.get(1) == 0 || key_g.get(1) == k-1)?false:true;

			g_dist.template get<0>(key) = info.LinId(key_g);

			// Count the point
			count++;

			++dom;
		}

		BOOST_REQUIRE_EQUAL(check,true);

		//! [Usage of a sub_grid iterator]

		}

		// Get the virtual cluster machine
		Vcluster & vcl = g_dist.getVC();

		// reduce
		vcl.sum(count);
		vcl.execute();

		// Check
		BOOST_REQUIRE_EQUAL(count,(size_t)(k-2)*(k-2));

		// check with a 1x1 square

		{

		grid_key_dx<2> one(k/2,k/2);
		grid_key_dx<2> one_end(k/2,k/2);

		count = 0;

		// get the sub-domain iterator
		auto dom = g_dist.getSubDomainIterator(one,one_end);

		while (dom.isNext())
		{
			auto key = dom.get();
			auto key_g = g_dist.getGKey(key);

			// key_g
			BOOST_REQUIRE_EQUAL(key_g.get(0),k/2);
			BOOST_REQUIRE_EQUAL(key_g.get(1),k/2);

			auto key_s_it = dom.getGKey(key);

			BOOST_REQUIRE_EQUAL(key_g.get(0),key_s_it.get(0));
			BOOST_REQUIRE_EQUAL(key_g.get(1),key_s_it.get(1));

			// Count the point
			count++;

			++dom;
		}

		// reduce
		vcl.sum(count);
		vcl.execute();

		BOOST_REQUIRE_EQUAL(count,1ul);
		}
	}
}

// Test decomposition grid iterator

void Test3D_decit(const Box<3,float> & domain, long int k)
{
	size_t k_bck = k;
	{
		Vcluster & v_cl = create_vcluster();

		if ( v_cl.getProcessingUnits() > 32 )
			return;

		long int big_step = k / 30;
		big_step = (big_step == 0)?1:big_step;
		long int small_step = 21;

		print_test( "Testing grid iterator from decomposition k<=",k);

		// 3D test
		for ( ; k >= 2 ; k-= (k > 2*big_step)?big_step:small_step )
		{
			BOOST_TEST_CHECKPOINT( "Testing grid iterator from decomposition k<=" << k );

			// grid size
			size_t sz[3];
			sz[0] = k;
			sz[1] = k;
			sz[2] = k;

			// factor
			float factor = pow(create_vcluster().getProcessingUnits()/2.0f,1.0f/3.0f);

			// Ghost
			Ghost<3,float> g(0.01 / factor);

			// Distributed grid with id decomposition
			grid_dist_id<3, float, Point_test<float>, CartDecomposition<3,float>> g_dist(sz,domain,g);

			// check the consistency of the decomposition
			bool val = g_dist.getDecomposition().check_consistency();
			BOOST_REQUIRE_EQUAL(val,true);

			// Grid sm
			grid_sm<3,void> info(sz);

			auto dom = g_dist.getDomainIterator();

			bool match = true;

			// create a grid iterator from the decomposition

			grid_dist_id_iterator_dec<CartDecomposition<3,float>> it_dec(g_dist.getDecomposition(),g_dist.getGridInfoVoid().getSize());

			while (dom.isNext())
			{
				auto key = dom.get();
				auto key_g = g_dist.getGKey(key);

				auto key_dec = it_dec.get();

				// Check if the two keys match
				match &= (key_dec == key_g);

				++dom;
				++it_dec;
			}

			BOOST_REQUIRE_EQUAL(match,true);
		}
	}

	k = k_bck;

	{
		Vcluster & v_cl = create_vcluster();

		if ( v_cl.getProcessingUnits() > 32 )
			return;

		long int big_step = k / 30;
		big_step = (big_step == 0)?1:big_step;
		long int small_step = 21;

		print_test( "Testing grid iterator from decomposition subset k<=",k);

		// 3D test
		for ( ; k >= 2 ; k-= (k > 2*big_step)?big_step:small_step )
		{
			BOOST_TEST_CHECKPOINT( "Testing grid iterator from decomposition k<=" << k );

			// grid size
			size_t sz[3];
			sz[0] = k;
			sz[1] = k;
			sz[2] = k;

			// factor
			float factor = pow(create_vcluster().getProcessingUnits()/2.0f,1.0f/3.0f);

			// Ghost
			Ghost<3,float> g(0.01 / factor);

			// Distributed grid with id decomposition
			grid_dist_id<3, float, Point_test<float>, CartDecomposition<3,float>> g_dist(sz,domain,g);

			// check the consistency of the decomposition
			bool val = g_dist.getDecomposition().check_consistency();
			BOOST_REQUIRE_EQUAL(val,true);

			// Grid sm
			grid_sm<3,void> info(sz);

			auto dom = g_dist.getSubDomainIterator({0,0,0},{(long int)sz[0]-2,(long int)sz[1]-2,(long int)sz[2]-2});

			bool match = true;

			// create a grid iterator from the decomposition

			grid_dist_id_iterator_dec<CartDecomposition<3,float>> it_dec(g_dist.getDecomposition(),sz,{0,0,0},{sz[0]-2,sz[1]-2,sz[2]-2});

			while (dom.isNext())
			{
				auto key = dom.get();
				auto key_g = g_dist.getGKey(key);

				auto key_dec = it_dec.get();

				// Check if the two keys match
				match &= (key_dec == key_g);

				++dom;
				++it_dec;
			}

			BOOST_REQUIRE_EQUAL(match,true);
		}
	}
}


// Test decomposition grid iterator

void Test3D_decskinit(const Box<3,float> & domain, long int k)
{
	{
		Vcluster & v_cl = create_vcluster();

		if ( v_cl.getProcessingUnits() > 32 )
			return;

		long int big_step = k / 30;
		big_step = (big_step == 0)?1:big_step;
		long int small_step = 21;

		print_test( "Testing grid skin iterator from decomposition k<=",k);

		// 3D test
		for ( ; k >= 2 ; k-= (k > 2*big_step)?big_step:small_step )
		{
			BOOST_TEST_CHECKPOINT( "Testing grid skin iterator from decomposition k<=" << k );

			// grid size
			size_t sz[3];
			sz[0] = k;
			sz[1] = k;
			sz[2] = k;

			if (k <= 9)
				continue;

			// factor
			float factor = pow(create_vcluster().getProcessingUnits()/2.0f,1.0f/3.0f);

			// Ghost
			Ghost<3,float> g(0.01 / factor);

			// Distributed grid with id decomposition
			grid_dist_id<3, float, Point_test<float>, CartDecomposition<3,float>> g_dist(sz,domain,g);

			// check the consistency of the decomposition
			bool val = g_dist.getDecomposition().check_consistency();
			BOOST_REQUIRE_EQUAL(val,true);

			// Grid sm
			grid_sm<3,void> info(sz);

			// create a grid skin iterator from the decomposition

			Box<3,size_t> A({3,3,3},{(size_t)k-3,(size_t)k-3,(size_t)k-3});
			Box<3,size_t> B = A;

			if (A.isValid() == false)
				continue;

			size_t bc[3] = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};
			grid_dist_id_iterator_dec_skin<CartDecomposition<3,float>> it_dec(g_dist.getDecomposition(),g_dist.getGridInfoVoid(),A,B,bc);

			size_t cnt = 0;

			bool tot_good = true;
			while (it_dec.isNext())
			{
				auto key_dec = it_dec.get();

				// one of the coordinate has to be or 3 or 8, none of
				// None of the coordinates must be bigger that

				bool eight_or_three = false;
				bool good = true;
				for (size_t i = 0; i < 3 ; i++)
				{
					if (key_dec.get(i) == 3 || key_dec.get(i) == k - 3)
						eight_or_three = true;

					if (key_dec.get(i) > k - 3 || key_dec.get(i) < 3 )
						good = false;
				}

				tot_good &= (eight_or_three) || good;

				cnt++;
				++it_dec;
			}

			create_vcluster().sum(cnt);
			create_vcluster().execute();

			BOOST_REQUIRE_EQUAL(cnt,(size_t)((k-5)*(k-5)*(k-5) - (k-7)*(k-7)*(k-7)));
			BOOST_REQUIRE_EQUAL(tot_good,true);
		}
	}
}

BOOST_AUTO_TEST_CASE( grid_dist_id_sub_iterator_test_use)
{
	// Domain
	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	long int k = 1024*1024*create_vcluster().getProcessingUnits();
	k = std::pow(k, 1/2.);

	Test2D_sub(domain,k);
}

BOOST_AUTO_TEST_CASE( grid_dist_id_decomposition_iterator )
{
	// Domain
	Box<3,float> domain3({0.0,0.0,0.0},{1.0,1.0,1.0});

	size_t k = 128*128*128*create_vcluster().getProcessingUnits();
	k = std::pow(k, 1/3.);
	Test3D_decit(domain3,k);
}

BOOST_AUTO_TEST_CASE( grid_dist_it_iterators_skin_test )
{
	// Domain
	Box<3,float> domain3({0.0,0.0,0.0},{1.0,1.0,1.0});

	size_t k = 128*128*128*create_vcluster().getProcessingUnits();
	k = std::pow(k, 1/3.);
	Test3D_decskinit(domain3,k);
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* SRC_GRID_ITERATORS_GRID_DIST_ID_ITERATORS_UNIT_TESTS_HPP_ */
