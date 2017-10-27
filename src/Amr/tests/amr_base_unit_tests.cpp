/*
 * amr_base_unit_test.cpp
 *
 *  Created on: Oct 5, 2017
 *      Author: i-bird
 */
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Grid/grid_dist_id.hpp"
#include "Point_test.hpp"
#include "Grid/tests/grid_dist_id_util_tests.hpp"

BOOST_AUTO_TEST_SUITE( amr_grid_dist_id_test )


BOOST_AUTO_TEST_CASE( grid_dist_id_amr )
{
	// Domain
	Box<2,float> domain2({0.0,0.0},{1.0,1.0});

	size_t sz[2] = {100,100};

	// Ghost
	Ghost<2,long int> g(1);

	// periodicity
	periodicity<2> pr = {{PERIODIC,PERIODIC}};

	openfpm::vector<Box<2,long int>> C_draw;
	C_draw.add(Box<2,long int>({20,20},{50,24}));
	C_draw.add(Box<2,long int>({51,20},{60,24}));
	C_draw.add(Box<2,long int>({61,20},{70,24}));
	C_draw.add(Box<2,long int>({20,25},{24,66}));
	C_draw.add(Box<2,long int>({15,67},{49,85}));
	C_draw.add(Box<2,long int>({50,76},{70,81}));
	C_draw.add(Box<2,long int>({30,25},{34,37}));
	C_draw.add(Box<2,long int>({50,66},{70,70}));

	size_t volume_key = 0;
	for (size_t i = 0 ; i < C_draw.size() ; i++)
	{
		volume_key += Box<2,long int>(C_draw.get(i)).getVolumeKey();
	}

	// Distributed grid with id decomposition
	grid_dist_id<2,float,Point_test<float>> g_dist(sz,domain2,g,pr,C_draw);

	// fill with gkey

	auto git = g_dist.getDomainIterator();
	grid_sm<2,void> gs(sz);

	size_t count = 0;

	while (git.isNext())
	{
		auto key = git.get();
		auto gkey = git.getGKey(key);

		g_dist.template get<0>(key) = gs.LinId(gkey);

		count++;

		++git;
	}

	Vcluster & vcl = create_vcluster();

	vcl.sum(count);
	vcl.execute();

	BOOST_REQUIRE_EQUAL(count,volume_key);

	g_dist.ghost_get<0>();

	// Check it is correct

	bool check = true;
	size_t check_count = 0;

	auto git2 = g_dist.getDomainGhostIterator();
	while (git2.isNext())
	{
		auto key = git2.get();
		auto gkey = git2.getGKey(key);

		float value = g_dist.template get<0>(key);

		// check if the point is inside or outside the domain

		for (size_t k = 0; k < C_draw.size() ; k++)
		{
			if (Box<2,long int>(C_draw.get(k)).isInside(gkey.toPoint()) == true)
			{
				check &= value == gs.LinId(gkey);

				// get the gdb_ext
				auto & gdb_ext = g_dist.getLocalGridsInfo();

				for (size_t s = 0 ; s < gdb_ext.size() ; s++)
				{
					Box<2,long int> bx = gdb_ext.get(s).Dbox;
					bx += gdb_ext.get(s).origin;
					if (bx.isInside(gkey.toPoint()))
					{
						check_count++;
						break;
					}
				}
				break;
			}
		}

		++git2;
	}

	vcl.sum(check_count);
	vcl.execute();

	BOOST_REQUIRE_EQUAL(check,true);
	BOOST_REQUIRE(check_count >= volume_key);
}

BOOST_AUTO_TEST_CASE( amr_grid_dist_id_iterator_test_use_2D)
{
	// Domain
	Box<2,float> domain({0.0,0.0},{1.0,1.0});

#ifdef TEST_COVERAGE_MODE
	long int k = 256*256*create_vcluster().getProcessingUnits();
#else
	long int k = 1024*1024*create_vcluster().getProcessingUnits();
#endif
	k = std::pow(k, 1/2.);

	long int big_step = k / 30;
	big_step = (big_step == 0)?1:big_step;
	long int small_step = 21;

	print_test( "AMR Testing 2D full grid k<=",k);

	// 2D test
	for ( ; k >= 2 ; k-= (k > 2*big_step)?big_step:small_step )
	{
		BOOST_TEST_CHECKPOINT( "AMR Testing 2D full grid k=" << k );

		//! [Create and access a distributed grid]

		// grid size
		size_t sz[2];
		sz[0] = k;
		sz[1] = k;

		// periodicity
		periodicity<2> pr = {{PERIODIC,PERIODIC}};

		float factor = pow(create_vcluster().getProcessingUnits()/2.0f,1.0f/2.0f);

		// Ghost
		Ghost<2,long int> g(1);

		openfpm::vector<Box<2,long int>> bx_def;
		bx_def.add(Box<2,long int>({0,0},{k-1,k-1}));

		// Distributed grid with id decomposition
		grid_dist_id<2, float, scalar<float>> g_dist(sz,domain,g,pr,bx_def);

		Test2D_core(g_dist,sz,k);
	}
}

BOOST_AUTO_TEST_SUITE_END()
