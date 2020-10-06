/*
 * grid_dist_id_HDF5_chckpnt_restart_test.hpp
 *
 *  Created on: Nov 9, 2016
 *      Author: Yaroslav Zaluzhnyi
 */

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "VCluster/py/py_runner.hpp"
#include "Grid/grid_dist_id.hpp"

void client_thread()
{
	system("python py_test/emulate_client.py");
}

BOOST_AUTO_TEST_SUITE( py_server_grid_test )

BOOST_AUTO_TEST_CASE( server_test )
{
	// Test grid periodic

	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	Vcluster<> & v_cl = create_vcluster();

	if ( v_cl.getProcessingUnits() > 32 )
	{return;}

	long int k = 13;

	BOOST_TEST_CHECKPOINT( "Testing grid periodic k<=" << k );

	// grid size
	size_t sz[3];
	sz[0] = k;
	sz[1] = k;
	sz[2] = k;

	// Ghost
	Ghost<3,long int> g(1);

	// periodicity
	periodicity<3> pr = {{PERIODIC,PERIODIC,PERIODIC}};

	// Distributed grid with id decomposition
	grid_dist_id<3, float, aggregate<float>> g_dist(sz,domain,g,pr);

	// launch a client thread

	auto lg = g_dist.getGridInfoVoid();
	auto it = g_dist.getDomainIterator();

	while(it.isNext())
	{
		auto key = it.get();
		auto gkey = it.getGKey(key);

		g_dist.template get<0>(key) = lg.LinId(gkey);

		++it;
	}

	do_breakpoint();

/*	for (int i = 1 ; i < 100; i++)
	{
		auto it = g_dist.getDomainIterator();

		while(it.isNext())
		{
			auto key = it.get();

			g_dist.template get<0>(key) = (float)g.LinId(key) / i;

			++it;
		}

		step_done();
	}*/
}

BOOST_AUTO_TEST_SUITE_END()


