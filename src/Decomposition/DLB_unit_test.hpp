/*
 * DLB_unit_test.hpp
 *
 *  Created on: Nov 24, 2015
 *      Author: Antonio Leo
 */

#ifndef DLB_UNIT_TEST_HPP_
#define DLB_UNIT_TEST_HPP_

#include "DLB.hpp"

BOOST_AUTO_TEST_SUITE (DLB_test)

BOOST_AUTO_TEST_CASE( DLB_test_use)
{
	// Vcluster
	Vcluster & vcl = *global_v_cluster;

	// Initialize the global VCluster
	init_global_v_cluster(&boost::unit_test::framework::master_test_suite().argc,&boost::unit_test::framework::master_test_suite().argv);

	// Init DLB tool
	DLB dlb(vcl);

	// Init dlb parameters
	dlb.setComputationCost(5);
	dlb.setSimulationStartTime(0);
	dlb.setSimulationEndTime(10);

	// Time value of the "unbalanced" process
	float t_high = 1;

	for(float t = dlb.getSimulationStartTime();
			t < dlb.getSimulationEndTime();
			t++)
	{
		dlb.startIteration(0);

		if(vcl.getProcessUnitID() == 0)
			dlb.endIteration(1);
		else{
			dlb.endIteration(t_high++);
		}


		bool rebalance = dlb.rebalanceNeeded();

		if(rebalance)
			t_high = 1;

		if(t == 6)
			BOOST_REQUIRE_EQUAL(rebalance,true);
	}

}

BOOST_AUTO_TEST_SUITE_END()

#endif /* DLB_UNIT_TEST_HPP_ */
