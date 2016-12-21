/*
 * unit_test_init_cleanup.hpp
 *
 *  Created on: Apr 17, 2015
 *      Author: Pietro Incardona
 */

#ifndef UNIT_TEST_INIT_CLEANUP_HPP_
#define UNIT_TEST_INIT_CLEANUP_HPP_

#include "VCluster/VCluster.hpp"

const char * test_dir;

struct ut_start
{
	//!
    ut_start()
    {
    	BOOST_TEST_MESSAGE("Initialize global VCluster");

    	openfpm_init(&boost::unit_test::framework::master_test_suite().argc,&boost::unit_test::framework::master_test_suite().argv);

#ifdef PERFORMANCE_TEST
    	test_dir = getenv("OPENFPM_PERFORMANCE_TEST_DIR");

    	if (test_dir == NULL)
    	{
    		std::cerr << "Error: " __FILE__ << ":" << __LINE__ << " in order to run the performance test you must set the environment variable $OPENFPM_PERFORMANCE_TEST_DIR to the test or an empty directory";
    		exit(1);
    	}
#endif
    }

    ~ut_start()
    {
    	BOOST_TEST_MESSAGE("Delete global VClster");
    	openfpm_finalize();
    }
};

//____________________________________________________________________________//

BOOST_GLOBAL_FIXTURE( ut_start );



#endif /* UNIT_TEST_INIT_CLEANUP_HPP_ */
