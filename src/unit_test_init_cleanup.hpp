/*
 * unit_test_init_cleanup.hpp
 *
 *  Created on: Apr 17, 2015
 *      Author: Pietro Incardona
 */

#ifndef UNIT_TEST_INIT_CLEANUP_HPP_
#define UNIT_TEST_INIT_CLEANUP_HPP_

#include "VCluster.hpp"

struct ut_start {
    ut_start()   { BOOST_TEST_MESSAGE("Initialize global VCluster"); openfpm_init(&boost::unit_test::framework::master_test_suite().argc,&boost::unit_test::framework::master_test_suite().argv); }
    ~ut_start()  { BOOST_TEST_MESSAGE("Delete global VClster"); openfpm_finalize(); }
};

//____________________________________________________________________________//

BOOST_GLOBAL_FIXTURE( ut_start );



#endif /* UNIT_TEST_INIT_CLEANUP_HPP_ */
