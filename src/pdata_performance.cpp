/*
 * performance.hpp
 *
 *  Created on: Mar 9, 2016
 *      Author: yaroslav
 */

#ifndef SRC_PDATA_PERFORMANCE_CPP_
#define SRC_PDATA_PERFORMANCE_CPP_

#include <iostream>
#include <mpi.h>
#include "config.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include "Vector/vector_dist.hpp"
#include "data_type/aggregate.hpp"
#include "Plot/GoogleChart.hpp"
#include "Point_test.hpp"
#include <sstream>

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

extern const char * test_dir;

// XML report
boost::property_tree::ptree pt;

#ifdef PERFORMANCE_TEST

BOOST_AUTO_TEST_SUITE( performance )

#include "Vector/performance/vector_dist_performance_util.hpp"
#include "Vector/performance/verlet_performance_tests.hpp"
#include "Vector/performance/cell_list_part_reorder.hpp"
#include "Vector/performance/cell_list_comp_reorder.hpp"


BOOST_AUTO_TEST_SUITE_END()

#endif

#endif /* SRC_PDATA_PERFORMANCE_CPP_ */
