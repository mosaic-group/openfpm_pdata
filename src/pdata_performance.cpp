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


#ifdef PERFORMANCE_TEST

#include "../../openfpm_numerics/src/interpolation/interpolation.hpp"
#include "Grid/grid_dist_id.hpp"
#include "Plot/GoogleChart.hpp"
#include "interpolation/mp4_kernel.hpp"
#include "Vector/performance/vector_dist_performance_util.hpp"
#include "util/performance/performance_util.hpp"

BOOST_AUTO_TEST_SUITE( performance )

#include "Vector/performance/verlet_performance_tests.hpp"
#include "Vector/performance/cell_list_part_reorder.hpp"
#include "Vector/performance/cell_list_comp_reorder.hpp"
#include "Vector/performance/vector_dist_gg_map_performance.hpp"
#include "Grid/performance/grid_dist_performance.hpp"

BOOST_AUTO_TEST_SUITE_END()

#endif

#endif /* SRC_PDATA_PERFORMANCE_CPP_ */
