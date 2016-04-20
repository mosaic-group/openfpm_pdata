#include <iostream>
#include "config.h"

#define NO_WARNING
#include "Graph/CartesianGraphFactory.hpp"

#define BOOST_DISABLE_ASSERTS


#define BOOST_TEST_MODULE "C++ test module for OpenFPM_pdata project"
#include <boost/test/included/unit_test.hpp>

#include "Grid/grid_dist_id.hpp"
#include "Point_test.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "memory/HeapMemory.hpp"
#include "Space/Shape/Box.hpp"
#include "util.hpp"

#include "unit_test_init_cleanup.hpp"
#include "Grid/staggered_grid_dist_unit_test.hpp"
#include "Decomposition/CartDecomposition_unit_test.hpp"
#include "Decomposition/ORB_unit_test.hpp"
#include "Graph/CartesianGraphFactory_unit_test.hpp"
#include "metis_util_unit_test.hpp"
#include "dec_optimizer_unit_test.hpp"
#include "Grid/grid_dist_id_unit_test.hpp"
#include "Vector/vector_dist_unit_test.hpp"
//#ifdef PERFORMANCE_TEST
#include "pdata_performance.hpp"
//#endif
//#include "Decomposition/nn_processor_unit_test.hpp"

