#include <iostream>

#define NO_WARNING
#include "Graph/CartesianGraphFactory.hpp"

#define BOOST_DISABLE_ASSERTS


#include "config.h"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

// initialization function:
bool init_unit_test()
{
  return true;
}

// entry point:
int main(int argc, char* argv[])
{
  return boost::unit_test::unit_test_main( &init_unit_test, argc, argv );
}

#include "debug.hpp"
#include "Grid/grid_dist_id.hpp"
#include "Point_test.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "Graph/dist_map_graph.hpp"
#include "memory/HeapMemory.hpp"
#include "Space/Shape/Box.hpp"
#include "util.hpp"

#include "unit_test_init_cleanup.hpp"
#include "Graph/CartesianGraphFactory_unit_test.hpp"
#include "Decomposition/CartDecomposition_unit_test.hpp"
#include "Decomposition/ORB_unit_test.hpp"
#include "Decomposition/Distribution/metis_util_unit_test.hpp"
#include "dec_optimizer_unit_test.hpp"
#include "Vector/vector_dist_unit_test.hpp"
#include "Decomposition/Distribution/Distribution_unit_tests.hpp"
#include "Grid/Iterators/grid_dist_id_iterators_unit_tests.hpp"
//#include "DLB/DLB_unit_test.hpp"
#include "Graph/dist_map_graph_unit_test.hpp"
#include "Graph/DistGraphFactory.hpp"
#include "Decomposition/nn_processor_unit_test.hpp"
#include "Grid/staggered_grid_dist_unit_test.hpp"
#include "Vector/vector_dist_MP_unit_tests.hpp"
//#include "antoniol_test_isolation.hpp"
