#include <iostream>
#include <thread>

size_t debug_tot_call = 0;

#define PRINT_STACKTRACE
#define CHECKFOR_POSNAN
#define CHECKFOR_POSINF
#define CHECKFOR_PROPNAN
#define CHECKFOR_PROPINF

#define NO_WARNING
#include "Graph/CartesianGraphFactory.hpp"

void timeout_cycle()
{
	// 6 seconds
	std::this_thread::sleep_for (std::chrono::seconds(900));

	std::cout << "Time Out" << std::endl;
	std::exit(1);
}


#define BOOST_DISABLE_ASSERTS


#include "config.h"
#undef VERSION

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

// initialization function:
bool init_unit_test()
{
//  std::thread to (timeout_cycle);
//  to.detach();
  return true;
}

// entry point
int main(int argc, char* argv[])
{
  return boost::unit_test::unit_test_main( &init_unit_test, argc, argv );
}

#include "Grid/grid_dist_id.hpp"
#include "Point_test.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "Graph/dist_map_graph.hpp"
#include "memory/HeapMemory.hpp"
#include "Space/Shape/Box.hpp"

#include "unit_test_init_cleanup.hpp"
#include "Graph/CartesianGraphFactory_unit_test.hpp"
#include "Decomposition/ORB_unit_test.hpp"
#include "Decomposition/Distribution/metis_util_unit_test.hpp"
#include "Decomposition/dec_optimizer_unit_test.hpp"
#include "Decomposition/Distribution/Distribution_unit_tests.hpp"
#include "Grid/Iterators/grid_dist_id_iterators_unit_tests.hpp"
//#include "DLB/DLB_unit_test.hpp"
#include "Graph/dist_map_graph_unit_test.hpp"
#include "Graph/DistGraphFactory.hpp"
#include "Vector/se_class3_vector_unit_tests.hpp"
#include "Vector/tests/vector_dist_dlb_test.hpp"
#include "Decomposition/Domain_NN_calculator_cart_unit_test.hpp"
//#include "antoniol_test_isolation.hpp"
