#include <iostream>
#include "config.h"
#include "Graph/CartesianGraphFactory.hpp"

#define BOOST_DISABLE_ASSERTS


#define BOOST_TEST_MODULE "C++ test module for OpenFPM_pdata project"
#include <boost/test/included/unit_test.hpp>

/*struct MPIFixture {
    MPIFixture() { MPI_Init(&boost::unit_test::framework::master_test_suite().argc,&boost::unit_test::framework::master_test_suite().argv);}
    ~MPIFixture() { MPI_Finalize(); }
};

BOOST_GLOBAL_FIXTURE(MPIFixture);*/

#include "Grid/grid_dist_id.hpp"
#include "Point_test.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "memory/HeapMemory.hpp"
#include "Space/Shape/Box.hpp"
#include "util.hpp"

#include "unit_test_init_cleanup.hpp"
#include "Decomposition/ORB_unit_test.hpp"
#include "Graph/CartesianGraphFactory_unit_test.hpp"
#include "metis_util_unit_test.hpp"
#include "dec_optimizer_unit_test.hpp"
#include "Grid/grid_dist_id_unit_test.hpp"
#include "Vector/vector_dist_unit_test.hpp"
#include "Decomposition/CartDecomposition_unit_test.hpp"

