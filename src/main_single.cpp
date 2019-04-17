#include <iostream>

size_t debug_tot_call = 0;

#define PRINT_STACKTRACE
#define CHECKFOR_POSNAN
#define CHECKFOR_POSINF
#define CHECKFOR_PROPNAN
#define CHECKFOR_PROPINF

#define BOOST_DISABLE_ASSERTS


#include "config.h"
#undef VERSION

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

#include "unit_test_init_cleanup.hpp"
