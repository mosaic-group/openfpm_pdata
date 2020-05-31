#define PRINT_RANK_TO_GPU


#include <hip/hip_runtime.h>
#include "initialize_wrapper.hpp"
#include "VCluster/VCluster.hpp"

void openfpm_init_wrapper(int * argc, char *** argv)
{
	openfpm_init(argc,argv);
}

void openfpm_finalize_wrapper()
{
	openfpm_finalize();
}
