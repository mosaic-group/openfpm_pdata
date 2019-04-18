#define PRINT_RANK_TO_GPU

#include "initialize_wrapper.hpp"
#include "VCluster/VCluster.hpp"

void openfpm_init_wrapper(int * argc, char *** argv)
{
	openfpm_init(argc,argv,init_options::in_situ_visualization);
}

void openfpm_finalize_wrapper()
{
	openfpm_finalize();
}
