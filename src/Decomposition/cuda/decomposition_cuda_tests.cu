#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "VCluster/VCluster.hpp"
#include "Decomposition/CartDecomposition.hpp"

#define SUB_UNIT_FACTOR 1024

BOOST_AUTO_TEST_SUITE( decomposition_to_gpu_test )


BOOST_AUTO_TEST_CASE( decomposition_to_gpu_test_use )
{
	auto & v_cl = create_vcluster();

	// Vcluster
	Vcluster & vcl = create_vcluster();

	CartDecomposition<3, float, CudaMemory, memory_traits_inte> dec(vcl);

	// Physical domain
	Box<3, float> box( { 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
	size_t div[3];

	// Get the number of processor and calculate the number of sub-domain
	// for each processor (SUB_UNIT_FACTOR=64)
	size_t n_proc = vcl.getProcessingUnits();
	size_t n_sub = n_proc * SUB_UNIT_FACTOR;

	// Set the number of sub-domains on each dimension (in a scalable way)
	for (int i = 0; i < 3; i++)
	{	div[i] = openfpm::math::round_big_2(pow(n_sub,1.0/3));}

	// Define ghost
	Ghost<3, float> g(0.01);

	// Boundary conditions
	size_t bc[] = { PERIODIC, PERIODIC, PERIODIC };

	// Decompose
	dec.setParameters(div,box,bc,g);
	dec.decompose();

	openfpm::vector_gpu<Point<3,float>> vg;
	vg.resize(10000);

	for (size_t i = 0 ; i < 10000 ; i++)
	{
		vg.template get<0>(i)[0] = (float)rand()/RAND_MAX;
		vg.template get<0>(i)[1] = (float)rand()/RAND_MAX;
		vg.template get<0>(i)[2] = (float)rand()/RAND_MAX;
	}

	vg.hostToDevice<0>();

	// process on GPU the processor ID for each particles

	auto ite = vg.getGPUIterator();

	openfpm::vector_gpu<aggregate<int,int>> proc_id_out;
	proc_id_out.resize(vg.size());

	process_id_proc_each_part<decltype(dec.toKernel()),decltype(vg.toKernel()),decltype(proc_id_out.toKernel())>
	<<<ite.wthr,ite.thr>>>
	(dec.toKernel(),vg.toKernel(),proc_id_out.toKernel(),v_cl.rank());

	proc_id_out.deviceToHost<0>();

	bool match = true;
	for (size_t i = 0 ; i < proc_id_out.size() ; i++)
	{
		Point<3,float> xp = vg.template get<0>(i);

		match &= proc_id_out.template get<0>(i) == dec.processorIDBC(xp);
	}
}

BOOST_AUTO_TEST_SUITE_END()
