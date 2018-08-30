#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "VCluster/VCluster.hpp"
#include "Vector/map_vector.hpp"
#include "Vector/cuda/vector_dist_cuda_funcs.cuh"
#include "Vector/util/vector_dist_funcs.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "util/cuda/scan_cuda.cuh"

#define SUB_UNIT_FACTOR 1024

BOOST_AUTO_TEST_SUITE( vector_dist_gpu_util_func_test )

BOOST_AUTO_TEST_CASE( decomposition_ie_ghost_gpu_test_use )
{
	auto & v_cl = create_vcluster();

	// Vcluster
	Vcluster<> & vcl = create_vcluster();

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
	Ghost<3, float> g(0.1);

	// Boundary conditions
	size_t bc[] = { PERIODIC, PERIODIC, PERIODIC };

	// Decompose
	dec.setParameters(div,box,bc,g);
	dec.decompose();

	// Get the local boxes

	int nsub = dec.getNSubDomain();
	int n_part = 10000 / nsub;

	openfpm::vector_gpu<Point<3,float>> vg;
	vg.resize(nsub*n_part);

	for (size_t k = 0 ; k < nsub ; k++)
	{
		SpaceBox<3,float> sp = dec.getSubDomain(k);

		for (size_t j = 0 ; j < n_part ; j++)
		{
			vg.template get<0>(k*n_part+j)[0] = (sp.getHigh(0) - sp.getLow(0))*((float)rand()/RAND_MAX) + sp.getLow(0);
			vg.template get<0>(k*n_part+j)[1] = (sp.getHigh(1) - sp.getLow(1))*((float)rand()/RAND_MAX) + sp.getLow(1);
			vg.template get<0>(k*n_part+j)[2] = (sp.getHigh(2) - sp.getLow(2))*((float)rand()/RAND_MAX) + sp.getLow(2);
		}
	}

	vg.hostToDevice<0>();

	// process on GPU the processor ID for each particles

	auto ite = vg.getGPUIterator();

	openfpm::vector_gpu<aggregate<unsigned int>> proc_id_out;
	proc_id_out.resize(vg.size()+1);
	proc_id_out.template get<0>(proc_id_out.size()-1) = 0;
	proc_id_out.template hostToDevice(proc_id_out.size()-1,proc_id_out.size()-1);

	num_proc_ghost_each_part<3,float,decltype(dec.toKernel()),decltype(vg.toKernel()),decltype(proc_id_out.toKernel())>
	<<<ite.wthr,ite.thr>>>
	(dec.toKernel(),vg.toKernel(),proc_id_out.toKernel());

	proc_id_out.deviceToHost<0>();

	bool match = true;
	for (size_t i = 0 ; i < vg.size() ; i++)
	{
		Point<3,float> xp = vg.template get<0>(i);

		match &= proc_id_out.template get<0>(i) == dec.ghost_processorID_N(xp);
	}

	BOOST_REQUIRE_EQUAL(match,true);

	////////////////////////// We now create the scan //////////////////////////////////////

    openfpm::vector<aggregate<unsigned int>,
                    CudaMemory,
                    typename memory_traits_inte<aggregate<unsigned int>>::type,
                    memory_traits_inte> starts;

    starts.resize(proc_id_out.size());

	// scan
	scan<unsigned int,unsigned int>(proc_id_out,starts);
	starts.deviceToHost<0>(starts.size()-1,starts.size()-1);

	size_t sz = starts.template get<0>(starts.size()-1);

	///////////////////////// We now test //////////////////////////

    openfpm::vector<aggregate<unsigned int,unsigned int>,
                    CudaMemory,
                    typename memory_traits_inte<aggregate<unsigned int,unsigned int>>::type,
                    memory_traits_inte> output;

    output.resize(sz);

	ite = vg.getGPUIterator();

	// we compute processor id for each particle
	proc_label_id_ghost<3,float,decltype(dec.toKernel()),decltype(vg.toKernel()),decltype(starts.toKernel()),decltype(output.toKernel())>
	<<<ite.wthr,ite.thr>>>
	(dec.toKernel(),vg.toKernel(),starts.toKernel(),output.toKernel());

	output.template deviceToHost<0,1>();

	for (size_t i = 0 ; i < output.size() ; i++)
	{
		std::cout << output.template get<0>(i) << "   " << output.template get<1>(i) << std::endl;
	}
}

BOOST_AUTO_TEST_CASE( decomposition_to_gpu_test_use )
{
	auto & v_cl = create_vcluster();

	// Vcluster
	Vcluster<> & vcl = create_vcluster();

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


BOOST_AUTO_TEST_CASE( vector_dist_gpu_find_buffer_offsets_test )
{
	openfpm::vector_gpu<aggregate<int,int>> vgp;
	openfpm::vector_gpu<aggregate<int,int>> offs;

	vgp.resize(200000);

	for (size_t k = 0 ; k < vgp.size() ; k++)
	{
		vgp.template get<0>(k) = k / 1000;
		vgp.template get<1>(k) = k / 1000;
	}

	offs.resize(220);

	CudaMemory mem;
	mem.allocate(sizeof(int));
	mem.fill(0);

	auto ite = vgp.getGPUIterator();
	vgp.hostToDevice<0>();

	find_buffer_offsets<decltype(vgp.toKernel()),decltype(offs.toKernel())><<<ite.wthr,ite.thr>>>(vgp.toKernel(),(int *)mem.getDevicePointer(),offs.toKernel());

	offs.template deviceToHost<0,1>();

	openfpm::vector<int> ofv;
	openfpm::vector<int> ofv2;

	for (size_t i = 0 ; i < ofv.size() ; i++)
	{
		ofv.add(offs.template get<0>(i));
		ofv2.add(offs.template get<1>(i));
	}

	ofv.sort();
	ofv2.sort();

	for (size_t i = 0 ; i < ofv.size() ; i++)
	{
		BOOST_REQUIRE_EQUAL(ofv.get(i),(i+1)*1000);
		BOOST_REQUIRE_EQUAL(ofv2.get(i),i);
	}
}

BOOST_AUTO_TEST_CASE(vector_dist_gpu_map_fill_send_buffer_test)
{
	openfpm::vector_gpu<aggregate<int,int>> m_opart;

    openfpm::vector<openfpm::vector<Point<3,float>,CudaMemory,typename memory_traits_inte<Point<3,float>>::type,memory_traits_inte,openfpm::grow_policy_identity>> m_pos;
    openfpm::vector<openfpm::vector<aggregate<float,float[2],float[3][3]>,CudaMemory,typename memory_traits_inte<aggregate<float,float[2],float[3][3]>>::type,memory_traits_inte,openfpm::grow_policy_identity>> m_prp;

    openfpm::vector_gpu<Point<3,float>> v_pos;
    openfpm::vector_gpu<aggregate<float,float[2],float[3][3]>> v_prp;

    unsigned int offset = 0;

    v_pos.resize(100000);
    v_prp.resize(v_pos.size());
    m_opart.resize(v_pos.size());

    for (size_t i = 0 ; i < v_pos.size() ; i++)
    {
    	v_pos.template get<0>(i)[0] = (float)rand()/RAND_MAX;
    	v_pos.template get<0>(i)[1] = (float)rand()/RAND_MAX;
    	v_pos.template get<0>(i)[2] = (float)rand()/RAND_MAX;

    	v_prp.template get<0>(i) = 5.0 + (float)rand()/RAND_MAX;
    	v_prp.template get<1>(i)[0] = 10.0 + (float)rand()/RAND_MAX;
    	v_prp.template get<1>(i)[1] = 11.0 + (float)rand()/RAND_MAX;
    	v_prp.template get<2>(i)[0][0] = 40.0 + (float)rand()/RAND_MAX;
    	v_prp.template get<2>(i)[0][1] = 50.0 + (float)rand()/RAND_MAX;
    	v_prp.template get<2>(i)[0][2] = 60.0 + (float)rand()/RAND_MAX;
    	v_prp.template get<2>(i)[1][0] = 70.0 + (float)rand()/RAND_MAX;
    	v_prp.template get<2>(i)[1][1] = 80.0 + (float)rand()/RAND_MAX;
    	v_prp.template get<2>(i)[1][2] = 150.0 + (float)rand()/RAND_MAX;
    	v_prp.template get<2>(i)[2][0] = 160.0 + (float)rand()/RAND_MAX;
    	v_prp.template get<2>(i)[2][1] = 170.0 + (float)rand()/RAND_MAX;
    	v_prp.template get<2>(i)[2][2] = 340.0 + (float)rand()/RAND_MAX;

    	int seg = i / 10000;
    	m_opart.template get<0>(i) = seg;
    	m_opart.template get<1>(i) = (9999 - i%10000) + seg * 10000;
    }

    m_pos.resize(10);
    m_prp.resize(10);

    for (size_t i = 0 ; i < m_pos.size() ; i++)
    {
    	m_pos.get(i).resize(10000);
    	m_prp.get(i).resize(10000);
    }

    v_pos.hostToDevice<0>();
    v_prp.hostToDevice<0,1,2>();

    m_opart.hostToDevice<0,1>();

    for (size_t i = 0 ; i < m_pos.size() ; i++)
    {
    	auto ite = m_pos.get(i).getGPUIterator();

		process_map_particles<decltype(m_opart.toKernel()),decltype(m_pos.get(i).toKernel()),decltype(m_prp.get(i).toKernel()),
																		   decltype(v_pos.toKernel()),decltype(v_prp.toKernel())>
						<<<ite.wthr,ite.thr>>>
						(m_opart.toKernel(),m_pos.get(i).toKernel(), m_prp.get(i).toKernel(),
											v_pos.toKernel(),v_prp.toKernel(),offset);

		m_pos.get(i).deviceToHost<0>();
		m_prp.get(i).deviceToHost<0,1,2>();

		bool match = true;

		for (size_t j = 0 ; j < m_pos.get(i).size() ; j++)
		{
			match &= (m_pos.get(i).template get<0>(j)[0] == v_pos.template get<0>(m_opart.template get<1>(offset+j))[0]);
			match &= (m_pos.get(i).template get<0>(j)[1] == v_pos.template get<0>(m_opart.template get<1>(offset+j))[1]);
			match &= (m_pos.get(i).template get<0>(j)[2] == v_pos.template get<0>(m_opart.template get<1>(offset+j))[2]);

			match &= (m_prp.get(i).template get<0>(j) == v_prp.template get<0>(m_opart.template get<1>(offset+j)));

			match &= (m_prp.get(i).template get<1>(j)[0] == v_prp.template get<1>(m_opart.template get<1>(offset+j))[0]);
			match &= (m_prp.get(i).template get<1>(j)[1] == v_prp.template get<1>(m_opart.template get<1>(offset+j))[1]);

			match &= (m_prp.get(i).template get<2>(j)[0][0] == v_prp.template get<2>(m_opart.template get<1>(offset+j))[0][0]);
			match &= (m_prp.get(i).template get<2>(j)[0][1] == v_prp.template get<2>(m_opart.template get<1>(offset+j))[0][1]);
			match &= (m_prp.get(i).template get<2>(j)[0][2] == v_prp.template get<2>(m_opart.template get<1>(offset+j))[0][2]);
			match &= (m_prp.get(i).template get<2>(j)[1][0] == v_prp.template get<2>(m_opart.template get<1>(offset+j))[1][0]);
			match &= (m_prp.get(i).template get<2>(j)[1][1] == v_prp.template get<2>(m_opart.template get<1>(offset+j))[1][1]);
			match &= (m_prp.get(i).template get<2>(j)[1][2] == v_prp.template get<2>(m_opart.template get<1>(offset+j))[1][2]);
			match &= (m_prp.get(i).template get<2>(j)[2][0] == v_prp.template get<2>(m_opart.template get<1>(offset+j))[2][0]);
			match &= (m_prp.get(i).template get<2>(j)[2][1] == v_prp.template get<2>(m_opart.template get<1>(offset+j))[2][1]);
			match &= (m_prp.get(i).template get<2>(j)[2][2] == v_prp.template get<2>(m_opart.template get<1>(offset+j))[2][2]);
		}

		BOOST_REQUIRE_EQUAL(match,true);

		offset += m_pos.get(i).size();
    }
}

BOOST_AUTO_TEST_SUITE_END()

