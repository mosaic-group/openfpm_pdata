#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "VCluster/VCluster.hpp"
#include "Vector/map_vector.hpp"
#include "Vector/cuda/vector_dist_cuda_funcs.cuh"
#include "Vector/util/vector_dist_funcs.hpp"

BOOST_AUTO_TEST_SUITE( vector_dist_gpu_util_func_test )

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

