#define BOOST_TEST_DYN_LINK

#define TEST1

#include <boost/test/unit_test.hpp>
#include "VCluster/VCluster.hpp"
#include "Vector/map_vector.hpp"
#include "Vector/cuda/vector_dist_cuda_funcs.cuh"
#include "Vector/util/vector_dist_funcs.hpp"
#include "Decomposition/CartDecomposition.hpp"
//#include "util/cuda/scan_cuda.cuh"
#include "Vector/vector_dist.hpp"
#include "util/cuda/scan_ofp.cuh"

#define SUB_UNIT_FACTOR 1024

BOOST_AUTO_TEST_SUITE( vector_dist_gpu_util_func_test )

BOOST_AUTO_TEST_CASE( vector_ghost_process_local_particles )
{
	typedef aggregate<float,float[3],float[3][3]> prop;

	openfpm::vector_gpu<prop> v_prp;
	v_prp.resize(10000);

	openfpm::vector_gpu<Point<3,float>> v_pos;
	v_pos.resize(10000);

	openfpm::vector_gpu<aggregate<unsigned int,unsigned int>> o_part_loc;

	for (size_t i = 0 ; i < v_prp.size() ; i++)
	{
		v_pos.template get<0>(i)[0] = (float)rand()/(float)RAND_MAX;
		v_pos.template get<0>(i)[1] = (float)rand()/(float)RAND_MAX;
		v_pos.template get<0>(i)[2] = (float)rand()/(float)RAND_MAX;

		v_prp.template get<0>(i) = i+12345;

		v_prp.template get<1>(i)[0] = i;
		v_prp.template get<1>(i)[1] = i+20000;
		v_prp.template get<1>(i)[2] = i+50000;

		v_prp.template get<2>(i)[0][0] = i+60000;
		v_prp.template get<2>(i)[0][1] = i+70000;
		v_prp.template get<2>(i)[0][2] = i+80000;
		v_prp.template get<2>(i)[1][0] = i+90000;
		v_prp.template get<2>(i)[1][1] = i+100000;
		v_prp.template get<2>(i)[1][2] = i+110000;
		v_prp.template get<2>(i)[2][0] = i+120000;
		v_prp.template get<2>(i)[2][1] = i+130000;
		v_prp.template get<2>(i)[2][1] = i+140000;
		v_prp.template get<2>(i)[2][2] = i+150000;
	}

	openfpm::vector_gpu<Box<3,float>> box_f_dev;
	openfpm::vector_gpu<aggregate<unsigned int>> box_f_sv;

	box_f_dev.resize(4);
	box_f_sv.resize(4);

	box_f_dev.template get<0>(0)[0] = 0.0;
	box_f_dev.template get<0>(0)[1] = 0.0;
	box_f_dev.template get<0>(0)[2] = 0.0;
	box_f_dev.template get<1>(0)[0] = 0.5;
	box_f_dev.template get<1>(0)[1] = 1.0;
	box_f_dev.template get<1>(0)[2] = 1.0;
	box_f_sv.template get<0>(0) = 0;

	box_f_dev.template get<0>(1)[0] = 0.0;
	box_f_dev.template get<0>(1)[1] = 0.0;
	box_f_dev.template get<0>(1)[2] = 0.0;
	box_f_dev.template get<1>(1)[0] = 0.3;
	box_f_dev.template get<1>(1)[1] = 1.0;
	box_f_dev.template get<1>(1)[2] = 1.0;
	box_f_sv.template get<0>(1) = 1;

	box_f_dev.template get<0>(2)[0] = 0.0;
	box_f_dev.template get<0>(2)[1] = 0.0;
	box_f_dev.template get<0>(2)[2] = 0.0;
	box_f_dev.template get<1>(2)[0] = 0.2;
	box_f_dev.template get<1>(2)[1] = 1.0;
	box_f_dev.template get<1>(2)[2] = 1.0;
	box_f_sv.template get<0>(2) = 2;

	box_f_dev.template get<0>(3)[0] = 0.0;
	box_f_dev.template get<0>(3)[1] = 0.0;
	box_f_dev.template get<0>(3)[2] = 0.0;
	box_f_dev.template get<1>(3)[0] = 0.1;
	box_f_dev.template get<1>(3)[1] = 1.0;
	box_f_dev.template get<1>(3)[2] = 1.0;
	box_f_sv.template get<0>(3) = 3;

	// Label the internal (assigned) particles
	auto ite = v_pos.getGPUIteratorTo(v_pos.size());

	o_part_loc.resize(v_pos.size()+1);
	o_part_loc.template get<0>(o_part_loc.size()-1) = 0;
	o_part_loc.template hostToDevice<0>(o_part_loc.size()-1,o_part_loc.size()-1);

	box_f_dev.hostToDevice<0,1>();
	box_f_sv.hostToDevice<0>();
	v_pos.hostToDevice<0>();
	v_prp.hostToDevice<0,1,2>();

	// label particle processor
	CUDA_LAUNCH_DIM3((num_shift_ghost_each_part<3,float,decltype(box_f_dev.toKernel()),decltype(box_f_sv.toKernel()),decltype(v_pos.toKernel()),decltype(o_part_loc.toKernel())>),
	ite.wthr,ite.thr,
	box_f_dev.toKernel(),box_f_sv.toKernel(),v_pos.toKernel(),o_part_loc.toKernel(),v_pos.size());

	o_part_loc.deviceToHost<0>();

	bool match = true;

	for (size_t i = 0 ; i < v_pos.size() ; i++)
	{
		if (v_pos.template get<0>(i)[0] >= 0.5)
		{match &= o_part_loc.template get<0>(i) == 0;}
		else if (v_pos.template get<0>(i)[0] >= 0.3)
		{match &= o_part_loc.template get<0>(i) == 1;}
		else if (v_pos.template get<0>(i)[0] >= 0.2)
		{match &= o_part_loc.template get<0>(i) == 2;}
		else if (v_pos.template get<0>(i)[0] >= 0.1)
		{match &= o_part_loc.template get<0>(i) == 3;}
		else
		{match &= o_part_loc.template get<0>(i) == 4;}
	}

	BOOST_REQUIRE_EQUAL(match,true);

	openfpm::vector_gpu<aggregate<unsigned int>> starts;
	starts.resize(o_part_loc.size());

	auto & v_cl = create_vcluster();
	openfpm::scan((unsigned int *)o_part_loc.template getDeviceBuffer<0>(), o_part_loc.size(), (unsigned int *)starts.template getDeviceBuffer<0>() , v_cl.getgpuContext());

	starts.deviceToHost<0>(starts.size()-1,starts.size()-1);
	size_t tot = starts.template get<0>(o_part_loc.size()-1);

	openfpm::vector<Point<3,float>,CudaMemory,memory_traits_inte> shifts;

	shifts.resize(4);

	shifts.template get<0>(0)[0] = 10.0;
	shifts.template get<0>(0)[1] = 0.0;
	shifts.template get<0>(0)[2] = 0.0;

	shifts.template get<0>(1)[0] = 20.0;
	shifts.template get<0>(1)[1] = 0.0;
	shifts.template get<0>(1)[2] = 0.0;

	shifts.template get<0>(2)[0] = 30.0;
	shifts.template get<0>(2)[1] = 0.0;
	shifts.template get<0>(2)[2] = 0.0;

	shifts.template get<0>(3)[0] = 40.0;
	shifts.template get<0>(3)[1] = 0.0;
	shifts.template get<0>(3)[2] = 0.0;

	size_t old = v_pos.size();
	v_pos.resize(v_pos.size() + tot);
	v_prp.resize(v_prp.size() + tot);

	shifts.template hostToDevice<0>();
	openfpm::vector_gpu<aggregate<unsigned int,unsigned int>> o_part_loc2;
	o_part_loc2.resize(tot);

	CUDA_LAUNCH_DIM3((shift_ghost_each_part<3,float,decltype(box_f_dev.toKernel()),decltype(box_f_sv.toKernel()),
			                     decltype(v_pos.toKernel()),decltype(v_prp.toKernel()),
			                     decltype(starts.toKernel()),decltype(shifts.toKernel()),
			                     decltype(o_part_loc2.toKernel())>),
	ite.wthr,ite.thr,
	box_f_dev.toKernel(),box_f_sv.toKernel(),
	 v_pos.toKernel(),v_prp.toKernel(),
	 starts.toKernel(),shifts.toKernel(),o_part_loc2.toKernel(),old,old);

	v_pos.deviceToHost<0>();
	o_part_loc2.deviceToHost<0,1>();
	v_prp.deviceToHost<0,1,2>();

	size_t base = old;
	size_t base_o = 0;
	for (size_t i = 0 ; i < old ; i++)
	{
		if (v_pos.template get<0>(i)[0] >= 0.5)
		{}
		else if (v_pos.template get<0>(i)[0] >= 0.3)
		{
			for (size_t j = 0 ; j < o_part_loc.template get<0>(i) ; j++)
			{
				match &= v_pos.template get<0>(base)[0] < 1.0 - (j+1.0)*10.0;
				match &= v_pos.template get<0>(base)[0] >= -(j+1.0)*10.0;

				match &= o_part_loc2.template get<0>(base_o) == i;
				match &= o_part_loc2.template get<1>(base_o) == j;

				////// We check the properties

				match &= v_prp.template get<0>(base) == v_prp.template get<0>(i);

				match &= v_prp.template get<1>(base)[0] == v_prp.template get<1>(i)[0];
				match &= v_prp.template get<1>(base)[1] == v_prp.template get<1>(i)[1];
				match &= v_prp.template get<1>(base)[2] == v_prp.template get<1>(i)[2];

				match &= v_prp.template get<2>(base)[0][0] == v_prp.template get<2>(i)[0][0];
				match &= v_prp.template get<2>(base)[0][1] == v_prp.template get<2>(i)[0][1];
				match &= v_prp.template get<2>(base)[0][2] == v_prp.template get<2>(i)[0][2];
				match &= v_prp.template get<2>(base)[1][0] == v_prp.template get<2>(i)[1][0];
				match &= v_prp.template get<2>(base)[1][1] == v_prp.template get<2>(i)[1][1];
				match &= v_prp.template get<2>(base)[1][2] == v_prp.template get<2>(i)[1][2];
				match &= v_prp.template get<2>(base)[2][0] == v_prp.template get<2>(i)[2][0];
				match &= v_prp.template get<2>(base)[2][1] == v_prp.template get<2>(i)[2][1];
				match &= v_prp.template get<2>(base)[2][2] == v_prp.template get<2>(i)[2][2];

				base++;
				base_o++;
			}
		}
		else if (v_pos.template get<0>(i)[0] >= 0.2)
		{
			for (size_t j = 0 ; j < o_part_loc.template get<0>(i) ; j++)
			{
				match &= v_pos.template get<0>(base)[0] < 1.0 - (j+1.0)*10.0;
				match &= v_pos.template get<0>(base)[0] >= -(j+1.0)*10.0;

				match &= o_part_loc2.template get<0>(base_o) == i;
				match &= o_part_loc2.template get<1>(base_o) == j;

				////// We check the properties

				match &= v_prp.template get<0>(base) == v_prp.template get<0>(i);

				match &= v_prp.template get<1>(base)[0] == v_prp.template get<1>(i)[0];
				match &= v_prp.template get<1>(base)[1] == v_prp.template get<1>(i)[1];
				match &= v_prp.template get<1>(base)[2] == v_prp.template get<1>(i)[2];

				match &= v_prp.template get<2>(base)[0][0] == v_prp.template get<2>(i)[0][0];
				match &= v_prp.template get<2>(base)[0][1] == v_prp.template get<2>(i)[0][1];
				match &= v_prp.template get<2>(base)[0][2] == v_prp.template get<2>(i)[0][2];
				match &= v_prp.template get<2>(base)[1][0] == v_prp.template get<2>(i)[1][0];
				match &= v_prp.template get<2>(base)[1][1] == v_prp.template get<2>(i)[1][1];
				match &= v_prp.template get<2>(base)[1][2] == v_prp.template get<2>(i)[1][2];
				match &= v_prp.template get<2>(base)[2][0] == v_prp.template get<2>(i)[2][0];
				match &= v_prp.template get<2>(base)[2][1] == v_prp.template get<2>(i)[2][1];
				match &= v_prp.template get<2>(base)[2][2] == v_prp.template get<2>(i)[2][2];


				base++;
				base_o++;
			}
		}
		else if (v_pos.template get<0>(i)[0] >= 0.1)
		{
			for (size_t j = 0 ; j < o_part_loc.template get<0>(i) ; j++)
			{
				match &= v_pos.template get<0>(base)[0] < 1.0 - (j+1.0)*10.0;
				match &= v_pos.template get<0>(base)[0] >= -(j+1.0)*10.0;

				match &= o_part_loc2.template get<0>(base_o) == i;
				match &= o_part_loc2.template get<1>(base_o) == j;

				////// We check the properties

				match &= v_prp.template get<0>(base) == v_prp.template get<0>(i);

				match &= v_prp.template get<1>(base)[0] == v_prp.template get<1>(i)[0];
				match &= v_prp.template get<1>(base)[1] == v_prp.template get<1>(i)[1];
				match &= v_prp.template get<1>(base)[2] == v_prp.template get<1>(i)[2];

				match &= v_prp.template get<2>(base)[0][0] == v_prp.template get<2>(i)[0][0];
				match &= v_prp.template get<2>(base)[0][1] == v_prp.template get<2>(i)[0][1];
				match &= v_prp.template get<2>(base)[0][2] == v_prp.template get<2>(i)[0][2];
				match &= v_prp.template get<2>(base)[1][0] == v_prp.template get<2>(i)[1][0];
				match &= v_prp.template get<2>(base)[1][1] == v_prp.template get<2>(i)[1][1];
				match &= v_prp.template get<2>(base)[1][2] == v_prp.template get<2>(i)[1][2];
				match &= v_prp.template get<2>(base)[2][0] == v_prp.template get<2>(i)[2][0];
				match &= v_prp.template get<2>(base)[2][1] == v_prp.template get<2>(i)[2][1];
				match &= v_prp.template get<2>(base)[2][2] == v_prp.template get<2>(i)[2][2];

				base++;
				base_o++;
			}
		}
		else
		{
			for (size_t j = 0 ; j < o_part_loc.template get<0>(i) ; j++)
			{
				match &= v_pos.template get<0>(base)[0] < 1.0 - (j+1.0)*10.0;
				match &= v_pos.template get<0>(base)[0] >= -(j+1.0)*10.0;

				match &= o_part_loc2.template get<0>(base_o) == i;
				match &= o_part_loc2.template get<1>(base_o) == j;

				////// We check the properties

				match &= v_prp.template get<0>(base) == v_prp.template get<0>(i);

				match &= v_prp.template get<1>(base)[0] == v_prp.template get<1>(i)[0];
				match &= v_prp.template get<1>(base)[1] == v_prp.template get<1>(i)[1];
				match &= v_prp.template get<1>(base)[2] == v_prp.template get<1>(i)[2];

				match &= v_prp.template get<2>(base)[0][0] == v_prp.template get<2>(i)[0][0];
				match &= v_prp.template get<2>(base)[0][1] == v_prp.template get<2>(i)[0][1];
				match &= v_prp.template get<2>(base)[0][2] == v_prp.template get<2>(i)[0][2];
				match &= v_prp.template get<2>(base)[1][0] == v_prp.template get<2>(i)[1][0];
				match &= v_prp.template get<2>(base)[1][1] == v_prp.template get<2>(i)[1][1];
				match &= v_prp.template get<2>(base)[1][2] == v_prp.template get<2>(i)[1][2];
				match &= v_prp.template get<2>(base)[2][0] == v_prp.template get<2>(i)[2][0];
				match &= v_prp.template get<2>(base)[2][1] == v_prp.template get<2>(i)[2][1];
				match &= v_prp.template get<2>(base)[2][2] == v_prp.template get<2>(i)[2][2];

				base++;
				base_o++;
			}
		}
	}

	BOOST_REQUIRE_EQUAL(match,true);

	////////////// Now we check that o_part_loc2 contain processble information

	openfpm::vector_gpu<Point<3,float>> v_pos2;
	openfpm::vector_gpu<prop> v_prp2;

	v_pos2.resize(old);
	v_prp2.resize(old);

	for (size_t i = 0 ; i < old ; i++)
	{
		v_pos2.template get<0>(i)[0] = v_pos.template get<0>(i)[0];
		v_pos2.template get<0>(i)[1] = v_pos.template get<0>(i)[1];
		v_pos2.template get<0>(i)[2] = v_pos.template get<0>(i)[2];

		v_prp2.template get<0>(i) = v_prp.template get<0>(i);

		v_prp2.template get<1>(i)[0] = v_prp.template get<1>(i)[0];
		v_prp2.template get<1>(i)[1] = v_prp.template get<1>(i)[1];
		v_prp2.template get<1>(i)[2] = v_prp.template get<1>(i)[2];

		v_prp2.template get<2>(i)[0][0] = v_prp.template get<2>(i)[0][0];
		v_prp2.template get<2>(i)[0][1] = v_prp.template get<2>(i)[0][1];
		v_prp2.template get<2>(i)[0][2] = v_prp.template get<2>(i)[0][2];
		v_prp2.template get<2>(i)[1][0] = v_prp.template get<2>(i)[1][0];
		v_prp2.template get<2>(i)[1][1] = v_prp.template get<2>(i)[1][1];
		v_prp2.template get<2>(i)[1][2] = v_prp.template get<2>(i)[1][2];
		v_prp2.template get<2>(i)[2][0] = v_prp.template get<2>(i)[2][0];
		v_prp2.template get<2>(i)[2][1] = v_prp.template get<2>(i)[2][1];
		v_prp2.template get<2>(i)[2][2] = v_prp.template get<2>(i)[2][2];
	}

	v_pos2.resize(v_pos.size());
	v_prp2.resize(v_prp.size());

	v_pos2.hostToDevice<0>();
	v_prp2.hostToDevice<0,1,2>();

	ite = o_part_loc2.getGPUIterator();

	CUDA_LAUNCH_DIM3((process_ghost_particles_local<true,3,decltype(o_part_loc2.toKernel()),decltype(v_pos2.toKernel()),decltype(v_prp2.toKernel()),decltype(shifts.toKernel())>),
	ite.wthr,ite.thr,
	o_part_loc2.toKernel(),v_pos2.toKernel(),v_prp2.toKernel(),shifts.toKernel(),old);

	v_pos2.template deviceToHost<0>();
	v_prp2.template deviceToHost<0,1,2>();

	for (size_t i = old ; i < v_pos.size() ; i++)
	{
		match &= v_pos.template get<0>(i)[0] == v_pos2.template get<0>(i)[0];
		match &= v_pos.template get<0>(i)[1] == v_pos2.template get<0>(i)[1];
		match &= v_pos.template get<0>(i)[2] == v_pos2.template get<0>(i)[2];

		match &= v_prp2.template get<0>(i) == v_prp.template get<0>(i);

		match &= v_prp2.template get<1>(i)[0] == v_prp.template get<1>(i)[0];
		match &= v_prp2.template get<1>(i)[1] == v_prp.template get<1>(i)[1];
		match &= v_prp2.template get<1>(i)[2] == v_prp.template get<1>(i)[2];

		match &= v_prp2.template get<2>(i)[0][0] == v_prp.template get<2>(i)[0][0];
		match &= v_prp2.template get<2>(i)[0][1] == v_prp.template get<2>(i)[0][1];
		match &= v_prp2.template get<2>(i)[0][2] == v_prp.template get<2>(i)[0][2];
		match &= v_prp2.template get<2>(i)[1][0] == v_prp.template get<2>(i)[1][0];
		match &= v_prp2.template get<2>(i)[1][1] == v_prp.template get<2>(i)[1][1];
		match &= v_prp2.template get<2>(i)[1][2] == v_prp.template get<2>(i)[1][2];
		match &= v_prp2.template get<2>(i)[2][0] == v_prp.template get<2>(i)[2][0];
		match &= v_prp2.template get<2>(i)[2][1] == v_prp.template get<2>(i)[2][1];
		match &= v_prp2.template get<2>(i)[2][2] == v_prp.template get<2>(i)[2][2];
	}

	BOOST_REQUIRE_EQUAL(match,true);
}

BOOST_AUTO_TEST_CASE( vector_ghost_fill_send_buffer_test )
{
	typedef aggregate<float,float[3],float[3][3]> prop;

	// Sending property object
	typedef object<typename object_creator<typename prop::type, 0,1,2>::type> prp_object;

	// send vector for each processor
	typedef openfpm::vector<prp_object,CudaMemory,memory_traits_inte> send_vector;

	openfpm::vector<send_vector> g_send_prp;

	auto & v_cl = create_vcluster();

	// Vcluster
	Vcluster<> & vcl = create_vcluster();

	openfpm::vector_gpu<prop> v_prp;
	v_prp.resize(10000);

	openfpm::vector_gpu<aggregate<unsigned int,unsigned int,unsigned int>> g_opart_device;

	for (size_t i = 0 ; i < v_prp.size() ; i++)
	{
		v_prp.template get<0>(i) = i+12345;

		v_prp.template get<1>(i)[0] = i;
		v_prp.template get<1>(i)[1] = i+20000;
		v_prp.template get<1>(i)[2] = i+50000;

		v_prp.template get<2>(i)[0][0] = i+60000;
		v_prp.template get<2>(i)[0][1] = i+70000;
		v_prp.template get<2>(i)[0][2] = i+80000;
		v_prp.template get<2>(i)[1][0] = i+90000;
		v_prp.template get<2>(i)[1][1] = i+100000;
		v_prp.template get<2>(i)[1][2] = i+110000;
		v_prp.template get<2>(i)[2][0] = i+120000;
		v_prp.template get<2>(i)[2][1] = i+130000;
		v_prp.template get<2>(i)[2][2] = i+140000;
	}

	v_prp.hostToDevice<0,1,2>();

	g_opart_device.resize(2*10000*3);

	for (size_t i = 0 ; i < 3 ; i++)
	{
		for (size_t j = 0 ; j < 10000 ; j++)
		{
			g_opart_device.template get<0>(i*2*10000 + j*2) = i;
			g_opart_device.template get<0>(i*2*10000 + j*2+1) = i;

			g_opart_device.template get<1>(i*2*10000 + j*2) = j;
			g_opart_device.template get<1>(i*2*10000 + j*2+1) = j;

			g_opart_device.template get<2>(i*2*10000 + j*2) = 0;
			g_opart_device.template get<2>(i*2*10000 + j*2+1) = 0;
		}
	}

	g_opart_device.hostToDevice<0,1,2>();

	g_send_prp.resize(3);

	bool match = true;
	size_t offset = 0;

	for (size_t i = 0 ; i < 3 ; i++)
	{
		g_send_prp.get(i).resize(2*10000);

		auto ite = g_send_prp.get(i).getGPUIterator();

		CUDA_LAUNCH_DIM3((process_ghost_particles_prp<decltype(g_opart_device.toKernel()),decltype(g_send_prp.get(i).toKernel()),decltype(v_prp.toKernel()),0,1,2>),
		ite.wthr,ite.thr,
		g_opart_device.toKernel(), g_send_prp.get(i).toKernel(),
		 v_prp.toKernel(),offset);

		offset += g_send_prp.get(i).size();

		///////////// TEST ////////////

		g_send_prp.get(i).deviceToHost<0,1,2>();

		for (size_t j = 0 ; j < 10000 ; j++)
		{
			match &= g_send_prp.get(i).template get<0>(2*j) == j+12345;

			match &= g_send_prp.get(i).template get<1>(2*j)[0] == j;
			match &= g_send_prp.get(i).template get<1>(2*j)[1] == j+20000;
			match &= g_send_prp.get(i).template get<1>(2*j)[2] == j+50000;

			match &= g_send_prp.get(i).template get<2>(2*j)[0][0] == j+60000;
			match &= g_send_prp.get(i).template get<2>(2*j)[0][1] == j+70000;
			match &= g_send_prp.get(i).template get<2>(2*j)[0][2] == j+80000;
			match &= g_send_prp.get(i).template get<2>(2*j)[1][0] == j+90000;
			match &= g_send_prp.get(i).template get<2>(2*j)[1][1] == j+100000;
			match &= g_send_prp.get(i).template get<2>(2*j)[1][2] == j+110000;
			match &= g_send_prp.get(i).template get<2>(2*j)[2][0] == j+120000;
			match &= g_send_prp.get(i).template get<2>(2*j)[2][1] == j+130000;
			match &= g_send_prp.get(i).template get<2>(2*j)[2][2] == j+140000;


			match = g_send_prp.get(i).template get<0>(2*j+1) == j+12345;

			match = g_send_prp.get(i).template get<1>(2*j+1)[0] == j;
			match = g_send_prp.get(i).template get<1>(2*j+1)[1] == j+20000;
			match = g_send_prp.get(i).template get<1>(2*j+1)[2] == j+50000;

			match = g_send_prp.get(i).template get<2>(2*j+1)[0][0] == j+60000;
			match = g_send_prp.get(i).template get<2>(2*j+1)[0][1] == j+70000;
			match = g_send_prp.get(i).template get<2>(2*j+1)[0][2] == j+80000;
			match = g_send_prp.get(i).template get<2>(2*j+1)[1][0] == j+90000;
			match = g_send_prp.get(i).template get<2>(2*j+1)[1][1] == j+100000;
			match = g_send_prp.get(i).template get<2>(2*j+1)[1][2] == j+110000;
			match = g_send_prp.get(i).template get<2>(2*j+1)[2][0] == j+120000;
			match = g_send_prp.get(i).template get<2>(2*j+1)[2][1] == j+130000;
			match = g_send_prp.get(i).template get<2>(2*j+1)[2][2] == j+140000;
		}
	}

	BOOST_REQUIRE_EQUAL(match,true);
}

BOOST_AUTO_TEST_CASE( decomposition_ie_ghost_gpu_test_use )
{
	auto & v_cl = create_vcluster();

	// Vcluster
	Vcluster<> & vcl = create_vcluster();

	typedef CartDecomposition<3, float, CudaMemory, memory_traits_inte> dec_type;

	dec_type dec(vcl);

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
			vg.template get<0>(k*n_part+j)[0] = (sp.getHigh(0) - sp.getLow(0))*((float)rand()/(float)RAND_MAX) + sp.getLow(0);
			vg.template get<0>(k*n_part+j)[1] = (sp.getHigh(1) - sp.getLow(1))*((float)rand()/(float)RAND_MAX) + sp.getLow(1);
			vg.template get<0>(k*n_part+j)[2] = (sp.getHigh(2) - sp.getLow(2))*((float)rand()/(float)RAND_MAX) + sp.getLow(2);
		}
	}

	vg.hostToDevice<0>();

	// process on GPU the processor ID for each particles

	auto ite = vg.getGPUIterator();

	openfpm::vector_gpu<aggregate<unsigned int>> proc_id_out;
	proc_id_out.resize(vg.size()+1);
	proc_id_out.template get<0>(proc_id_out.size()-1) = 0;
	proc_id_out.template hostToDevice(proc_id_out.size()-1,proc_id_out.size()-1);

	CUDA_LAUNCH_DIM3((num_proc_ghost_each_part<3,float,decltype(dec.toKernel()),decltype(vg.toKernel()),decltype(proc_id_out.toKernel())>),
	ite.wthr,ite.thr,
	dec.toKernel(),vg.toKernel(),proc_id_out.toKernel());

/*	proc_id_out.deviceToHost<0>();

	bool match = true;
	for (size_t i = 0 ; i < vg.size() ; i++)
	{
		Point<3,float> xp = vg.template get<0>(i);

		match &= proc_id_out.template get<0>(i) == dec.ghost_processorID_N(xp);

		const openfpm::vector<std::pair<size_t, size_t>> & vp_id = dec.template ghost_processorID_pair<typename dec_type::lc_processor_id, typename dec_type::shift_id>(xp, UNIQUE);

		match &= proc_id_out.template get<0>(i) == vp_id.size();
	}


	BOOST_REQUIRE_EQUAL(match,true);

	////////////////////////// We now create the scan //////////////////////////////////////

    openfpm::vector<aggregate<unsigned int>,
                    CudaMemory,
                    memory_traits_inte> starts;

    starts.resize(proc_id_out.size());

    openfpm::scan((unsigned int *)proc_id_out.template getDeviceBuffer<0>(),proc_id_out.size(),(unsigned int *)starts.template getDeviceBuffer<0>(),v_cl.getgpuContext());

	starts.deviceToHost<0>(starts.size()-1,starts.size()-1);

	size_t sz = starts.template get<0>(starts.size()-1);

	///////////////////////// we collect the processor and shift id //////////////////////////

    openfpm::vector<aggregate<unsigned int,long unsigned int>,
                    CudaMemory,
                    memory_traits_inte> output;

    output.resize(sz);

	ite = vg.getGPUIterator();

	// we compute processor id for each particle
	CUDA_LAUNCH_DIM3((proc_label_id_ghost<3,float,decltype(dec.toKernel()),decltype(vg.toKernel()),decltype(starts.toKernel()),decltype(output.toKernel())>),
	ite.wthr,ite.thr,
	dec.toKernel(),vg.toKernel(),starts.toKernel(),output.toKernel());

	output.template deviceToHost<0,1>();

	//////////////////// TESTING //////////////////////////

	starts.deviceToHost<0>();

	match = true;

	for (size_t i = 0 ; i < starts.size() - 1 ; i++)
	{
		size_t base = starts.template get<0>(i);
		size_t sz = starts.template get<0>(i+1) - base;

		if (sz != 0)
		{
			size_t pid = output.template get<1>(base) & 0xFFFFFFFF;
			Point<3,float> xp = vg.template get<0>(pid);

			openfpm::vector<proc_box_id> tmp_sort1;
			openfpm::vector<proc_box_id> tmp_sort2;

			const openfpm::vector<std::pair<size_t, size_t>> & vp_id = dec.template ghost_processorID_pair<typename dec_type::lc_processor_id, typename dec_type::shift_id>(xp, UNIQUE);

			tmp_sort1.resize(vp_id.size());

			for (size_t j = 0 ; j < vp_id.size() ; j++)
			{
				tmp_sort1.get(j).proc_id = dec.IDtoProc(vp_id.get(j).first);
				tmp_sort1.get(j).box_id = 0;
				tmp_sort1.get(j).shift_id = vp_id.get(j).second;
			}

			tmp_sort1.sort();

			tmp_sort2.resize(sz);
			for (size_t j = 0 ; j < sz ; j++)
			{
				tmp_sort2.get(j).proc_id = output.template get<0>(base+j);
				tmp_sort2.get(j).box_id = 0;
				tmp_sort2.get(j).shift_id = output.template get<1>(base+j) >> 32;
			}

			tmp_sort2.sort();

			match &= tmp_sort1.size() == tmp_sort2.size();

			for (size_t j = 0 ; j < tmp_sort1.size() ; j++)
			{
				match &= tmp_sort1.get(j).proc_id == tmp_sort2.get(j).proc_id;
				match &= tmp_sort1.get(j).shift_id == tmp_sort2.get(j).shift_id;
			}
		}
	}

	BOOST_REQUIRE_EQUAL(match,true);*/
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
		vg.template get<0>(i)[0] = (float)rand()/(float)RAND_MAX;
		vg.template get<0>(i)[1] = (float)rand()/(float)RAND_MAX;
		vg.template get<0>(i)[2] = (float)rand()/(float)RAND_MAX;
	}

	vg.hostToDevice<0>();

	// process on GPU the processor ID for each particles

	auto ite = vg.getGPUIterator();

	openfpm::vector_gpu<aggregate<int,int,int>> proc_id_out;
	proc_id_out.resize(vg.size());

	openfpm::vector_gpu<aggregate<int,int,int>> dev_counter;
	dev_counter.resize(v_cl.size());
	dev_counter.fill<0>(0);
	dev_counter.fill<1>(0);
	dev_counter.fill<2>(0);

	CUDA_LAUNCH_DIM3((process_id_proc_each_part<3,float,decltype(dec.toKernel()),decltype(vg.toKernel()),decltype(proc_id_out.toKernel()),decltype(dev_counter.toKernel())>),
	ite.wthr,ite.thr,
	dec.toKernel(),vg.toKernel(),proc_id_out.toKernel(),dev_counter.toKernel(),v_cl.rank());


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
	vgp.hostToDevice<0,1>();

	CUDA_LAUNCH((find_buffer_offsets<1,decltype(vgp.toKernel()),decltype(offs.toKernel())>),ite,vgp.toKernel(),(int *)mem.getDevicePointer(),offs.toKernel());

	offs.template deviceToHost<0,1>();

	mem.deviceToHost();
	int n_ele = *(int *)mem.getPointer();
	BOOST_REQUIRE_EQUAL(n_ele,199);

	openfpm::vector<int> ofv;
	openfpm::vector<int> ofv2;

	for (size_t i = 0 ; i < n_ele ; i++)
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

BOOST_AUTO_TEST_CASE(vector_dist_reorder_lbl)
{
	openfpm::vector_gpu<aggregate<int,int,int>> lbl_p;
	openfpm::vector_gpu<aggregate<int>> starts;

	lbl_p.resize(100);
	starts.resize(10);

	for (int i = 0 ; i < 10 ; i++) // <------ particle id
	{
		for (int j = 0 ; j < 10 ; j++) // <----- processor
		{
			lbl_p.template get<2>(i*10+j) = i;
			lbl_p.template get<1>(i*10+j) = j;
		}
		starts.template get<0>(i) = (i*10);
	}

	// move lbl and starts to gpu
	starts.template hostToDevice<0>();
	lbl_p.template hostToDevice<1,2>();

	auto ite = lbl_p.getGPUIterator();

	CUDA_LAUNCH_DIM3((reorder_lbl<decltype(lbl_p.toKernel()),decltype(starts.toKernel())>),ite.wthr,ite.thr,lbl_p.toKernel(),starts.toKernel());

	starts.template deviceToHost<0>();
	lbl_p.template deviceToHost<0,1,2>();

	for (int i = 0 ; i < 10 ; i++) // <------ particle id
	{
		for (int j = 0 ; j < 10 ; j++) // <----- processor
		{
			BOOST_REQUIRE_EQUAL(lbl_p.template get<0>(j*10+i),i*10+j);
		}
	}
}

BOOST_AUTO_TEST_CASE(vector_dist_merge_sort)
{
	openfpm::vector_gpu<aggregate<float[3],float[3],float[3]>> v_prp;
	openfpm::vector_gpu<Point<3,float>> v_pos;

	openfpm::vector_gpu<aggregate<float[3],float[3],float[3]>> v_prp_out;
	openfpm::vector_gpu<Point<3,float>> v_pos_out;

	openfpm::vector_gpu<aggregate<int>> ns_to_s;

	v_prp.resize(10000);
	v_pos.resize(10000);
	v_prp_out.resize(10000);
	v_pos_out.resize(10000);
	ns_to_s.resize(10000);

	for (int i = 0 ; i < 10000 ; i++) // <------ particle id
	{
		v_pos_out.template get<0>(i)[0] = i;
		v_pos_out.template get<0>(i)[1] = i+10000;
		v_pos_out.template get<0>(i)[2] = i+20000;

		v_pos.template get<0>(i)[0] = 0;
		v_pos.template get<0>(i)[1] = 0;
		v_pos.template get<0>(i)[2] = 0;

		v_prp_out.template get<0>(i)[0] = i+60123;
		v_prp_out.template get<0>(i)[1] = i+73543;
		v_prp_out.template get<0>(i)[2] = i+82432;

		v_prp_out.template get<1>(i)[0] = i+80123;
		v_prp_out.template get<1>(i)[1] = i+93543;
		v_prp_out.template get<1>(i)[2] = i+102432;

		v_prp_out.template get<2>(i)[0] = i+110123;
		v_prp_out.template get<2>(i)[1] = i+123543;
		v_prp_out.template get<2>(i)[2] = i+132432;

		v_prp.template get<0>(i)[0] = 0;
		v_prp.template get<0>(i)[1] = 0;
		v_prp.template get<0>(i)[2] = 0;

		v_prp.template get<1>(i)[0] = 0;
		v_prp.template get<1>(i)[1] = 0;
		v_prp.template get<1>(i)[2] = 0;

		v_prp.template get<2>(i)[0] = 0;
		v_prp.template get<2>(i)[1] = 0;
		v_prp.template get<2>(i)[2] = 0;

		ns_to_s.template get<0>(i) = 10000-i-1;
	}

	v_prp.template hostToDevice<0,1,2>();
	v_prp_out.template hostToDevice<0,1,2>();
	v_pos.template hostToDevice<0>();
	v_pos_out.template hostToDevice<0>();
	ns_to_s.template hostToDevice<0>();

	auto ite = v_pos.getGPUIterator();

	CUDA_LAUNCH_DIM3((merge_sort_part<false,decltype(v_pos.toKernel()),decltype(v_prp.toKernel()),decltype(ns_to_s.toKernel()),0>),ite.wthr,ite.thr,v_pos.toKernel(),v_prp.toKernel(),
																								 v_pos_out.toKernel(),v_prp_out.toKernel(),
																								 ns_to_s.toKernel());

	v_prp.template deviceToHost<0,1,2>();

	bool match = true;
	for (int i = 0 ; i < 10000 ; i++) // <------ particle id
	{
		match &= v_prp_out.template get<0>(10000-i-1)[0] == v_prp.template get<0>(i)[0];
		match &= v_prp_out.template get<0>(10000-i-1)[1] == v_prp.template get<0>(i)[1];
		match &= v_prp_out.template get<0>(10000-i-1)[2] == v_prp.template get<0>(i)[2];

		match &= v_prp.template get<1>(10000-i-1)[0] == 0;
		match &= v_prp.template get<1>(10000-i-1)[1] == 0;
		match &= v_prp.template get<1>(10000-i-1)[2] == 0;

		match &= v_prp.template get<2>(10000-i-1)[0] == 0;
		match &= v_prp.template get<2>(10000-i-1)[1] == 0;
		match &= v_prp.template get<2>(10000-i-1)[2] == 0;
	}

	BOOST_REQUIRE_EQUAL(match,true);

	CUDA_LAUNCH_DIM3((merge_sort_part<false,decltype(v_pos.toKernel()),decltype(v_prp.toKernel()),decltype(ns_to_s.toKernel()),1,2>),ite.wthr,ite.thr,v_pos.toKernel(),v_prp.toKernel(),
																								 v_pos_out.toKernel(),v_prp_out.toKernel(),
																								 ns_to_s.toKernel());

	v_prp.template deviceToHost<0,1,2>();
	v_pos.template deviceToHost<0>();

	for (int i = 0 ; i < 10000 ; i++) // <------ particle id
	{
		match &= v_prp_out.template get<0>(10000-i-1)[0] == v_prp.template get<0>(i)[0];
		match &= v_prp_out.template get<0>(10000-i-1)[1] == v_prp.template get<0>(i)[1];
		match &= v_prp_out.template get<0>(10000-i-1)[2] == v_prp.template get<0>(i)[2];

		match &= v_prp_out.template get<1>(10000-i-1)[0] == v_prp.template get<1>(i)[0];
		match &= v_prp_out.template get<1>(10000-i-1)[1] == v_prp.template get<1>(i)[1];
		match &= v_prp_out.template get<1>(10000-i-1)[2] == v_prp.template get<1>(i)[2];

		match &= v_prp_out.template get<2>(10000-i-1)[0] == v_prp.template get<2>(i)[0];
		match &= v_prp_out.template get<2>(10000-i-1)[1] == v_prp.template get<2>(i)[1];
		match &= v_prp_out.template get<2>(10000-i-1)[2] == v_prp.template get<2>(i)[2];

		match &= v_pos.template get<0>(10000-i-1)[0] == 0;
		match &= v_pos.template get<0>(10000-i-1)[1] == 0;
		match &= v_pos.template get<0>(10000-i-1)[2] == 0;
	}

	BOOST_REQUIRE_EQUAL(match,true);

	CUDA_LAUNCH_DIM3((merge_sort_part<true,decltype(v_pos.toKernel()),decltype(v_prp.toKernel()),decltype(ns_to_s.toKernel())>),ite.wthr,ite.thr,v_pos.toKernel(),v_prp.toKernel(),
																								 v_pos_out.toKernel(),v_prp_out.toKernel(),
																								 ns_to_s.toKernel());

	v_prp.template deviceToHost<0,1,2>();
	v_pos.template deviceToHost<0>();

	for (int i = 0 ; i < 10000 ; i++) // <------ particle id
	{


		match &= v_prp_out.template get<0>(10000-i-1)[0] == v_prp.template get<0>(i)[0];
		match &= v_prp_out.template get<0>(10000-i-1)[1] == v_prp.template get<0>(i)[1];
		match &= v_prp_out.template get<0>(10000-i-1)[2] == v_prp.template get<0>(i)[2];

		match &= v_prp_out.template get<1>(10000-i-1)[0] == v_prp.template get<1>(i)[0];
		match &= v_prp_out.template get<1>(10000-i-1)[1] == v_prp.template get<1>(i)[1];
		match &= v_prp_out.template get<1>(10000-i-1)[2] == v_prp.template get<1>(i)[2];

		match &= v_prp_out.template get<2>(10000-i-1)[0] == v_prp.template get<2>(i)[0];
		match &= v_prp_out.template get<2>(10000-i-1)[1] == v_prp.template get<2>(i)[1];
		match &= v_prp_out.template get<2>(10000-i-1)[2] == v_prp.template get<2>(i)[2];


		match &= v_pos_out.template get<0>(10000-i-1)[0] == v_pos.template get<0>(i)[0];
		match &= v_pos_out.template get<0>(10000-i-1)[1] == v_pos.template get<0>(i)[1];
		match &= v_pos_out.template get<0>(10000-i-1)[2] == v_pos.template get<0>(i)[2];
	}

	BOOST_REQUIRE_EQUAL(match,true);
}

BOOST_AUTO_TEST_CASE(vector_dist_gpu_map_fill_send_buffer_test)
{
	openfpm::vector_gpu<aggregate<int,int>> m_opart;

    openfpm::vector<openfpm::vector<Point<3,float>,CudaMemory,memory_traits_inte,openfpm::grow_policy_identity>> m_pos;
    openfpm::vector<openfpm::vector<aggregate<float,float[2],float[3][3]>,CudaMemory,memory_traits_inte,openfpm::grow_policy_identity>> m_prp;

    openfpm::vector_gpu<Point<3,float>> v_pos;
    openfpm::vector_gpu<aggregate<float,float[2],float[3][3]>> v_prp;

    unsigned int offset = 0;

    v_pos.resize(100000);
    v_prp.resize(v_pos.size());
    m_opart.resize(v_pos.size());

    for (size_t i = 0 ; i < v_pos.size() ; i++)
    {
    	v_pos.template get<0>(i)[0] = (float)rand()/(float)RAND_MAX;
    	v_pos.template get<0>(i)[1] = (float)rand()/(float)RAND_MAX;
    	v_pos.template get<0>(i)[2] = (float)rand()/(float)RAND_MAX;

    	v_prp.template get<0>(i) = 5.0 + (float)rand()/(float)RAND_MAX;
    	v_prp.template get<1>(i)[0] = 10.0 + (float)rand()/(float)RAND_MAX;
    	v_prp.template get<1>(i)[1] = 11.0 + (float)rand()/(float)RAND_MAX;
    	v_prp.template get<2>(i)[0][0] = 40.0 + (float)rand()/(float)RAND_MAX;
    	v_prp.template get<2>(i)[0][1] = 50.0 + (float)rand()/(float)RAND_MAX;
    	v_prp.template get<2>(i)[0][2] = 60.0 + (float)rand()/(float)RAND_MAX;
    	v_prp.template get<2>(i)[1][0] = 70.0 + (float)rand()/(float)RAND_MAX;
    	v_prp.template get<2>(i)[1][1] = 80.0 + (float)rand()/(float)RAND_MAX;
    	v_prp.template get<2>(i)[1][2] = 150.0 + (float)rand()/(float)RAND_MAX;
    	v_prp.template get<2>(i)[2][0] = 160.0 + (float)rand()/(float)RAND_MAX;
    	v_prp.template get<2>(i)[2][1] = 170.0 + (float)rand()/(float)RAND_MAX;
    	v_prp.template get<2>(i)[2][2] = 340.0 + (float)rand()/(float)RAND_MAX;

    	int seg = i / 10000;
    	m_opart.template get<1>(i) = seg;
    	m_opart.template get<0>(i) = (9999 - i%10000) + seg * 10000;
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

		CUDA_LAUNCH_DIM3((process_map_particles<decltype(m_opart.toKernel()),decltype(m_pos.get(i).toKernel()),decltype(m_prp.get(i).toKernel()),
																		   decltype(v_pos.toKernel()),decltype(v_prp.toKernel())>),
						ite.wthr,ite.thr,
						m_opart.toKernel(),m_pos.get(i).toKernel(), m_prp.get(i).toKernel(),
											v_pos.toKernel(),v_prp.toKernel(),offset);

		m_pos.get(i).deviceToHost<0>();
		m_prp.get(i).deviceToHost<0,1,2>();

		bool match = true;

		for (size_t j = 0 ; j < m_pos.get(i).size() ; j++)
		{
			match &= (m_pos.get(i).template get<0>(j)[0] == v_pos.template get<0>(m_opart.template get<0>(offset+j))[0]);
			match &= (m_pos.get(i).template get<0>(j)[1] == v_pos.template get<0>(m_opart.template get<0>(offset+j))[1]);
			match &= (m_pos.get(i).template get<0>(j)[2] == v_pos.template get<0>(m_opart.template get<0>(offset+j))[2]);

			match &= (m_prp.get(i).template get<0>(j) == v_prp.template get<0>(m_opart.template get<0>(offset+j)));

			match &= (m_prp.get(i).template get<1>(j)[0] == v_prp.template get<1>(m_opart.template get<0>(offset+j))[0]);
			match &= (m_prp.get(i).template get<1>(j)[1] == v_prp.template get<1>(m_opart.template get<0>(offset+j))[1]);

			match &= (m_prp.get(i).template get<2>(j)[0][0] == v_prp.template get<2>(m_opart.template get<0>(offset+j))[0][0]);
			match &= (m_prp.get(i).template get<2>(j)[0][1] == v_prp.template get<2>(m_opart.template get<0>(offset+j))[0][1]);
			match &= (m_prp.get(i).template get<2>(j)[0][2] == v_prp.template get<2>(m_opart.template get<0>(offset+j))[0][2]);
			match &= (m_prp.get(i).template get<2>(j)[1][0] == v_prp.template get<2>(m_opart.template get<0>(offset+j))[1][0]);
			match &= (m_prp.get(i).template get<2>(j)[1][1] == v_prp.template get<2>(m_opart.template get<0>(offset+j))[1][1]);
			match &= (m_prp.get(i).template get<2>(j)[1][2] == v_prp.template get<2>(m_opart.template get<0>(offset+j))[1][2]);
			match &= (m_prp.get(i).template get<2>(j)[2][0] == v_prp.template get<2>(m_opart.template get<0>(offset+j))[2][0]);
			match &= (m_prp.get(i).template get<2>(j)[2][1] == v_prp.template get<2>(m_opart.template get<0>(offset+j))[2][1]);
			match &= (m_prp.get(i).template get<2>(j)[2][2] == v_prp.template get<2>(m_opart.template get<0>(offset+j))[2][2]);
		}

		BOOST_REQUIRE_EQUAL(match,true);

		offset += m_pos.get(i).size();
    }
}

template<unsigned int prp>
void vector_dist_remove_marked_type()
{
	auto & v_cl = create_vcluster();

	if (v_cl.size() > 16)
	{return;}

	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	// set the ghost based on the radius cut off (make just a little bit smaller than the spacing)
	Ghost<3,float> g(0.1);

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	vector_dist_gpu<3,float,aggregate<float,float,int,int>> vd(50000*v_cl.size(),domain,bc,g);

	// Fill the position

	auto it = vd.getDomainIterator();

	while(it.isNext())
	{
		auto p = it.get();

		vd.getPos(p)[0] = (float)rand() / (float)RAND_MAX;
		vd.getPos(p)[1] = (float)rand() / (float)RAND_MAX;
		vd.getPos(p)[2] = (float)rand() / (float)RAND_MAX;

		++it;
	}

	vd.map();
	vd.template ghost_get<>();

	it = vd.getDomainIterator();

	float fc = 1.0;
	float dc = 1.0;
	int ic = 1;
	int sc = 1;

	while(it.isNext())
	{
		auto p = it.get();

		vd.template getProp<0>(p) = fc;
		vd.template getProp<1>(p) = dc;
		vd.template getProp<2>(p) = ic;
		vd.template getProp<3>(p) = sc;

		vd.template getProp<prp>(p) = (ic % 3 == 0);

		fc += 1.0;
		dc += 1.0;
		ic += 1;
		sc += 1;

		++it;
	}

	size_t sz = vd.size_local() - vd.size_local()/3;

	vd.template hostToDeviceProp<0,1,2,3>();

	remove_marked<prp>(vd);

	BOOST_REQUIRE_EQUAL(vd.size_local(),sz);

	vd.template deviceToHostProp<0,1,2,3>();

	auto it2 = vd.getDomainIterator();

	// There should not be number divisible by 3

	bool test = true;

	while(it2.isNext())
	{
		auto p = it2.get();

		if (prp != 0)
		{test &= ((int)vd.template getProp<0>(p) % 3 != 0);}

		if (prp != 1)
		{test &= ((int)vd.template getProp<1>(p) % 3 != 0);}

		if (prp != 2)
		{test &= ((int)vd.template getProp<2>(p) % 3 != 0);}

		if (prp != 3)
		{test &= ((int)vd.template getProp<3>(p) % 3 != 0);}

		if (test == false)
		{
			if (prp != 0)
			{std::cout << (int)vd.template getProp<0>(p) << std::endl;}

			if (prp != 1)
			{std::cout << (int)vd.template getProp<1>(p) << std::endl;}

			if (prp != 2)
			{std::cout << (int)vd.template getProp<2>(p) << std::endl;}

			if (prp != 3)
			{std::cout << (int)vd.template getProp<3>(p) << std::endl;}

			break;
		}

		++it2;
	}

	BOOST_REQUIRE_EQUAL(test,true);


	// We test where we do not remove anything

	size_t size_old = vd.size_local();

	// Because remove_marked is destructive we have to reset the property
	vd.getPropVector().template fill<prp>(0);

	remove_marked<prp>(vd);

	BOOST_REQUIRE_EQUAL(vd.size_local(),size_old);

	// Now we try to remove all
	vd.getPropVector().template fill<prp>(1);

	remove_marked<prp>(vd);

	BOOST_REQUIRE_EQUAL(vd.size_local(),0);
}

BOOST_AUTO_TEST_CASE(vector_dist_remove_marked)
{
	vector_dist_remove_marked_type<0>();
	vector_dist_remove_marked_type<1>();
	vector_dist_remove_marked_type<2>();
	vector_dist_remove_marked_type<3>();
}


BOOST_AUTO_TEST_CASE( vector_dist_particle_NN_MP_iteration_gpu )
{
	typedef  aggregate<size_t,size_t,size_t> part_prop;

	Vcluster<> & v_cl = create_vcluster();

	if (v_cl.getProcessingUnits() > 24)
	{return;}

	float L = 1000.0;

    // set the seed
	// create the random generator engine
    std::default_random_engine eg;
    eg.seed(v_cl.rank()*4533);
    std::uniform_real_distribution<float> ud(-L,L);

    long int k = 4096 * v_cl.getProcessingUnits();

	long int big_step = k / 4;
	big_step = (big_step == 0)?1:big_step;

	BOOST_TEST_CHECKPOINT( "Testing 3D periodic vector symmetric cell-list k=" << k );

	Box<3,float> box({-L,-L,-L},{L,L,L});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	float r_cut = 100.0;

	// ghost
	Ghost<3,float> ghost(r_cut);

	// Distributed vector
	vector_dist_gpu<3,float,part_prop> vd(k,box,bc,ghost,BIND_DEC_TO_GHOST);

/*	size_t start = vd.init_size_accum(k);

	auto it = vd.getIterator();

	while (it.isNext())
	{
		auto key = it.get();

		vd.getPosWrite(key)[0] = ud(eg);
		vd.getPosWrite(key)[1] = ud(eg);
		vd.getPosWrite(key)[2] = ud(eg);

		// Fill some properties randomly

		vd.template getPropWrite<0>(key) = 0;
		vd.template getPropWrite<1>(key) = 0;
		vd.template getPropWrite<2>(key) = key.getKey() + start;

		++it;
	}

	vd.map();

	// sync the ghost
	vd.template ghost_get<0,2>();

	auto NN = vd.getCellList(r_cut);
	auto p_it = vd.getDomainIterator();

	while (p_it.isNext())
	{
		auto p = p_it.get();

		Point<3,float> xp = vd.getPosRead(p);

		auto Np = NN.getNNIterator(NN.getCell(xp));

		while (Np.isNext())
		{
			auto q = Np.get();

			if (p.getKey() == q)
			{
				++Np;
				continue;
			}

			// repulsive

			Point<3,float> xq = vd.getPosRead(q);
			Point<3,float> f = (xp - xq);

			float distance = f.norm();

			// Particle should be inside 2 * r_cut range

			if (distance < r_cut )
			{
				vd.template getPropWrite<0>(p)++;
			}

			++Np;
		}

		++p_it;
	}

	// We now divide the particles on 4 phases

	openfpm::vector<vector_dist_gpu<3,float,part_prop>> phases;
	phases.add( vector_dist_gpu<3,float,part_prop>(vd.getDecomposition(),0));
	phases.add( vector_dist_gpu<3,float,part_prop>(phases.get(0).getDecomposition(),0));
	phases.add( vector_dist_gpu<3,float,part_prop>(phases.get(0).getDecomposition(),0));
	phases.add( vector_dist_gpu<3,float,part_prop>(phases.get(0).getDecomposition(),0));

	auto it2 = vd.getDomainIterator();

	while (it2.isNext())
	{
		auto p = it2.get();

		if (p.getKey() % 4 == 0)
		{
			phases.get(0).add();
			phases.get(0).getLastPos()[0] = vd.getPos(p)[0];
			phases.get(0).getLastPos()[1] = vd.getPos(p)[1];
			phases.get(0).getLastPos()[2] = vd.getPos(p)[2];

			phases.get(0).getLastProp<1>() = 0;

			phases.get(0).template getLastProp<2>() = vd.template getProp<2>(p);
		}
		else if (p.getKey() % 4 == 1)
		{
			phases.get(1).add();
			phases.get(1).getLastPos()[0] = vd.getPos(p)[0];
			phases.get(1).getLastPos()[1] = vd.getPos(p)[1];
			phases.get(1).getLastPos()[2] = vd.getPos(p)[2];

			phases.get(1).getLastProp<1>() = 0;

			phases.get(1).template getLastProp<2>() = vd.template getProp<2>(p);
		}
		else if (p.getKey() % 4 == 2)
		{
			phases.get(2).add();
			phases.get(2).getLastPos()[0] = vd.getPos(p)[0];
			phases.get(2).getLastPos()[1] = vd.getPos(p)[1];
			phases.get(2).getLastPos()[2] = vd.getPos(p)[2];

			phases.get(2).getLastProp<1>() = 0;

			phases.get(2).template getLastProp<2>() = vd.template getProp<2>(p);
		}
		else
		{
			phases.get(3).add();
			phases.get(3).getLastPos()[0] = vd.getPos(p)[0];
			phases.get(3).getLastPos()[1] = vd.getPos(p)[1];
			phases.get(3).getLastPos()[2] = vd.getPos(p)[2];

			phases.get(3).getLastProp<1>() = 0;

			phases.get(3).template getLastProp<2>() = vd.template getProp<2>(p);
		}

		++it2;
	}

	// now we synchronize the ghosts

	for (size_t i = 0 ; i < phases.size() ; i++)
	{
		phases.get(i).template ghost_get<0,1,2>();
	}

	typedef decltype(phases.get(0).getCellListSym(r_cut)) cell_list_type;

	openfpm::vector<cell_list_type> NN_ptr;

	for (size_t i = 0 ; i < phases.size() ; i++)
	{
		NN_ptr.add(phases.get(i).getCellListSym(r_cut));
	}

	// We now interact all the phases

	for (size_t i = 0 ; i < phases.size() ; i++)
	{
		for (size_t j = 0 ; j < phases.size() ; j++)
		{
			auto p_it2 = phases.get(i).getDomainIterator();

			while (p_it2.isNext())
			{
				auto p = p_it2.get();

				Point<3,float> xp = phases.get(i).getPosRead(p);

				auto Np = NN_ptr.get(j).getNNIteratorSymMP<NO_CHECK>(NN_ptr.get(j).getCell(xp),p.getKey(),phases.get(i).getPosVector(),phases.get(j).getPosVector());

				while (Np.isNext())
				{
					auto q = Np.get();

					if (p.getKey() == q && i == j)
					{
						++Np;
						continue;
					}

					// repulsive

					Point<3,float> xq = phases.get(j).getPosRead(q);
					Point<3,float> f = (xp - xq);

					float distance = f.norm();

					// Particle should be inside r_cut range

					if (distance < r_cut )
					{
						phases.get(i).template getPropWrite<1>(p)++;
						phases.get(j).template getPropWrite<1>(q)++;
					}

					++Np;
				}

				++p_it2;
			}
		}
	}

	for (size_t i = 0 ; i < phases.size() ; i++)
	{
		phases.get(i).template ghost_put<add_,1>();
	}

	auto p_it3 = vd.getDomainIterator();

	bool ret = true;
	while (p_it3.isNext())
	{
		auto p = p_it3.get();

		int ph;

		if (p.getKey() % 4 == 0)
		{ph = 0;}
		else if (p.getKey() % 4 == 1)
		{ph = 1;}
		else if (p.getKey() % 4 == 2)
		{ph = 2;}
		else
		{ph = 3;}

		size_t pah = p.getKey()/4;
		ret &= phases.get(ph).template getPropRead<1>(pah) == vd.template getPropRead<0>(p);

		if (ret == false)
		{
			std::cout << "Error on particle: " << vd.template getPropRead<2>(p) << "   " << v_cl.rank() << std::endl;

			std::cout << "phase " << ph << " particle " << pah << "   " <<  phases.get(ph).template getPropRead<1>(pah) << "  A  " << vd.template getPropRead<0>(p) << std::endl;

			break;
		}

		++p_it3;
	}

	BOOST_REQUIRE_EQUAL(ret,true);*/
}

BOOST_AUTO_TEST_SUITE_END()

