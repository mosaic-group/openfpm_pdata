#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "VCluster/VCluster.hpp"
#include "Vector/map_vector.hpp"
#include "Vector/cuda/vector_dist_cuda_funcs.cuh"
#include "Vector/util/vector_dist_funcs.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "util/cuda/scan_cuda.cuh"
#include "util/cuda/moderngpu/kernel_scan.hxx"

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
		v_pos.template get<0>(i)[0] = (float)rand()/RAND_MAX;
		v_pos.template get<0>(i)[1] = (float)rand()/RAND_MAX;
		v_pos.template get<0>(i)[2] = (float)rand()/RAND_MAX;

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
	num_shift_ghost_each_part<3,float,decltype(box_f_dev.toKernel()),decltype(v_pos.toKernel()),decltype(o_part_loc.toKernel())>
	<<<ite.wthr,ite.thr>>>
	(box_f_dev.toKernel(),v_pos.toKernel(),o_part_loc.toKernel());

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
	mgpu::scan((unsigned int *)o_part_loc.template getDeviceBuffer<0>(), o_part_loc.size(), (unsigned int *)starts.template getDeviceBuffer<0>() , v_cl.getmgpuContext());

	starts.deviceToHost<0>(starts.size()-1,starts.size()-1);
	size_t tot = starts.template get<0>(o_part_loc.size()-1);

	openfpm::vector<Point<3,float>,CudaMemory,typename memory_traits_inte<Point<3,float>>::type,memory_traits_inte> shifts;

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

	shift_ghost_each_part<3,float,decltype(box_f_dev.toKernel()),decltype(box_f_sv.toKernel()),
			                     decltype(v_pos.toKernel()),decltype(v_prp.toKernel()),
			                     decltype(starts.toKernel()),decltype(shifts.toKernel()),
			                     decltype(o_part_loc2.toKernel())>
	<<<ite.wthr,ite.thr>>>
	(box_f_dev.toKernel(),box_f_sv.toKernel(),
	 v_pos.toKernel(),v_prp.toKernel(),
	 starts.toKernel(),shifts.toKernel(),o_part_loc2.toKernel(),old);

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

	process_ghost_particles_local<true,3,decltype(o_part_loc2.toKernel()),decltype(v_pos2.toKernel()),decltype(v_prp2.toKernel()),decltype(shifts.toKernel())>
	<<<ite.wthr,ite.thr>>>
	(o_part_loc2.toKernel(),v_pos2.toKernel(),v_prp2.toKernel(),shifts.toKernel(),old);

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
	typedef openfpm::vector<prp_object,CudaMemory,typename memory_traits_inte<prp_object>::type,memory_traits_inte> send_vector;

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

		process_ghost_particles_prp<decltype(g_opart_device.toKernel()),decltype(g_send_prp.get(i).toKernel()),decltype(v_prp.toKernel()),0,1,2>
		<<<ite.wthr,ite.thr>>>
		(g_opart_device.toKernel(), g_send_prp.get(i).toKernel(),
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

		const openfpm::vector<std::pair<size_t, size_t>> & vp_id = dec.template ghost_processorID_pair<typename dec_type::lc_processor_id, typename dec_type::shift_id>(xp, UNIQUE);

		match &= proc_id_out.template get<0>(i) == vp_id.size();
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

	///////////////////////// we collect the processor and shift id //////////////////////////

    openfpm::vector<aggregate<unsigned int,long unsigned int>,
                    CudaMemory,
                    typename memory_traits_inte<aggregate<unsigned int,long unsigned int>>::type,
                    memory_traits_inte> output;

    output.resize(sz);

	ite = vg.getGPUIterator();

	// we compute processor id for each particle
	proc_label_id_ghost<3,float,decltype(dec.toKernel()),decltype(vg.toKernel()),decltype(starts.toKernel()),decltype(output.toKernel())>
	<<<ite.wthr,ite.thr>>>
	(dec.toKernel(),vg.toKernel(),starts.toKernel(),output.toKernel());

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

	BOOST_REQUIRE_EQUAL(match,true);
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

	find_buffer_offsets<1,decltype(vgp.toKernel()),decltype(offs.toKernel())><<<ite.wthr,ite.thr>>>(vgp.toKernel(),(int *)mem.getDevicePointer(),offs.toKernel());

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

