#define BOOST_TEST_DYN_LINK
#include "config.h"
#include <boost/test/unit_test.hpp>
#include "VCluster/VCluster.hpp"
#include <Vector/vector_dist.hpp>
#include "Vector/tests/vector_dist_util_unit_tests.hpp"

#define SUB_UNIT_FACTOR 1024

template<unsigned int dim , typename vector_dist_type>
__global__ void move_parts_gpu_test(vector_dist_type vd)
{
	auto p = GET_PARTICLE(vd);

#pragma unroll
	for (int i = 0 ; i < dim ; i++)
	{
		vd.getPos(p)[i] += 0.05;
	}
}

BOOST_AUTO_TEST_SUITE( vector_dist_gpu_test )

void print_test(std::string test, size_t sz)
{
	if (create_vcluster().getProcessUnitID() == 0)
		std::cout << test << " " << sz << "\n";
}


__global__  void initialize_props(vector_dist_ker<3, float, aggregate<float, float [3], float[3]>> vd)
{
	auto p = GET_PARTICLE(vd);

	vd.template getProp<0>(p) = vd.getPos(p)[0] + vd.getPos(p)[1] + vd.getPos(p)[2];

	vd.template getProp<1>(p)[0] = vd.getPos(p)[0] + vd.getPos(p)[1];
	vd.template getProp<1>(p)[1] = vd.getPos(p)[0] + vd.getPos(p)[2];
	vd.template getProp<1>(p)[2] = vd.getPos(p)[1] + vd.getPos(p)[2];
}

template<typename T,typename CellList_type>
__global__  void calculate_force(vector_dist_ker<3, T, aggregate<T, T[3], T [3]>> vd,
		                         vector_dist_ker<3, T, aggregate<T, T[3], T [3]>> vd_sort,
		                         CellList_type cl,
		                         int rank)
{
	auto p = GET_PARTICLE(vd);

	Point<3,T> xp = vd.getPos(p);

    auto it = cl.getNNIterator(cl.getCell(xp));

    Point<3,T> force1({0.0,0.0,0.0});
    Point<3,T> force2({0.0,0.0,0.0});

    while (it.isNext())
    {
    	auto q1 = it.get_sort();
    	auto q2 = it.get();

    	if (q2 == p) {++it; continue;}

    	Point<3,T> xq_1 = vd_sort.getPos(q1);
    	Point<3,T> xq_2 = vd.getPos(q2);

    	Point<3,T> r1 = xq_1 - xp;
    	Point<3,T> r2 = xq_2 - xp;

    	// Normalize

    	if (r1.norm() > 1e-6)
    	{
    		r1 /= r1.norm();
    		force1 += vd_sort.template getProp<0>(q1)*r1;
    	}
    	if (r2.norm() > 1e-6)
    	{
    		r2 /= r2.norm();
    		force2 += vd.template getProp<0>(q2)*r2;
    	}

    	++it;
    }

    vd.template getProp<1>(p)[0] = force1.get(0);
    vd.template getProp<1>(p)[1] = force1.get(1);
    vd.template getProp<1>(p)[2] = force1.get(2);

    vd.template getProp<2>(p)[0] = force2.get(0);
    vd.template getProp<2>(p)[1] = force2.get(1);
    vd.template getProp<2>(p)[2] = force2.get(2);
}

template<typename T, typename CellList_type>
__global__  void calculate_force_full_sort(vector_dist_ker<3, T, aggregate<T, T[3], T [3]>> vd,
		                         	 	   CellList_type cl, int rank)
{
	unsigned int p;
	GET_PARTICLE_SORT(p,cl);

	Point<3,T> xp = vd.getPos(p);

    auto it = cl.getNNIterator(cl.getCell(xp));

    Point<3,T> force1({0.0,0.0,0.0});

    while (it.isNext())
    {
    	auto q1 = it.get_sort();

    	if (q1 == p) {++it; continue;}

    	Point<3,T> xq_1 = vd.getPos(q1);

    	Point<3,T> r1 = xq_1 - xp;

    	// Normalize

    	if (r1.norm() > 1e-6)
    	{
    		r1 /= r1.norm();

    		force1 += vd.template getProp<0>(q1)*r1;
    	}

    	++it;
    }

    vd.template getProp<1>(p)[0] = force1.get(0);
    vd.template getProp<1>(p)[1] = force1.get(1);
    vd.template getProp<1>(p)[2] = force1.get(2);
}

template<typename CellList_type, typename vector_type>
bool check_force(CellList_type & NN_cpu, vector_type & vd)
{
	typedef typename vector_type::stype St;

	auto it6 = vd.getDomainIterator();

	bool match = true;

	while (it6.isNext())
	{
		auto p = it6.get();

		Point<3,St> xp = vd.getPos(p);

		// Calculate on CPU

		Point<3,St> force({0.0,0.0,0.0});

		auto NNc = NN_cpu.getNNIterator(NN_cpu.getCell(xp));

		while (NNc.isNext())
		{
			auto q = NNc.get();

	    	if (q == p.getKey()) {++NNc; continue;}

	    	Point<3,St> xq_2 = vd.getPos(q);
	    	Point<3,St> r2 = xq_2 - xp;

	    	// Normalize

	    	if (r2.norm() > 1e-6)
	    	{
	    		r2 /= r2.norm();
	    		force += vd.template getProp<0>(q)*r2;
	    	}

			++NNc;
		}

		match &= fabs(vd.template getProp<1>(p)[0] - vd.template getProp<2>(p)[0]) < 0.0003;
		match &= fabs(vd.template getProp<1>(p)[1] - vd.template getProp<2>(p)[1]) < 0.0003;
		match &= fabs(vd.template getProp<1>(p)[2] - vd.template getProp<2>(p)[2]) < 0.0003;

		match &= fabs(vd.template getProp<1>(p)[0] - force.get(0)) < 0.0003;
		match &= fabs(vd.template getProp<1>(p)[1] - force.get(1)) < 0.0003;
		match &= fabs(vd.template getProp<1>(p)[2] - force.get(2)) < 0.0003;

		if (match == false)
		{
			std::cout << "ERROR: " << vd.template getProp<1>(p)[0]  << "   " << vd.template getProp<2>(p)[0] << std::endl;
	                std::cout << "ERROR: " << vd.template getProp<1>(p)[1]  << "   " << vd.template getProp<2>(p)[1] << std::endl;
	                std::cout << "ERROR: " << vd.template getProp<1>(p)[2]  << "   " << vd.template getProp<2>(p)[2] << std::endl;

	                std::cout << p.getKey() << " ERROR2: " << vd.template getProp<1>(p)[0] << "   " <<  force.get(0) << std::endl;
	                std::cout << p.getKey() << " ERROR2: " << vd.template getProp<1>(p)[1] << "   " <<  force.get(1) << std::endl;
	                std::cout << p.getKey() << " ERROR2: " << vd.template getProp<1>(p)[2] << "   " <<  force.get(2) << std::endl;


			break;
		}

		++it6;
	}

	return match;
}

BOOST_AUTO_TEST_CASE( vector_dist_gpu_ghost_get )
{
	auto & v_cl = create_vcluster();

	if (v_cl.size() > 16)
	{return;}

	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	// set the ghost based on the radius cut off (make just a little bit smaller than the spacing)
	Ghost<3,float> g(0.1);

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	vector_dist_gpu<3,float,aggregate<float,float[3],float[3]>> vd(1000,domain,bc,g);

	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		vd.getPos(p)[0] = (float)rand() / (float)RAND_MAX;
		vd.getPos(p)[1] = (float)rand() / (float)RAND_MAX;
		vd.getPos(p)[2] = (float)rand() / (float)RAND_MAX;

		vd.template getProp<0>(p) = vd.getPos(p)[0] + vd.getPos(p)[1] + vd.getPos(p)[2];

		vd.template getProp<1>(p)[0] = vd.getPos(p)[0] + vd.getPos(p)[1];
		vd.template getProp<1>(p)[1] = vd.getPos(p)[0] + vd.getPos(p)[2];
		vd.template getProp<1>(p)[2] = vd.getPos(p)[1] + vd.getPos(p)[2];

		vd.template getProp<2>(p)[0] = vd.getPos(p)[0] + 3.0*vd.getPos(p)[1];
		vd.template getProp<2>(p)[1] = vd.getPos(p)[0] + 3.0*vd.getPos(p)[2];
		vd.template getProp<2>(p)[2] = vd.getPos(p)[1] + 3.0*vd.getPos(p)[2];


		++it;
	}

	// Ok we redistribute the particles (CPU based)
	vd.map();

	vd.template ghost_get<0,1,2>();

	// Now we check the the ghost contain the correct information

	bool check = true;

	auto itg = vd.getDomainAndGhostIterator();

	while (itg.isNext())
	{
		auto p = itg.get();

		check &= (vd.template getProp<0>(p) == vd.getPos(p)[0] + vd.getPos(p)[1] + vd.getPos(p)[2]);

		check &= (vd.template getProp<1>(p)[0] == vd.getPos(p)[0] + vd.getPos(p)[1]);
		check &= (vd.template getProp<1>(p)[1] == vd.getPos(p)[0] + vd.getPos(p)[2]);
		check &= (vd.template getProp<1>(p)[2] == vd.getPos(p)[1] + vd.getPos(p)[2]);

		check &= (vd.template getProp<2>(p)[0] == vd.getPos(p)[0] + 3.0*vd.getPos(p)[1]);
		check &= (vd.template getProp<2>(p)[1] == vd.getPos(p)[0] + 3.0*vd.getPos(p)[2]);
		check &= (vd.template getProp<2>(p)[2] == vd.getPos(p)[1] + 3.0*vd.getPos(p)[2]);

		++itg;
	}

	size_t tot_s = vd.size_local_with_ghost();

	v_cl.sum(tot_s);
	v_cl.execute();

	// We check that we check something
	BOOST_REQUIRE(tot_s > 1000);
}

template<typename vector_type, typename CellList_type, typename CellList_type_cpu>
void check_cell_list_cpu_and_gpu(vector_type & vd, CellList_type & NN, CellList_type_cpu & NN_cpu)
{
	const auto it5 = vd.getDomainIteratorGPU(32);

	CUDA_LAUNCH((calculate_force<typename vector_type::stype,decltype(NN.toKernel())>),it5,vd.toKernel(),vd.toKernel_sorted(),NN.toKernel(),create_vcluster().rank());

	vd.template deviceToHostProp<1,2>();

	bool test = check_force(NN_cpu,vd);
	BOOST_REQUIRE_EQUAL(test,true);

	// We reset the property 1 on device

	auto rst = vd.getDomainIterator();

	while (rst.isNext())
	{
		auto p = rst.get();

		vd.template getProp<1>(p)[0] = 0.0;
		vd.template getProp<1>(p)[1] = 0.0;
		vd.template getProp<1>(p)[2] = 0.0;

		++rst;
	}

	vd.template hostToDeviceProp<1>();

	// We do exactly the same test as before, but now we completely use the sorted version

	CUDA_LAUNCH((calculate_force_full_sort<typename vector_type::stype,decltype(NN.toKernel())>),it5,vd.toKernel_sorted(),NN.toKernel(),create_vcluster().rank());

	vd.template merge_sort<1>(NN);
	vd.template deviceToHostProp<1>();

	test = check_force(NN_cpu,vd);
	BOOST_REQUIRE_EQUAL(test,true);
}

template<typename CellList_type>
void vector_dist_gpu_test_impl()
{
	auto & v_cl = create_vcluster();

	if (v_cl.size() > 16)
	{return;}

	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	// set the ghost based on the radius cut off (make just a little bit smaller than the spacing)
	Ghost<3,float> g(0.1);

	// Boundary conditions
	size_t bc[3]={NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

	vector_dist_gpu<3,float,aggregate<float,float[3],float[3]>> vd(10000,domain,bc,g);

	srand(55067*create_vcluster().rank());

	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		int x = rand();
		int y = rand();
		int z = rand();

		vd.getPos(p)[0] = (float)x / (float)RAND_MAX;
		vd.getPos(p)[1] = (float)y / (float)RAND_MAX;
		vd.getPos(p)[2] = (float)z / (float)RAND_MAX;

		Point<3,float> xp = vd.getPos(p);

		++it;
	}

	// Ok we redistribute the particles (CPU based)
	vd.map();

	size_t size_l = vd.size_local();

	v_cl.sum(size_l);
	v_cl.execute();

	BOOST_REQUIRE_EQUAL(size_l,10000);


	auto & ct = vd.getDecomposition();

	bool noOut = true;
	size_t cnt = 0;

	auto it2 = vd.getDomainIterator();

	while (it2.isNext())
	{
		auto p = it2.get();

		noOut &= ct.isLocal(vd.getPos(p));

		cnt++;
		++it2;
	}

	BOOST_REQUIRE_EQUAL(noOut,true);
	BOOST_REQUIRE_EQUAL(cnt,vd.size_local());

	// now we offload all the properties

	const auto it3 = vd.getDomainIteratorGPU();

	// offload to device
	vd.hostToDevicePos();

	CUDA_LAUNCH_DIM3(initialize_props,it3.wthr,it3.thr,vd.toKernel());

	// now we check what we initialized

	vd.deviceToHostProp<0,1>();

	auto it4 = vd.getDomainIterator();

	while (it4.isNext())
	{
		auto p = it4.get();

		BOOST_REQUIRE_CLOSE(vd.template getProp<0>(p),vd.getPos(p)[0] + vd.getPos(p)[1] + vd.getPos(p)[2],0.01);

		BOOST_REQUIRE_CLOSE(vd.template getProp<1>(p)[0],vd.getPos(p)[0] + vd.getPos(p)[1],0.01);
		BOOST_REQUIRE_CLOSE(vd.template getProp<1>(p)[1],vd.getPos(p)[0] + vd.getPos(p)[2],0.01);
		BOOST_REQUIRE_CLOSE(vd.template getProp<1>(p)[2],vd.getPos(p)[1] + vd.getPos(p)[2],0.01);

		++it4;
	}

	// here we do a ghost_get
	vd.ghost_get<0>();

	// Double ghost get to check crashes
	vd.ghost_get<0>();

	// we re-offload what we received
	vd.hostToDevicePos();
	vd.template hostToDeviceProp<0>();

	auto NN = vd.template getCellListGPU<CellList_type>(0.1);
	auto NN_cpu = vd.getCellList(0.1);
	check_cell_list_cpu_and_gpu(vd,NN,NN_cpu);

	auto NN_up = vd.template getCellListGPU<CellList_type>(0.1);

	NN_up.clear();
	vd.updateCellList(NN_up);
	check_cell_list_cpu_and_gpu(vd,NN_up,NN_cpu);
}

template<typename CellList_type>
void vector_dist_gpu_make_sort_test_impl()
{
	auto & v_cl = create_vcluster();

	if (v_cl.size() > 16)
	{return;}

	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	// set the ghost based on the radius cut off (make just a little bit smaller than the spacing)
	Ghost<3,float> g(0.1);

	// Boundary conditions
	size_t bc[3]={NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

	vector_dist_gpu<3,float,aggregate<float,float[3],float[3]>> vd(10000,domain,bc,g);

	srand(55067*create_vcluster().rank());

	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		int x = rand();
		int y = rand();
		int z = rand();

		vd.getPos(p)[0] = (float)x / (float)RAND_MAX;
		vd.getPos(p)[1] = (float)y / (float)RAND_MAX;
		vd.getPos(p)[2] = (float)z / (float)RAND_MAX;

		++it;
	}

	vd.hostToDevicePos();

	// Ok we redistribute the particles
	vd.map(RUN_ON_DEVICE);

	auto it3 = vd.getDomainIteratorGPU();

	CUDA_LAUNCH_DIM3(initialize_props,it3.wthr,it3.thr,vd.toKernel());

	// Here we check make sort does not mess-up particles we use a Cell-List to check that
	// the two cell-list constructed are identical

	vd.deviceToHostPos();

	auto NN_cpu1 = vd.getCellList(0.1);
	auto NN = vd.template getCellListGPU<CellList_type>(0.1);
	vd.make_sort(NN);

	vd.deviceToHostPos();

	auto NN_cpu2 = vd.getCellList(0.1);

	// here we compare the two cell-lists

	bool match = true;
	for (size_t i = 0 ; i < NN_cpu1.getNCells() ; i++)
	{
		match &= NN_cpu1.getNelements(i) == NN_cpu2.getNelements(i);
	}

	BOOST_REQUIRE_EQUAL(match,true);

	// In this second step we check that we can use make_sort_from to check we can sort partifcles even
	// when ghost are filled

	// Here we get do a make sort
	NN = vd.template getCellListGPU<CellList_type>(0.1);
	vd.make_sort(NN);

	openfpm::vector_gpu<aggregate<float,float[3],float[3]>> tmp_prp = vd.getPropVector();
	openfpm::vector_gpu<Point<3,float>> tmp_pos = vd.getPosVector();

	vd.deviceToHostPos();
	tmp_pos.template deviceToHost<0>();

	// here we do a ghost_get
	vd.ghost_get<0>(RUN_ON_DEVICE);

	// Here we get do a make sort
	NN = vd.template getCellListGPU<CellList_type>(0.1);

	CUDA_CHECK()

	vd.make_sort_from(NN);

	// Check

	tmp_pos.deviceToHost<0>();
	vd.deviceToHostPos();

	match = true;
	for (size_t i = 0 ; i < vd.size_local() ; i++)
	{
		Point<3,float> p1 = vd.getPos(i);
		Point<3,float> p2 = tmp_pos.template get<0>(i);

		// They must be in the same cell
		auto c1 = NN.getCell(p1);
		auto c2 = NN.getCell(p1);

		match &= c1 == c2;
	}

	BOOST_REQUIRE_EQUAL(match,true);
}


BOOST_AUTO_TEST_CASE(vector_dist_gpu_make_sort_sparse)
{
	vector_dist_gpu_make_sort_test_impl<CELLLIST_GPU_SPARSE<3,float>>();
}

BOOST_AUTO_TEST_CASE(vector_dist_gpu_make_sort)
{
	vector_dist_gpu_make_sort_test_impl<CellList_gpu<3,float,CudaMemory,shift_only<3, float>>>();
}

BOOST_AUTO_TEST_CASE( vector_dist_gpu_test)
{
	vector_dist_gpu_test_impl<CellList_gpu<3,float,CudaMemory,shift_only<3, float>>>();
}

BOOST_AUTO_TEST_CASE( vector_dist_gpu_test_sparse)
{
	vector_dist_gpu_test_impl<CELLLIST_GPU_SPARSE<3,float>>();
}

template<typename St>
void vdist_calc_gpu_test()
{
	auto & v_cl = create_vcluster();

	if (v_cl.size() > 16)
	{return;}

	Box<3,St> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	// set the ghost based on the radius cut off (make just a little bit smaller than the spacing)
	Ghost<3,St> g(0.1);

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	//! [Create a gpu vector]

	vector_dist_gpu<3,St,aggregate<St,St[3],St[3]>> vd(1000,domain,bc,g);

	//! [Create a gpu vector]

	//! [Fill gpu vector and move to GPU]

	srand(v_cl.rank()*10000);
	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		vd.getPos(p)[0] = (St)rand() / (float)RAND_MAX;
		vd.getPos(p)[1] = (St)rand() / (float)RAND_MAX;
		vd.getPos(p)[2] = (St)rand() / (float)RAND_MAX;

		vd.template getProp<0>(p) = vd.getPos(p)[0] + vd.getPos(p)[1] + vd.getPos(p)[2];

		vd.template getProp<1>(p)[0] = vd.getPos(p)[0];
		vd.template getProp<1>(p)[1] = vd.getPos(p)[1];
		vd.template getProp<1>(p)[2] = vd.getPos(p)[2];

		vd.template getProp<2>(p)[0] = vd.getPos(p)[0] + vd.getPos(p)[1];
		vd.template getProp<2>(p)[1] = vd.getPos(p)[0] + vd.getPos(p)[2];
		vd.template getProp<2>(p)[2] = vd.getPos(p)[1] + vd.getPos(p)[2];

		++it;
	}

	// move on device
	vd.hostToDevicePos();
	vd.template hostToDeviceProp<0,1,2>();

	// Ok we redistribute the particles (GPU based)
	vd.map(RUN_ON_DEVICE);

	//! [Fill gpu vector and move to GPU]

	vd.deviceToHostPos();
	vd.template deviceToHostProp<0,1,2>();

	// Reset the host part

	auto it3 = vd.getDomainIterator();

	while (it3.isNext())
	{
		auto p = it3.get();

		vd.getPos(p)[0] = 1.0;
		vd.getPos(p)[1] = 1.0;
		vd.getPos(p)[2] = 1.0;

		vd.template getProp<0>(p) = 0.0;

		vd.template getProp<0>(p) = 0.0;
		vd.template getProp<0>(p) = 0.0;
		vd.template getProp<0>(p) = 0.0;

		vd.template getProp<0>(p) = 0.0;
		vd.template getProp<0>(p) = 0.0;
		vd.template getProp<0>(p) = 0.0;

		++it3;
	}

	// we move from Device to CPU

	vd.deviceToHostPos();
	vd.template deviceToHostProp<0,1,2>();

	// Check

	auto it2 = vd.getDomainIterator();

	bool match = true;
	while (it2.isNext())
	{
		auto p = it2.get();

		match &= vd.template getProp<0>(p) == vd.getPos(p)[0] + vd.getPos(p)[1] + vd.getPos(p)[2];

		match &= vd.template getProp<1>(p)[0] == vd.getPos(p)[0];
		match &= vd.template getProp<1>(p)[1] == vd.getPos(p)[1];
		match &= vd.template getProp<1>(p)[2] == vd.getPos(p)[2];

		match &= vd.template getProp<2>(p)[0] == vd.getPos(p)[0] + vd.getPos(p)[1];
		match &= vd.template getProp<2>(p)[1] == vd.getPos(p)[0] + vd.getPos(p)[2];
		match &= vd.template getProp<2>(p)[2] == vd.getPos(p)[1] + vd.getPos(p)[2];

		++it2;
	}

	BOOST_REQUIRE_EQUAL(match,true);

	// count local particles

	size_t l_cnt = 0;
	size_t nl_cnt = 0;
	size_t n_out = 0;

	// Domain + ghost box
	Box<3,St> dom_ext = domain;
	dom_ext.enlarge(g);

	auto it5 = vd.getDomainIterator();
	count_local_n_local<3>(vd,it5,bc,domain,dom_ext,l_cnt,nl_cnt,n_out);

	BOOST_REQUIRE_EQUAL(n_out,0);
	BOOST_REQUIRE_EQUAL(l_cnt,vd.size_local());

	// we do 10 gpu steps (using a cpu vector to check that map and ghost get work as expented)

	for (size_t i = 0 ; i < 10 ; i++)
	{
		vd.map(RUN_ON_DEVICE);

		vd.deviceToHostPos();
		vd.template deviceToHostProp<0,1,2>();

		// To test we copy on a cpu distributed vector and we do a map

		vector_dist<3,St,aggregate<St,St[3],St[3]>> vd_cpu(vd.getDecomposition().template duplicate_convert<HeapMemory,memory_traits_lin>(),0);

		auto itc = vd.getDomainIterator();

		while (itc.isNext())
		{
			auto p = itc.get();

			vd_cpu.add();

			vd_cpu.getLastPos()[0] = vd.getPos(p)[0];
			vd_cpu.getLastPos()[1] = vd.getPos(p)[1];
			vd_cpu.getLastPos()[2] = vd.getPos(p)[2];

			vd_cpu.template getLastProp<0>() = vd.template getProp<0>(p);

			vd_cpu.template getLastProp<1>()[0] = vd.template getProp<1>(p)[0];
			vd_cpu.template getLastProp<1>()[1] = vd.template getProp<1>(p)[1];
			vd_cpu.template getLastProp<1>()[2] = vd.template getProp<1>(p)[2];

			vd_cpu.template getLastProp<2>()[0] = vd.template getProp<2>(p)[0];
			vd_cpu.template getLastProp<2>()[1] = vd.template getProp<2>(p)[1];
			vd_cpu.template getLastProp<2>()[2] = vd.template getProp<2>(p)[2];

			++itc;
		}

		vd_cpu.template ghost_get<0,1,2>();

		//! [Fill the ghost on GPU]

		vd.template ghost_get<0,1,2>(RUN_ON_DEVICE);

		//! [Fill the ghost on GPU]

		vd.deviceToHostPos();
		vd.template deviceToHostProp<0,1,2>();

		match = true;

		// Particle on the gpu ghost and cpu ghost are not ordered in the same way so we have to reorder

		struct part
		{
			Point<3,St> xp;


			St prp0;
			St prp1[3];
			St prp2[3];

			bool operator<(const part & tmp) const
			{
				if (xp.get(0) < tmp.xp.get(0))
				{return true;}
				else if (xp.get(0) > tmp.xp.get(0))
				{return false;}

				if (xp.get(1) < tmp.xp.get(1))
				{return true;}
				else if (xp.get(1) > tmp.xp.get(1))
				{return false;}

				if (xp.get(2) < tmp.xp.get(2))
				{return true;}
				else if (xp.get(2) > tmp.xp.get(2))
				{return false;}

				return false;
			}
		};

		openfpm::vector<part> cpu_sort;
		openfpm::vector<part> gpu_sort;

		cpu_sort.resize(vd_cpu.size_local_with_ghost() - vd_cpu.size_local());
		gpu_sort.resize(vd.size_local_with_ghost() - vd.size_local());

		BOOST_REQUIRE_EQUAL(cpu_sort.size(),gpu_sort.size());

		size_t cnt = 0;

		auto itc2 = vd.getGhostIterator();
		while (itc2.isNext())
		{
			auto p = itc2.get();

			cpu_sort.get(cnt).xp.get(0) = vd_cpu.getPos(p)[0];
			gpu_sort.get(cnt).xp.get(0) = vd.getPos(p)[0];
			cpu_sort.get(cnt).xp.get(1) = vd_cpu.getPos(p)[1];
			gpu_sort.get(cnt).xp.get(1) = vd.getPos(p)[1];
			cpu_sort.get(cnt).xp.get(2) = vd_cpu.getPos(p)[2];
			gpu_sort.get(cnt).xp.get(2) = vd.getPos(p)[2];

			cpu_sort.get(cnt).prp0 = vd_cpu.template getProp<0>(p);
			gpu_sort.get(cnt).prp0 = vd.template getProp<0>(p);

			cpu_sort.get(cnt).prp1[0] = vd_cpu.template getProp<1>(p)[0];
			gpu_sort.get(cnt).prp1[0] = vd.template getProp<1>(p)[0];
			cpu_sort.get(cnt).prp1[1] = vd_cpu.template getProp<1>(p)[1];
			gpu_sort.get(cnt).prp1[1] = vd.template getProp<1>(p)[1];
			cpu_sort.get(cnt).prp1[2] = vd_cpu.template getProp<1>(p)[2];
			gpu_sort.get(cnt).prp1[2] = vd.template getProp<1>(p)[2];

			cpu_sort.get(cnt).prp2[0] = vd_cpu.template getProp<2>(p)[0];
			gpu_sort.get(cnt).prp2[0] = vd.template getProp<2>(p)[0];
			cpu_sort.get(cnt).prp2[1] = vd_cpu.template getProp<2>(p)[1];
			gpu_sort.get(cnt).prp2[1] = vd.template getProp<2>(p)[1];
			cpu_sort.get(cnt).prp2[2] = vd_cpu.template getProp<2>(p)[2];
			gpu_sort.get(cnt).prp2[2] = vd.template getProp<2>(p)[2];

			++cnt;
			++itc2;
		}

		cpu_sort.sort();
		gpu_sort.sort();

		for (size_t i = 0 ; i < cpu_sort.size() ; i++)
		{
			match &= cpu_sort.get(i).xp.get(0) == gpu_sort.get(i).xp.get(0);
			match &= cpu_sort.get(i).xp.get(1) == gpu_sort.get(i).xp.get(1);
			match &= cpu_sort.get(i).xp.get(2) == gpu_sort.get(i).xp.get(2);

			match &= cpu_sort.get(i).prp0 == gpu_sort.get(i).prp0;
			match &= cpu_sort.get(i).prp1[0] == gpu_sort.get(i).prp1[0];
			match &= cpu_sort.get(i).prp1[1] == gpu_sort.get(i).prp1[1];
			match &= cpu_sort.get(i).prp1[2] == gpu_sort.get(i).prp1[2];

			match &= cpu_sort.get(i).prp2[0] == gpu_sort.get(i).prp2[0];
			match &= cpu_sort.get(i).prp2[1] == gpu_sort.get(i).prp2[1];
			match &= cpu_sort.get(i).prp2[2] == gpu_sort.get(i).prp2[2];
		}

		BOOST_REQUIRE_EQUAL(match,true);

		// move particles on gpu

		auto ite = vd.getDomainIteratorGPU();
		CUDA_LAUNCH_DIM3((move_parts_gpu_test<3,decltype(vd.toKernel())>),ite.wthr,ite.thr,vd.toKernel());
	}
}

BOOST_AUTO_TEST_CASE( vector_dist_map_on_gpu_test)
{
	vdist_calc_gpu_test<float>();
	vdist_calc_gpu_test<double>();
}

BOOST_AUTO_TEST_CASE(vector_dist_reduce)
{
	auto & v_cl = create_vcluster();

	if (v_cl.size() > 16)
	{return;}

	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	// set the ghost based on the radius cut off (make just a little bit smaller than the spacing)
	Ghost<3,float> g(0.1);

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	vector_dist_gpu<3,float,aggregate<float,double,int,size_t>> vd(5000*v_cl.size(),domain,bc,g);

	auto it = vd.getDomainIterator();

	float fc = 1.0;
	double dc = 1.0;
	int ic = 1.0;
	size_t sc = 1.0;

	while(it.isNext())
	{
		auto p = it.get();

		vd.template getProp<0>(p) = fc;
		vd.template getProp<1>(p) = dc;
		vd.template getProp<2>(p) = ic;
		vd.template getProp<3>(p) = sc;

		fc += 1.0;
		dc += 1.0;
		ic += 1;
		sc += 1;

		++it;
	}

	vd.template hostToDeviceProp<0,1,2,3>();

	float redf = reduce_local<0,_add_>(vd);
	double redd = reduce_local<1,_add_>(vd);
	int redi = reduce_local<2,_add_>(vd);
	size_t reds = reduce_local<3,_add_>(vd);

	BOOST_REQUIRE_EQUAL(redf,(vd.size_local()+1.0)*(vd.size_local())/2.0);
	BOOST_REQUIRE_EQUAL(redd,(vd.size_local()+1.0)*(vd.size_local())/2.0);
	BOOST_REQUIRE_EQUAL(redi,(vd.size_local()+1)*(vd.size_local())/2);
	BOOST_REQUIRE_EQUAL(reds,(vd.size_local()+1)*(vd.size_local())/2);

	float redf2 = reduce_local<0,_max_>(vd);
	double redd2 = reduce_local<1,_max_>(vd);
	int redi2 = reduce_local<2,_max_>(vd);
	size_t reds2 = reduce_local<3,_max_>(vd);

	BOOST_REQUIRE_EQUAL(redf2,vd.size_local());
	BOOST_REQUIRE_EQUAL(redd2,vd.size_local());
	BOOST_REQUIRE_EQUAL(redi2,vd.size_local());
	BOOST_REQUIRE_EQUAL(reds2,vd.size_local());
}

template<typename CellList_type>
void vector_dist_dlb_on_cuda_impl(size_t k,double r_cut)
{
	std::random_device r;

    std::seed_seq seed2{/*r() +*/ create_vcluster().rank(),
    					/*r() +*/ create_vcluster().rank(),
    					/*r() +*/ create_vcluster().rank(),
    					/*r() +*/ create_vcluster().rank(),
    					/*r() +*/ create_vcluster().rank(),
    					/*r() +*/ create_vcluster().rank(),
    					/*r() +*/ create_vcluster().rank(),
    					/*r() +*/ create_vcluster().rank()};
    std::mt19937 e2(seed2);

	typedef vector_dist_gpu<3,double,aggregate<double,double[3],double[3]>> vector_type;

	Vcluster<> & v_cl = create_vcluster();

	if (v_cl.getProcessingUnits() > 8)
		return;

	std::uniform_real_distribution<double> unif(0.0,0.3);

	Box<3,double> domain({0.0,0.0,0.0},{1.0,1.0,1.0});
	Ghost<3,double> g(0.1);
	size_t bc[3] = {PERIODIC,PERIODIC,PERIODIC};

	vector_type vd(0,domain,bc,g,DEC_GRAN(2048));

	// Only processor 0 initialy add particles on a corner of a domain

	if (v_cl.getProcessUnitID() == 0)
	{
		for(size_t i = 0 ; i < k ; i++)
		{
			vd.add();

			vd.getLastPos()[0] = unif(e2);
			vd.getLastPos()[1] = unif(e2);
			vd.getLastPos()[2] = unif(e2);
		}
	}

	// Move to GPU
	vd.hostToDevicePos();
	vd.template hostToDeviceProp<0>();

	vd.map(RUN_ON_DEVICE);
	vd.template ghost_get<>(RUN_ON_DEVICE);

	// now move to CPU

	vd.deviceToHostPos();
	vd.template deviceToHostProp<0>();

	// Get the neighborhood of each particles

	auto VV = vd.getVerlet(r_cut);

	// store the number of neighborhood for each particles

	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		vd.template getProp<0>(p) = VV.getNNPart(p.getKey());

		++it;
	}

	// Move to GPU
	vd.template hostToDeviceProp<0>();

	ModelSquare md;
	md.factor = 10;
	vd.addComputationCosts(md);
	vd.getDecomposition().decompose();
	vd.map(RUN_ON_DEVICE);

	vd.deviceToHostPos();
	// Move info to CPU for addComputationcosts

	vd.addComputationCosts(md);

	openfpm::vector<size_t> loads;
	size_t load = vd.getDecomposition().getDistribution().getProcessorLoad();
	v_cl.allGather(load,loads);
	v_cl.execute();

	for (size_t i = 0 ; i < loads.size() ; i++)
	{
		double load_f = load;
		double load_fc = loads.get(i);

		BOOST_REQUIRE_CLOSE(load_f,load_fc,7.0);
	}

	BOOST_REQUIRE(vd.size_local() != 0);

	Point<3,double> v({1.0,1.0,1.0});

	for (size_t i = 0 ; i < 25 ; i++)
	{
		// move particles to CPU and move the particles by 0.1

		vd.deviceToHostPos();

		auto it = vd.getDomainIterator();

		while (it.isNext())
		{
			auto p = it.get();

			vd.getPos(p)[0] += v.get(0) * 0.09;
			vd.getPos(p)[1] += v.get(1) * 0.09;
			vd.getPos(p)[2] += v.get(2) * 0.09;

			++it;
		}

		//Back to GPU
		vd.hostToDevicePos();
		vd.map(RUN_ON_DEVICE);
		vd.template ghost_get<0>(RUN_ON_DEVICE);

		vd.deviceToHostPos();
		vd.template deviceToHostProp<0,1,2>();

		// Check calc forces
		auto NN_gpu = vd.template getCellListGPU<CellList_type>(r_cut);
		auto NN_cpu = vd.getCellList(r_cut);
		check_cell_list_cpu_and_gpu(vd,NN_gpu,NN_cpu);

		auto VV2 = vd.getVerlet(r_cut);

		auto it2 = vd.getDomainIterator();

		bool match = true;
		while (it2.isNext())
		{
			auto p = it2.get();

			match &= vd.template getProp<0>(p) == VV2.getNNPart(p.getKey());

			++it2;
		}

		BOOST_REQUIRE_EQUAL(match,true);

		ModelSquare md;
		vd.addComputationCosts(md);
		vd.getDecomposition().redecompose(200);
		vd.map(RUN_ON_DEVICE);

		BOOST_REQUIRE(vd.size_local() != 0);

		vd.template ghost_get<0>(RUN_ON_DEVICE);
		vd.deviceToHostPos();
		vd.template deviceToHostProp<0>();

		vd.addComputationCosts(md);

		openfpm::vector<size_t> loads;
		size_t load = vd.getDecomposition().getDistribution().getProcessorLoad();
		v_cl.allGather(load,loads);
		v_cl.execute();

		for (size_t i = 0 ; i < loads.size() ; i++)
		{
			double load_f = load;
			double load_fc = loads.get(i);

#ifdef ENABLE_ASAN
			BOOST_REQUIRE_CLOSE(load_f,load_fc,30.0);
#else
			BOOST_REQUIRE_CLOSE(load_f,load_fc,10.0);
#endif
		}
	}
}

template<typename CellList_type>
void vector_dist_dlb_on_cuda_impl_async(size_t k,double r_cut)
{
	std::random_device r;

    std::seed_seq seed2{r() + create_vcluster().rank(),
    					r() + create_vcluster().rank(),
    					r() + create_vcluster().rank(),
    					r() + create_vcluster().rank(),
    					r() + create_vcluster().rank(),
    					r() + create_vcluster().rank(),
    					r() + create_vcluster().rank(),
    					r() + create_vcluster().rank()};
    std::mt19937 e2(seed2);

	typedef vector_dist_gpu<3,double,aggregate<double,double[3],double[3]>> vector_type;

	Vcluster<> & v_cl = create_vcluster();

	if (v_cl.getProcessingUnits() > 8)
		return;

	std::uniform_real_distribution<double> unif(0.0,0.3);

	Box<3,double> domain({0.0,0.0,0.0},{1.0,1.0,1.0});
	Ghost<3,double> g(0.1);
	size_t bc[3] = {PERIODIC,PERIODIC,PERIODIC};

	vector_type vd(0,domain,bc,g,DEC_GRAN(2048));

	// Only processor 0 initialy add particles on a corner of a domain

	if (v_cl.getProcessUnitID() == 0)
	{
		for(size_t i = 0 ; i < k ; i++)
		{
			vd.add();

			vd.getLastPos()[0] = unif(e2);
			vd.getLastPos()[1] = unif(e2);
			vd.getLastPos()[2] = unif(e2);
		}
	}

	// Move to GPU
	vd.hostToDevicePos();
	vd.template hostToDeviceProp<0>();

	vd.map(RUN_ON_DEVICE);
	vd.template Ighost_get<>(RUN_ON_DEVICE);
	vd.template ghost_wait<>(RUN_ON_DEVICE);

	// now move to CPU

	vd.deviceToHostPos();
	vd.template deviceToHostProp<0>();

	// Get the neighborhood of each particles

	auto VV = vd.getVerlet(r_cut);

	// store the number of neighborhood for each particles

	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		vd.template getProp<0>(p) = VV.getNNPart(p.getKey());

		++it;
	}

	// Move to GPU
	vd.template hostToDeviceProp<0>();

	ModelSquare md;
	md.factor = 10;
	vd.addComputationCosts(md);
	vd.getDecomposition().decompose();
	vd.map(RUN_ON_DEVICE);

	vd.deviceToHostPos();
	// Move info to CPU for addComputationcosts

	vd.addComputationCosts(md);

	openfpm::vector<size_t> loads;
	size_t load = vd.getDecomposition().getDistribution().getProcessorLoad();
	v_cl.allGather(load,loads);
	v_cl.execute();

	for (size_t i = 0 ; i < loads.size() ; i++)
	{
		double load_f = load;
		double load_fc = loads.get(i);

		BOOST_REQUIRE_CLOSE(load_f,load_fc,7.0);
	}

	BOOST_REQUIRE(vd.size_local() != 0);

	Point<3,double> v({1.0,1.0,1.0});

	for (size_t i = 0 ; i < 25 ; i++)
	{
		// move particles to CPU and move the particles by 0.1

		vd.deviceToHostPos();

		auto it = vd.getDomainIterator();

		while (it.isNext())
		{
			auto p = it.get();

			vd.getPos(p)[0] += v.get(0) * 0.09;
			vd.getPos(p)[1] += v.get(1) * 0.09;
			vd.getPos(p)[2] += v.get(2) * 0.09;

			++it;
		}

		// Back to GPU
		vd.hostToDevicePos();
		vd.map(RUN_ON_DEVICE);
		vd.template Ighost_get<0>(RUN_ON_DEVICE);
		vd.template ghost_wait<0>(RUN_ON_DEVICE);
		vd.deviceToHostPos();
		vd.template deviceToHostProp<0,1,2>();

		// Check calc forces
		auto NN_gpu = vd.template getCellListGPU<CellList_type>(r_cut);
		auto NN_cpu = vd.getCellList(r_cut);
		check_cell_list_cpu_and_gpu(vd,NN_gpu,NN_cpu);

		auto VV2 = vd.getVerlet(r_cut);

		auto it2 = vd.getDomainIterator();

		bool match = true;
		while (it2.isNext())
		{
			auto p = it2.get();

			match &= vd.template getProp<0>(p) == VV2.getNNPart(p.getKey());

			++it2;
		}

		BOOST_REQUIRE_EQUAL(match,true);

		ModelSquare md;
		vd.addComputationCosts(md);
		vd.getDecomposition().redecompose(200);
		vd.map(RUN_ON_DEVICE);

		BOOST_REQUIRE(vd.size_local() != 0);

//		vd.template ghost_get<0>(RUN_ON_DEVICE);
//		vd.template ghost_get<0>(RUN_ON_DEVICE);
		vd.template Ighost_get<0>(RUN_ON_DEVICE);
		vd.template ghost_wait<0>(RUN_ON_DEVICE);
		vd.deviceToHostPos();
		vd.template deviceToHostProp<0>();

		vd.addComputationCosts(md);

		openfpm::vector<size_t> loads;
		size_t load = vd.getDecomposition().getDistribution().getProcessorLoad();
		v_cl.allGather(load,loads);
		v_cl.execute();

		for (size_t i = 0 ; i < loads.size() ; i++)
		{
			double load_f = load;
			double load_fc = loads.get(i);

			BOOST_REQUIRE_CLOSE(load_f,load_fc,10.0);
		}
	}
}

BOOST_AUTO_TEST_CASE(vector_dist_dlb_on_cuda_async)
{
	vector_dist_dlb_on_cuda_impl_async<CellList_gpu<3,double,CudaMemory,shift_only<3,double>,unsigned int,int,false>>(50000,0.01);
}

BOOST_AUTO_TEST_CASE(vector_dist_dlb_on_cuda)
{
	vector_dist_dlb_on_cuda_impl<CellList_gpu<3,double,CudaMemory,shift_only<3,double>,unsigned int,int,false>>(50000,0.01);
}

BOOST_AUTO_TEST_CASE(vector_dist_dlb_on_cuda_sparse)
{
	vector_dist_dlb_on_cuda_impl<CELLLIST_GPU_SPARSE<3,double>>(50000,0.01);
}

BOOST_AUTO_TEST_CASE(vector_dist_dlb_on_cuda2)
{
	if (create_vcluster().size() <= 3)
	{return;};

	#ifndef CUDA_ON_CPU
	vector_dist_dlb_on_cuda_impl<CellList_gpu<3,double,CudaMemory,shift_only<3,double>,unsigned int,int,false>>(1000000,0.01);
	#endif
}

BOOST_AUTO_TEST_CASE(vector_dist_dlb_on_cuda3)
{
	if (create_vcluster().size() < 8)
	{return;}

	#ifndef CUDA_ON_CPU
	vector_dist_dlb_on_cuda_impl<CellList_gpu<3,double,CudaMemory,shift_only<3,double>,unsigned int,int,false>>(15000000,0.005);
	#endif
}


BOOST_AUTO_TEST_CASE(vector_dist_keep_prop_on_cuda)
{
	typedef vector_dist_gpu<3,double,aggregate<double,double[3],double[3][3]>> vector_type;

	Vcluster<> & v_cl = create_vcluster();

	if (v_cl.getProcessingUnits() > 8)
		return;

	Box<3,double> domain({0.0,0.0,0.0},{1.0,1.0,1.0});
	Ghost<3,double> g(0.1);
	size_t bc[3] = {PERIODIC,PERIODIC,PERIODIC};

	vector_type vd(0,domain,bc,g,DEC_GRAN(2048));

	// Only processor 0 initialy add particles on a corner of a domain

	if (v_cl.getProcessUnitID() == 0)
	{
		for(size_t i = 0 ; i < 50000 ; i++)
		{
			vd.add();

			vd.getLastPos()[0] = ((double)rand())/RAND_MAX * 0.3;
			vd.getLastPos()[1] = ((double)rand())/RAND_MAX * 0.3;
			vd.getLastPos()[2] = ((double)rand())/RAND_MAX * 0.3;
		}
	}

	// Move to GPU
	vd.hostToDevicePos();
	vd.template hostToDeviceProp<0>();

	vd.map(RUN_ON_DEVICE);
	vd.template ghost_get<>(RUN_ON_DEVICE);

	// now move to CPU

	vd.deviceToHostPos();
	vd.template deviceToHostProp<0>();


	// store the number of neighborhood for each particles

	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		vd.template getProp<0>(p) = 0.0;

		vd.template getProp<1>(p)[0] = 1000.0;
		vd.template getProp<1>(p)[1] = 2000.0;
		vd.template getProp<1>(p)[2] = 3000.0;

		vd.template getProp<2>(p)[0][0] = 6000,0;
		vd.template getProp<2>(p)[0][1] = 7000.0;
		vd.template getProp<2>(p)[0][2] = 8000.0;
		vd.template getProp<2>(p)[1][0] = 9000.0;
		vd.template getProp<2>(p)[1][1] = 10000.0;
		vd.template getProp<2>(p)[1][2] = 11000.0;
		vd.template getProp<2>(p)[2][0] = 12000.0;
		vd.template getProp<2>(p)[2][1] = 13000.0;
		vd.template getProp<2>(p)[2][2] = 14000.0;

		++it;
	}

	// Move to GPU
	vd.template hostToDeviceProp<0,1,2>();

	ModelSquare md;
	md.factor = 10;
	vd.addComputationCosts(md);
	vd.getDecomposition().decompose();
	vd.map(RUN_ON_DEVICE);

	vd.deviceToHostPos();
	// Move info to CPU for addComputationcosts

	vd.addComputationCosts(md);

	openfpm::vector<size_t> loads;
	size_t load = vd.getDecomposition().getDistribution().getProcessorLoad();
	v_cl.allGather(load,loads);
	v_cl.execute();

	for (size_t i = 0 ; i < loads.size() ; i++)
	{
		double load_f = load;
		double load_fc = loads.get(i);

		BOOST_REQUIRE_CLOSE(load_f,load_fc,7.0);
	}

	BOOST_REQUIRE(vd.size_local() != 0);

	Point<3,double> v({1.0,1.0,1.0});

	int base = 0;

	for (size_t i = 0 ; i < 25 ; i++)
	{
		if (i % 2 == 0)
		{
			// move particles to CPU and move the particles by 0.1

			vd.deviceToHostPos();

			auto it = vd.getDomainIterator();

			while (it.isNext())
			{
				auto p = it.get();

				vd.getPos(p)[0] += v.get(0) * 0.09;
				vd.getPos(p)[1] += v.get(1) * 0.09;
				vd.getPos(p)[2] += v.get(2) * 0.09;

				++it;
			}

			//Back to GPU
			vd.hostToDevicePos();
			vd.map(RUN_ON_DEVICE);
			vd.template ghost_get<>(RUN_ON_DEVICE);
			vd.deviceToHostPos();
			vd.template deviceToHostProp<0,1,2>();

			ModelSquare md;
			vd.addComputationCosts(md);
			vd.getDecomposition().redecompose(200);
			vd.map(RUN_ON_DEVICE);

			BOOST_REQUIRE(vd.size_local() != 0);

			vd.template ghost_get<0>(RUN_ON_DEVICE);
			vd.deviceToHostPos();
			vd.template deviceToHostProp<0,1,2>();

			vd.addComputationCosts(md);

			openfpm::vector<size_t> loads;
			size_t load = vd.getDecomposition().getDistribution().getProcessorLoad();
			v_cl.allGather(load,loads);
			v_cl.execute();

			for (size_t i = 0 ; i < loads.size() ; i++)
			{
				double load_f = load;
				double load_fc = loads.get(i);

				BOOST_REQUIRE_CLOSE(load_f,load_fc,10.0);
			}
		}
		else
		{
			vd.template deviceToHostProp<0,1,2>();

			auto it2 = vd.getDomainIterator();

			bool match = true;
			while (it2.isNext())
			{
				auto p = it2.get();

				vd.template getProp<0>(p) += 1;

				vd.template getProp<1>(p)[0] += 1.0;
				vd.template getProp<1>(p)[1] += 1.0;
				vd.template getProp<1>(p)[2] += 1.0;

				vd.template getProp<2>(p)[0][0] += 1.0;
				vd.template getProp<2>(p)[0][1] += 1.0;
				vd.template getProp<2>(p)[0][2] += 1.0;
				vd.template getProp<2>(p)[1][0] += 1.0;
				vd.template getProp<2>(p)[1][1] += 1.0;
				vd.template getProp<2>(p)[1][2] += 1.0;
				vd.template getProp<2>(p)[2][0] += 1.0;
				vd.template getProp<2>(p)[2][1] += 1.0;
				vd.template getProp<2>(p)[2][2] += 1.0;

				++it2;
			}

			vd.template hostToDeviceProp<0,1,2>();

			++base;

			vd.template ghost_get<0,1,2>(RUN_ON_DEVICE | KEEP_PROPERTIES);
			vd.template deviceToHostProp<0,1,2>();

			// Check that the ghost contain the correct information

			auto itg = vd.getGhostIterator();

			while (itg.isNext())
			{
				auto p = itg.get();

				match &= vd.template getProp<0>(p) == base;

				match &= vd.template getProp<1>(p)[0] == base + 1000.0;
				match &= vd.template getProp<1>(p)[1] == base + 2000.0;
				match &= vd.template getProp<1>(p)[2] == base + 3000.0;

				match &= vd.template getProp<2>(p)[0][0] == base + 6000.0;
				match &= vd.template getProp<2>(p)[0][1] == base + 7000.0;
				match &= vd.template getProp<2>(p)[0][2] == base + 8000.0;
				match &= vd.template getProp<2>(p)[1][0] == base + 9000.0;
				match &= vd.template getProp<2>(p)[1][1] == base + 10000.0;
				match &= vd.template getProp<2>(p)[1][2] == base + 11000.0;
				match &= vd.template getProp<2>(p)[2][0] == base + 12000.0;
				match &= vd.template getProp<2>(p)[2][1] == base + 13000.0;
				match &= vd.template getProp<2>(p)[2][2] == base + 14000.0;

				++itg;
			}

			BOOST_REQUIRE_EQUAL(match,true);
		}
	}
}

struct type_is_one
{
	__device__ static bool check(int c)
	{
		return c == 1;
	}
};

BOOST_AUTO_TEST_CASE(vector_dist_get_index_set)
{
	Box<3,double> domain({0.0,0.0,0.0},{1.0,1.0,1.0});
	Ghost<3,double> g(0.1);
	size_t bc[3] = {PERIODIC,PERIODIC,PERIODIC};

	if (create_vcluster().size() >= 16)
	{return;}

	Vcluster<> & v_cl = create_vcluster();

	vector_dist_gpu<3,double,aggregate<int,double>> vdg(10000,domain,bc,g,DEC_GRAN(128));

	auto it = vdg.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		vdg.getPos(p)[0] = (double)rand() / RAND_MAX;
		vdg.getPos(p)[1] = (double)rand() / RAND_MAX;
		vdg.getPos(p)[2] = (double)rand() / RAND_MAX;

		vdg.template getProp<0>(p) = (int)((double)rand() / RAND_MAX / 0.5);

		vdg.template getProp<1>(p) = (double)rand() / RAND_MAX;

		++it;
	}

	vdg.map();

	vdg.hostToDeviceProp<0,1>();
	vdg.hostToDevicePos();

	auto cl = vdg.getCellListGPU(0.1);

	// than we get a cell-list to force reorder

	openfpm::vector_gpu<aggregate<unsigned int>> ids;

	get_indexes_by_type<0,type_is_one>(vdg.getPropVectorSort(),ids,vdg.size_local(),v_cl.getgpuContext());

	// test

	ids.template deviceToHost<0>();

	auto & vs = vdg.getPropVectorSort();
	vs.template deviceToHost<0>();

	bool match = true;

	for (int i = 0 ; i < ids.size() ; i++)
	{
		if (vs.template get<0>(ids.template get<0>(i)) != 1)
		{match = false;}
	}

	BOOST_REQUIRE_EQUAL(match,true);
}

BOOST_AUTO_TEST_CASE(vector_dist_compare_host_device)
{
	Box<3,double> domain({0.0,0.0,0.0},{1.0,1.0,1.0});
	Ghost<3,double> g(0.1);
	size_t bc[3] = {PERIODIC,PERIODIC,PERIODIC};

	if (create_vcluster().size() >= 16)
	{return;}

	vector_dist_gpu<3,double,aggregate<double,double[3],double[3][3]>> vdg(10000,domain,bc,g,DEC_GRAN(128));

	auto it = vdg.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		vdg.getPos(p)[0] = (double)rand() / RAND_MAX;
		vdg.getPos(p)[1] = (double)rand() / RAND_MAX;
		vdg.getPos(p)[2] = (double)rand() / RAND_MAX;

		vdg.template getProp<0>(p) = (double)rand() / RAND_MAX;

		vdg.template getProp<1>(p)[0] = (double)rand() / RAND_MAX;
		vdg.template getProp<1>(p)[1] = (double)rand() / RAND_MAX;
		vdg.template getProp<1>(p)[2] = (double)rand() / RAND_MAX;

		vdg.template getProp<2>(p)[0][0] = (double)rand() / RAND_MAX;
		vdg.template getProp<2>(p)[0][1] = (double)rand() / RAND_MAX;
		vdg.template getProp<2>(p)[0][2] = (double)rand() / RAND_MAX;
		vdg.template getProp<2>(p)[1][0] = (double)rand() / RAND_MAX;
		vdg.template getProp<2>(p)[1][1] = (double)rand() / RAND_MAX;
		vdg.template getProp<2>(p)[1][2] = (double)rand() / RAND_MAX;
		vdg.template getProp<2>(p)[2][0] = (double)rand() / RAND_MAX;
		vdg.template getProp<2>(p)[2][1] = (double)rand() / RAND_MAX;
		vdg.template getProp<2>(p)[2][2] = (double)rand() / RAND_MAX;

		++it;
	}

	vdg.map();

	vdg.hostToDeviceProp<0,1,2>();
	vdg.hostToDevicePos();

	bool test = vdg.compareHostAndDevicePos(0.00001,0.00000001);
	BOOST_REQUIRE_EQUAL(test,true);

	vdg.getPos(100)[0] = 0.99999999;

	test = vdg.compareHostAndDevicePos(0.00001,0.00000001);
	BOOST_REQUIRE_EQUAL(test,false);

	vdg.hostToDevicePos();
	vdg.getPos(100)[0] = 0.99999999;

	test = vdg.compareHostAndDevicePos(0.00001,0.00000001);
	BOOST_REQUIRE_EQUAL(test,true);

	////////////////////////////////////////////////// PROP VECTOR

	test = vdg.compareHostAndDeviceProp<1>(0.00001,0.00000001);
	BOOST_REQUIRE_EQUAL(test,true);

	vdg.getProp<1>(103)[0] = 0.99999999;

	test = vdg.compareHostAndDeviceProp<1>(0.00001,0.00000001);
	BOOST_REQUIRE_EQUAL(test,false);

	vdg.hostToDeviceProp<1>();
	vdg.getProp<1>(103)[0] = 0.99999999;

	test = vdg.compareHostAndDeviceProp<1>(0.00001,0.00000001);
	BOOST_REQUIRE_EQUAL(test,true);

	////////////////////////////////////////////////// PROP scalar


	test = vdg.compareHostAndDeviceProp<0>(0.00001,0.00000001);
	BOOST_REQUIRE_EQUAL(test,true);

	vdg.getProp<0>(105) = 0.99999999;

	test = vdg.compareHostAndDeviceProp<0>(0.00001,0.00000001);
	BOOST_REQUIRE_EQUAL(test,false);

	vdg.hostToDeviceProp<0>();
	vdg.getProp<0>(105) = 0.99999999;

	test = vdg.compareHostAndDeviceProp<0>(0.00001,0.00000001);
	BOOST_REQUIRE_EQUAL(test,true);


	////////////////////////////////////////////////// PROP scalar


	test = vdg.compareHostAndDeviceProp<2>(0.00001,0.00000001);
	BOOST_REQUIRE_EQUAL(test,true);

	vdg.getProp<2>(108)[1][2] = 0.99999999;

	test = vdg.compareHostAndDeviceProp<2>(0.00001,0.00000001);
	BOOST_REQUIRE_EQUAL(test,false);

	vdg.hostToDeviceProp<2>();
	vdg.getProp<2>(108)[1][2] = 0.99999999;

	test = vdg.compareHostAndDeviceProp<2>(0.00001,0.00000001);
	BOOST_REQUIRE_EQUAL(test,true);
}

template<typename vector_dist_type>
__global__ void assign_to_ghost(vector_dist_type vds)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;

	if (i >= vds.size())	{return;}

	vds.template getProp<0>(i) = 1000.0 + i;

	vds.template getProp<1>(i)[0] = 2000.0 + i;
	vds.template getProp<1>(i)[1] = 3000.0 + i;
	vds.template getProp<1>(i)[2] = 4000.0 + i;

	vds.template getProp<2>(i)[0][0] = 12000.0 + i;
	vds.template getProp<2>(i)[0][1] = 13000.0 + i;
	vds.template getProp<2>(i)[0][2] = 14000.0 + i;
	vds.template getProp<2>(i)[1][0] = 22000.0 + i;
	vds.template getProp<2>(i)[1][1] = 23000.0 + i;
	vds.template getProp<2>(i)[1][2] = 24000.0 + i;
	vds.template getProp<2>(i)[2][0] = 32000.0 + i;
	vds.template getProp<2>(i)[2][1] = 33000.0 + i;
	vds.template getProp<2>(i)[2][2] = 34000.0 + i;

}

BOOST_AUTO_TEST_CASE(vector_dist_domain_and_ghost_test)
{
        Box<3,double> domain({0.0,0.0,0.0},{1.0,1.0,1.0});
        Ghost<3,double> g(0.1);
        size_t bc[3] = {PERIODIC,PERIODIC,PERIODIC};

        if (create_vcluster().size() >= 16)
        {return;}

        vector_dist_gpu<3,double,aggregate<double,double[3],double[3][3]>> vdg(10000,domain,bc,g,DEC_GRAN(128));

        auto ite = vdg.getDomainAndGhostIteratorGPU();

        CUDA_LAUNCH(assign_to_ghost,ite,vdg.toKernel());

        vdg.template deviceToHostProp<0,1,2>();


        auto it = vdg.getDomainAndGhostIterator();

        bool check = true;

        while (it.isNext())
        {
                auto k = it.get();

                check &= vdg.template getProp<0>(k) == 1000.0 + k.getKey();

                check &= vdg.template getProp<1>(k)[0] == 2000.0 + k.getKey();
                check &= vdg.template getProp<1>(k)[1] == 3000.0 + k.getKey();
                check &= vdg.template getProp<1>(k)[2] == 4000.0 + k.getKey();

                check &= vdg.template getProp<2>(k)[0][0] == 12000.0 + k.getKey();
                check &= vdg.template getProp<2>(k)[0][1] == 13000.0 + k.getKey();
                check &= vdg.template getProp<2>(k)[0][2] == 14000.0 + k.getKey();
                check &= vdg.template getProp<2>(k)[1][0] == 22000.0 + k.getKey();
                check &= vdg.template getProp<2>(k)[1][1] == 23000.0 + k.getKey();
                check &= vdg.template getProp<2>(k)[1][2] == 24000.0 + k.getKey();
                check &= vdg.template getProp<2>(k)[2][0] == 32000.0 + k.getKey();
                check &= vdg.template getProp<2>(k)[2][1] == 33000.0 + k.getKey();
                check &= vdg.template getProp<2>(k)[2][2] == 34000.0 + k.getKey();

                ++it;
        }


        BOOST_REQUIRE_EQUAL(check,true);
}

template<typename vT>
__global__ void launch_overflow(vT vs, vT vs2)
{
	vs2.template getProp<1>(57)[0];
}

BOOST_AUTO_TEST_CASE(vector_dist_overflow_se_class1)
{
	Box<3,double> domain({0.0,0.0,0.0},{1.0,1.0,1.0});
	Ghost<3,double> g(0.1);
	size_t bc[3] = {PERIODIC,PERIODIC,PERIODIC};

	if (create_vcluster().size() >= 16)
	{return;}

	std::cout << "****** TEST ERROR MESSAGE BEGIN ********" << std::endl;

	vector_dist_gpu<3,double,aggregate<double,double[3],double[3][3]>> vdg(0,domain,bc,g,DEC_GRAN(128));
	vector_dist_gpu<3,double,aggregate<double,double[3],double[3][3]>> vdg2(0,domain,bc,g,DEC_GRAN(128));


	vdg.setCapacity(100);

	ite_gpu<1> ite;

	ite.wthr.x = 1;
	ite.wthr.y = 1;
	ite.wthr.z = 1;
	ite.thr.x = 1;
	ite.thr.y = 1;
	ite.thr.z = 1;

	try
	{
		CUDA_LAUNCH(launch_overflow,ite,vdg.toKernel(),vdg2.toKernel());
	}
	catch(...)
	{
		std::cout << "SE_CLASS1 Catch" << std::endl;
	};

	std::cout << "****** TEST ERROR MESSAGE END ********" << std::endl;
}



BOOST_AUTO_TEST_CASE( vector_dist_ghost_put_gpu )
{

// This example require float atomic, until C++20 we are out of luck
#ifndef CUDIFY_USE_OPENMP

	Vcluster<> & v_cl = create_vcluster();

	long int k = 25*25*25*create_vcluster().getProcessingUnits();
	k = std::pow(k, 1/3.);

	if (v_cl.getProcessingUnits() > 48)
		return;

	print_test("Testing 3D periodic ghost put GPU k=",k);
	BOOST_TEST_CHECKPOINT( "Testing 3D periodic ghost put k=" << k );

	long int big_step = k / 30;
	big_step = (big_step == 0)?1:big_step;
	long int small_step = 21;

	// 3D test
	for ( ; k >= 2 ; k-= (k > 2*big_step)?big_step:small_step )
	{
		float r_cut = 1.3 / k;
		float r_g = 1.5 / k;

		Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

		// Boundary conditions
		size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

		// ghost
		Ghost<3,float> ghost(r_g);

		typedef  aggregate<float,float,float> part_prop;

		// Distributed vector
		vector_dist_gpu<3,float, part_prop > vd(0,box,bc,ghost);

		auto it = vd.getGridIterator({(size_t)k,(size_t)k,(size_t)k});

		while (it.isNext())
		{
			auto key = it.get();

			vd.add();

			vd.getLastPosWrite()[0] = key.get(0)*it.getSpacing(0);
			vd.getLastPosWrite()[1] = key.get(1)*it.getSpacing(1);
			vd.getLastPosWrite()[2] = key.get(2)*it.getSpacing(2);

			// Fill some properties randomly

			vd.getLastPropWrite<0>() = 0.0;

			vd.getLastPropWrite<2>() = 0.0;

			++it;
		}

		vd.map();

		vd.hostToDevicePos();
		vd.template hostToDeviceProp<0,2>();
		// sync the ghost
		vd.ghost_get<0,2>(RUN_ON_DEVICE);
		vd.template deviceToHostProp<0,2>();
		vd.deviceToHostPos();

		{
			auto NN = vd.getCellList(r_cut);
			float a = 1.0f*k*k;

			// run trough all the particles + ghost

			auto it2 = vd.getDomainIterator();

			while (it2.isNext())
			{
				// particle p
				auto p = it2.get();
				Point<3,float> xp = vd.getPos(p);

				// Get an iterator over the neighborhood particles of p
				auto Np = NN.getNNIterator<NO_CHECK>(NN.getCell(xp));

				// For each neighborhood particle ...
				while (Np.isNext())
				{
					auto q = Np.get();
					Point<3,float> xq = vd.getPosRead(q);

					float dist = xp.distance(xq);

					if (dist < r_cut)
					{
						vd.getPropWrite<0>(q) += a*(-dist*dist+r_cut*r_cut);
						vd.getPropWrite<2>(q) += a*(-dist*dist+r_cut*r_cut) / 2;
					}

					++Np;
				}

				++it2;
			}

			vd.hostToDevicePos();
			vd.template hostToDeviceProp<0,2>();
			vd.template ghost_put<add_atomic_,0,2>(RUN_ON_DEVICE);
			vd.template deviceToHostProp<0,2>();
			vd.deviceToHostPos();

			bool ret = true;
			auto it3 = vd.getDomainIterator();

			float constant = vd.getProp<0>(it3.get());
			float constanta = vd.getProp<2>(it3.get());
			float eps = 0.001;

			while (it3.isNext())
			{
				float constant2 = vd.getProp<0>(it3.get());
				float constant3 = vd.getProp<2>(it3.get());
				if (fabs(constant - constant2)/constant > eps || fabs(constanta - constant3)/constanta > eps)
				{
					Point<3,float> p = vd.getPosRead(it3.get());

					std::cout << p.toString() << "    " <<  constant2 << "/" << constant << "/" << constant3 << "    " << v_cl.getProcessUnitID() << std::endl;
					ret = false;
					break;
				}

				++it3;
			}
			BOOST_REQUIRE_EQUAL(ret,true);
		}

		auto itp = vd.getDomainAndGhostIterator();
		while (itp.isNext())
		{
			auto key = itp.get();

			vd.getPropWrite<0>(key) = 0.0;
			vd.getPropWrite<2>(key) = 0.0;

			++itp;
		}

		{
			auto NN = vd.getCellList(r_cut);
			float a = 1.0f*k*k;

			// run trough all the particles + ghost

			auto it2 = vd.getDomainIterator();

			while (it2.isNext())
			{
				// particle p
				auto p = it2.get();
				Point<3,float> xp = vd.getPosRead(p);

				// Get an iterator over the neighborhood particles of p
				auto Np = NN.getNNIterator<NO_CHECK>(NN.getCell(xp));

				// For each neighborhood particle ...
				while (Np.isNext())
				{
					auto q = Np.get();
					Point<3,float> xq = vd.getPosRead(q);

					float dist = xp.distance(xq);

					if (dist < r_cut)
					{
						vd.getPropWrite<0>(q) += a*(-dist*dist+r_cut*r_cut);
						vd.getPropWrite<2>(q) += a*(-dist*dist+r_cut*r_cut);
					}

					++Np;
				}

				++it2;
			}

			vd.hostToDevicePos();
			vd.template hostToDeviceProp<0,2>();
			vd.template ghost_put<add_atomic_,0>(RUN_ON_DEVICE);
			vd.template ghost_put<add_atomic_,2>(RUN_ON_DEVICE);
			vd.template deviceToHostProp<0,2>();
			vd.deviceToHostPos();

			bool ret = true;
			auto it3 = vd.getDomainIterator();

			float constant = vd.getPropRead<0>(it3.get());
			float constanta = vd.getPropRead<2>(it3.get());
			float eps = 0.001;

			while (it3.isNext())
			{
				float constant2 = vd.getPropRead<0>(it3.get());
				float constant3 = vd.getPropRead<0>(it3.get());
				if (fabs(constant - constant2)/constant > eps || fabs(constanta - constant3)/constanta > eps)
				{
					Point<3,float> p = vd.getPosRead(it3.get());

					std::cout << p.toString() << "    " <<  constant2 << "/" << constant << "/" << constant3 << "    " << v_cl.getProcessUnitID() << std::endl;
					ret = false;
					break;
				}

				++it3;
			}
			BOOST_REQUIRE_EQUAL(ret,true);
		}
	}

#endif
}

BOOST_AUTO_TEST_SUITE_END()
