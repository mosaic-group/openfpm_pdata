#define BOOST_TEST_DYN_LINK
#include "config.h"
#include <boost/test/unit_test.hpp>
#include "VCluster/VCluster.hpp"
#include <Vector/vector_dist.hpp>
#include "Vector/tests/vector_dist_util_unit_tests.hpp"

#define SUB_UNIT_FACTOR 1024

template<unsigned int dim , typename vector_dist_type>
__global__ void move_parts_gpu_test(vector_dist_type vecDist)
{
	auto p = GET_PARTICLE(vecDist);

#pragma unroll
	for (int i = 0 ; i < dim ; i++)
	{
		vecDist.getPos(p)[i] += 0.05;
	}
}

BOOST_AUTO_TEST_SUITE( vector_dist_gpu_test )

void print_test(std::string test, size_t sz)
{
	if (create_vcluster().getProcessUnitID() == 0)
		std::cout << test << " " << sz << "\n";
}


__global__  void initialize_props(vector_dist_ker<3, float, aggregate<float, float [3], float[3]>> vecDist)
{
	auto p = GET_PARTICLE(vecDist);

	vecDist.template getProp<0>(p) = vecDist.getPos(p)[0] + vecDist.getPos(p)[1] + vecDist.getPos(p)[2];

	vecDist.template getProp<1>(p)[0] = vecDist.getPos(p)[0] + vecDist.getPos(p)[1];
	vecDist.template getProp<1>(p)[1] = vecDist.getPos(p)[0] + vecDist.getPos(p)[2];
	vecDist.template getProp<1>(p)[2] = vecDist.getPos(p)[1] + vecDist.getPos(p)[2];

	vecDist.template getProp<2>(p)[0] = vecDist.getPos(p)[0] + vecDist.getPos(p)[1];
	vecDist.template getProp<2>(p)[1] = vecDist.getPos(p)[0] + vecDist.getPos(p)[2];
	vecDist.template getProp<2>(p)[2] = vecDist.getPos(p)[1] + vecDist.getPos(p)[2];
}

template<typename T,typename CellList_type>
__global__  void calculate_force(
	vector_dist_ker<3, T, aggregate<T, T[3], T [3]>> vecDist,
	CellList_type cellList,
	int rank)
{
	size_t p = GET_PARTICLE(vecDist);

	Point<3,T> xp = vecDist.getPos(p);

	auto it = cellList.getNNIterator(cellList.getCell(xp));

	Point<3,T> force({0.0,0.0,0.0});

	while (it.isNext())
	{
		auto q = it.get();

		if (q == p) {++it; continue;}

		Point<3,T> xq = vecDist.getPos(q);
		Point<3,T> r = xq - xp;

		if (r.norm() > 1e-6)
		{
			r /= r.norm();
			force += vecDist.template getProp<0>(q)*r;
		}

		++it;
	}

	vecDist.template getProp<1>(p)[0] = force.get(0);
	vecDist.template getProp<1>(p)[1] = force.get(1);
	vecDist.template getProp<1>(p)[2] = force.get(2);
}

template<typename T,typename CellList_type>
__global__  void calculate_force_sort(
	vector_dist_ker<3, T, aggregate<T, T[3], T [3]>> vecDistSort,
	CellList_type cellList,
	int rank)
{
	size_t p; GET_PARTICLE_SORT(p, cellList);

	Point<3,T> xp = vecDistSort.getPos(p);
	Point<3,T> force({0.0,0.0,0.0});

	auto it = cellList.getNNIterator(cellList.getCell(xp));

	while (it.isNext())
	{
		auto q = it.get_sort();

		if (q == p) {++it; continue;}

		Point<3,T> xq = vecDistSort.getPos(q);
		Point<3,T> r = xq - xp;

		// Normalize

		if (r.norm() > 1e-6)
		{
			r /= r.norm();
			force += vecDistSort.template getProp<0>(q)*r;
		}

		++it;
	}

	vecDistSort.template getProp<1>(p)[0] = force.get(0);
	vecDistSort.template getProp<1>(p)[1] = force.get(1);
	vecDistSort.template getProp<1>(p)[2] = force.get(2);
}

template<typename CellList_type, typename vector_type>
bool check_force(CellList_type & cellList, vector_type & vecDist)
{
	typedef typename vector_type::stype St;

	auto it6 = vecDist.getDomainIterator();

	bool match = true;

	while (it6.isNext())
	{
		auto p = it6.get();

		Point<3,St> xp = vecDist.getPos(p);

		// Calculate on CPU

		Point<3,St> force({0.0,0.0,0.0});

		auto NNc = cellList.getNNIterator(cellList.getCell(xp));

		while (NNc.isNext())
		{
			auto q = NNc.get();

			if (q == p.getKey()) {++NNc; continue;}

			Point<3,St> xq_2 = vecDist.getPos(q);
			Point<3,St> r2 = xq_2 - xp;

			// Normalize

			if (r2.norm() > 1e-6)
			{
				r2 /= r2.norm();
				force += vecDist.template getProp<0>(q)*r2;
			}

			++NNc;
		}

		match &= fabs(vecDist.template getProp<1>(p)[0] - force.get(0)) < 0.0003;
		match &= fabs(vecDist.template getProp<1>(p)[1] - force.get(1)) < 0.0003;
		match &= fabs(vecDist.template getProp<1>(p)[2] - force.get(2)) < 0.0003;

		if (match == false)
		{
			std::cout << p.getKey() << " ERROR: " << vecDist.template getProp<1>(p)[0] << "   " <<  force.get(0) << std::endl;
			// std::cout << p.getKey() << " ERROR: " << vecDist.template getProp<1>(p)[1] << "   " <<  force.get(1) << std::endl;
			// std::cout << p.getKey() << " ERROR: " << vecDist.template getProp<1>(p)[2] << "   " <<  force.get(2) << std::endl;

			// break;
		}

		++it6;
	}

	return match;
}

BOOST_AUTO_TEST_CASE( vector_dist_gpu_ghost_get )
{
	auto & vCluster = create_vcluster();

	if (vCluster.size() > 16)
	{return;}

	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	// set the ghost based on the radius cut off (make just a little bit smaller than the spacing)
	Ghost<3,float> g(0.1);

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	vector_dist_gpu<3,float,aggregate<float,float[3],float[3]>> vecDist(1000,domain,bc,g);

	auto it = vecDist.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		vecDist.getPos(p)[0] = (float)rand() / (float)RAND_MAX;
		vecDist.getPos(p)[1] = (float)rand() / (float)RAND_MAX;
		vecDist.getPos(p)[2] = (float)rand() / (float)RAND_MAX;

		vecDist.template getProp<0>(p) = vecDist.getPos(p)[0] + vecDist.getPos(p)[1] + vecDist.getPos(p)[2];

		vecDist.template getProp<1>(p)[0] = vecDist.getPos(p)[0] + vecDist.getPos(p)[1];
		vecDist.template getProp<1>(p)[1] = vecDist.getPos(p)[0] + vecDist.getPos(p)[2];
		vecDist.template getProp<1>(p)[2] = vecDist.getPos(p)[1] + vecDist.getPos(p)[2];

		vecDist.template getProp<2>(p)[0] = vecDist.getPos(p)[0] + 3.0*vecDist.getPos(p)[1];
		vecDist.template getProp<2>(p)[1] = vecDist.getPos(p)[0] + 3.0*vecDist.getPos(p)[2];
		vecDist.template getProp<2>(p)[2] = vecDist.getPos(p)[1] + 3.0*vecDist.getPos(p)[2];


		++it;
	}

	// Ok we redistribute the particles (CPU based)
	vecDist.map();

	vecDist.template ghost_get<0,1,2>();

	// Now we check the the ghost contain the correct information

	bool check = true;

	auto itg = vecDist.getDomainAndGhostIterator();

	while (itg.isNext())
	{
		auto p = itg.get();

		check &= (vecDist.template getProp<0>(p) == vecDist.getPos(p)[0] + vecDist.getPos(p)[1] + vecDist.getPos(p)[2]);

		check &= (vecDist.template getProp<1>(p)[0] == vecDist.getPos(p)[0] + vecDist.getPos(p)[1]);
		check &= (vecDist.template getProp<1>(p)[1] == vecDist.getPos(p)[0] + vecDist.getPos(p)[2]);
		check &= (vecDist.template getProp<1>(p)[2] == vecDist.getPos(p)[1] + vecDist.getPos(p)[2]);

		check &= (vecDist.template getProp<2>(p)[0] == vecDist.getPos(p)[0] + 3.0*vecDist.getPos(p)[1]);
		check &= (vecDist.template getProp<2>(p)[1] == vecDist.getPos(p)[0] + 3.0*vecDist.getPos(p)[2]);
		check &= (vecDist.template getProp<2>(p)[2] == vecDist.getPos(p)[1] + 3.0*vecDist.getPos(p)[2]);

		++itg;
	}

	size_t tot_s = vecDist.size_local_with_ghost();

	vCluster.sum(tot_s);
	vCluster.execute();

	// We check that we check something
	BOOST_REQUIRE(tot_s > 1000);
}

template<typename vector_type, typename CellList_type, typename CellList_type_cpu>
void compareCellListCpuGpu(vector_type & vecDist, CellList_type & NN, CellList_type_cpu & cellList)
{
	const auto it5 = vecDist.getDomainIteratorGPU(32);

	CUDA_LAUNCH((calculate_force<typename vector_type::stype,decltype(NN.toKernel())>),
		it5,
		vecDist.toKernel(),
		NN.toKernel(),
		(int)create_vcluster().rank()
	);

	vecDist.template deviceToHostProp<1>();

	bool test = check_force(cellList,vecDist);
	BOOST_REQUIRE_EQUAL(test,true);
}

template<typename vector_type, typename CellList_type, typename CellList_type_cpu>
void compareCellListCpuGpuSorted(vector_type & vecDistSort, CellList_type & cellListGPU, CellList_type_cpu & cellList)
{
	const auto it5 = vecDistSort.getDomainIteratorGPU(32);

	CUDA_LAUNCH((calculate_force_sort<typename vector_type::stype,decltype(cellListGPU.toKernel())>),
		it5,
		vecDistSort.toKernel(),
		cellListGPU.toKernel(),
		(int)create_vcluster().rank()
	);

	vecDistSort.template restoreOrder<1>(cellListGPU);

	vecDistSort.template deviceToHostProp<1>();
	vecDistSort.deviceToHostPos();

	bool test = check_force(cellList,vecDistSort);
	BOOST_REQUIRE_EQUAL(test,true);
}

template<typename CellList_type, bool sorted>
void vector_dist_gpu_test_impl()
{
	auto & vCluster = create_vcluster();

	if (vCluster.size() > 16)
	{return;}

	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	// set the ghost based on the radius cut off (make just a little bit smaller than the spacing)
	Ghost<3,float> g(0.1);

	// Boundary conditions
	size_t bc[3]={NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

	vector_dist_gpu<3,float,aggregate<float,float[3],float[3]>> vecDist(10000,domain,bc,g);

	srand(55067*create_vcluster().rank());

	auto it = vecDist.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		int x = rand();
		int y = rand();
		int z = rand();

		vecDist.getPos(p)[0] = (float)x / (float)RAND_MAX;
		vecDist.getPos(p)[1] = (float)y / (float)RAND_MAX;
		vecDist.getPos(p)[2] = (float)z / (float)RAND_MAX;

		Point<3,float> xp = vecDist.getPos(p);

		++it;
	}

	// Ok we redistribute the particles (CPU based)
	vecDist.map();

	size_t size_l = vecDist.size_local();

	vCluster.sum(size_l);
	vCluster.execute();

	BOOST_REQUIRE_EQUAL(size_l,10000);

	auto & dec = vecDist.getDecomposition();

	bool noOut = true;
	size_t cnt = 0;

	auto it2 = vecDist.getDomainIterator();

	while (it2.isNext())
	{
		auto p = it2.get();

		noOut &= dec.isLocal(vecDist.getPos(p));

		cnt++;
		++it2;
	}

	BOOST_REQUIRE_EQUAL(noOut,true);
	BOOST_REQUIRE_EQUAL(cnt,vecDist.size_local());

	// now we offload all the properties

	const auto it3 = vecDist.getDomainIteratorGPU();

	// offload to device
	vecDist.hostToDevicePos();

	CUDA_LAUNCH_DIM3(initialize_props,it3.wthr,it3.thr,vecDist.toKernel());

	// now we check what we initialized

	vecDist.deviceToHostProp<0,1>();

	auto it4 = vecDist.getDomainIterator();

	while (it4.isNext())
	{
		auto p = it4.get();

		BOOST_REQUIRE_CLOSE(vecDist.template getProp<0>(p),vecDist.getPos(p)[0] + vecDist.getPos(p)[1] + vecDist.getPos(p)[2],0.01);

		BOOST_REQUIRE_CLOSE(vecDist.template getProp<1>(p)[0],vecDist.getPos(p)[0] + vecDist.getPos(p)[1],0.01);
		BOOST_REQUIRE_CLOSE(vecDist.template getProp<1>(p)[1],vecDist.getPos(p)[0] + vecDist.getPos(p)[2],0.01);
		BOOST_REQUIRE_CLOSE(vecDist.template getProp<1>(p)[2],vecDist.getPos(p)[1] + vecDist.getPos(p)[2],0.01);

		++it4;
	}

	// here we do a ghost_get
	vecDist.ghost_get<0>();

	// Double ghost get to check crashes
	vecDist.ghost_get<0>();

	// we re-offload what we received
	vecDist.hostToDevicePos();
	vecDist.template hostToDeviceProp<0>();

	size_t opt = (sorted)? CL_GPU_REORDER : 0;
	auto cellList = vecDist.getCellList(0.1);
	auto cellListGPU = vecDist.template getCellListGPU<CellList_type>(0.1, opt | CL_NON_SYMMETRIC);

	if (sorted)
	{
		vecDist.hostToDevicePos();
		vecDist.template hostToDeviceProp<0>();

		vecDist.template updateCellListGPU<0>(cellListGPU);
		compareCellListCpuGpuSorted(vecDist,cellListGPU,cellList);
	}

	else
	{
		vecDist.updateCellListGPU(cellListGPU);
		compareCellListCpuGpu(vecDist,cellListGPU,cellList);
	}
}

template<typename CellList_type>
void vector_dist_gpu_make_sort_test_impl()
{
	auto & vCluster = create_vcluster();

	if (vCluster.size() > 16)
	{return;}

	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	// set the ghost based on the radius cut off (make just a little bit smaller than the spacing)
	Ghost<3,float> g(0.1);

	// Boundary conditions
	size_t bc[3]={NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

	vector_dist_gpu<3,float,aggregate<float,float[3],float[3]>> vecDist(10000,domain,bc,g);

	srand(55067*create_vcluster().rank());

	auto it = vecDist.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		int x = rand();
		int y = rand();
		int z = rand();

		vecDist.getPos(p)[0] = (float)x / (float)RAND_MAX;
		vecDist.getPos(p)[1] = (float)y / (float)RAND_MAX;
		vecDist.getPos(p)[2] = (float)z / (float)RAND_MAX;

		++it;
	}

	vecDist.hostToDevicePos();

	// redistribute the particles
	vecDist.map(RUN_ON_DEVICE);

	auto it3 = vecDist.getDomainIteratorGPU();

	CUDA_LAUNCH_DIM3(initialize_props,it3.wthr,it3.thr,vecDist.toKernel());

	vecDist.template deviceToHostProp<0,1,2>();

	// Here we check CL_GPU_REORDER correctly restores particle order
	// we use a Cell-List to check that the two cell-list constructed are identical

	vecDist.deviceToHostPos();

	openfpm::vector_gpu<Point<3,float>> tmpPos = vecDist.getPosVector();
	openfpm::vector_gpu<aggregate<float,float[3],float[3]>> tmpPrp = vecDist.getPropVector();

	auto NN_cpu1 = vecDist.getCellList(0.1);
	auto NN = vecDist.template getCellListGPU<CellList_type>(0.1, CL_NON_SYMMETRIC | CL_GPU_REORDER);

	vecDist.hostToDevicePos();
	vecDist.template hostToDeviceProp<0,1,2>();

	vecDist.template updateCellListGPU<0,1,2>(NN);
	vecDist.template restoreOrder<0,1,2>(NN);

	vecDist.deviceToHostPos();
	vecDist.template deviceToHostProp<0,1,2>();

	vecDist.hostToDevicePos();
	vecDist.template hostToDeviceProp<0>();

	vecDist.template updateCellListGPU<0>(NN);
	vecDist.template restoreOrder<0>(NN);

	vecDist.deviceToHostPos();
	vecDist.template deviceToHostProp<0>();

	vecDist.hostToDevicePos();
	vecDist.template hostToDeviceProp<1>();

	vecDist.template updateCellListGPU<1>(NN);
	vecDist.template restoreOrder<1>(NN);

	vecDist.deviceToHostPos();
	vecDist.template deviceToHostProp<1>();

	vecDist.hostToDevicePos();
	vecDist.template hostToDeviceProp<2>();

	vecDist.template updateCellListGPU<2>(NN);
	vecDist.template restoreOrder<2>(NN);

	vecDist.deviceToHostPos();
	vecDist.template deviceToHostProp<2>();

	auto NN_cpu2 = vecDist.getCellList(0.1);

	// compare two cell-lists
	bool match = true;
	for (size_t i = 0 ; i < NN_cpu1.getNCells() ; i++)
	{
		match &= NN_cpu1.getNelements(i) == NN_cpu2.getNelements(i);
	}
	BOOST_REQUIRE_EQUAL(match,true);

	match = true;
	for (size_t i = 0 ; i < vecDist.size_local() ; i++)
	{
		Point<3,float> p1 = vecDist.getPos(i);
		Point<3,float> p2 = tmpPos.template get<0>(i);

		// They must be in the same cell
		auto c1 = NN.getCell(p1);
		auto c2 = NN.getCell(p2);

		match &= c1 == c2;
	}
	BOOST_REQUIRE_EQUAL(match,true);

	match = true;
	for (size_t i = 0 ; i < vecDist.size_local() ; i++)
	{
		for (int j = 0; j < 3; ++j)
			match &= (vecDist.getPos(i)[j] - tmpPos.template get<0>(i)[j]) < 0.0003;

		match &= (vecDist.getProp<0>(i) - tmpPrp.template get<0>(i)) < 0.0003;

		for (int j = 0; j < 3; ++j)
			match &= (vecDist.getProp<1>(i)[j] - tmpPrp.template get<1>(i)[j]) < 0.0003;

		for (int j = 0; j < 3; ++j)
			match &= (vecDist.getProp<2>(i)[j] - tmpPrp.template get<2>(i)[j]) < 0.0003;
	}

	BOOST_REQUIRE_EQUAL(match,true);

	CUDA_CHECK();
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
	vector_dist_gpu_test_impl<CellList_gpu<3,float,CudaMemory,shift_only<3, float>>, false>();
}

BOOST_AUTO_TEST_CASE( vector_dist_gpu_test_sorted)
{
	vector_dist_gpu_test_impl<CellList_gpu<3,float,CudaMemory,shift_only<3, float>>, true>();
}

BOOST_AUTO_TEST_CASE( vector_dist_gpu_test_sparse)
{
	vector_dist_gpu_test_impl<CELLLIST_GPU_SPARSE<3,float>, false>();
}

BOOST_AUTO_TEST_CASE( vector_dist_gpu_test_sparse_sorted)
{
	vector_dist_gpu_test_impl<CELLLIST_GPU_SPARSE<3,float>, true>();
}

template<typename St>
void vdist_calc_gpu_test()
{
	auto & vCluster = create_vcluster();

	if (vCluster.size() > 16)
	{return;}

	Box<3,St> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	// set the ghost based on the radius cut off (make just a little bit smaller than the spacing)
	Ghost<3,St> g(0.1);

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	//! [Create a gpu vector]

	vector_dist_gpu<3,St,aggregate<St,St[3],St[3]>> vecDist(1000,domain,bc,g);

	//! [Create a gpu vector]

	//! [Fill gpu vector and move to GPU]

	srand(vCluster.rank()*10000);
	auto it = vecDist.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		vecDist.getPos(p)[0] = (St)rand() / (float)RAND_MAX;
		vecDist.getPos(p)[1] = (St)rand() / (float)RAND_MAX;
		vecDist.getPos(p)[2] = (St)rand() / (float)RAND_MAX;

		vecDist.template getProp<0>(p) = vecDist.getPos(p)[0] + vecDist.getPos(p)[1] + vecDist.getPos(p)[2];

		vecDist.template getProp<1>(p)[0] = vecDist.getPos(p)[0];
		vecDist.template getProp<1>(p)[1] = vecDist.getPos(p)[1];
		vecDist.template getProp<1>(p)[2] = vecDist.getPos(p)[2];

		vecDist.template getProp<2>(p)[0] = vecDist.getPos(p)[0] + vecDist.getPos(p)[1];
		vecDist.template getProp<2>(p)[1] = vecDist.getPos(p)[0] + vecDist.getPos(p)[2];
		vecDist.template getProp<2>(p)[2] = vecDist.getPos(p)[1] + vecDist.getPos(p)[2];

		++it;
	}

	// move on device
	vecDist.hostToDevicePos();
	vecDist.template hostToDeviceProp<0,1,2>();

	// Ok we redistribute the particles (GPU based)
	vecDist.map(RUN_ON_DEVICE);

	//! [Fill gpu vector and move to GPU]

	vecDist.deviceToHostPos();
	vecDist.template deviceToHostProp<0,1,2>();

	// Reset the host part

	auto it3 = vecDist.getDomainIterator();

	while (it3.isNext())
	{
		auto p = it3.get();

		vecDist.getPos(p)[0] = 1.0;
		vecDist.getPos(p)[1] = 1.0;
		vecDist.getPos(p)[2] = 1.0;

		vecDist.template getProp<0>(p) = 0.0;

		vecDist.template getProp<0>(p) = 0.0;
		vecDist.template getProp<0>(p) = 0.0;
		vecDist.template getProp<0>(p) = 0.0;

		vecDist.template getProp<0>(p) = 0.0;
		vecDist.template getProp<0>(p) = 0.0;
		vecDist.template getProp<0>(p) = 0.0;

		++it3;
	}

	// we move from Device to CPU

	vecDist.deviceToHostPos();
	vecDist.template deviceToHostProp<0,1,2>();

	// Check

	auto it2 = vecDist.getDomainIterator();

	bool match = true;
	while (it2.isNext())
	{
		auto p = it2.get();

		match &= vecDist.template getProp<0>(p) == vecDist.getPos(p)[0] + vecDist.getPos(p)[1] + vecDist.getPos(p)[2];

		match &= vecDist.template getProp<1>(p)[0] == vecDist.getPos(p)[0];
		match &= vecDist.template getProp<1>(p)[1] == vecDist.getPos(p)[1];
		match &= vecDist.template getProp<1>(p)[2] == vecDist.getPos(p)[2];

		match &= vecDist.template getProp<2>(p)[0] == vecDist.getPos(p)[0] + vecDist.getPos(p)[1];
		match &= vecDist.template getProp<2>(p)[1] == vecDist.getPos(p)[0] + vecDist.getPos(p)[2];
		match &= vecDist.template getProp<2>(p)[2] == vecDist.getPos(p)[1] + vecDist.getPos(p)[2];

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

	auto it5 = vecDist.getDomainIterator();
	count_local_n_local<3>(vecDist,it5,bc,domain,dom_ext,l_cnt,nl_cnt,n_out);

	BOOST_REQUIRE_EQUAL(n_out,0);
	BOOST_REQUIRE_EQUAL(l_cnt,vecDist.size_local());

	// we do 10 gpu steps (using a cpu vector to check that map and ghost get work as expented)

	for (size_t i = 0 ; i < 10 ; i++)
	{
		vecDist.map(RUN_ON_DEVICE);

		vecDist.deviceToHostPos();
		vecDist.template deviceToHostProp<0,1,2>();

		// To test we copy on a cpu distributed vector and we do a map

		vector_dist<3,St,aggregate<St,St[3],St[3]>> vd_cpu(vecDist.getDecomposition().template duplicate_convert<HeapMemory,memory_traits_lin>(),0);

		auto itc = vecDist.getDomainIterator();

		while (itc.isNext())
		{
			auto p = itc.get();

			vd_cpu.add();

			vd_cpu.getLastPos()[0] = vecDist.getPos(p)[0];
			vd_cpu.getLastPos()[1] = vecDist.getPos(p)[1];
			vd_cpu.getLastPos()[2] = vecDist.getPos(p)[2];

			vd_cpu.template getLastProp<0>() = vecDist.template getProp<0>(p);

			vd_cpu.template getLastProp<1>()[0] = vecDist.template getProp<1>(p)[0];
			vd_cpu.template getLastProp<1>()[1] = vecDist.template getProp<1>(p)[1];
			vd_cpu.template getLastProp<1>()[2] = vecDist.template getProp<1>(p)[2];

			vd_cpu.template getLastProp<2>()[0] = vecDist.template getProp<2>(p)[0];
			vd_cpu.template getLastProp<2>()[1] = vecDist.template getProp<2>(p)[1];
			vd_cpu.template getLastProp<2>()[2] = vecDist.template getProp<2>(p)[2];

			++itc;
		}

		vd_cpu.template ghost_get<0,1,2>();

		//! [Fill the ghost on GPU]

		vecDist.template ghost_get<0,1,2>(RUN_ON_DEVICE);

		//! [Fill the ghost on GPU]

		vecDist.deviceToHostPos();
		vecDist.template deviceToHostProp<0,1,2>();

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
		gpu_sort.resize(vecDist.size_local_with_ghost() - vecDist.size_local());

		BOOST_REQUIRE_EQUAL(cpu_sort.size(),gpu_sort.size());

		size_t cnt = 0;

		auto itc2 = vecDist.getGhostIterator();
		while (itc2.isNext())
		{
			auto p = itc2.get();

			cpu_sort.get(cnt).xp.get(0) = vd_cpu.getPos(p)[0];
			gpu_sort.get(cnt).xp.get(0) = vecDist.getPos(p)[0];
			cpu_sort.get(cnt).xp.get(1) = vd_cpu.getPos(p)[1];
			gpu_sort.get(cnt).xp.get(1) = vecDist.getPos(p)[1];
			cpu_sort.get(cnt).xp.get(2) = vd_cpu.getPos(p)[2];
			gpu_sort.get(cnt).xp.get(2) = vecDist.getPos(p)[2];

			cpu_sort.get(cnt).prp0 = vd_cpu.template getProp<0>(p);
			gpu_sort.get(cnt).prp0 = vecDist.template getProp<0>(p);

			cpu_sort.get(cnt).prp1[0] = vd_cpu.template getProp<1>(p)[0];
			gpu_sort.get(cnt).prp1[0] = vecDist.template getProp<1>(p)[0];
			cpu_sort.get(cnt).prp1[1] = vd_cpu.template getProp<1>(p)[1];
			gpu_sort.get(cnt).prp1[1] = vecDist.template getProp<1>(p)[1];
			cpu_sort.get(cnt).prp1[2] = vd_cpu.template getProp<1>(p)[2];
			gpu_sort.get(cnt).prp1[2] = vecDist.template getProp<1>(p)[2];

			cpu_sort.get(cnt).prp2[0] = vd_cpu.template getProp<2>(p)[0];
			gpu_sort.get(cnt).prp2[0] = vecDist.template getProp<2>(p)[0];
			cpu_sort.get(cnt).prp2[1] = vd_cpu.template getProp<2>(p)[1];
			gpu_sort.get(cnt).prp2[1] = vecDist.template getProp<2>(p)[1];
			cpu_sort.get(cnt).prp2[2] = vd_cpu.template getProp<2>(p)[2];
			gpu_sort.get(cnt).prp2[2] = vecDist.template getProp<2>(p)[2];

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

		auto ite = vecDist.getDomainIteratorGPU();
		CUDA_LAUNCH_DIM3((move_parts_gpu_test<3,decltype(vecDist.toKernel())>),ite.wthr,ite.thr,vecDist.toKernel());
	}
}

BOOST_AUTO_TEST_CASE( vector_dist_map_on_gpu_test)
{
	vdist_calc_gpu_test<float>();
	vdist_calc_gpu_test<double>();
}

BOOST_AUTO_TEST_CASE(vector_dist_reduce)
{
	auto & vCluster = create_vcluster();

	if (vCluster.size() > 16)
	{return;}

	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	// set the ghost based on the radius cut off (make just a little bit smaller than the spacing)
	Ghost<3,float> g(0.1);

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	vector_dist_gpu<3,float,aggregate<float,double,int,size_t>> vecDist(5000*vCluster.size(),domain,bc,g);

	auto it = vecDist.getDomainIterator();

	float fc = 1.0;
	double dc = 1.0;
	int ic = 1.0;
	size_t sc = 1.0;

	while(it.isNext())
	{
		auto p = it.get();

		vecDist.template getProp<0>(p) = fc;
		vecDist.template getProp<1>(p) = dc;
		vecDist.template getProp<2>(p) = ic;
		vecDist.template getProp<3>(p) = sc;

		fc += 1.0;
		dc += 1.0;
		ic += 1;
		sc += 1;

		++it;
	}

	vecDist.template hostToDeviceProp<0,1,2,3>();

	float redf = reduce_local<0,_add_>(vecDist);
	double redd = reduce_local<1,_add_>(vecDist);
	int redi = reduce_local<2,_add_>(vecDist);
	size_t reds = reduce_local<3,_add_>(vecDist);

	BOOST_REQUIRE_EQUAL(redf,(vecDist.size_local()+1.0)*(vecDist.size_local())/2.0);
	BOOST_REQUIRE_EQUAL(redd,(vecDist.size_local()+1.0)*(vecDist.size_local())/2.0);
	BOOST_REQUIRE_EQUAL(redi,(vecDist.size_local()+1)*(vecDist.size_local())/2);
	BOOST_REQUIRE_EQUAL(reds,(vecDist.size_local()+1)*(vecDist.size_local())/2);

	float redf2 = reduce_local<0,_max_>(vecDist);
	double redd2 = reduce_local<1,_max_>(vecDist);
	int redi2 = reduce_local<2,_max_>(vecDist);
	size_t reds2 = reduce_local<3,_max_>(vecDist);

	BOOST_REQUIRE_EQUAL(redf2,vecDist.size_local());
	BOOST_REQUIRE_EQUAL(redd2,vecDist.size_local());
	BOOST_REQUIRE_EQUAL(redi2,vecDist.size_local());
	BOOST_REQUIRE_EQUAL(reds2,vecDist.size_local());
}

template<typename CellList_type, bool sorted>
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

	Vcluster<> & vCluster = create_vcluster();

	if (vCluster.getProcessingUnits() > 8)
		return;

	std::uniform_real_distribution<double> unif(0.0,0.3);

	Box<3,double> domain({0.0,0.0,0.0},{1.0,1.0,1.0});
	Ghost<3,double> g(0.1);
	size_t bc[3] = {PERIODIC,PERIODIC,PERIODIC};

	vector_type vecDist(0,domain,bc,g,DEC_GRAN(2048));

	// Only processor 0 initialy add particles on a corner of a domain

	if (vCluster.getProcessUnitID() == 0)
	{
		for(size_t i = 0 ; i < k ; i++)
		{
			vecDist.add();

			vecDist.getLastPos()[0] = unif(e2);
			vecDist.getLastPos()[1] = unif(e2);
			vecDist.getLastPos()[2] = unif(e2);
		}
	}

	// Move to GPU
	vecDist.hostToDevicePos();
	vecDist.template hostToDeviceProp<0>();

	vecDist.map(RUN_ON_DEVICE);
	vecDist.template ghost_get<>(RUN_ON_DEVICE);

	// now move to CPU

	vecDist.deviceToHostPos();
	vecDist.template deviceToHostProp<0>();

	// Get the neighborhood of each particles

	auto VV = vecDist.getVerlet(r_cut);

	// store the number of neighborhood for each particles

	auto it = vecDist.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		vecDist.template getProp<0>(p) = VV.getNNPart(p.getKey());

		++it;
	}

	// Move to GPU
	vecDist.template hostToDeviceProp<0>();

	ModelSquare md;
	md.factor = 10;
	vecDist.addComputationCosts(md);
	vecDist.getDecomposition().decompose();
	vecDist.map(RUN_ON_DEVICE);

	vecDist.deviceToHostPos();
	// Move info to CPU for addComputationcosts

	vecDist.addComputationCosts(md);

	openfpm::vector<size_t> loads;
	size_t load = vecDist.getDecomposition().getDistribution().getProcessorLoad();
	vCluster.allGather(load,loads);
	vCluster.execute();

	for (size_t i = 0 ; i < loads.size() ; i++)
	{
		double load_f = load;
		double load_fc = loads.get(i);

		BOOST_REQUIRE_CLOSE(load_f,load_fc,7.0);
	}

	BOOST_REQUIRE(vecDist.size_local() != 0);

	Point<3,double> v({1.0,1.0,1.0});

	for (size_t i = 0 ; i < 25 ; i++)
	{
		// move particles to CPU and move the particles by 0.1

		vecDist.deviceToHostPos();

		auto it = vecDist.getDomainIterator();

		while (it.isNext())
		{
			auto p = it.get();

			vecDist.getPos(p)[0] += v.get(0) * 0.09;
			vecDist.getPos(p)[1] += v.get(1) * 0.09;
			vecDist.getPos(p)[2] += v.get(2) * 0.09;

			++it;
		}

		//Back to GPU
		vecDist.hostToDevicePos();
		vecDist.map(RUN_ON_DEVICE);
		vecDist.template ghost_get<0>(RUN_ON_DEVICE);

		vecDist.deviceToHostPos();
		vecDist.template deviceToHostProp<0,1,2>();

		// Check calc forces
		size_t opt = (sorted)? CL_GPU_REORDER : 0;
		auto cellList = vecDist.getCellList(r_cut);
		auto cellListGPU = vecDist.template getCellListGPU<CellList_type>(r_cut, opt | CL_NON_SYMMETRIC);

		if (sorted)
		{
			vecDist.hostToDevicePos();
			vecDist.template hostToDeviceProp<0>();

			vecDist.template updateCellListGPU<0>(cellListGPU);
			compareCellListCpuGpuSorted(vecDist,cellListGPU,cellList);
		}

		else
		{
			vecDist.updateCellListGPU(cellListGPU);
			compareCellListCpuGpu(vecDist,cellListGPU,cellList);
		}

		auto VV2 = vecDist.getVerlet(r_cut);

		auto it2 = vecDist.getDomainIterator();

		bool match = true;
		while (it2.isNext())
		{
			auto p = it2.get();

			match &= vecDist.template getProp<0>(p) == VV2.getNNPart(p.getKey());

			++it2;
		}

		BOOST_REQUIRE_EQUAL(match,true);

		ModelSquare md;
		vecDist.addComputationCosts(md);
		vecDist.getDecomposition().redecompose(200);
		vecDist.map(RUN_ON_DEVICE);

		BOOST_REQUIRE(vecDist.size_local() != 0);

		vecDist.template ghost_get<0>(RUN_ON_DEVICE);
		vecDist.deviceToHostPos();
		vecDist.template deviceToHostProp<0>();

		vecDist.addComputationCosts(md);

		openfpm::vector<size_t> loads;
		size_t load = vecDist.getDecomposition().getDistribution().getProcessorLoad();
		vCluster.allGather(load,loads);
		vCluster.execute();

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

template<typename CellList_type, bool sorted>
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

	Vcluster<> & vCluster = create_vcluster();

	if (vCluster.getProcessingUnits() > 8)
		return;

	std::uniform_real_distribution<double> unif(0.0,0.3);

	Box<3,double> domain({0.0,0.0,0.0},{1.0,1.0,1.0});
	Ghost<3,double> g(0.1);
	size_t bc[3] = {PERIODIC,PERIODIC,PERIODIC};

	vector_type vecDist(0,domain,bc,g,DEC_GRAN(2048));

	// Only processor 0 initialy add particles on a corner of a domain

	if (vCluster.getProcessUnitID() == 0)
	{
		for(size_t i = 0 ; i < k ; i++)
		{
			vecDist.add();

			vecDist.getLastPos()[0] = unif(e2);
			vecDist.getLastPos()[1] = unif(e2);
			vecDist.getLastPos()[2] = unif(e2);
		}
	}

	// Move to GPU
	vecDist.hostToDevicePos();
	vecDist.template hostToDeviceProp<0>();

	vecDist.map(RUN_ON_DEVICE);
	vecDist.template Ighost_get<>(RUN_ON_DEVICE);
	vecDist.template ghost_wait<>(RUN_ON_DEVICE);

	// now move to CPU

	vecDist.deviceToHostPos();
	vecDist.template deviceToHostProp<0>();

	// Get the neighborhood of each particles

	auto VV = vecDist.getVerlet(r_cut);

	// store the number of neighborhood for each particles

	auto it = vecDist.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		vecDist.template getProp<0>(p) = VV.getNNPart(p.getKey());

		++it;
	}

	// Move to GPU
	vecDist.template hostToDeviceProp<0>();

	ModelSquare md;
	md.factor = 10;
	vecDist.addComputationCosts(md);
	vecDist.getDecomposition().decompose();
	vecDist.map(RUN_ON_DEVICE);

	vecDist.deviceToHostPos();
	// Move info to CPU for addComputationcosts

	vecDist.addComputationCosts(md);

	openfpm::vector<size_t> loads;
	size_t load = vecDist.getDecomposition().getDistribution().getProcessorLoad();
	vCluster.allGather(load,loads);
	vCluster.execute();

	for (size_t i = 0 ; i < loads.size() ; i++)
	{
		double load_f = load;
		double load_fc = loads.get(i);

		BOOST_REQUIRE_CLOSE(load_f,load_fc,7.0);
	}

	BOOST_REQUIRE(vecDist.size_local() != 0);

	Point<3,double> v({1.0,1.0,1.0});

	for (size_t i = 0 ; i < 25 ; i++)
	{
		// move particles to CPU and move the particles by 0.1

		vecDist.deviceToHostPos();

		auto it = vecDist.getDomainIterator();

		while (it.isNext())
		{
			auto p = it.get();

			vecDist.getPos(p)[0] += v.get(0) * 0.09;
			vecDist.getPos(p)[1] += v.get(1) * 0.09;
			vecDist.getPos(p)[2] += v.get(2) * 0.09;

			++it;
		}

		// Back to GPU
		vecDist.hostToDevicePos();
		vecDist.map(RUN_ON_DEVICE);
		vecDist.template Ighost_get<0>(RUN_ON_DEVICE);
		vecDist.template ghost_wait<0>(RUN_ON_DEVICE);
		vecDist.deviceToHostPos();
		vecDist.template deviceToHostProp<0,1,2>();

		// Check calc forces
		size_t opt = (sorted)? CL_GPU_REORDER : 0;
		auto cellList = vecDist.getCellList(r_cut);
		auto cellListGPU = vecDist.template getCellListGPU<CellList_type>(r_cut, opt | CL_NON_SYMMETRIC);

		if (sorted)
		{
			vecDist.hostToDevicePos();
			vecDist.template hostToDeviceProp<0>();

			vecDist.template updateCellListGPU<0>(cellListGPU);
			compareCellListCpuGpuSorted(vecDist,cellListGPU,cellList);
		}

		else
		{
			vecDist.updateCellListGPU(cellListGPU);
			compareCellListCpuGpu(vecDist,cellListGPU,cellList);
		}

		auto VV2 = vecDist.getVerlet(r_cut);

		auto it2 = vecDist.getDomainIterator();

		bool match = true;
		while (it2.isNext())
		{
			auto p = it2.get();

			match &= vecDist.template getProp<0>(p) == VV2.getNNPart(p.getKey());

			++it2;
		}

		BOOST_REQUIRE_EQUAL(match,true);

		ModelSquare md;
		vecDist.addComputationCosts(md);
		vecDist.getDecomposition().redecompose(200);
		vecDist.map(RUN_ON_DEVICE);

		BOOST_REQUIRE(vecDist.size_local() != 0);

		vecDist.template Ighost_get<0>(RUN_ON_DEVICE);
		vecDist.template ghost_wait<0>(RUN_ON_DEVICE);
		vecDist.deviceToHostPos();
		vecDist.template deviceToHostProp<0>();

		vecDist.addComputationCosts(md);

		openfpm::vector<size_t> loads;
		size_t load = vecDist.getDecomposition().getDistribution().getProcessorLoad();
		vCluster.allGather(load,loads);
		vCluster.execute();

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
	vector_dist_dlb_on_cuda_impl_async<CellList_gpu<3,double,CudaMemory,shift_only<3,double>,false>, false>(50000,0.01);
}

BOOST_AUTO_TEST_CASE(vector_dist_dlb_on_cuda_async_sorted)
{
	vector_dist_dlb_on_cuda_impl_async<CellList_gpu<3,double,CudaMemory,shift_only<3,double>,false>, true>(50000,0.01);
}

BOOST_AUTO_TEST_CASE(vector_dist_dlb_on_cuda)
{
	vector_dist_dlb_on_cuda_impl<CellList_gpu<3,double,CudaMemory,shift_only<3,double>,false>, false>(50000,0.01);
}

BOOST_AUTO_TEST_CASE(vector_dist_dlb_on_cuda_sorted)
{
	vector_dist_dlb_on_cuda_impl<CellList_gpu<3,double,CudaMemory,shift_only<3,double>,false>, true>(50000,0.01);
}

BOOST_AUTO_TEST_CASE(vector_dist_dlb_on_cuda_sparse)
{
	vector_dist_dlb_on_cuda_impl<CELLLIST_GPU_SPARSE<3,double>, false>(50000,0.01);
}

BOOST_AUTO_TEST_CASE(vector_dist_dlb_on_cuda_sparse_sorted)
{
	vector_dist_dlb_on_cuda_impl<CELLLIST_GPU_SPARSE<3,double>, true>(50000,0.01);
}

BOOST_AUTO_TEST_CASE(vector_dist_dlb_on_cuda2)
{
	if (create_vcluster().size() <= 3)
	{return;};

	#ifndef CUDA_ON_CPU
	vector_dist_dlb_on_cuda_impl<CellList_gpu<3,double,CudaMemory,shift_only<3,double>,false>, false>(1000000,0.01);
	#endif
}

BOOST_AUTO_TEST_CASE(vector_dist_dlb_on_cuda3)
{
	if (create_vcluster().size() < 8)
	{return;}

	#ifndef CUDA_ON_CPU
	vector_dist_dlb_on_cuda_impl<CellList_gpu<3,double,CudaMemory,shift_only<3,double>,false>, false>(15000000,0.005);
	#endif
}


BOOST_AUTO_TEST_CASE(vector_dist_keep_prop_on_cuda)
{
	typedef vector_dist_gpu<3,double,aggregate<double,double[3],double[3][3]>> vector_type;

	Vcluster<> & vCluster = create_vcluster();

	if (vCluster.getProcessingUnits() > 8)
		return;

	Box<3,double> domain({0.0,0.0,0.0},{1.0,1.0,1.0});
	Ghost<3,double> g(0.1);
	size_t bc[3] = {PERIODIC,PERIODIC,PERIODIC};

	vector_type vecDist(0,domain,bc,g,DEC_GRAN(2048));

	// Only processor 0 initialy add particles on a corner of a domain

	if (vCluster.getProcessUnitID() == 0)
	{
		for(size_t i = 0 ; i < 50000 ; i++)
		{
			vecDist.add();

			vecDist.getLastPos()[0] = ((double)rand())/RAND_MAX * 0.3;
			vecDist.getLastPos()[1] = ((double)rand())/RAND_MAX * 0.3;
			vecDist.getLastPos()[2] = ((double)rand())/RAND_MAX * 0.3;
		}
	}

	// Move to GPU
	vecDist.hostToDevicePos();
	vecDist.template hostToDeviceProp<0>();

	vecDist.map(RUN_ON_DEVICE);
	vecDist.template ghost_get<>(RUN_ON_DEVICE);

	// now move to CPU

	vecDist.deviceToHostPos();
	vecDist.template deviceToHostProp<0>();


	// store the number of neighborhood for each particles

	auto it = vecDist.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		vecDist.template getProp<0>(p) = 0.0;

		vecDist.template getProp<1>(p)[0] = 1000.0;
		vecDist.template getProp<1>(p)[1] = 2000.0;
		vecDist.template getProp<1>(p)[2] = 3000.0;

		vecDist.template getProp<2>(p)[0][0] = 6000,0;
		vecDist.template getProp<2>(p)[0][1] = 7000.0;
		vecDist.template getProp<2>(p)[0][2] = 8000.0;
		vecDist.template getProp<2>(p)[1][0] = 9000.0;
		vecDist.template getProp<2>(p)[1][1] = 10000.0;
		vecDist.template getProp<2>(p)[1][2] = 11000.0;
		vecDist.template getProp<2>(p)[2][0] = 12000.0;
		vecDist.template getProp<2>(p)[2][1] = 13000.0;
		vecDist.template getProp<2>(p)[2][2] = 14000.0;

		++it;
	}

	// Move to GPU
	vecDist.template hostToDeviceProp<0,1,2>();

	ModelSquare md;
	md.factor = 10;
	vecDist.addComputationCosts(md);
	vecDist.getDecomposition().decompose();
	vecDist.map(RUN_ON_DEVICE);

	vecDist.deviceToHostPos();
	// Move info to CPU for addComputationcosts

	vecDist.addComputationCosts(md);

	openfpm::vector<size_t> loads;
	size_t load = vecDist.getDecomposition().getDistribution().getProcessorLoad();
	vCluster.allGather(load,loads);
	vCluster.execute();

	for (size_t i = 0 ; i < loads.size() ; i++)
	{
		double load_f = load;
		double load_fc = loads.get(i);

		BOOST_REQUIRE_CLOSE(load_f,load_fc,7.0);
	}

	BOOST_REQUIRE(vecDist.size_local() != 0);

	Point<3,double> v({1.0,1.0,1.0});

	int base = 0;

	for (size_t i = 0 ; i < 25 ; i++)
	{
		if (i % 2 == 0)
		{
			// move particles to CPU and move the particles by 0.1

			vecDist.deviceToHostPos();

			auto it = vecDist.getDomainIterator();

			while (it.isNext())
			{
				auto p = it.get();

				vecDist.getPos(p)[0] += v.get(0) * 0.09;
				vecDist.getPos(p)[1] += v.get(1) * 0.09;
				vecDist.getPos(p)[2] += v.get(2) * 0.09;

				++it;
			}

			//Back to GPU
			vecDist.hostToDevicePos();
			vecDist.map(RUN_ON_DEVICE);
			vecDist.template ghost_get<>(RUN_ON_DEVICE);
			vecDist.deviceToHostPos();
			vecDist.template deviceToHostProp<0,1,2>();

			ModelSquare md;
			vecDist.addComputationCosts(md);
			vecDist.getDecomposition().redecompose(200);
			vecDist.map(RUN_ON_DEVICE);

			BOOST_REQUIRE(vecDist.size_local() != 0);

			vecDist.template ghost_get<0>(RUN_ON_DEVICE);
			vecDist.deviceToHostPos();
			vecDist.template deviceToHostProp<0,1,2>();

			vecDist.addComputationCosts(md);

			openfpm::vector<size_t> loads;
			size_t load = vecDist.getDecomposition().getDistribution().getProcessorLoad();
			vCluster.allGather(load,loads);
			vCluster.execute();

			for (size_t i = 0 ; i < loads.size() ; i++)
			{
				double load_f = load;
				double load_fc = loads.get(i);

				BOOST_REQUIRE_CLOSE(load_f,load_fc,10.0);
			}
		}
		else
		{
			vecDist.template deviceToHostProp<0,1,2>();

			auto it2 = vecDist.getDomainIterator();

			bool match = true;
			while (it2.isNext())
			{
				auto p = it2.get();

				vecDist.template getProp<0>(p) += 1;

				vecDist.template getProp<1>(p)[0] += 1.0;
				vecDist.template getProp<1>(p)[1] += 1.0;
				vecDist.template getProp<1>(p)[2] += 1.0;

				vecDist.template getProp<2>(p)[0][0] += 1.0;
				vecDist.template getProp<2>(p)[0][1] += 1.0;
				vecDist.template getProp<2>(p)[0][2] += 1.0;
				vecDist.template getProp<2>(p)[1][0] += 1.0;
				vecDist.template getProp<2>(p)[1][1] += 1.0;
				vecDist.template getProp<2>(p)[1][2] += 1.0;
				vecDist.template getProp<2>(p)[2][0] += 1.0;
				vecDist.template getProp<2>(p)[2][1] += 1.0;
				vecDist.template getProp<2>(p)[2][2] += 1.0;

				++it2;
			}

			vecDist.template hostToDeviceProp<0,1,2>();

			++base;

			vecDist.template ghost_get<0,1,2>(RUN_ON_DEVICE | KEEP_PROPERTIES);
			vecDist.template deviceToHostProp<0,1,2>();

			// Check that the ghost contain the correct information

			auto itg = vecDist.getGhostIterator();

			while (itg.isNext())
			{
				auto p = itg.get();

				match &= vecDist.template getProp<0>(p) == base;

				match &= vecDist.template getProp<1>(p)[0] == base + 1000.0;
				match &= vecDist.template getProp<1>(p)[1] == base + 2000.0;
				match &= vecDist.template getProp<1>(p)[2] == base + 3000.0;

				match &= vecDist.template getProp<2>(p)[0][0] == base + 6000.0;
				match &= vecDist.template getProp<2>(p)[0][1] == base + 7000.0;
				match &= vecDist.template getProp<2>(p)[0][2] == base + 8000.0;
				match &= vecDist.template getProp<2>(p)[1][0] == base + 9000.0;
				match &= vecDist.template getProp<2>(p)[1][1] == base + 10000.0;
				match &= vecDist.template getProp<2>(p)[1][2] == base + 11000.0;
				match &= vecDist.template getProp<2>(p)[2][0] == base + 12000.0;
				match &= vecDist.template getProp<2>(p)[2][1] == base + 13000.0;
				match &= vecDist.template getProp<2>(p)[2][2] == base + 14000.0;

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

	Vcluster<> & vCluster = create_vcluster();

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
	vdg.updateCellListGPU(cl);

	// than we get a cell-list to force reorder

	openfpm::vector_gpu<aggregate<unsigned int>> ids;

	get_indexes_by_type<0,type_is_one>(vdg.getPropVector(),ids,vdg.size_local(),vCluster.getGpuContext());

	// test

	ids.template deviceToHost<0>();

	auto & vs = vdg.getPropVector();
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

	Vcluster<> & vCluster = create_vcluster();

	long int k = 25*25*25*create_vcluster().getProcessingUnits();
	k = std::pow(k, 1/3.);

	if (vCluster.getProcessingUnits() > 48)
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
		vector_dist_gpu<3,float, part_prop > vecDist(0,box,bc,ghost);

		auto it = vecDist.getGridIterator({(size_t)k,(size_t)k,(size_t)k});

		while (it.isNext())
		{
			auto key = it.get();

			vecDist.add();

			vecDist.getLastPosWrite()[0] = key.get(0)*it.getSpacing(0);
			vecDist.getLastPosWrite()[1] = key.get(1)*it.getSpacing(1);
			vecDist.getLastPosWrite()[2] = key.get(2)*it.getSpacing(2);

			// Fill some properties randomly

			vecDist.getLastPropWrite<0>() = 0.0;

			vecDist.getLastPropWrite<2>() = 0.0;

			++it;
		}

		vecDist.map();

		vecDist.hostToDevicePos();
		vecDist.template hostToDeviceProp<0,2>();
		// sync the ghost
		vecDist.ghost_get<0,2>(RUN_ON_DEVICE);
		vecDist.template deviceToHostProp<0,2>();
		vecDist.deviceToHostPos();

		{
			auto NN = vecDist.getCellList(r_cut);
			float a = 1.0f*k*k;

			// run trough all the particles + ghost

			auto it2 = vecDist.getDomainIterator();

			while (it2.isNext())
			{
				// particle p
				auto p = it2.get();
				Point<3,float> xp = vecDist.getPos(p);

				// Get an iterator over the neighborhood particles of p
				auto Np = NN.getNNIterator(NN.getCell(xp));

				// For each neighborhood particle ...
				while (Np.isNext())
				{
					auto q = Np.get();
					Point<3,float> xq = vecDist.getPosRead(q);

					float dist = xp.distance(xq);

					if (dist < r_cut)
					{
						vecDist.getPropWrite<0>(q) += a*(-dist*dist+r_cut*r_cut);
						vecDist.getPropWrite<2>(q) += a*(-dist*dist+r_cut*r_cut) / 2;
					}

					++Np;
				}

				++it2;
			}

			vecDist.hostToDevicePos();
			vecDist.template hostToDeviceProp<0,2>();
			vecDist.template ghost_put<add_atomic_,0,2>(RUN_ON_DEVICE);
			vecDist.template deviceToHostProp<0,2>();
			vecDist.deviceToHostPos();

			bool ret = true;
			auto it3 = vecDist.getDomainIterator();

			float constant = vecDist.getProp<0>(it3.get());
			float constanta = vecDist.getProp<2>(it3.get());
			float eps = 0.001;

			while (it3.isNext())
			{
				float constant2 = vecDist.getProp<0>(it3.get());
				float constant3 = vecDist.getProp<2>(it3.get());
				if (fabs(constant - constant2)/constant > eps || fabs(constanta - constant3)/constanta > eps)
				{
					Point<3,float> p = vecDist.getPosRead(it3.get());

					std::cout << p.toString() << "    " <<  constant2 << "/" << constant << "/" << constant3 << "    " << vCluster.getProcessUnitID() << std::endl;
					ret = false;
					break;
				}

				++it3;
			}
			BOOST_REQUIRE_EQUAL(ret,true);
		}

		auto itp = vecDist.getDomainAndGhostIterator();
		while (itp.isNext())
		{
			auto key = itp.get();

			vecDist.getPropWrite<0>(key) = 0.0;
			vecDist.getPropWrite<2>(key) = 0.0;

			++itp;
		}

		{
			auto NN = vecDist.getCellList(r_cut);
			float a = 1.0f*k*k;

			// run trough all the particles + ghost

			auto it2 = vecDist.getDomainIterator();

			while (it2.isNext())
			{
				// particle p
				auto p = it2.get();
				Point<3,float> xp = vecDist.getPosRead(p);

				// Get an iterator over the neighborhood particles of p
				auto Np = NN.getNNIterator(NN.getCell(xp));

				// For each neighborhood particle ...
				while (Np.isNext())
				{
					auto q = Np.get();
					Point<3,float> xq = vecDist.getPosRead(q);

					float dist = xp.distance(xq);

					if (dist < r_cut)
					{
						vecDist.getPropWrite<0>(q) += a*(-dist*dist+r_cut*r_cut);
						vecDist.getPropWrite<2>(q) += a*(-dist*dist+r_cut*r_cut);
					}

					++Np;
				}

				++it2;
			}

			vecDist.hostToDevicePos();
			vecDist.template hostToDeviceProp<0,2>();
			vecDist.template ghost_put<add_atomic_,0>(RUN_ON_DEVICE);
			vecDist.template ghost_put<add_atomic_,2>(RUN_ON_DEVICE);
			vecDist.template deviceToHostProp<0,2>();
			vecDist.deviceToHostPos();

			bool ret = true;
			auto it3 = vecDist.getDomainIterator();

			float constant = vecDist.getPropRead<0>(it3.get());
			float constanta = vecDist.getPropRead<2>(it3.get());
			float eps = 0.001;

			while (it3.isNext())
			{
				float constant2 = vecDist.getPropRead<0>(it3.get());
				float constant3 = vecDist.getPropRead<0>(it3.get());
				if (fabs(constant - constant2)/constant > eps || fabs(constanta - constant3)/constanta > eps)
				{
					Point<3,float> p = vecDist.getPosRead(it3.get());

					std::cout << p.toString() << "    " <<  constant2 << "/" << constant << "/" << constant3 << "    " << vCluster.getProcessUnitID() << std::endl;
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
