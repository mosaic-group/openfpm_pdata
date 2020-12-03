/*
 * SparseGridGpu single GPU benchmarks using the distributed interface.
 * 
 * Tommaso Bianucci, 17 Nov 2020
 */

#define BOOST_TEST_DYN_LINK
#define DISABLE_MPI_WRITTERS
// Generic perf-test stuff
#include <boost/test/unit_test.hpp>
// #include "performancePlots.hpp"
#include <iostream>
extern std::string suiteURI;
// extern report_sparse_grid_tests report_sparsegrid_funcs;q
// extern std::set<std::string> testSet;
// Specific include
#include "Grid/grid_dist_id.hpp"

void print_test_name()
{
	std::cout << std::endl; // Empty line
	std::string fullName = boost::unit_test::framework::current_test_case().full_name();
	std::cout << ": Test=" << fullName << std::endl;
}

void print_results_insert(const size_t gridSize,
					const size_t pitch,
					const float occupancy,
					const size_t numElements,
					const float elapsedTime,
					const float insertionRate)
{
	std::cout << ":: gridSize=" << gridSize
		<< "; pitch=" << pitch
		<< "; occupancy=" << ( (float) numElements ) / gridSize
		<< "; numElements=" << numElements
		<< "; time=" << elapsedTime
		<< "; insertion rate=" << insertionRate << " GElem/s"
		<< std::endl;
}
void print_results_stencil(const size_t gridSize,
					const size_t pitch,
					const float occupancy,
					const size_t numElements,
					const float elapsedTime,
					const float processingRate,
					const float gflops)
{
	std::cout << ":: gridSize=" << gridSize
		<< "; pitch=" << pitch
		<< "; occupancy=" << ( (float) numElements ) / gridSize
		<< "; numElements=" << numElements
		<< "; time=" << elapsedTime
		<< "; processing rate=" << processingRate << " GElem/s"
		<< "; throughput=" << gflops << " GFlops/s"
		<< std::endl;
}

template <typename SgridT, typename BoxT, typename cT>
float insertConcentricSpheres2D(SgridT &gdist, 
									BoxT box, 
									cT c, 
									const unsigned int pitch, 
									const float occupancy,
									const unsigned int iterations=1)
{
	// unsigned int r1 = 10, r2 = 10+r1;
	unsigned int r1 = (unsigned int) round(pitch*(1-occupancy)), r2 = pitch;

	typedef typename GetAddBlockType<SgridT>::type InsertBlockT;

	auto c0=c[0], c1=c[1];

	// Start timing
	timer ts;
	cudaDeviceSynchronize();
	ts.start();

	for (unsigned int it=0; it<iterations; ++it)
	{
		gdist.addPoints(box.getKP1(),box.getKP2(),
				        [c0,c1,r1,r2] __device__ (int i, int j)
				        {
				        	// This controls if to insert in position (i,j)
				        	// The interval for insertion is [r1,r2) (half-closed)
				        	float r = sqrtf(
				        					(i-c0)*(i-c0)
				        					+ (j-c1)*(j-c1)
				        				  ); // 2-norm
				        	r = fmodf(r, r2);
							return r>=r1;
				        },
				        [it] __device__ (InsertBlockT & data, int i, int j)
				        {
				        	data.template get<0>() = it*(i + j);
				        }
				        );

		gdist.template flush<smax_<0>>(flush_type::FLUSH_ON_DEVICE);
		// gdist.template ghost_get<0,1>(RUN_ON_DEVICE);
	}

	// Stop timing
	cudaDeviceSynchronize();
	ts.stop();

	gdist.template deviceToHost<0>();

    return ts.getwct();
}
template <typename SgridT, typename BoxT, typename cT>
float insertConcentricSpheres3D(SgridT &gdist, 
									BoxT box, 
									cT c, 
									const unsigned int pitch, 
									const float occupancy,
									const unsigned int iterations=1)
{
	// unsigned int r1 = 10, r2 = 10+r1;
	unsigned int r1 = (unsigned int) round(pitch*(1-occupancy)), r2 = pitch;

	typedef typename GetAddBlockType<SgridT>::type InsertBlockT;

	auto c0=c[0], c1=c[1], c2=c[2];

	// Start timing
	timer ts;
	cudaDeviceSynchronize();
	ts.start();

	for (unsigned int it=0; it<iterations; ++it)
	{
		gdist.addPoints(box.getKP1(),box.getKP2(),
				        [c0,c1,c2,r1,r2] __device__ (int i, int j, int k)
				        {
				        	// This controls if to insert in position (i,j,k)
				        	// The interval for insertion is [r1,r2) (half-closed)
				        	float r = sqrtf(
				        					(i-c0)*(i-c0)
				        					+ (j-c1)*(j-c1)
				        					+ (k-c2)*(k-c2)
				        				  ); // 2-norm
				        	r = fmodf(r, r2);
							return r>=r1;
				        },
				        [it] __device__ (InsertBlockT & data, int i, int j, int k)
				        {
				        	data.template get<0>() = it*(i + j + k);
				        }
				        );

		gdist.template flush<smax_<0>>(flush_type::FLUSH_ON_DEVICE);
		// gdist.template ghost_get<0,1>(RUN_ON_DEVICE);
	}

	// Stop timing
	cudaDeviceSynchronize();
	ts.stop();

	gdist.template deviceToHost<0>();

    return ts.getwct();
}

template <typename SgridT, typename BoxT>
float insertFullGrid2D(SgridT &gdist, 
						BoxT box,
						const unsigned int iterations=1)
{
	typedef typename GetAddBlockType<SgridT>::type InsertBlockT;

	// Start timing
	timer ts;
	cudaDeviceSynchronize();
	ts.start();

	for (unsigned int it=0; it<iterations; ++it)
	{
		gdist.addPoints(box.getKP1(),box.getKP2(),
				        [] __device__ (int i, int j)
				        {
				        	// This controls if to insert in position (i,j,k)
							return true;
				        },
				        [it] __device__ (InsertBlockT & data, int i, int j)
				        {
				        	data.template get<0>() = it*(i + j);
				        }
				        );

		gdist.template flush<smax_<0>>(flush_type::FLUSH_ON_DEVICE);
		// gdist.template ghost_get<0,1>(RUN_ON_DEVICE);
	}

	// Stop timing
	cudaDeviceSynchronize();
	ts.stop();

	gdist.template deviceToHost<0>();

    return ts.getwct();
}
template <typename SgridT, typename BoxT>
float insertFullGrid3D(SgridT &gdist, 
						BoxT box,
						const unsigned int iterations=1)
{
	typedef typename GetAddBlockType<SgridT>::type InsertBlockT;

	// Start timing
	timer ts;
	cudaDeviceSynchronize();
	ts.start();

	for (unsigned int it=0; it<iterations; ++it)
	{
		gdist.addPoints(box.getKP1(),box.getKP2(),
				        [] __device__ (int i, int j, int k)
				        {
				        	// This controls if to insert in position (i,j,k)
							return true;
				        },
				        [it] __device__ (InsertBlockT & data, int i, int j, int k)
				        {
				        	data.template get<0>() = it*(i + j + k);
				        }
				        );

		gdist.template flush<smax_<0>>(flush_type::FLUSH_ON_DEVICE);
		// gdist.template ghost_get<0,1>(RUN_ON_DEVICE);
	}

	// Stop timing
	cudaDeviceSynchronize();
	ts.stop();

	gdist.template deviceToHost<0>();

    return ts.getwct();
}

BOOST_AUTO_TEST_SUITE(performance)
BOOST_AUTO_TEST_SUITE(SparseGridGpu_dist_single)

BOOST_AUTO_TEST_SUITE(dim_2D)

template <unsigned int be, unsigned int tbe, unsigned int iterations=100>
void test_2D_insert_spheres(const float occupancy = 0.5, const unsigned int pitch = 32)
{
	// Get local SGG below
	// auto & sg = sgrid.get_loc_grid(0)

	const size_t sz[2] = {10*1024,10*1024};
	periodicity<2> bc = {PERIODIC,PERIODIC};
	Ghost<2,long int> g(1);
	Box<2,float> domain({0.0,0.0},{1.0,1.0});
	sgrid_dist_id_gpu_z_cb<2,float,aggregate<float>, be,tbe*tbe> gdist(sz,domain,g,bc);
	gdist.template setBackgroundValue<0>(666);

	// Box<2,size_t> box({1,1},{sz[0]-1,sz[1]-1});
	Box<2,size_t> box({0,0},{sz[0]-1,sz[1]-1});

	// Insert the concentric spheres on GPU
	// const float occupancy = 0.5;
	// const unsigned int pitch = 32;
	// const unsigned int iterations = 10;
    size_t c[3] = { sz[0]/2, sz[1]/2, 0 };

    auto elapsedTime = insertConcentricSpheres2D(gdist, box, c, pitch, occupancy, iterations);

    auto numElements = gdist.size_local_inserted();
    size_t gridSize = sz[0]*sz[1];
    auto occupancyEmp = ((float) numElements) / gridSize;
    auto insertionRate = 1e-9*numElements*iterations/elapsedTime; // In MElem/s

    print_results_insert(gridSize, pitch, occupancyEmp, numElements, elapsedTime, insertionRate);
}
BOOST_AUTO_TEST_CASE(insert_spheres_2x2_2x2)
{
	print_test_name();
	test_2D_insert_spheres<2,2,50>(0.05, 32);
	test_2D_insert_spheres<2,2,50>(0.1, 32);
	test_2D_insert_spheres<2,2,50>(0.25, 32);
	test_2D_insert_spheres<2,2,50>(0.5, 32);
	test_2D_insert_spheres<2,2,50>(0.75, 32);
	test_2D_insert_spheres<2,2,50>(0.9, 32);
	test_2D_insert_spheres<2,2,50>(1.0, 32);
	test_2D_insert_spheres<2,2,50>(0.05, 64);
	test_2D_insert_spheres<2,2,50>(0.1, 64);
	test_2D_insert_spheres<2,2,50>(0.25, 64);
	test_2D_insert_spheres<2,2,50>(0.5, 64);
	test_2D_insert_spheres<2,2,50>(0.75, 64);
	test_2D_insert_spheres<2,2,50>(0.9, 64);
	test_2D_insert_spheres<2,2,50>(1.0, 64);
}
BOOST_AUTO_TEST_CASE(insert_spheres_4x4_4x4)
{
	print_test_name();
	test_2D_insert_spheres<4,4,512>(0.05, 32);
	test_2D_insert_spheres<4,4,512>(0.1, 32);
	test_2D_insert_spheres<4,4,512>(0.25, 32);
	test_2D_insert_spheres<4,4,512>(0.5, 32);
	test_2D_insert_spheres<4,4,512>(0.75, 32);
	test_2D_insert_spheres<4,4,512>(0.9, 32);
	test_2D_insert_spheres<4,4,512>(1.0, 32);
	test_2D_insert_spheres<4,4,512>(0.05, 64);
	test_2D_insert_spheres<4,4,512>(0.1, 64);
	test_2D_insert_spheres<4,4,512>(0.25, 64);
	test_2D_insert_spheres<4,4,512>(0.5, 64);
	test_2D_insert_spheres<4,4,512>(0.75, 64);
	test_2D_insert_spheres<4,4,512>(0.9, 64);
	test_2D_insert_spheres<4,4,512>(1.0, 64);
}
BOOST_AUTO_TEST_CASE(insert_spheres_8x8_8x8)
{
	print_test_name();
	test_2D_insert_spheres<8,8,1000>(0.05, 32);
	test_2D_insert_spheres<8,8,1000>(0.1, 32);
	test_2D_insert_spheres<8,8,1000>(0.25, 32);
	test_2D_insert_spheres<8,8,1000>(0.5, 32);
	test_2D_insert_spheres<8,8,1000>(0.75, 32);
	test_2D_insert_spheres<8,8,1000>(0.9, 32);
	test_2D_insert_spheres<8,8,1000>(1.0, 32);
	test_2D_insert_spheres<8,8,1000>(0.05, 64);
	test_2D_insert_spheres<8,8,1000>(0.1, 64);
	test_2D_insert_spheres<8,8,1000>(0.25, 64);
	test_2D_insert_spheres<8,8,1000>(0.5, 64);
	test_2D_insert_spheres<8,8,1000>(0.75, 64);
	test_2D_insert_spheres<8,8,1000>(0.9, 64);
	test_2D_insert_spheres<8,8,1000>(1.0, 64);
}
BOOST_AUTO_TEST_CASE(insert_spheres_16x16_16x16)
{
	print_test_name();
	test_2D_insert_spheres<16,16,1000>(0.05, 32);
	test_2D_insert_spheres<16,16,1000>(0.1, 32);
	test_2D_insert_spheres<16,16,1000>(0.25, 32);
	test_2D_insert_spheres<16,16,1000>(0.5, 32);
	test_2D_insert_spheres<16,16,1000>(0.75, 32);
	test_2D_insert_spheres<16,16,1000>(0.9, 32);
	test_2D_insert_spheres<16,16,1000>(1.0, 32);
	test_2D_insert_spheres<16,16,1000>(0.05, 64);
	test_2D_insert_spheres<16,16,1000>(0.1, 64);
	test_2D_insert_spheres<16,16,1000>(0.25, 64);
	test_2D_insert_spheres<16,16,1000>(0.5, 64);
	test_2D_insert_spheres<16,16,1000>(0.75, 64);
	test_2D_insert_spheres<16,16,1000>(0.9, 64);
	test_2D_insert_spheres<16,16,1000>(1.0, 64);
}
BOOST_AUTO_TEST_CASE(insert_spheres_32x32_32x32)
{
	print_test_name();
	test_2D_insert_spheres<32,32,1000>(0.05, 32);
	test_2D_insert_spheres<32,32,1000>(0.1, 32);
	test_2D_insert_spheres<32,32,1000>(0.25, 32);
	test_2D_insert_spheres<32,32,1000>(0.5, 32);
	test_2D_insert_spheres<32,32,1000>(0.75, 32);
	test_2D_insert_spheres<32,32,1000>(0.9, 32);
	test_2D_insert_spheres<32,32,1000>(1.0, 32);
	test_2D_insert_spheres<32,32,1000>(0.05, 64);
	test_2D_insert_spheres<32,32,1000>(0.1, 64);
	test_2D_insert_spheres<32,32,1000>(0.25, 64);
	test_2D_insert_spheres<32,32,1000>(0.5, 64);
	test_2D_insert_spheres<32,32,1000>(0.75, 64);
	test_2D_insert_spheres<32,32,1000>(0.9, 64);
	test_2D_insert_spheres<32,32,1000>(1.0, 64);
}

template <unsigned int be, unsigned int tbe, unsigned int iterations=100>
void test_2D_insert_full()
{
	// Get local SGG below
	// auto & sg = sgrid.get_loc_grid(0)

	const size_t sz[2] = {10*1024,10*1024};
	periodicity<2> bc = {PERIODIC,PERIODIC};
	Ghost<2,long int> g(1);
	Box<2,float> domain({0.0,0.0},{1.0,1.0});
	sgrid_dist_id_gpu_z_cb<2,float,aggregate<float>, be,tbe*tbe> gdist(sz,domain,g,bc);
	gdist.template setBackgroundValue<0>(666);

	// Box<2,size_t> box({1,1},{sz[0]-1,sz[1]-1});
	Box<2,size_t> box({0,0},{sz[0]-1,sz[1]-1});

	// const unsigned int iterations = 50;

	// Insert full on GPU
    auto elapsedTime = insertFullGrid2D(gdist, box, iterations);

    auto numElements = gdist.size_local_inserted();
    size_t gridSize = sz[0]*sz[1];
    auto occupancyEmp = ((float) numElements) / gridSize;
    auto insertionRate = 1e-9*numElements*iterations/elapsedTime; // In MElem/s

    print_results_insert(gridSize, 1U, occupancyEmp, numElements, elapsedTime, insertionRate);
}
BOOST_AUTO_TEST_CASE(insert_full_2x2_2x2)
{
	print_test_name();
	test_2D_insert_full<2,2,50>();
}
BOOST_AUTO_TEST_CASE(insert_full_4x4_4x4)
{
	print_test_name();
	test_2D_insert_full<4,4,512>();
}
BOOST_AUTO_TEST_CASE(insert_full_8x8_8x8)
{
	print_test_name();
	test_2D_insert_full<8,8,1000>();
}
BOOST_AUTO_TEST_CASE(insert_full_16x16_16x16)
{
	print_test_name();
	test_2D_insert_full<16,16,1000>();
}
BOOST_AUTO_TEST_CASE(insert_full_32x32_32x32)
{
	print_test_name();
	test_2D_insert_full<32,32,1000>();
}

template <unsigned int be, unsigned int tbe, unsigned int iterations=100>
void test_2D_stencil_spheres(const float occupancy = 0.5, const unsigned int pitch = 32)
{
	const size_t sz[2] = {10*1024,10*1024};
	periodicity<2> bc = {PERIODIC,PERIODIC};
	Ghost<2,long int> g(1);
	Box<2,float> domain({0.0,0.0},{1.0,1.0});
	sgrid_dist_id_gpu_z_cb<2,float,aggregate<float, float>, be, tbe*tbe> gdist(sz,domain,g,bc);
	gdist.template setBackgroundValue<0>(666);

	// Box<2,size_t> box({1,1},{sz[0]-1,sz[1]-1});
	Box<2,size_t> box({0,0},{sz[0]-1,sz[1]-1});

	// Insert the concentric spheres on GPU
	// const float occupancy = 0.1;
	// const unsigned int pitch = 32;
	// const unsigned int pitch = 10;
	// const unsigned int pitch = 100;
    size_t c[3] = { sz[0]/2, sz[1]/2, 0};

    auto elapsedTime_insert = insertConcentricSpheres2D(gdist, box, c, pitch, occupancy);
    // gdist.tagBoundaries();
    // gdist.findNeighbours();

    // Convolve a stencil
    // GetCpBlockType<GridType, property, stencilSize>
    // typedef typename GetCpBlockType<decltype(gdist),0,1>::type CpBlockType;

    // const unsigned int iterations = 1000;

    timer ts;
    cudaDeviceSynchronize();
    ts.start();

    for (unsigned int it=0; it<iterations; ++it)
    {
		gdist.template conv_cross<0,1,1>(
			{2,2,0},
			{(int)sz[0]-3,(int)sz[1]-3,0},
				[] __device__ (float & u, 
								cross_stencil<2,float> & cs)
			{
				return u 
						+ (cs.xm[0] + cs.xp[0] 
						+ cs.xm[1] + cs.xp[1] 
						- 4.0*u)*0.1;
			}
		);
    	cudaDeviceSynchronize(); // We don't want an overlapping-kernel mess here!
		gdist.template conv_cross<1,0,1>(
			{2,2,0},
			{(int)sz[0]-3,(int)sz[1]-3,0},
				[] __device__ (float & u, 
								cross_stencil<2,float> & cs)
			{
				return u 
						+ (cs.xm[0] + cs.xp[0] 
						+ cs.xm[1] + cs.xp[1] 
						- 4.0*u)*0.1;
			}
		);
    	cudaDeviceSynchronize(); // We don't want an overlapping-kernel mess here!
	}

	// cudaDeviceSynchronize();
	ts.stop();

	gdist.template deviceToHost<0,1>();

	float elapsedTime = ts.getwct();

    auto numElements = gdist.size_local_inserted();
    size_t gridSize = sz[0]*sz[1];
    auto occupancyEmp = ((float) numElements) / gridSize;
    auto processingRate = 1e-9*numElements*2*iterations/elapsedTime; // In GElem/s
    auto gflops = 7*processingRate; // GFlops/s

	print_results_stencil(gridSize, pitch, occupancyEmp, numElements, elapsedTime, processingRate, gflops);
}
BOOST_AUTO_TEST_CASE(stencil_spheres_2x2_2x2)
{
	print_test_name();
	test_2D_stencil_spheres<2,2,100>(0.05, 32);
	test_2D_stencil_spheres<2,2,100>(0.1, 32);
	test_2D_stencil_spheres<2,2,100>(0.25, 32);
	test_2D_stencil_spheres<2,2,100>(0.5, 32);
	test_2D_stencil_spheres<2,2,100>(0.75, 32);
	test_2D_stencil_spheres<2,2,100>(0.9, 32);
	test_2D_stencil_spheres<2,2,100>(1.0, 32);
	test_2D_stencil_spheres<2,2,100>(0.05, 64);
	test_2D_stencil_spheres<2,2,100>(0.1, 64);
	test_2D_stencil_spheres<2,2,100>(0.25, 64);
	test_2D_stencil_spheres<2,2,100>(0.5, 64);
	test_2D_stencil_spheres<2,2,100>(0.75, 64);
	test_2D_stencil_spheres<2,2,100>(0.9, 64);
	test_2D_stencil_spheres<2,2,100>(1.0, 64);
}
BOOST_AUTO_TEST_CASE(stencil_spheres_4x4_4x4)
{
	print_test_name();
	test_2D_stencil_spheres<4,4,512>(0.05, 32);
	test_2D_stencil_spheres<4,4,512>(0.1, 32);
	test_2D_stencil_spheres<4,4,512>(0.25, 32);
	test_2D_stencil_spheres<4,4,512>(0.5, 32);
	test_2D_stencil_spheres<4,4,512>(0.75, 32);
	test_2D_stencil_spheres<4,4,512>(0.9, 32);
	test_2D_stencil_spheres<4,4,512>(1.0, 32);
	test_2D_stencil_spheres<4,4,512>(0.05, 64);
	test_2D_stencil_spheres<4,4,512>(0.1, 64);
	test_2D_stencil_spheres<4,4,512>(0.25, 64);
	test_2D_stencil_spheres<4,4,512>(0.5, 64);
	test_2D_stencil_spheres<4,4,512>(0.75, 64);
	test_2D_stencil_spheres<4,4,512>(0.9, 64);
	test_2D_stencil_spheres<4,4,512>(1.0, 64);
}
BOOST_AUTO_TEST_CASE(stencil_spheres_8x8_8x8)
{
	print_test_name();
	test_2D_stencil_spheres<8,8,1000>(0.05, 32);
	test_2D_stencil_spheres<8,8,1000>(0.1, 32);
	test_2D_stencil_spheres<8,8,1000>(0.25, 32);
	test_2D_stencil_spheres<8,8,1000>(0.5, 32);
	test_2D_stencil_spheres<8,8,1000>(0.75, 32);
	test_2D_stencil_spheres<8,8,1000>(0.9, 32);
	test_2D_stencil_spheres<8,8,1000>(1.0, 32);
	test_2D_stencil_spheres<8,8,1000>(0.05, 64);
	test_2D_stencil_spheres<8,8,1000>(0.1, 64);
	test_2D_stencil_spheres<8,8,1000>(0.25, 64);
	test_2D_stencil_spheres<8,8,1000>(0.5, 64);
	test_2D_stencil_spheres<8,8,1000>(0.75, 64);
	test_2D_stencil_spheres<8,8,1000>(0.9, 64);
	test_2D_stencil_spheres<8,8,1000>(1.0, 64);
}
BOOST_AUTO_TEST_CASE(stencil_spheres_16x16_16x16)
{
	print_test_name();
	test_2D_stencil_spheres<16,16,1000>(0.05, 32);
	test_2D_stencil_spheres<16,16,1000>(0.1, 32);
	test_2D_stencil_spheres<16,16,1000>(0.25, 32);
	test_2D_stencil_spheres<16,16,1000>(0.5, 32);
	test_2D_stencil_spheres<16,16,1000>(0.75, 32);
	test_2D_stencil_spheres<16,16,1000>(0.9, 32);
	test_2D_stencil_spheres<16,16,1000>(1.0, 32);
	test_2D_stencil_spheres<16,16,1000>(0.05, 64);
	test_2D_stencil_spheres<16,16,1000>(0.1, 64);
	test_2D_stencil_spheres<16,16,1000>(0.25, 64);
	test_2D_stencil_spheres<16,16,1000>(0.5, 64);
	test_2D_stencil_spheres<16,16,1000>(0.75, 64);
	test_2D_stencil_spheres<16,16,1000>(0.9, 64);
	test_2D_stencil_spheres<16,16,1000>(1.0, 64);
}
BOOST_AUTO_TEST_CASE(stencil_spheres_32x32_32x32)
{
	print_test_name();
	test_2D_stencil_spheres<32,32,1000>(0.05, 32);
	test_2D_stencil_spheres<32,32,1000>(0.1, 32);
	test_2D_stencil_spheres<32,32,1000>(0.25, 32);
	test_2D_stencil_spheres<32,32,1000>(0.5, 32);
	test_2D_stencil_spheres<32,32,1000>(0.75, 32);
	test_2D_stencil_spheres<32,32,1000>(0.9, 32);
	test_2D_stencil_spheres<32,32,1000>(1.0, 32);
	test_2D_stencil_spheres<32,32,1000>(0.05, 64);
	test_2D_stencil_spheres<32,32,1000>(0.1, 64);
	test_2D_stencil_spheres<32,32,1000>(0.25, 64);
	test_2D_stencil_spheres<32,32,1000>(0.5, 64);
	test_2D_stencil_spheres<32,32,1000>(0.75, 64);
	test_2D_stencil_spheres<32,32,1000>(0.9, 64);
	test_2D_stencil_spheres<32,32,1000>(1.0, 64);
}

template <unsigned int be, unsigned int tbe, unsigned int iterations=100>
void test_2D_stencil_full()
{
	const size_t sz[2] = {10*1024,10*1024};
	periodicity<2> bc = {PERIODIC,PERIODIC};
	Ghost<2,long int> g(1);
	Box<2,float> domain({0.0,0.0},{1.0,1.0});
	sgrid_dist_id_gpu_z_cb<2,float,aggregate<float, float>, be, tbe*tbe> gdist(sz,domain,g,bc);
	gdist.template setBackgroundValue<0>(666);

	// Box<2,size_t> box({1,1},{sz[0]-1,sz[1]-1});
	Box<2,size_t> box({0,0},{sz[0]-1,sz[1]-1});

	// Insert full on GPU
    auto elapsedTime_insert = insertFullGrid2D(gdist, box);
    // gdist.tagBoundaries();
    // gdist.findNeighbours();

    // Convolve a stencil
    // GetCpBlockType<GridType, property, stencilSize>
    // typedef typename GetCpBlockType<decltype(gdist),0,1>::type CpBlockType;

    // const unsigned int iterations = 1000;

    timer ts;
    cudaDeviceSynchronize();
    ts.start();

    for (unsigned int it=0; it<iterations; ++it)
    {
		gdist.template conv_cross<0,1,1>(
			{2,2,0},
			{(int)sz[0]-3,(int)sz[1]-3,0},
				[] __device__ (float & u, 
								cross_stencil<2,float> & cs)
			{
				return u 
						+ (cs.xm[0] + cs.xp[0] 
						+ cs.xm[1] + cs.xp[1] 
						- 4.0*u)*0.1;
			}
		);
    	cudaDeviceSynchronize(); // We don't want an overlapping-kernel mess here!
		gdist.template conv_cross<1,0,1>(
			{2,2,0},
			{(int)sz[0]-3,(int)sz[1]-3,0},
				[] __device__ (float & u, 
								cross_stencil<2,float> & cs)
			{
				return u 
						+ (cs.xm[0] + cs.xp[0] 
						+ cs.xm[1] + cs.xp[1] 
						- 4.0*u)*0.1;
			}
		);
    	cudaDeviceSynchronize(); // We don't want an overlapping-kernel mess here!
	}

	// cudaDeviceSynchronize();
	ts.stop();

	gdist.template deviceToHost<0,1>();

	float elapsedTime = ts.getwct();

    auto numElements = gdist.size_local_inserted();
    size_t gridSize = sz[0]*sz[1];
    auto occupancyEmp = ((float) numElements) / gridSize;
    auto processingRate = 1e-9*numElements*2*iterations/elapsedTime; // In GElem/s
    auto gflops = 7*processingRate; // GFlops/s

	print_results_stencil(gridSize, 1U, occupancyEmp, numElements, elapsedTime, processingRate, gflops);
}
BOOST_AUTO_TEST_CASE(stencil_full_2x2_2x2)
{
	print_test_name();
	test_2D_stencil_full<2,2,100>();
}
BOOST_AUTO_TEST_CASE(stencil_full_4x4_4x4)
{
	print_test_name();
	test_2D_stencil_full<4,4,512>();
}
BOOST_AUTO_TEST_CASE(stencil_full_8x8_8x8)
{
	print_test_name();
	test_2D_stencil_full<8,8,1000>();
}
BOOST_AUTO_TEST_CASE(stencil_full_16x16_16x16)
{
	print_test_name();
	test_2D_stencil_full<16,16,1000>();
}
BOOST_AUTO_TEST_CASE(stencil_full_32x32_32x32)
{
	print_test_name();
	test_2D_stencil_full<32,32,1000>();
}

BOOST_AUTO_TEST_SUITE_END() //dim_2D

BOOST_AUTO_TEST_SUITE(dim_3D)

template <unsigned int be, unsigned int tbe, unsigned int iterations=100>
void test_3D_insert_spheres(const float occupancy = 0.5, const unsigned int pitch = 32)
{
	//const size_t sz[3] = {2*512,2*512,512};
    const size_t sz[3] = {2*256,2*256,256};
	periodicity<3> bc = {PERIODIC,PERIODIC,PERIODIC};
	Ghost<3,long int> g(1);
	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});
	sgrid_dist_id_gpu_z_cb<3,float,aggregate<float>, be,tbe*tbe*tbe> gdist(sz,domain,g,bc);
	gdist.template setBackgroundValue<0>(666);

	// Box<3,size_t> box({1,1,1},{sz[0]-1,sz[1]-1,sz[2]-1});
	Box<3,size_t> box({0,0,0},{sz[0]-1,sz[1]-1,sz[2]-1});

	// Insert the concentric spheres on GPU
	// const unsigned int iterations = 100;
    size_t c[3] = { sz[0]/2, sz[1]/2, sz[2]/2 };

    auto elapsedTime = insertConcentricSpheres3D(gdist, box, c, pitch, occupancy, iterations);

    auto numElements = gdist.size_local_inserted();
    size_t gridSize = sz[0]*sz[1]*sz[2];
    auto occupancyEmp = ((float) numElements) / gridSize;
    auto insertionRate = 1e-9*numElements*iterations/elapsedTime; // In MElem/s

    print_results_insert(gridSize, pitch, occupancyEmp, numElements, elapsedTime, insertionRate);
}
BOOST_AUTO_TEST_CASE(insert_spheres_2x2x2_2x2x2)
{
	print_test_name();
	test_3D_insert_spheres<2,2,50>(0.05, 32);
	test_3D_insert_spheres<2,2,50>(0.1, 32);
	test_3D_insert_spheres<2,2,50>(0.25, 32);
	test_3D_insert_spheres<2,2,50>(0.5, 32);
	test_3D_insert_spheres<2,2,50>(0.75, 32);
	test_3D_insert_spheres<2,2,50>(0.9, 32);
	test_3D_insert_spheres<2,2,50>(1.0, 32);
	test_3D_insert_spheres<2,2,50>(0.05, 64);
	test_3D_insert_spheres<2,2,50>(0.1, 64);
	test_3D_insert_spheres<2,2,50>(0.25, 64);
	test_3D_insert_spheres<2,2,50>(0.5, 64);
	test_3D_insert_spheres<2,2,50>(0.75, 64);
	test_3D_insert_spheres<2,2,50>(0.9, 64);
	test_3D_insert_spheres<2,2,50>(1.0, 64);
}
BOOST_AUTO_TEST_CASE(insert_spheres_4x4x4_4x4x4)
{
	print_test_name();
	test_3D_insert_spheres<4,4,1000>(0.05, 32);
	test_3D_insert_spheres<4,4,1000>(0.1, 32);
	test_3D_insert_spheres<4,4,1000>(0.25, 32);
	test_3D_insert_spheres<4,4,1000>(0.5, 32);
	test_3D_insert_spheres<4,4,1000>(0.75, 32);
	test_3D_insert_spheres<4,4,1000>(0.9, 32);
	test_3D_insert_spheres<4,4,1000>(1.0, 32);
	test_3D_insert_spheres<4,4,1000>(0.05, 64);
	test_3D_insert_spheres<4,4,1000>(0.1, 64);
	test_3D_insert_spheres<4,4,1000>(0.25, 64);
	test_3D_insert_spheres<4,4,1000>(0.5, 64);
	test_3D_insert_spheres<4,4,1000>(0.75, 64);
	test_3D_insert_spheres<4,4,1000>(0.9, 64);
	test_3D_insert_spheres<4,4,1000>(1.0, 64);
}
BOOST_AUTO_TEST_CASE(insert_spheres_8x8x8_8x8x8)
{
	print_test_name();
	test_3D_insert_spheres<8,8,1000>(0.05, 32);
	test_3D_insert_spheres<8,8,1000>(0.1, 32);
	test_3D_insert_spheres<8,8,1000>(0.25, 32);
	test_3D_insert_spheres<8,8,1000>(0.5, 32);
	test_3D_insert_spheres<8,8,1000>(0.75, 32);
	test_3D_insert_spheres<8,8,1000>(0.9, 32);
	test_3D_insert_spheres<8,8,1000>(1.0, 32);
	test_3D_insert_spheres<8,8,1000>(0.05, 64);
	test_3D_insert_spheres<8,8,1000>(0.1, 64);
	test_3D_insert_spheres<8,8,1000>(0.25, 64);
	test_3D_insert_spheres<8,8,1000>(0.5, 64);
	test_3D_insert_spheres<8,8,1000>(0.75, 64);
	test_3D_insert_spheres<8,8,1000>(0.9, 64);
	test_3D_insert_spheres<8,8,1000>(1.0, 64);
}
// BOOST_AUTO_TEST_CASE(insert_spheres_16x16x16_16x16x16)
// {
// 	print_test_name();
// 	test_3D_insert_spheres<16,16,100>(0.05, 32);
// 	test_3D_insert_spheres<16,16,100>(0.1, 32);
// 	test_3D_insert_spheres<16,16,100>(0.25, 32);
// 	test_3D_insert_spheres<16,16,100>(0.5, 32);
// 	test_3D_insert_spheres<16,16,100>(0.75, 32);
// 	test_3D_insert_spheres<16,16,100>(0.9, 32);
// 	test_3D_insert_spheres<16,16,100>(1.0, 32);
// 	test_3D_insert_spheres<16,16,100>(0.05, 64);
// 	test_3D_insert_spheres<16,16,100>(0.1, 64);
// 	test_3D_insert_spheres<16,16,100>(0.25, 64);
// 	test_3D_insert_spheres<16,16,100>(0.5, 64);
// 	test_3D_insert_spheres<16,16,100>(0.75, 64);
// 	test_3D_insert_spheres<16,16,100>(0.9, 64);
// 	test_3D_insert_spheres<16,16,100>(1.0, 64);
// }
// BOOST_AUTO_TEST_CASE(insert_spheres_32x32x32_32x32x32)
// {
// 	print_test_name();
// 	test_3D_insert_spheres<32,32,10>(0.05, 32);
// 	test_3D_insert_spheres<32,32,10>(0.1, 32);
// 	test_3D_insert_spheres<32,32,10>(0.25, 32);
// 	test_3D_insert_spheres<32,32,10>(0.5, 32);
// 	test_3D_insert_spheres<32,32,10>(0.75, 32);
// 	test_3D_insert_spheres<32,32,10>(0.9, 32);
// 	test_3D_insert_spheres<32,32,10>(1.0, 32);
// 	test_3D_insert_spheres<32,32,10>(0.05, 64);
// 	test_3D_insert_spheres<32,32,10>(0.1, 64);
// 	test_3D_insert_spheres<32,32,10>(0.25, 64);
// 	test_3D_insert_spheres<32,32,10>(0.5, 64);
// 	test_3D_insert_spheres<32,32,10>(0.75, 64);
// 	test_3D_insert_spheres<32,32,10>(0.9, 64);
// 	test_3D_insert_spheres<32,32,10>(1.0, 64);
// }

template <unsigned int be, unsigned int tbe, unsigned int iterations=100>
void test_3D_insert_full()
{
	//const size_t sz[3] = {2*512,2*512,512};
    const size_t sz[3] = {2*256,2*256,256};
	periodicity<3> bc = {PERIODIC,PERIODIC,PERIODIC};
	Ghost<3,long int> g(1);
	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});
	sgrid_dist_id_gpu_z_cb<3,float,aggregate<float>, be,tbe*tbe*tbe> gdist(sz,domain,g,bc);
	gdist.template setBackgroundValue<0>(666);

	// Box<3,size_t> box({1,1,1},{sz[0]-1,sz[1]-1,sz[2]-1});
	Box<3,size_t> box({0,0,0},{sz[0]-1,sz[1]-1,sz[2]-1});

	// const unsigned int iterations = 100;

	// Insert full on GPU
    auto elapsedTime = insertFullGrid3D(gdist, box, iterations);

    auto numElements = gdist.size_local_inserted();
    size_t gridSize = sz[0]*sz[1]*sz[2];
    auto occupancyEmp = ((float) numElements) / gridSize;
    auto insertionRate = 1e-9*numElements*iterations/elapsedTime; // In MElem/s

    print_results_insert(gridSize, 1U, occupancyEmp, numElements, elapsedTime, insertionRate);
}
BOOST_AUTO_TEST_CASE(insert_full_2x2x2_2x2x2)
{
	print_test_name();
	test_3D_insert_full<2,2,10>();
}
BOOST_AUTO_TEST_CASE(insert_full_4x4x4_4x4x4)
{
	print_test_name();
	test_3D_insert_full<4,4,100>();
}
BOOST_AUTO_TEST_CASE(insert_full_8x8x8_8x8x8)
{
	print_test_name();
	test_3D_insert_full<8,8,100>();
}
// BOOST_AUTO_TEST_CASE(insert_full_16x16x16_16x16x16)
// {
// 	print_test_name();
// 	test_3D_insert_full<16,16,100>();
// }
// BOOST_AUTO_TEST_CASE(insert_full_32x32x32_32x32x32)
// {
// 	print_test_name();
// 	test_3D_insert_full<32,32,10>();
// }

template <unsigned int be, unsigned int tbe, unsigned int iterations=512>
void test_3D_stencil_spheres(const float occupancy = 0.5, const unsigned int pitch = 32)
{
	//const size_t sz[3] = {2*512,2*512,512};
    const size_t sz[3] = {2*256,2*256,256};
	periodicity<3> bc = {PERIODIC,PERIODIC,PERIODIC};
	Ghost<3,long int> g(1);
	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});
	sgrid_dist_id_gpu_z_cb<3,float,aggregate<float,float>, be, tbe*tbe*tbe> gdist(sz,domain,g,bc);
	gdist.template setBackgroundValue<0>(666);

	// Box<3,size_t> box({1,1,1},{sz[0]-1,sz[1]-1,sz[2]-1});
	Box<3,size_t> box({0,0,0},{sz[0]-1,sz[1]-1,sz[2]-1});

	// Insert the concentric spheres on GPU
	// const float occupancy = 0.1;
	// const unsigned int pitch = 32;
	// const unsigned int pitch = 10;
	// const unsigned int pitch = 100;
    size_t c[3] = { sz[0]/2, sz[1]/2, sz[2]/2 };

    auto elapsedTime_insert = insertConcentricSpheres3D(gdist, box, c, pitch, occupancy);
    // gdist.tagBoundaries();
    // gdist.findNeighbours();

    // Convolve a stencil
    // GetCpBlockType<GridType, property, stencilSize>
    // typedef typename GetCpBlockType<decltype(gdist),0,1>::type CpBlockType;

    // const unsigned int iterations = 512;

    timer ts;
    cudaDeviceSynchronize();
    ts.start();

    for (unsigned int it=0; it<iterations; ++it)
    {
		gdist.template conv_cross<0,1,1>(
			{2,2,2},
			{(int)sz[0]-3,(int)sz[1]-3,(int)sz[2]-3},
				[] __device__ (float & u, 
								cross_stencil<3,float> & cs)
			{
				return u 
						+ (cs.xm[0] + cs.xp[0] 
						+ cs.xm[1] + cs.xp[1] 
						+ cs.xm[2] + cs.xp[2] 
						- 6.0*u)*0.1;
			}
		);
    	cudaDeviceSynchronize(); // We don't want an overlapping-kernel mess here!
		gdist.template conv_cross<1,0,1>(
			{2,2,2},
			{(int)sz[0]-3,(int)sz[1]-3,(int)sz[2]-3},
				[] __device__ (float & u, 
								cross_stencil<3,float> & cs)
			{
				return u 
						+ (cs.xm[0] + cs.xp[0] 
						+ cs.xm[1] + cs.xp[1] 
						+ cs.xm[2] + cs.xp[2] 
						- 6.0*u)*0.1;
			}
		);
    	cudaDeviceSynchronize(); // We don't want an overlapping-kernel mess here!
	}

	// cudaDeviceSynchronize();
	ts.stop();

	gdist.template deviceToHost<0,1>();

	float elapsedTime = ts.getwct();

    auto numElements = gdist.size_local_inserted();
    size_t gridSize = sz[0]*sz[1]*sz[2];
    auto occupancyEmp = ((float) numElements) / gridSize;
    auto processingRate = 1e-9*numElements*2*iterations/elapsedTime; // In GElem/s
    auto gflops = 9*processingRate; // GFlops/s

	print_results_stencil(gridSize, pitch, occupancyEmp, numElements, elapsedTime, processingRate, gflops);
}
BOOST_AUTO_TEST_CASE(stencil_spheres_2x2x2_2x2x2)
{
	print_test_name();
	test_3D_stencil_spheres<2,2,512>(0.05, 32);
	test_3D_stencil_spheres<2,2,512>(0.1, 32);
	test_3D_stencil_spheres<2,2,512>(0.25, 32);
	test_3D_stencil_spheres<2,2,512>(0.5, 32);
	test_3D_stencil_spheres<2,2,512>(0.75, 32);
	test_3D_stencil_spheres<2,2,512>(0.9, 32);
	test_3D_stencil_spheres<2,2,512>(1.0, 32);
	test_3D_stencil_spheres<2,2,512>(0.05, 64);
	test_3D_stencil_spheres<2,2,512>(0.1, 64);
	test_3D_stencil_spheres<2,2,512>(0.25, 64);
	test_3D_stencil_spheres<2,2,512>(0.5, 64);
	test_3D_stencil_spheres<2,2,512>(0.75, 64);
	test_3D_stencil_spheres<2,2,512>(0.9, 64);
	test_3D_stencil_spheres<2,2,512>(1.0, 64);
}
BOOST_AUTO_TEST_CASE(stencil_spheres_4x4x4_4x4x4)
{
	print_test_name();
	test_3D_stencil_spheres<4,4,512>(0.05, 32);
	test_3D_stencil_spheres<4,4,512>(0.1, 32);
	test_3D_stencil_spheres<4,4,512>(0.25, 32);
	test_3D_stencil_spheres<4,4,512>(0.5, 32);
	test_3D_stencil_spheres<4,4,512>(0.75, 32);
	test_3D_stencil_spheres<4,4,512>(0.9, 32);
	test_3D_stencil_spheres<4,4,512>(1.0, 32);
	test_3D_stencil_spheres<4,4,512>(0.05, 64);
	test_3D_stencil_spheres<4,4,512>(0.1, 64);
	test_3D_stencil_spheres<4,4,512>(0.25, 64);
	test_3D_stencil_spheres<4,4,512>(0.5, 64);
	test_3D_stencil_spheres<4,4,512>(0.75, 64);
	test_3D_stencil_spheres<4,4,512>(0.9, 64);
	test_3D_stencil_spheres<4,4,512>(1.0, 64);
}
BOOST_AUTO_TEST_CASE(stencil_spheres_8x8x8_8x8x8)
{
	print_test_name();
	test_3D_stencil_spheres<8,8,512>(0.05, 32);
	test_3D_stencil_spheres<8,8,512>(0.1, 32);
	test_3D_stencil_spheres<8,8,512>(0.25, 32);
	test_3D_stencil_spheres<8,8,512>(0.5, 32);
	test_3D_stencil_spheres<8,8,512>(0.75, 32);
	test_3D_stencil_spheres<8,8,512>(0.9, 32);
	test_3D_stencil_spheres<8,8,512>(1.0, 32);
	test_3D_stencil_spheres<8,8,512>(0.05, 64);
	test_3D_stencil_spheres<8,8,512>(0.1, 64);
	test_3D_stencil_spheres<8,8,512>(0.25, 64);
	test_3D_stencil_spheres<8,8,512>(0.5, 64);
	test_3D_stencil_spheres<8,8,512>(0.75, 64);
	test_3D_stencil_spheres<8,8,512>(0.9, 64);
	test_3D_stencil_spheres<8,8,512>(1.0, 64);
}
// BOOST_AUTO_TEST_CASE(stencil_spheres_16x16x16_16x16x16)
// {
// 	print_test_name();
// 	test_3D_stencil_spheres<16,16,512>(0.05, 32);
// 	test_3D_stencil_spheres<16,16,512>(0.1, 32);
// 	test_3D_stencil_spheres<16,16,512>(0.25, 32);
// 	test_3D_stencil_spheres<16,16,512>(0.5, 32);
// 	test_3D_stencil_spheres<16,16,512>(0.75, 32);
// 	test_3D_stencil_spheres<16,16,512>(0.9, 32);
// 	test_3D_stencil_spheres<16,16,512>(1.0, 32);
// 	test_3D_stencil_spheres<16,16,512>(0.05, 64);
// 	test_3D_stencil_spheres<16,16,512>(0.1, 64);
// 	test_3D_stencil_spheres<16,16,512>(0.25, 64);
// 	test_3D_stencil_spheres<16,16,512>(0.5, 64);
// 	test_3D_stencil_spheres<16,16,512>(0.75, 64);
// 	test_3D_stencil_spheres<16,16,512>(0.9, 64);
// 	test_3D_stencil_spheres<16,16,512>(1.0, 64);
// }
// BOOST_AUTO_TEST_CASE(stencil_spheres_32x32x32_32x32x32)
// {
// 	print_test_name();
// 	test_3D_stencil_spheres<32,32,512>(0.05, 32);
// 	test_3D_stencil_spheres<32,32,512>(0.1, 32);
// 	test_3D_stencil_spheres<32,32,512>(0.25, 32);
// 	test_3D_stencil_spheres<32,32,512>(0.5, 32);
// 	test_3D_stencil_spheres<32,32,512>(0.75, 32);
// 	test_3D_stencil_spheres<32,32,512>(0.9, 32);
// 	test_3D_stencil_spheres<32,32,512>(1.0, 32);
// 	test_3D_stencil_spheres<32,32,512>(0.05, 64);
// 	test_3D_stencil_spheres<32,32,512>(0.1, 64);
// 	test_3D_stencil_spheres<32,32,512>(0.25, 64);
// 	test_3D_stencil_spheres<32,32,512>(0.5, 64);
// 	test_3D_stencil_spheres<32,32,512>(0.75, 64);
// 	test_3D_stencil_spheres<32,32,512>(0.9, 64);
// 	test_3D_stencil_spheres<32,32,512>(1.0, 64);
// }

template <unsigned int be, unsigned int tbe, unsigned int iterations=512>
void test_3D_stencil_full()
{
	//const size_t sz[3] = {2*512,2*512,512};
    const size_t sz[3] = {2*256,2*256,256};
	periodicity<3> bc = {PERIODIC,PERIODIC,PERIODIC};
	Ghost<3,long int> g(1);
	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});
	sgrid_dist_id_gpu_z_cb<3,float,aggregate<float,float>, be, tbe*tbe*tbe> gdist(sz,domain,g,bc);
	gdist.template setBackgroundValue<0>(666);

	// Box<3,size_t> box({1,1,1},{sz[0]-1,sz[1]-1,sz[2]-1});
	Box<3,size_t> box({0,0,0},{sz[0]-1,sz[1]-1,sz[2]-1});

	// Insert full on GPU
    auto elapsedTime_insert = insertFullGrid3D(gdist, box);
    // gdist.tagBoundaries();
    // gdist.findNeighbours();

    // Convolve a stencil
    // GetCpBlockType<GridType, property, stencilSize>
    // typedef typename GetCpBlockType<decltype(gdist),0,1>::type CpBlockType;

    // const unsigned int iterations = 512;

    timer ts;
    cudaDeviceSynchronize();
    ts.start();

    for (unsigned int it=0; it<iterations; ++it)
    {
		gdist.template conv_cross<0,1,1>(
			{2,2,2},
			{(int)sz[0]-3,(int)sz[1]-3,(int)sz[2]-3},
				[] __device__ (float & u, 
								cross_stencil<3,float> & cs)
			{
				return u 
						+ (cs.xm[0] + cs.xp[0] 
						+ cs.xm[1] + cs.xp[1] 
						+ cs.xm[2] + cs.xp[2] 
						- 6.0*u)*0.1;
			}
		);
    	cudaDeviceSynchronize(); // We don't want an overlapping-kernel mess here!
		gdist.template conv_cross<1,0,1>(
			{2,2,2},
			{(int)sz[0]-3,(int)sz[1]-3,(int)sz[2]-3},
				[] __device__ (float & u, 
								cross_stencil<3,float> & cs)
			{
				return u 
						+ (cs.xm[0] + cs.xp[0] 
						+ cs.xm[1] + cs.xp[1] 
						+ cs.xm[2] + cs.xp[2] 
						- 6.0*u)*0.1;
			}
		);
    	cudaDeviceSynchronize(); // We don't want an overlapping-kernel mess here!
	}

	// cudaDeviceSynchronize();
	ts.stop();

	gdist.template deviceToHost<0,1>();

	float elapsedTime = ts.getwct();

    auto numElements = gdist.size_local_inserted();
    size_t gridSize = sz[0]*sz[1]*sz[2];
    auto occupancyEmp = ((float) numElements) / gridSize;
    auto processingRate = 1e-9*numElements*2*iterations/elapsedTime; // In GElem/s
    auto gflops = 9*processingRate; // GFlops/s

	print_results_stencil(gridSize, 1U, occupancyEmp, numElements, elapsedTime, processingRate, gflops);
}
BOOST_AUTO_TEST_CASE(stencil_full_2x2x2_2x2x2)
{
	print_test_name();
	test_3D_stencil_full<2,2,512>();
}
BOOST_AUTO_TEST_CASE(stencil_full_4x4x4_4x4x4)
{
	print_test_name();
	test_3D_stencil_full<4,4,512>();
}
BOOST_AUTO_TEST_CASE(stencil_full_8x8x8_8x8x8)
{
	print_test_name();
	test_3D_stencil_full<8,8,512>();
}
// BOOST_AUTO_TEST_CASE(stencil_full_16x16x16_16x16x16)
// {
// 	print_test_name();
// 	test_3D_stencil_full<16,16,512>();
// }
// BOOST_AUTO_TEST_CASE(stencil_full_32x32x32_32x32x32)
// {
// 	print_test_name();
// 	test_3D_stencil_full<32,32,512>();
// }

BOOST_AUTO_TEST_SUITE_END() //dim_3D
BOOST_AUTO_TEST_SUITE_END() //SparseGridGpu_dist_single
BOOST_AUTO_TEST_SUITE_END() //performance

//eof
