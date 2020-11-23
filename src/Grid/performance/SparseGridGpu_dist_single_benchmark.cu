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

// Work plan:
// - Insert benchmarks
// - Stencil benchmarks
// - Stencil insert benchmarks

template <typename SgridT, typename BoxT, typename cT>
float insertConcentricSpheres(SgridT &gdist, 
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
				        				  );
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
float insertFullGrid(SgridT &gdist, 
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

BOOST_AUTO_TEST_CASE(insert_spheres)
{
	size_t sz[3] = {2*500,2*500,500};
	periodicity<3> bc = {PERIODIC,PERIODIC,PERIODIC};
	Ghost<3,long int> g(1);
	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});
	sgrid_dist_id_gpu<3,float,aggregate<float>> gdist(sz,domain,g,bc);
	gdist.template setBackgroundValue<0>(666);

	// Box<3,size_t> box({1,1,1},{sz[0]-1,sz[1]-1,sz[2]-1});
	Box<3,size_t> box({0,0,0},{sz[0]-1,sz[1]-1,sz[2]-1});

	// Insert the concentric spheres on GPU
	const float occupancy = 0.5;
	const unsigned int pitch = 32;
	const unsigned int iterations = 100;
    unsigned int c[3] = { sz[0]/2, sz[1]/2, sz[2]/2 };

    auto elapsedTime = insertConcentricSpheres(gdist, box, c, pitch, occupancy, iterations);

    auto numElements = gdist.size_local_inserted();
    size_t gridSize = sz[0]*sz[1]*sz[2];
    auto insertionRate = 1e-9*numElements*iterations/elapsedTime; // In MElem/s

	std::cout << "::: numElements: " << numElements << std::endl;
	std::cout << "::: gridSize: " << gridSize << std::endl;
	std::cout << "::: occupancy: " << ( (float) numElements ) / gridSize << std::endl;
	std::cout << "::: Time: " << elapsedTime << std::endl;
	std::cout << "::: Insertion rate: " << insertionRate << " GElem/s" << std::endl;
}

BOOST_AUTO_TEST_CASE(insert_full)
{
	size_t sz[3] = {2*500,2*500,500};
	periodicity<3> bc = {PERIODIC,PERIODIC,PERIODIC};
	Ghost<3,long int> g(1);
	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});
	sgrid_dist_id_gpu<3,float,aggregate<float>> gdist(sz,domain,g,bc);
	gdist.template setBackgroundValue<0>(666);

	// Box<3,size_t> box({1,1,1},{sz[0]-1,sz[1]-1,sz[2]-1});
	Box<3,size_t> box({0,0,0},{sz[0]-1,sz[1]-1,sz[2]-1});

	const unsigned int iterations = 100;

	// Insert full on GPU
    auto elapsedTime = insertFullGrid(gdist, box, iterations);

    auto numElements = gdist.size_local_inserted();
    size_t gridSize = sz[0]*sz[1]*sz[2];
    auto insertionRate = 1e-9*numElements*iterations/elapsedTime; // In MElem/s

	std::cout << "::: numElements: " << numElements << std::endl;
	std::cout << "::: gridSize: " << gridSize << std::endl;
	std::cout << "::: occupancy: " << ( (float) numElements ) / gridSize << std::endl;
	std::cout << "::: Time: " << elapsedTime << std::endl;
	std::cout << "::: Insertion rate: " << insertionRate << " GElem/s" << std::endl;
}

BOOST_AUTO_TEST_CASE(stencil_full)
{
	size_t sz[3] = {2*500,2*500,500};
	periodicity<3> bc = {PERIODIC,PERIODIC,PERIODIC};
	Ghost<3,long int> g(1);
	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});
	sgrid_dist_id_gpu<3,float,aggregate<float,float>> gdist(sz,domain,g,bc);
	gdist.template setBackgroundValue<0>(666);

	// Box<3,size_t> box({1,1,1},{sz[0]-1,sz[1]-1,sz[2]-1});
	Box<3,size_t> box({0,0,0},{sz[0]-1,sz[1]-1,sz[2]-1});

	// Insert full on GPU
    auto elapsedTime_insert = insertFullGrid(gdist, box);

    // Convolve a stencil
    // GetCpBlockType<GridType, property, stencilSize>
    // typedef typename GetCpBlockType<decltype(gdist),0,1>::type CpBlockType;

    const unsigned int iterations = 100;

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
	}

	cudaDeviceSynchronize();
	ts.stop();

	gdist.deviceToHost<0,1>();

	float elapsedTime = ts.getwct();

    auto numElements = gdist.size_local_inserted();
    size_t gridSize = sz[0]*sz[1]*sz[2];
    auto processingRate = 1e-9*numElements*2*iterations/elapsedTime; // In GElem/s
    auto gflops = 9*processingRate; // GFlops/s

	std::cout << "::: numElements: " << numElements << std::endl;
	std::cout << "::: gridSize: " << gridSize << std::endl;
	std::cout << "::: occupancy: " << ( (float) numElements ) / gridSize << std::endl;
	std::cout << "::: Time: " << elapsedTime << std::endl;
	std::cout << "::: Processing rate: " << processingRate << " GElem/s" << std::endl;
	std::cout << "::: Throughput: " << gflops << " GFlops/s" << std::endl;

}

BOOST_AUTO_TEST_SUITE_END() //SparseGridGpu_dist_single
BOOST_AUTO_TEST_SUITE_END() //performance

//eof
