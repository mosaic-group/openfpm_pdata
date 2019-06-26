#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include "Grid/grid_dist_id.hpp"

BOOST_AUTO_TEST_SUITE( sgrid_gpu_test_suite )

template<unsigned int p>
struct insert_kernel
{
	template<typename SparseGridGpu_type>
	__device__ void operator()(SparseGridGpu_type & sg, ite_gpu<SparseGridGpu_type::d> & ite, float c)
	{
	    sg.init();

	    const auto bDimX = blockDim.x;
	    const auto bDimY = blockDim.y;
	    const auto bDimZ = blockDim.z;
	    const auto bIdX = blockIdx.x;
	    const auto bIdY = blockIdx.y;
	    const auto bIdZ = blockIdx.z;
	    const auto tIdX = threadIdx.x;
	    const auto tIdY = threadIdx.y;
	    const auto tIdZ = threadIdx.z;
	    int x = bIdX * bDimX + tIdX;
	    int y = bIdY * bDimY + tIdY;
	    int z = bIdZ * bDimZ + tIdZ;

	    if (x+ite.start.get(0) > ite.stop.get(0))
	    {return;}
	    if (SparseGridGpu_type::d >= 2 && y+ite.start.get(1) > ite.stop.get(1))
	    {return;}
	    if (SparseGridGpu_type::d >= 3 && z+ite.start.get(1) > ite.stop.get(2))
	    {return;}

	    grid_key_dx<SparseGridGpu_type::d, size_t> coord({x+ite.start.get(0), y+ite.start.get(1), z+ite.start.get(2)});

	//    size_t pos = sg.getLinId(coord);
	//    printf("insertValues: bDim=(%d,%d), bId=(%d,%d), tId=(%d,%d) : "
	//           "pos=%ld, coord={%d,%d}, value=%d\n",
	//           bDimX, bDimY,
	//           bIdX, bIdY,
	//           tIdX, tIdY,
	//           pos,
	//           x, y,
	//           x); //debug

	    sg.template insert<p>(coord) = c;

	    __syncthreads();

	    sg.flush_block_insert();

	    // Compiler avoid warning
	    y++;
	    z++;
	}
};

template<unsigned int p>
struct stencil_kernel
{
	template<typename SparseGridGpu_type>
	__device__ void operator()(SparseGridGpu_type & sg, ite_gpu<SparseGridGpu_type::d> & ite, float c)
	{
		// TODO
	}
};

BOOST_AUTO_TEST_CASE( sgrid_gpu_test_base )
{
	size_t sz[2] = {17,17};
	periodicity<2> bc = {PERIODIC,PERIODIC};

	Ghost<2,long int> g(1);

	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	sgrid_dist_id_gpu<2,float,aggregate<float>> gdist(sz,domain,g,bc);

	gdist.template setBackgroundValue<0>(666);

	/////// CPU insert

/*	auto it = gdist.getGridIterator(box.getKP1(),box.getKP2());

	while (it.isNext())
	{
		auto p = it.get_dist();

		gdist.template insert<0>(p) = 1.0;

		++it;
	}

	gdist.template flush<>();

	Box<2,size_t> box2({0,0},{15,15});
	auto it2 = gdist.getGridIterator(box2.getKP1(),box2.getKP2());

	while (it2.isNext())
	{
		auto p = it2.get_dist();

		std::cout << gdist.template get<0>(p) << std::endl;

		++it2;
	}*/

	/////// host to device

	/////// GPU insert + flush

	Box<2,size_t> box({1,1},{1,1});
	auto it = gdist.getGridIterator(box.getKP1(),box.getKP2());

	/////// GPU Run kernel

	gdist.setInsertBuffer(128);

	float c = 5.0;

	gdist.template iterateGridGPU<insert_kernel<0>>(it,c);
	gdist.template flush<smax_<0>>(flush_type::FLUSH_ON_DEVICE);

	gdist.template deviceToHost<0>();

	{
		Box<2,size_t> box2({0,0},{15,15});

		auto it = gdist.getGridIterator(box2.getKP1(),box2.getKP2());

		while (it.isNext())
		{
			auto p = it.get_dist();
			auto p2 = it.get();

			if (p2.get(0) == box.getLow(0) && p2.get(1) == box.getLow(1))
			{
				BOOST_REQUIRE_EQUAL(gdist.template get<0>(p), 5.0);
			}
			else
			{
				BOOST_REQUIRE_EQUAL(gdist.template get<0>(p), 666.0);
			}

			++it;
		}
	}

	//

	c = 3.0;

	Box<2,size_t> box3({3,3},{11,11});

	auto it3 = gdist.getGridIterator(box3.getKP1(),box3.getKP2());
	gdist.setInsertBuffer(128);

	gdist.template iterateGridGPU<insert_kernel<0>>(it3,c);
	gdist.template flush<smax_<0>>(flush_type::FLUSH_ON_DEVICE);
	gdist.template deviceToHost<0>();

	{
		Box<2,size_t> box2({0,0},{15,15});

		auto it = gdist.getGridIterator(box2.getKP1(),box2.getKP2());

		while (it.isNext())
		{
			auto p = it.get_dist();
			auto p2 = it.get();

			Point<2,size_t> p2_ = p2.toPoint();

			if (box.isInside(p2_))
			{
				BOOST_REQUIRE_EQUAL(gdist.template get<0>(p), 5.0);
			}
			else if (box3.isInside(p2_))
			{
				BOOST_REQUIRE_EQUAL(gdist.template get<0>(p), 3.0);
			}
			else
			{
				BOOST_REQUIRE_EQUAL(gdist.template get<0>(p), 666.0);
			}

			++it;
		}
	}

	////////////////////////////////////

	gdist.setInsertBuffer(128);
	gdist.template iterateGPU<stencil_kernel<0>>();


}


BOOST_AUTO_TEST_SUITE_END()
