#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include "Grid/grid_dist_id.hpp"

BOOST_AUTO_TEST_SUITE( sgrid_gpu_test_suite )

template<unsigned int p>
struct insert_kernel2D
{
	template<typename SparseGridGpu_type, typename ite_type>
	__device__ void operator()(SparseGridGpu_type & sg, ite_type & ite, float c)
	{
		GRID_ID_2_GLOBAL(ite);

	    sg.init();

	    sg.template insert<p>(key) = c + keyg.get(0) + keyg.get(1);

	    __syncthreads();

	    sg.flush_block_insert();
	}
};

template<unsigned int p>
struct insert_kernel3D
{
	template<typename SparseGridGpu_type, typename ite_type>
	__device__ void operator()(SparseGridGpu_type & sg, ite_type & ite, float c)
	{
		GRID_ID_3_GLOBAL(ite);

	    sg.init();

	    sg.template insert<p>(key) = c + keyg.get(0) + keyg.get(1) + keyg.get(2);

	    __syncthreads();

	    sg.flush_block_insert();
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

	/////// GPU insert + flush

	Box<2,size_t> box({1,1},{1,1});
	auto it = gdist.getGridIterator(box.getKP1(),box.getKP2());

	/////// GPU Run kernel

	gdist.setInsertBuffer(1);

	float c = 5.0;

	gdist.template iterateGridGPU<insert_kernel2D<0>>(it,c);
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
				BOOST_REQUIRE_EQUAL(gdist.template get<0>(p), 7.0);
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

	gdist.template iterateGridGPU<insert_kernel2D<0>>(it3,c);
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
				BOOST_REQUIRE_EQUAL(gdist.template get<0>(p), 7.0);
			}
			else if (box3.isInside(p2_))
			{
				float tst = c + p2.get(0) + p2.get(1);

				BOOST_REQUIRE_EQUAL(gdist.template get<0>(p), tst);
			}
			else
			{
				BOOST_REQUIRE_EQUAL(gdist.template get<0>(p), 666.0);
			}

			++it;
		}
	}
}


BOOST_AUTO_TEST_CASE( sgrid_gpu_test_output )
{
	auto & v_cl = create_vcluster();

	if (v_cl.size() > 3){return;}

	size_t sz[2] = {17,17};
	periodicity<2> bc = {PERIODIC,PERIODIC};

	Ghost<2,long int> g(1);

	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	sgrid_dist_id_gpu<2,float,aggregate<float>> gdist(sz,domain,g,bc);

	gdist.template setBackgroundValue<0>(666);

	/////// GPU insert + flush

	Box<2,size_t> box({1,1},{15,15});
	auto it = gdist.getGridIterator(box.getKP1(),box.getKP2());

	/////// GPU Run kernel

	gdist.setInsertBuffer(128);

	float c = 5.0;

	gdist.template iterateGridGPU<insert_kernel2D<0>>(it,c);
	gdist.template flush<smax_<0>>(flush_type::FLUSH_ON_DEVICE);

	gdist.template deviceToHost<0>();

	gdist.write("sgrid_gpu_output");

	std::string file_test("sgrid_gpu_output_" + std::to_string(v_cl.size()) + "_" + std::to_string(v_cl.rank())  + ".vtk");
	std::string file("sgrid_gpu_output_" + std::to_string(v_cl.rank()) + ".vtk");

	bool test = compare(file,"test_data/" + file_test);

	BOOST_REQUIRE_EQUAL(true,test);
}

BOOST_AUTO_TEST_CASE( sgrid_gpu_test_ghost_get )
{
	size_t sz[2] = {17,17};
	periodicity<2> bc = {PERIODIC,PERIODIC};

	Ghost<2,long int> g(1);

	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	sgrid_dist_id_gpu<2,float,aggregate<float>> gdist(sz,domain,g,bc);

	gdist.template setBackgroundValue<0>(666);

	/////// GPU insert + flush

	Box<2,size_t> box({1,1},{15,15});
	auto it = gdist.getGridIterator(box.getKP1(),box.getKP2());

	/////// GPU Run kernel

	gdist.setInsertBuffer(225);

	float c = 5.0;

	gdist.template iterateGridGPU<insert_kernel2D<0>>(it,c);
	gdist.template flush<smax_<0>>(flush_type::FLUSH_ON_DEVICE);

	gdist.template deviceToHost<0>();
//	gdist.write("broken");

//	gdist.template ghost_get<0>(RUN_ON_DEVICE);
}



BOOST_AUTO_TEST_SUITE_END()
