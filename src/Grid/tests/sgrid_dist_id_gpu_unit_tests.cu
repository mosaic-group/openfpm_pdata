#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include "Grid/grid_dist_id.hpp"

BOOST_AUTO_TEST_SUITE( sgrid_gpu_test_suite )

template<unsigned int p>
struct insert_kernel
{
	template<typename SparseGridGpu_type, typename ite_type>
	__device__ void operator()(SparseGridGpu_type & sg, ite_type & ite, float c)
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
	    if (SparseGridGpu_type::d >= 3 && z+ite.start.get(2) > ite.stop.get(2))
	    {return;}

	    grid_key_dx<SparseGridGpu_type::d, size_t> coord;
	    grid_key_dx<SparseGridGpu_type::d, size_t> coord_glob;

	    if (SparseGridGpu_type::d >= 2)
	    {
	    	coord.set_d(0,x+ite.start.get(0));
	    	coord_glob.set_d(0,x+ite.start.get(0)+ite.origin.get(0));
	    	coord.set_d(1,y+ite.start.get(1));
	    	coord_glob.set_d(1,y+ite.start.get(1)+ite.origin.get(1));
	    }
	    else if (SparseGridGpu_type::d >= 3)
	    {
		    coord.set_d(0,x+ite.start.get(0));
		    coord_glob.set_d(0,x+ite.start.get(0)+ite.origin.get(0));
		    coord.set_d(1,y+ite.start.get(1));
		    coord_glob.set_d(1,y+ite.start.get(1)+ite.origin.get(1));
		    coord.set_d(2,z+ite.start.get(2));
		    coord_glob.set_d(2,z+ite.start.get(2)+ite.origin.get(2));
	    }


	    if (SparseGridGpu_type::d >= 2)
	    {sg.template insert<p>(coord) = c + coord_glob.get(0) + coord_glob.get(1);}
	    else
	    {sg.template insert<p>(coord) = c + coord_glob.get(0) + coord_glob.get(1) + coord_glob.get(2);}

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

	/////// GPU insert + flush

	Box<2,size_t> box({1,1},{1,1});
	auto it = gdist.getGridIterator(box.getKP1(),box.getKP2());

	/////// GPU Run kernel

	gdist.setInsertBuffer(1);

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
				BOOST_REQUIRE_EQUAL(gdist.template get<0>(p), 7.0);
			}
			else
			{
				if (gdist.template get<0>(p) != 666.0)
				{
					float f = gdist.template get<0>(p);
					std::cout << "ERROR: " << gdist.template get<0>(p) << std::endl;
				}

				BOOST_REQUIRE_EQUAL(gdist.template get<0>(p), 666.0);
			}

			++it;
		}
	}

	return;

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

	gdist.template iterateGridGPU<insert_kernel<0>>(it,c);
	gdist.template flush<smax_<0>>(flush_type::FLUSH_ON_DEVICE);

	gdist.template deviceToHost<0>();

	gdist.write("sgrid_gpu_output");

	std::string file_test("sgrid_gpu_output_" + std::to_string(v_cl.size()) + "_" + std::to_string(v_cl.rank())  + ".vtk");
	std::string file("sgrid_gpu_output_" + std::to_string(v_cl.rank()) + ".vtk");

	bool test = compare(file,"test_data/" + file_test);

	BOOST_REQUIRE_EQUAL(true,test);
}

BOOST_AUTO_TEST_SUITE_END()
