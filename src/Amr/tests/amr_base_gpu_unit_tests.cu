/*
 * amr_base_gpu_unit_tests.cu
 *
 *  Created on: Aug 28, 2019
 *      Author: i-bird
 */

/*
 * amr_base_unit_test.cpp
 *
 *  Created on: Oct 5, 2017
 *      Author: i-bird
 */

#include <hip/hip_runtime.h>
#include "config.h"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Amr/grid_dist_amr.hpp"
#include "Point_test.hpp"
#include "Grid/tests/grid_dist_id_util_tests.hpp"

struct amr_launch_sparse
{
	template<typename grid_type, typename ite_type>
	__device__ void operator()(grid_type & grid, ite_type itg, float spacing, Point<3,float> center)
	{
		GRID_ID_3_GLOBAL(itg);

	    __shared__ bool is_block_empty;

	    if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0)
	    {is_block_empty = true;}

	    grid.init();

	    int offset = 0;
	    grid_key_dx<3,int> blk;
	    bool out = grid.getInsertBlockOffset(itg,key,blk,offset);

	    auto blockId = grid.getBlockLinId(blk);

	    const float x = keyg.get(0)*spacing - center.get(0);
	    const float y = keyg.get(1)*spacing - center.get(1);
	    const float z = keyg.get(2)*spacing - center.get(2);

	    float radius = sqrt((float) (x*x + y*y + z*z));

	    bool is_active = radius < 0.4 && radius > 0.3;

	    if (is_active == true)
	    {is_block_empty = false;}

	    __syncthreads();

	    if (is_block_empty == false)
	    {
	        auto ec = grid.insertBlock(blockId);

	        if ( is_active == true)
	        {
	            ec.template get<0>()[offset] = x+y+z;
	            ec.template get<grid_type::pMask>()[offset] = 1;
	        }
	    }

	    __syncthreads();

	    grid.flush_block_insert();
	}
};



BOOST_AUTO_TEST_SUITE( amr_grid_dist_id_test )


BOOST_AUTO_TEST_CASE( grid_dist_id_amr_gpu )
{
	// Domain
	Box<3,float> domain3({0.0,0.0,0.0},{1.0,1.0,1.0});


	Ghost<3,long int> g(1);
	sgrid_dist_amr_gpu<3,float,aggregate<float>> amr_g(domain3,g);

	size_t g_sz[3] = {4,4,4};

	size_t n_lvl = 6;

	amr_g.initLevels(n_lvl,g_sz);

	for (size_t i = 0 ; i < amr_g.getNLvl() ; i++)
	{
		// Fill the AMR with something

		size_t count = 0;

		auto it = amr_g.getGridIteratorGPU(i);
		it.setGPUInsertBuffer(1);

		Point<3,float> center({0.5,0.5,0.5});

		it.launch(amr_launch_sparse(),it.getSpacing(0),center);
		amr_g.getDistGrid(i).template flush<smax_<0>>(FLUSH_ON_DEVICE);

		amr_g.getDistGrid(i).template deviceToHost<0>();

		auto it2 = amr_g.getDistGrid(i).getDomainIterator();

		while (it2.isNext())
		{
			auto key = it2.get();
			auto keyg = it2.getGKey(key);

			count++;

			++it2;
		}

		auto & v_cl = create_vcluster();

		v_cl.sum(count);
		v_cl.execute();

		switch(i)
		{
		case 0:
			BOOST_REQUIRE_EQUAL(count,0);
			break;
		case 1:
			BOOST_REQUIRE_EQUAL(count,30);
			break;
		case 2:
			BOOST_REQUIRE_EQUAL(count,282);
			break;
		case 3:
			BOOST_REQUIRE_EQUAL(count,2192);
			break;
		case 4:
			BOOST_REQUIRE_EQUAL(count,16890);
			break;
		case 5:
			BOOST_REQUIRE_EQUAL(count,136992);
			break;
		}
	}

	// Iterate across all the levels initialized
/*	auto it = amr_g.getDomainIterator();

	size_t count = 0;

	while (it.isNext())
	{
		count++;

		++it;
	}

	Vcluster<> & v_cl = create_vcluster();

	v_cl.sum(count);
	v_cl.execute();

	BOOST_REQUIRE_EQUAL(count,correct_result);

	auto itc = amr_g.getDomainIteratorCells();

	size_t count_c = 0;

	while (itc.isNext())
	{
		count_c++;

		++itc;
	}

	v_cl.sum(count_c);
	v_cl.execute();

	auto it_level = amr_g.getDomainIteratorCells(3);

	while (it_level.isNext())
	{
		auto key = it_level.get();

		amr_g.template get<0>(3,key);

		++it_level;
	}

	BOOST_REQUIRE_EQUAL(count_c,correct_result_cell);*/
}

BOOST_AUTO_TEST_SUITE_END()
