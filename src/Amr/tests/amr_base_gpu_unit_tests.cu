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

BOOST_AUTO_TEST_CASE( grid_dist_id_amr_gpu_link_test )
{
	auto & v_cl = create_vcluster();

	// Domain
	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	Ghost<2,long int> g(1);
	sgrid_dist_amr_gpu<2,float,aggregate<float>> amr_g(domain,g);

	size_t g_sz[2] = {17,17};

	size_t n_lvl = 3;

	amr_g.initLevels(n_lvl,g_sz);

	grid_key_dx<2> start({5,5});
	grid_key_dx<2> start_lvl_dw({9,9});
	grid_key_dx<2> stop_lvl_dw({12,12});
	grid_key_dx<2> start_lvl_dw2({19,19});
	grid_key_dx<2> stop_lvl_dw2({23,23});

	auto it = amr_g.getGridIterator(0,start,start);
	auto it2 = amr_g.getGridIterator(1,start_lvl_dw,stop_lvl_dw);
	auto it3 = amr_g.getGridIterator(2,start_lvl_dw2,stop_lvl_dw2);
//	it.setGPUInsertBuffer(4);

	auto & lvl_0 = amr_g.getDistGrid(0);
	auto & lvl_1 = amr_g.getDistGrid(1);
	auto & lvl_2 = amr_g.getDistGrid(2);

	// Add points in level 0

	while (it.isNext())
	{
		auto key = it.get_dist();

		lvl_0.template insertFlush<0>(key) = 1.0;

		++it;
	}

	while (it2.isNext())
	{
		auto key = it2.get_dist();

		lvl_1.template insertFlush<0>(key) = 2.0;

		++it2;
	}

	while (it3.isNext())
	{
		auto key = it3.get_dist();

		lvl_2.template insertFlush<0>(key) = 3.0;

		++it3;
	}

	amr_g.hostToDevice<0>();
	amr_g.tagBoundaries<NNStar<2>>();
	amr_g.construct_level_connections();

	/////////////////////////////////////////////////////////////

	auto & lvl_zero_d = amr_g.getDistGrid(0);
	auto & lvl_one_d = amr_g.getDistGrid(1);
	auto & lvl_two_d = amr_g.getDistGrid(2);

	// For each local grid

	for (int i = 0 ; i < lvl_zero_d.getN_loc_grid() ; i++)
	{

		// Check
		auto & lvl_zero = lvl_zero_d.get_loc_grid(i);
		auto & lvl_one = lvl_one_d.get_loc_grid(i);
		auto & lvl_two = lvl_two_d.get_loc_grid(i);

		auto & offs_dw_link = lvl_zero.getDownLinksOffsets();
		auto & dw_links = lvl_zero.getDownLinks();

		BOOST_REQUIRE_EQUAL(offs_dw_link.size(),1);
		BOOST_REQUIRE_EQUAL(dw_links.size(),4);

		auto & indexL0 = lvl_zero.private_get_blockMap().getIndexBuffer();
		auto & indexL1 = lvl_one.private_get_blockMap().getIndexBuffer();
		auto & indexL2 = lvl_two.private_get_blockMap().getIndexBuffer();

		auto & dataL0 = lvl_zero.private_get_blockMap().getDataBuffer();
		auto & dataL1 = lvl_one.private_get_blockMap().getDataBuffer();
		auto & dataL2 = lvl_two.private_get_blockMap().getDataBuffer();

		dw_links.template deviceToHost<0,1>();

		BOOST_REQUIRE_EQUAL(dataL1.template get<0>(dw_links.template get<0>(0))[dw_links.template get<1>(0)],2);
		BOOST_REQUIRE_EQUAL(dataL1.template get<0>(dw_links.template get<0>(1))[dw_links.template get<1>(1)],2);
		BOOST_REQUIRE_EQUAL(dataL1.template get<0>(dw_links.template get<0>(2))[dw_links.template get<1>(2)],2);
		BOOST_REQUIRE_EQUAL(dataL1.template get<0>(dw_links.template get<0>(3))[dw_links.template get<1>(3)],2);

		auto & offs_dw_link_1 = lvl_one.getDownLinksOffsets();
		auto & dw_links_1 = lvl_one.getDownLinks();

		BOOST_REQUIRE_EQUAL(offs_dw_link_1.size(),12);
		BOOST_REQUIRE_EQUAL(dw_links_1.size(),9);

		dw_links_1.template deviceToHost<0,1>();

		BOOST_REQUIRE_EQUAL(dataL2.template get<0>(dw_links_1.template get<0>(0))[dw_links_1.template get<1>(0)],3);
		BOOST_REQUIRE_EQUAL(dataL2.template get<0>(dw_links_1.template get<0>(1))[dw_links_1.template get<1>(1)],3);
		BOOST_REQUIRE_EQUAL(dataL2.template get<0>(dw_links_1.template get<0>(2))[dw_links_1.template get<1>(2)],3);
		BOOST_REQUIRE_EQUAL(dataL2.template get<0>(dw_links_1.template get<0>(3))[dw_links_1.template get<1>(3)],3);
		BOOST_REQUIRE_EQUAL(dataL2.template get<0>(dw_links_1.template get<0>(4))[dw_links_1.template get<1>(4)],3);
		BOOST_REQUIRE_EQUAL(dataL2.template get<0>(dw_links_1.template get<0>(5))[dw_links_1.template get<1>(5)],3);
		BOOST_REQUIRE_EQUAL(dataL2.template get<0>(dw_links_1.template get<0>(6))[dw_links_1.template get<1>(6)],3);
		BOOST_REQUIRE_EQUAL(dataL2.template get<0>(dw_links_1.template get<0>(7))[dw_links_1.template get<1>(7)],3);
		BOOST_REQUIRE_EQUAL(dataL2.template get<0>(dw_links_1.template get<0>(8))[dw_links_1.template get<1>(8)],3);
	}

/*	grid_key_dx<2> k({8,8});
	grid_key_dx<2> k2({16,16});

	lvl_zero.insertFlush<0>(k) = 1.0;
	lvl_one.insertFlush<0>(k2) = 5.0;

	lvl_one.template hostToDevice<0>();
	lvl_one.tagBoundaries(v_cl.getmgpuContext());*/

	/////////////////////////////////////////////////////////////


}

BOOST_AUTO_TEST_SUITE_END()
