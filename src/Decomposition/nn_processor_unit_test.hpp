/*
 * nn_processor_unit_test.hpp
 *
 *  Created on: Dec 16, 2015
 *      Author: i-bird
 */

#ifndef SRC_DECOMPOSITION_NN_PROCESSOR_UNIT_TEST_HPP_
#define SRC_DECOMPOSITION_NN_PROCESSOR_UNIT_TEST_HPP_

#include "VCluster.hpp"

void create_decomposition2x2(openfpm::vector<openfpm::vector<long unsigned int>> & box_nn_processor, openfpm::vector<SpaceBox<2,float>> & sub_domains)
{
	Vcluster & v_cl = create_vcluster();

	box_nn_processor.add();

	if (v_cl.getProcessUnitID() == 0)
	{
		box_nn_processor.get(0).add(1);
		box_nn_processor.get(0).add(2);
		box_nn_processor.get(0).add(3);

		sub_domains.add(Box<2,float>({0.0,0.0},{0.5,0.5}));
	}
	else if (v_cl.getProcessUnitID() == 1)
	{
		box_nn_processor.get(0).add(0);
		box_nn_processor.get(0).add(2);
		box_nn_processor.get(0).add(3);

		sub_domains.add(Box<2,float>({0.5,0.0},{1.0,0.5}));
	}
	else if (v_cl.getProcessUnitID() == 2)
	{
		box_nn_processor.get(0).add(1);
		box_nn_processor.get(0).add(0);
		box_nn_processor.get(0).add(3);

		sub_domains.add(Box<2,float>({0.0,0.5},{0.5,1.0}));
	}
	else if (v_cl.getProcessUnitID() == 3)
	{
		box_nn_processor.get(0).add(1);
		box_nn_processor.get(0).add(2);
		box_nn_processor.get(0).add(0);

		sub_domains.add(Box<2,float>({0.5,0.5},{1.0,1.0}));
	}
}

BOOST_AUTO_TEST_SUITE( nn_processor_test )

BOOST_AUTO_TEST_CASE( nn_processor_np_test)
{
	Vcluster & v_cl = create_vcluster();

	/*!
	 *
	 * We test this situation
	 *
	 * \verbatim
		+-------+-------+
		|       |       |
		|   0   |   1   |
		|       |       |
		|       |       |
		+---------------+
		|       |       |
		|   2   |   3   |
		|       |       |
		|       |       |
		+-------+-------+

	 * \endverbatim
	 *
	 *
	 */

	if (v_cl.getProcessingUnits() != 4)
		return;

	openfpm::vector<openfpm::vector<long unsigned int>> box_nn_processor;
	openfpm::vector<SpaceBox<2,float>> sub_domains;

	create_decomposition2x2(box_nn_processor,sub_domains);

	nn_prcs<2,float> nnp(v_cl);
	nnp.create(box_nn_processor, sub_domains);

	BOOST_REQUIRE_EQUAL(nnp.getNNProcessors(),3ul);

	if (v_cl.getProcessUnitID() == 0)
	{
		BOOST_REQUIRE_EQUAL(nnp.getNRealSubdomains(1),1);
		BOOST_REQUIRE_EQUAL(nnp.getNRealSubdomains(2),1);
		BOOST_REQUIRE_EQUAL(nnp.getNRealSubdomains(3),1);

		const openfpm::vector< ::Box<2,float> > & nsubs1 = nnp.getNearSubdomains(1);
		const openfpm::vector< ::Box<2,float> > & nsubs2 = nnp.getNearSubdomains(2);
		const openfpm::vector< ::Box<2,float> > & nsubs3 = nnp.getNearSubdomains(3);

		SpaceBox<2,float> b1_a = nsubs1.get(0);
		SpaceBox<2,float> b2_a = nsubs2.get(0);
		SpaceBox<2,float> b3_a = nsubs3.get(0);


		SpaceBox<2,float> b1_b = Box<2,float>({0.5,0.0},{1.0,0.5});
		SpaceBox<2,float> b2_b = Box<2,float>({0.0,0.5},{0.5,1.0});
		SpaceBox<2,float> b3_b = Box<2,float>({0.5,0.5},{1.0,1.0});

		bool ret1 = b1_a == b1_b;
		bool ret2 = b2_a == b2_b;
		bool ret3 = b3_a == b3_b;

		BOOST_REQUIRE_EQUAL(ret1,true);
		BOOST_REQUIRE_EQUAL(ret2,true);
		BOOST_REQUIRE_EQUAL(ret3,true);
	}
	else if (v_cl.getProcessUnitID() == 1)
	{
		BOOST_REQUIRE_EQUAL(nnp.getNRealSubdomains(0),1);
		BOOST_REQUIRE_EQUAL(nnp.getNRealSubdomains(2),1);
		BOOST_REQUIRE_EQUAL(nnp.getNRealSubdomains(3),1);

		const openfpm::vector< ::Box<2,float> > & nsubs1 = nnp.getNearSubdomains(0);
		const openfpm::vector< ::Box<2,float> > & nsubs2 = nnp.getNearSubdomains(2);
		const openfpm::vector< ::Box<2,float> > & nsubs3 = nnp.getNearSubdomains(3);

		SpaceBox<2,float> b1_a = nsubs1.get(0);
		SpaceBox<2,float> b2_a = nsubs2.get(0);
		SpaceBox<2,float> b3_a = nsubs3.get(0);


		SpaceBox<2,float> b1_b = Box<2,float>({0.0,0.0},{0.5,0.5});
		SpaceBox<2,float> b2_b = Box<2,float>({0.0,0.5},{0.5,1.0});
		SpaceBox<2,float> b3_b = Box<2,float>({0.5,0.5},{1.0,1.0});

		bool ret1 = b1_a == b1_b;
		bool ret2 = b2_a == b2_b;
		bool ret3 = b3_a == b3_b;

		BOOST_REQUIRE_EQUAL(ret1,true);
		BOOST_REQUIRE_EQUAL(ret2,true);
		BOOST_REQUIRE_EQUAL(ret3,true);
	}
	else if (v_cl.getProcessUnitID() == 2)
	{
		BOOST_REQUIRE_EQUAL(nnp.getNRealSubdomains(1),1);
		BOOST_REQUIRE_EQUAL(nnp.getNRealSubdomains(0),1);
		BOOST_REQUIRE_EQUAL(nnp.getNRealSubdomains(3),1);

		const openfpm::vector< ::Box<2,float> > & nsubs1 = nnp.getNearSubdomains(1);
		const openfpm::vector< ::Box<2,float> > & nsubs2 = nnp.getNearSubdomains(0);
		const openfpm::vector< ::Box<2,float> > & nsubs3 = nnp.getNearSubdomains(3);

		SpaceBox<2,float> b1_a = nsubs1.get(0);
		SpaceBox<2,float> b2_a = nsubs2.get(0);
		SpaceBox<2,float> b3_a = nsubs3.get(0);

		SpaceBox<2,float> b1_b = Box<2,float>({0.5,0.0},{1.0,0.5});
		SpaceBox<2,float> b2_b = Box<2,float>({0.0,0.0},{0.5,0.5});
		SpaceBox<2,float> b3_b = Box<2,float>({0.5,0.5},{1.0,1.0});

		bool ret1 = b1_a == b1_b;
		bool ret2 = b2_a == b2_b;
		bool ret3 = b3_a == b3_b;

		BOOST_REQUIRE_EQUAL(ret1,true);
		BOOST_REQUIRE_EQUAL(ret2,true);
		BOOST_REQUIRE_EQUAL(ret3,true);
	}
	else if (v_cl.getProcessUnitID() == 3)
	{
		BOOST_REQUIRE_EQUAL(nnp.getNRealSubdomains(0),1);
		BOOST_REQUIRE_EQUAL(nnp.getNRealSubdomains(1),1);
		BOOST_REQUIRE_EQUAL(nnp.getNRealSubdomains(2),1);

		const openfpm::vector< ::Box<2,float> > & nsubs1 = nnp.getNearSubdomains(0);
		const openfpm::vector< ::Box<2,float> > & nsubs2 = nnp.getNearSubdomains(1);
		const openfpm::vector< ::Box<2,float> > & nsubs3 = nnp.getNearSubdomains(2);

		SpaceBox<2,float> b1_a = nsubs1.get(0);
		SpaceBox<2,float> b2_a = nsubs2.get(0);
		SpaceBox<2,float> b3_a = nsubs3.get(0);

		SpaceBox<2,float> b1_b = Box<2,float>({0.0,0.0},{0.5,0.5});
		SpaceBox<2,float> b2_b = Box<2,float>({0.5,0.0},{1.0,0.5});
		SpaceBox<2,float> b3_b = Box<2,float>({0.0,0.5},{0.5,1.0});

		bool ret1 = b1_a == b1_b;
		bool ret2 = b2_a == b2_b;
		bool ret3 = b3_a == b3_b;

		BOOST_REQUIRE_EQUAL(ret1,true);
		BOOST_REQUIRE_EQUAL(ret2,true);
		BOOST_REQUIRE_EQUAL(ret3,true);
	}
}

BOOST_AUTO_TEST_CASE( nn_processor_box_periodic_test)
{
	// Vcluster
	Vcluster & v_cl = create_vcluster();

	/*!
	 *
	 * We test this situation
	 *
	 * \verbatim
		+-------+-------+
		|       |       |
		|   0   |   1   |
		|       |       |
		|       |       |
		+---------------+
		|       |       |
		|   2   |   3   |
		|       |       |
		|       |       |
		+-------+-------+

	 * \endverbatim
	 *
	 *
	 */

	if (v_cl.getProcessingUnits() != 4)
		return;

	Box<2,float> domain({0.0,0.0},{1.0,1.0});
	const size_t bc[2] = {PERIODIC,PERIODIC};

	Ghost<2,float> ghost(0.01);

	openfpm::vector<openfpm::vector<long unsigned int>> box_nn_processor;
	openfpm::vector<SpaceBox<2,float>> sub_domains;

	create_decomposition2x2(box_nn_processor,sub_domains);

	//////////////

	nn_prcs<2,float> nnp(v_cl);
	nnp.create(box_nn_processor, sub_domains);

	// check that nn_processor contain the correct boxes

	nnp.applyBC(domain,ghost,bc);

	if (v_cl.getProcessUnitID() == 0)
	{
		BOOST_REQUIRE_EQUAL(nnp.getNearSubdomains(1).size(),4);
		BOOST_REQUIRE_EQUAL(nnp.getNearSubdomains(2).size(),4);
		BOOST_REQUIRE_EQUAL(nnp.getNearSubdomains(3).size(),4);

		openfpm::vector<Box<2,float>> bv;

		bv.add(Box<2,float>({0.5,0},{1.0,0.5}));
		bv.add(Box<2,float>({-0.5,0.0},{0.0,0.5}));
		bv.add(Box<2,float>({0.5,1.0},{1.0,1.5}));
		bv.add(Box<2,float>({-0.5,1.0},{0.0,1.5}));

		bool ret = nnp.getNearSubdomains(1) == bv;
		BOOST_REQUIRE_EQUAL(ret,true);

		bv.clear();

		bv.add(Box<2,float>({0.0,0.5},{0.5,1.0}));
		bv.add(Box<2,float>({1.0,0.5},{1.5,1.0}));
		bv.add(Box<2,float>({0.0,-0.5},{0.5,0.0}));
		bv.add(Box<2,float>({1.0,-0.5},{1.5,0.0}));

		ret = nnp.getNearSubdomains(2) == bv;
		BOOST_REQUIRE_EQUAL(ret,true);

		bv.clear();

		bv.add(Box<2,float>({0.5,0.5},{1.0,1.0}));
		bv.add(Box<2,float>({-0.5,0.5},{0.0,1.0}));
		bv.add(Box<2,float>({0.5,-0.5},{1.0,0.0}));
		bv.add(Box<2,float>({-0.5,-0.5},{0.0,0.0}));

		ret = nnp.getNearSubdomains(3) == bv;
		BOOST_REQUIRE_EQUAL(ret,true);
	}
	else if (v_cl.getProcessUnitID() == 1)
	{
		BOOST_REQUIRE_EQUAL(nnp.getNearSubdomains(0).size(),4);
		BOOST_REQUIRE_EQUAL(nnp.getNearSubdomains(2).size(),4);
		BOOST_REQUIRE_EQUAL(nnp.getNearSubdomains(3).size(),4);

		openfpm::vector<Box<2,float>> bv;

		bv.add(Box<2,float>({0.0,0},{0.5,0.5}));
		bv.add(Box<2,float>({1.0,0.0},{1.5,0.5}));
		bv.add(Box<2,float>({0.0,1.0},{0.5,1.5}));
		bv.add(Box<2,float>({1.0,1.0},{1.5,1.5}));

		bool ret = nnp.getNearSubdomains(0) == bv;
		BOOST_REQUIRE_EQUAL(ret,true);

		bv.clear();

		bv.add(Box<2,float>({0.0,0.5},{0.5,1.0}));
		bv.add(Box<2,float>({1.0,0.5},{1.5,1.0}));
		bv.add(Box<2,float>({0.0,-0.5},{0.5,0.0}));
		bv.add(Box<2,float>({1.0,-0.5},{1.5,0.0}));

		ret = nnp.getNearSubdomains(2) == bv;
		BOOST_REQUIRE_EQUAL(ret,true);

		bv.clear();

		bv.add(Box<2,float>({0.5,0.5},{1.0,1.0}));
		bv.add(Box<2,float>({-0.5,0.5},{0.0,1.0}));
		bv.add(Box<2,float>({0.5,-0.5},{1.0,0.0}));
		bv.add(Box<2,float>({-0.5,-0.5},{0.0,0.0}));

		ret = nnp.getNearSubdomains(3) == bv;
		BOOST_REQUIRE_EQUAL(ret,true);
	}
	else if (v_cl.getProcessUnitID() == 2)
	{
		BOOST_REQUIRE_EQUAL(nnp.getNearSubdomains(0).size(),4);
		BOOST_REQUIRE_EQUAL(nnp.getNearSubdomains(1).size(),4);
		BOOST_REQUIRE_EQUAL(nnp.getNearSubdomains(3).size(),4);

		openfpm::vector<Box<2,float>> bv;

		bv.add(Box<2,float>({0.0,0},{0.5,0.5}));
		bv.add(Box<2,float>({1.0,0.0},{1.5,0.5}));
		bv.add(Box<2,float>({0.0,1.0},{0.5,1.5}));
		bv.add(Box<2,float>({1.0,1.0},{1.5,1.5}));

		bool ret = nnp.getNearSubdomains(0) == bv;
		BOOST_REQUIRE_EQUAL(ret,true);

		bv.clear();

		bv.add(Box<2,float>({0.5,0},{1.0,0.5}));
		bv.add(Box<2,float>({-0.5,0.0},{0.0,0.5}));
		bv.add(Box<2,float>({0.5,1.0},{1.0,1.5}));
		bv.add(Box<2,float>({-0.5,1.0},{0.0,1.5}));

		ret = nnp.getNearSubdomains(1) == bv;
		BOOST_REQUIRE_EQUAL(ret,true);

		bv.clear();

		bv.add(Box<2,float>({0.5,0.5},{1.0,1.0}));
		bv.add(Box<2,float>({-0.5,0.5},{0.0,1.0}));
		bv.add(Box<2,float>({0.5,-0.5},{1.0,0.0}));
		bv.add(Box<2,float>({-0.5,-0.5},{0.0,0.0}));

		ret = nnp.getNearSubdomains(3) == bv;
		BOOST_REQUIRE_EQUAL(ret,true);
	}
	else if (v_cl.getProcessUnitID() == 3)
	{
		BOOST_REQUIRE_EQUAL(nnp.getNearSubdomains(0).size(),4);
		BOOST_REQUIRE_EQUAL(nnp.getNearSubdomains(1).size(),4);
		BOOST_REQUIRE_EQUAL(nnp.getNearSubdomains(2).size(),4);

		openfpm::vector<Box<2,float>> bv;

		bv.add(Box<2,float>({0.0,0},{0.5,0.5}));
		bv.add(Box<2,float>({1.0,0.0},{1.5,0.5}));
		bv.add(Box<2,float>({0.0,1.0},{0.5,1.5}));
		bv.add(Box<2,float>({1.0,1.0},{1.5,1.5}));

		bool ret = nnp.getNearSubdomains(0) == bv;
		BOOST_REQUIRE_EQUAL(ret,true);

		bv.clear();

		bv.add(Box<2,float>({0.5,0},{1.0,0.5}));
		bv.add(Box<2,float>({-0.5,0.0},{0.0,0.5}));
		bv.add(Box<2,float>({0.5,1.0},{1.0,1.5}));
		bv.add(Box<2,float>({-0.5,1.0},{0.0,1.5}));

		ret = nnp.getNearSubdomains(1) == bv;
		BOOST_REQUIRE_EQUAL(ret,true);

		bv.clear();

		bv.add(Box<2,float>({0.0,0.5},{0.5,1.0}));
		bv.add(Box<2,float>({1.0,0.5},{1.5,1.0}));
		bv.add(Box<2,float>({0.0,-0.5},{0.5,0.0}));
		bv.add(Box<2,float>({1.0,-0.5},{1.5,0.0}));

/*		for (size_t i = 0 ; i < nnp.getNearSubdomains(2).size() ; i++)
		{
			Box<2,float> b = nnp.getNearSubdomains(2).get(i);
			std::cout << "BOX: " << b.toString() << std::endl;
		}*/

		ret = nnp.getNearSubdomains(2) == bv;
		BOOST_REQUIRE_EQUAL(ret,true);
	}
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* SRC_DECOMPOSITION_NN_PROCESSOR_UNIT_TEST_HPP_ */
