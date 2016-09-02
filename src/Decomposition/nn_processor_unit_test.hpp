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

	BOOST_REQUIRE_EQUAL(nnp.getNNProcessors(),3);

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
	constexpr unsigned int dim = 3;
	typedef float T;

	Box<dim,T> domain({0.0,0.0,0.0},{1.0,1.0,1.0});
	Point<dim,T> middle({0.5,0.5,0.5});

	const size_t bc[dim] = {PERIODIC,PERIODIC,PERIODIC};

	// Vcluster
	Vcluster & v_cl = create_vcluster();

	Ghost<dim,T> ghost(0.01);

	//////////////

	nn_prcs<dim,T> nnp(v_cl);

/*	std::unordered_map<size_t, N_box<dim,T>> & nnp_sub = nnp.get_nn_processor_subdomains();
	openfpm::vector<size_t> & nnp_np = nnp.get_nn_processors();

	// we add the boxes

	size_t tot_n = 0;
	HyperCube<dim> hyp;

	for (long int i = dim-1 ; i >= 0 ; i--)
	{
		std::vector<comb<dim>> cmbs = hyp.getCombinations_R(i);

		for (size_t j = 0 ; j < cmbs.size() ; j++)
		{
			// Create a fake processor number
			size_t prc = i;

			Point<dim,T> p1 = (middle * toPoint<dim,T>::convert(cmbs[j]) + middle)* 1.0/1.1;
			Point<dim,T> p2 = p1 + Point<dim,T>({0.1,0.1,0.1}) * 1.0/1.1;

			Box<dim,T> bx(p1,p2);
			nnp_sub[prc+1].id = prc;
			nnp_sub[prc+1].bx.add(bx);

			tot_n++;
		}
	}

	for (size_t i = 0; i < dim; i++)
	{
		nnp_np.add(i+1);
	}*/

	// check that nn_processor contain the correct boxes

	nnp.applyBC(domain,ghost,bc);

	if (v_cl.getProcessUnitID() == 0)
	{

	}
//	BOOST_REQUIRE_EQUAL(nnp.getNearSubdomains(nnp.IDtoProc(2)).size(),12ul);
//	BOOST_REQUIRE_EQUAL(nnp.getNearSubdomains(nnp.IDtoProc(0)).size(),8ul*8ul);
//	BOOST_REQUIRE_EQUAL(nnp.getNearSubdomains(nnp.IDtoProc(1)).size(),12ul*4ul);
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* SRC_DECOMPOSITION_NN_PROCESSOR_UNIT_TEST_HPP_ */
