/*
 * vector_dist_dlb_test.hpp
 *
 *  Created on: Feb 21, 2017
 *      Author: i-bird
 */

#ifndef SRC_VECTOR_VECTOR_DIST_DLB_TEST_HPP_
#define SRC_VECTOR_VECTOR_DIST_DLB_TEST_HPP_

BOOST_AUTO_TEST_SUITE( vector_dist_dlb_test )

BOOST_AUTO_TEST_CASE( vector_dist_dlb_test_part )
{
	Vcluster & v_cl = create_vcluster();

	if (v_cl.getProcessingUnits() > 8)
		return;

	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});
	Ghost<3,float> g(0.1);
	size_t bc[3] = {PERIODIC,PERIODIC,PERIODIC};

	vector_dist<3,float,aggregate<float>> vd(0,domain,bc,g,DEC_GRAN(2048));

	// Only processor 0 initialy add particles on a corner of a domain

	if (v_cl.getProcessUnitID() == 0)
	{
		for(size_t i = 0 ; i < 10000 ; i++)
		{
			vd.add();

			vd.getLastPos()[0] = ((float)rand())/RAND_MAX * 0.3;
			vd.getLastPos()[1] = ((float)rand())/RAND_MAX * 0.3;
			vd.getLastPos()[2] = ((float)rand())/RAND_MAX * 0.3;
		}
	}



	vd.map();
	vd.ghost_get<>();

	ModelSquare md;
	md.factor = 10;
	vd.addComputationCosts(md);
	vd.getDecomposition().decompose();
	vd.map();

	vd.addComputationCosts(md);

	openfpm::vector<size_t> loads;
	size_t load = vd.getDecomposition().getDistribution().getProcessorLoad();
	v_cl.allGather(load,loads);
	v_cl.execute();

	for (size_t i = 0 ; i < loads.size() ; i++)
	{
		float load_f = load;
		float load_fc = loads.get(i);

		BOOST_REQUIRE_CLOSE(load_f,load_fc,7.0);
	}

	BOOST_REQUIRE(vd.size_local() != 0);

	Point<3,float> v({1.0,1.0,1.0});

	for (size_t i = 0 ; i < 25 ; i++)
	{
		// move the particles by 0.1

		auto it = vd.getDomainIterator();

		while (it.isNext())
		{
			auto p = it.get();

			vd.getPos(p)[0] += v.get(0) * 0.09;
			vd.getPos(p)[1] += v.get(1) * 0.09;
			vd.getPos(p)[2] += v.get(2) * 0.09;

			++it;
		}
		vd.map();

		ModelSquare md;
		vd.addComputationCosts(md);
		vd.getDecomposition().redecompose(200);
		vd.map();

		BOOST_REQUIRE(vd.size_local() != 0);

		vd.ghost_get<>();

		vd.addComputationCosts(md);

		openfpm::vector<size_t> loads;
		size_t load = vd.getDecomposition().getDistribution().getProcessorLoad();
		v_cl.allGather(load,loads);
		v_cl.execute();

		for (size_t i = 0 ; i < loads.size() ; i++)
		{
			float load_f = load;
			float load_fc = loads.get(i);

			BOOST_REQUIRE_CLOSE(load_f,load_fc,7.0);
		}
	}
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* SRC_VECTOR_VECTOR_DIST_DLB_TEST_HPP_ */
