/*
 * vector_dist_dlb_test.hpp
 *
 *  Created on: Feb 21, 2017
 *      Author: i-bird
 */

#ifndef SRC_VECTOR_VECTOR_DIST_DLB_TEST_HPP_
#define SRC_VECTOR_VECTOR_DIST_DLB_TEST_HPP_

#include "DLB/LB_Model.hpp"
#include "Vector/vector_dist.hpp"

BOOST_AUTO_TEST_SUITE( vector_dist_dlb_test )

template<typename vector_type>
void mp_test_template(vector_type & vd0, vector_type & vd1, vector_type & vd2, vector_type & vd3)
{
	Vcluster<> & v_cl = create_vcluster();

	// Only processor 0 initialy add particles on a corner of a domain

	if (v_cl.getProcessUnitID() == 0)
	{
		for(size_t i = 0 ; i < 50000 ; i++)
		{
			vd0.add();
			vd1.add();
			vd2.add();
			vd3.add();

			// coverty[dont_call]
			vd0.getLastPos()[0] = ((float)rand())/RAND_MAX * 0.3;
			// coverty[dont_call]
			vd0.getLastPos()[1] = ((float)rand())/RAND_MAX * 0.3;
			// coverty[dont_call]
			vd0.getLastPos()[2] = ((float)rand())/RAND_MAX * 0.3;

			// coverty[dont_call]
			vd1.getLastPos()[0] = ((float)rand())/RAND_MAX * 0.3 + 0.1;
			// coverty[dont_call]
			vd1.getLastPos()[1] = ((float)rand())/RAND_MAX * 0.3 + 0.1;
			// coverty[dont_call]
			vd1.getLastPos()[2] = ((float)rand())/RAND_MAX * 0.3 + 0.1;

			// coverty[dont_call]
			vd2.getLastPos()[0] = ((float)rand())/RAND_MAX * 0.3 + 0.2;
			// coverty[dont_call]
			vd2.getLastPos()[1] = ((float)rand())/RAND_MAX * 0.3 + 0.2;
			// coverty[dont_call]
			vd2.getLastPos()[2] = ((float)rand())/RAND_MAX * 0.3 + 0.2;

			// coverty[dont_call]
			vd3.getLastPos()[0] = ((float)rand())/RAND_MAX * 0.3 + 0.3;
			// coverty[dont_call]
			vd3.getLastPos()[1] = ((float)rand())/RAND_MAX * 0.3 + 0.3;
			// coverty[dont_call]
			vd3.getLastPos()[2] = ((float)rand())/RAND_MAX * 0.3 + 0.3;
		}
	}

	vd0.map();
	vd0.template ghost_get<>();
	vd1.map();
	vd1.template ghost_get<>();
	vd2.map();
	vd2.template ghost_get<>();
	vd3.map();
	vd3.template ghost_get<>();

	ModelSquare md;
	md.factor = 1;
	vd0.initializeComputationCosts();
	vd0.addComputationCosts(vd0,md);
	vd0.addComputationCosts(vd1,md);
	vd0.addComputationCosts(vd2,md);
	vd0.addComputationCosts(vd3,md);
	vd0.finalizeComputationCosts();

	vd0.getDecomposition().decompose();

	// Copy the decomposition back to the other
	vd1.getDecomposition() = vd0.getDecomposition();
	vd2.getDecomposition() = vd0.getDecomposition();
	vd3.getDecomposition() = vd0.getDecomposition();

	vd0.map();
	vd1.map();
	vd2.map();
	vd3.map();

	vd0.initializeComputationCosts();
	vd0.addComputationCosts(vd0,md);
	vd0.addComputationCosts(vd1,md);
	vd0.addComputationCosts(vd2,md);
	vd0.addComputationCosts(vd3,md);
	vd0.finalizeComputationCosts();

	openfpm::vector<size_t> loads;
	size_t load = vd0.getDecomposition().getDistribution().getProcessorLoad();
	v_cl.allGather(load,loads);
	v_cl.execute();

	for (size_t i = 0 ; i < loads.size() ; i++)
	{
		float load_f = load;
		float load_fc = loads.get(i);

		BOOST_REQUIRE_CLOSE(load_f,load_fc,7.0);
	}

	Point<3,float> v({1.0,1.0,1.0});

	for (size_t i = 0 ; i < 25 ; i++)
	{
		// move the particles by 0.1

		{
		auto it = vd0.getDomainIterator();

		while (it.isNext())
		{
			auto p = it.get();

			vd0.getPos(p)[0] += v.get(0) * 0.09;
			vd0.getPos(p)[1] += v.get(1) * 0.09;
			vd0.getPos(p)[2] += v.get(2) * 0.09;

			++it;
		}
		}

		{
		auto it = vd1.getDomainIterator();
		while (it.isNext())
		{
			auto p = it.get();

			vd1.getPos(p)[0] += v.get(0) * 0.06;
			vd1.getPos(p)[1] += v.get(1) * 0.06;
			vd1.getPos(p)[2] += v.get(2) * 0.06;

			++it;
		}
		}

		{
		auto it = vd2.getDomainIterator();
		while (it.isNext())
		{
			auto p = it.get();

			vd2.getPos(p)[0] += v.get(0) * 0.06;
			vd2.getPos(p)[1] += v.get(1) * 0.06;
			vd2.getPos(p)[2] += v.get(2) * 0.06;

			++it;
		}
		}

		{
		auto it = vd3.getDomainIterator();
		while (it.isNext())
		{
			auto p = it.get();

			vd3.getPos(p)[0] += v.get(0) * 0.06;
			vd3.getPos(p)[1] += v.get(1) * 0.06;
			vd3.getPos(p)[2] += v.get(2) * 0.06;

			++it;
		}
		}

		vd0.map();
		vd1.map();
		vd2.map();
		vd3.map();

		ModelSquare md;
		vd0.initializeComputationCosts();
		vd0.addComputationCosts(vd0,md);
		vd0.addComputationCosts(vd1,md);
		vd0.addComputationCosts(vd2,md);
		vd0.addComputationCosts(vd3,md);
		vd0.finalizeComputationCosts();

		vd0.getDecomposition().redecompose(200);

		// Copy the decomposition back to the other
		vd1.getDecomposition() = vd0.getDecomposition();
		vd2.getDecomposition() = vd0.getDecomposition();
		vd3.getDecomposition() = vd0.getDecomposition();
		vd0.map();
		vd1.map();
		vd2.map();
		vd3.map();

		vd0.template ghost_get<>();
		vd1.template ghost_get<>();
		vd2.template ghost_get<>();
		vd3.template ghost_get<>();

		vd0.initializeComputationCosts();
		vd0.addComputationCosts(vd0,md);
		vd0.addComputationCosts(vd1,md);
		vd0.addComputationCosts(vd2,md);
		vd0.addComputationCosts(vd3,md);
		vd0.finalizeComputationCosts();

		openfpm::vector<size_t> loads;
		size_t load = vd0.getDecomposition().getDistribution().getProcessorLoad();
		v_cl.allGather(load,loads);
		v_cl.execute();

		for (size_t i = 0 ; i < loads.size() ; i++)
		{
			float load_f = load;
			float load_fc = loads.get(i);

			BOOST_REQUIRE_CLOSE(load_f,load_fc,10.0);
		}
	}
}

template<typename vector_type> void test_dlb_vector()
{
	Vcluster<> & v_cl = create_vcluster();

	if (v_cl.getProcessingUnits() > 8)
		return;

	Box<3,double> domain({0.0,0.0,0.0},{1.0,1.0,1.0});
	Ghost<3,double> g(0.1);
	size_t bc[3] = {PERIODIC,PERIODIC,PERIODIC};

	vector_type vd(0,domain,bc,g,DEC_GRAN(2048));

	// Only processor 0 initialy add particles on a corner of a domain

	if (v_cl.getProcessUnitID() == 0)
	{
		for(size_t i = 0 ; i < 50000 ; i++)
		{
			vd.add();

			vd.getLastPos()[0] = ((double)rand())/RAND_MAX * 0.3;
			vd.getLastPos()[1] = ((double)rand())/RAND_MAX * 0.3;
			vd.getLastPos()[2] = ((double)rand())/RAND_MAX * 0.3;
		}
	}

	vd.map();

	vd.template ghost_get<>();

	// Get the neighborhood of each particles

	auto VV = vd.getVerlet(0.01);

	// store the number of neighborhood for each particles

	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		vd.template getProp<0>(p) = VV.getNNPart(p.getKey());

		++it;
	}

	ModelSquare md;
	md.factor = 10;
	vd.map();

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
		double load_f = load;
		double load_fc = loads.get(i);

		BOOST_REQUIRE_CLOSE(load_f,load_fc,7.0);
	}

	BOOST_REQUIRE(vd.size_local() != 0);

	Point<3,double> v({1.0,1.0,1.0});

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

		vd.template ghost_get<>();

		auto VV2 = vd.getVerlet(0.01);

		auto it2 = vd.getDomainIterator();

		bool match = true;
		while (it2.isNext())
		{
			auto p = it2.get();

			match &= vd.template getProp<0>(p) == VV2.getNNPart(p.getKey());

			++it2;
		}

		BOOST_REQUIRE_EQUAL(match,true);

		ModelSquare md;
		vd.addComputationCosts(md);
		vd.getDecomposition().redecompose(200);
		vd.map();

		BOOST_REQUIRE(vd.size_local() != 0);

		vd.template ghost_get<>();

		vd.addComputationCosts(md);

		openfpm::vector<size_t> loads;
		size_t load = vd.getDecomposition().getDistribution().getProcessorLoad();
		v_cl.allGather(load,loads);
		v_cl.execute();

		for (size_t i = 0 ; i < loads.size() ; i++)
		{
			double load_f = load;
			double load_fc = loads.get(i);

			BOOST_REQUIRE_CLOSE(load_f,load_fc,10.0);
		}
	}
}


template<typename vector_type> void test_dlb_multi_phase_vector()
{
	Vcluster<> & v_cl = create_vcluster();

	if (v_cl.getProcessingUnits() > 8)
		return;

	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});
	Ghost<3,float> g(0.1);
	size_t bc[3] = {PERIODIC,PERIODIC,PERIODIC};

	vector_type vd0(0,domain,bc,g,DEC_GRAN(2048));
	vector_type vd1(0,domain,bc,g,DEC_GRAN(2048));
	vector_type vd2(0,domain,bc,g,DEC_GRAN(2048));
	vector_type vd3(0,domain,bc,g,DEC_GRAN(2048));

	mp_test_template(vd0,vd1,vd2,vd3);
}



template<typename vector_type> void test_dlb_multi_phase_v_vector()
{
	Vcluster<> & v_cl = create_vcluster();

	if (v_cl.getProcessingUnits() > 8)
		return;

	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});
	Ghost<3,float> g(0.1);
	size_t bc[3] = {PERIODIC,PERIODIC,PERIODIC};

	openfpm::vector<vector_type> v_phases;
	{
		vector_type vd0(0,domain,bc,g,DEC_GRAN(2048));
		v_phases.add(vd0);
		v_phases.add(vector_type(vd0.getDecomposition(),0));
		v_phases.add(vector_type(vd0.getDecomposition(),0));
		v_phases.add(vector_type(vd0.getDecomposition(),0));
	}

	auto & vd0 = v_phases.get(0);
	auto & vd1 = v_phases.get(1);
	auto & vd2 = v_phases.get(2);
	auto & vd3 = v_phases.get(3);

	mp_test_template(vd0,vd1,vd2,vd3);
}

BOOST_AUTO_TEST_CASE( vector_dist_dlb_test_part )
{
	test_dlb_vector<vector_dist<3,double,aggregate<double>>>();
}

BOOST_AUTO_TEST_CASE( vector_dist_dlb_multi_phase_test_part )
{
	test_dlb_multi_phase_vector<vector_dist<3,float,aggregate<float>>>();
}

BOOST_AUTO_TEST_CASE( vector_dist_dlb_multi_phase_v_test_part )
{
	test_dlb_multi_phase_v_vector<vector_dist<3,float,aggregate<float>>>();
}

BOOST_AUTO_TEST_CASE( vector_dist_dlb_metis_test_part )
{
	test_dlb_vector<vector_dist<3,
	                            double,
								aggregate<double>,
	                            CartDecomposition<3,double,HeapMemory,memory_traits_lin,MetisDistribution<3,double>>>>();
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* SRC_VECTOR_VECTOR_DIST_DLB_TEST_HPP_ */
