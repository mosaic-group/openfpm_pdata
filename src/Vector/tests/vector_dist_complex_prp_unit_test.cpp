/*
 * vector_dist_complex_prp_unit_test.hpp
 *
 *  Created on: Sep 18, 2016
 *      Author: i-bird
 */
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Vector/vector_dist.hpp"
#include "Vector/performance/vector_dist_performance_util.hpp"
#include "vector_dist_util_unit_tests.hpp"

extern void print_test_v(std::string test, size_t sz);
extern long int decrement(long int k, long int step);

BOOST_AUTO_TEST_CASE( vector_dist_periodic_complex_prp_test_use_3d )
{
	Vcluster & v_cl = create_vcluster();

	if (v_cl.getProcessingUnits() > 48)
		return;

    // set the seed
	// create the random generator engine
	std::srand(v_cl.getProcessUnitID());
    std::default_random_engine eg;
    std::uniform_real_distribution<float> ud(0.0f, 1.0f);

#ifdef TEST_COVERAGE_MODE
    long int k = 24288 * v_cl.getProcessingUnits();
#else
    long int k = 124288 * v_cl.getProcessingUnits();
#endif

	long int big_step = k / 4;
	big_step = (big_step == 0)?1:big_step;

	print_test_v( "Testing 3D vector full complex prp vector k<=",k);

	// 3D test
	for ( ; k >= 2 ; k-= decrement(k,big_step) )
	{
		BOOST_TEST_CHECKPOINT( "Testing 3D full complex prp vector k=" << k );

		Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

		// Boundary conditions
		size_t bc[3]={NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

		// factor
		float factor = pow(create_vcluster().getProcessingUnits()/2.0f,1.0f/3.0f);

		// ghost
		Ghost<3,float> ghost(0.05 / factor);

		// ghost2 (a little bigger because of round off error)
		Ghost<3,float> ghost2(0.05001 / factor);

		// Distributed vector
		vector_dist<3,float, aggregate<openfpm::vector<float>,grid_cpu<3,aggregate<double,double[3]>> > > vd(k,box,bc,ghost);

		auto it = vd.getIterator();

		while (it.isNext())
		{
			auto key = it.get();

			vd.getPos(key)[0] = ud(eg);
			vd.getPos(key)[1] = ud(eg);
			vd.getPos(key)[2] = ud(eg);

			// initalize vector property and grid

			vd.getProp<0>(key).add(vd.getPos(key)[0]);
			vd.getProp<0>(key).add(vd.getPos(key)[1]);
			vd.getProp<0>(key).add(vd.getPos(key)[2]);

			vd.getProp<0>(key).add(vd.getPos(key)[0]+vd.getPos(key)[1]);
			vd.getProp<0>(key).add(vd.getPos(key)[1]+vd.getPos(key)[2]);
			vd.getProp<0>(key).add(vd.getPos(key)[0]+vd.getPos(key)[2]);

			// Grid
			size_t sz[] = {3,3,3};
			vd.getProp<1>(key).resize(sz);

			vd.getProp<1>(key).get<0>({0,0,0}) = vd.getPos(key)[0];
			vd.getProp<1>(key).get<0>({0,0,1}) = vd.getPos(key)[1];
			vd.getProp<1>(key).get<0>({0,0,2}) = vd.getPos(key)[2];

			vd.getProp<1>(key).get<0>({0,1,0}) = 2.0*vd.getPos(key)[0];
			vd.getProp<1>(key).get<0>({0,1,1}) = 2.0*vd.getPos(key)[1];
			vd.getProp<1>(key).get<0>({0,1,2}) = 2.0*vd.getPos(key)[2];

			vd.getProp<1>(key).get<0>({0,2,0}) = 3.0*vd.getPos(key)[0];
			vd.getProp<1>(key).get<0>({0,2,1}) = 3.0*vd.getPos(key)[1];
			vd.getProp<1>(key).get<0>({0,2,2}) = 3.0*vd.getPos(key)[2];

			vd.getProp<1>(key).get<0>({1,0,0}) = 4.0*vd.getPos(key)[0];
			vd.getProp<1>(key).get<0>({1,0,1}) = 4.0*vd.getPos(key)[1];
			vd.getProp<1>(key).get<0>({1,0,2}) = 4.0*vd.getPos(key)[2];

			vd.getProp<1>(key).get<0>({1,1,0}) = 5.0*vd.getPos(key)[0];
			vd.getProp<1>(key).get<0>({1,1,1}) = 5.0*vd.getPos(key)[1];
			vd.getProp<1>(key).get<0>({1,1,2}) = 5.0*vd.getPos(key)[2];

			vd.getProp<1>(key).get<0>({1,2,0}) = 6.0*vd.getPos(key)[0];
			vd.getProp<1>(key).get<0>({1,2,1}) = 6.0*vd.getPos(key)[1];
			vd.getProp<1>(key).get<0>({1,2,2}) = 6.0*vd.getPos(key)[2];

			vd.getProp<1>(key).get<0>({2,0,0}) = 7.0*vd.getPos(key)[0];
			vd.getProp<1>(key).get<0>({2,0,1}) = 7.0*vd.getPos(key)[1];
			vd.getProp<1>(key).get<0>({2,0,2}) = 7.0*vd.getPos(key)[2];

			vd.getProp<1>(key).get<0>({2,1,0}) = 8.0*vd.getPos(key)[0];
			vd.getProp<1>(key).get<0>({2,1,1}) = 8.0*vd.getPos(key)[1];
			vd.getProp<1>(key).get<0>({2,1,2}) = 8.0*vd.getPos(key)[2];

			vd.getProp<1>(key).get<0>({2,2,0}) = 9.0*vd.getPos(key)[0];
			vd.getProp<1>(key).get<0>({2,2,1}) = 9.0*vd.getPos(key)[1];
			vd.getProp<1>(key).get<0>({2,2,2}) = 9.0*vd.getPos(key)[2];

			++it;
		}

		vd.map();

		// sync the ghost
		vd.ghost_get<0,1>();

		// Domain + ghost
		Box<3,float> dom_ext = box;
		dom_ext.enlarge(ghost2);

		// Iterate on all particles domain + ghost
		size_t l_cnt = 0;
		size_t nl_cnt = 0;
		size_t n_out = 0;

		auto it2 = vd.getIterator();
		count_local_n_local<3,vector_dist<3,float, aggregate<openfpm::vector<float>,grid_cpu<3,aggregate<double,double[3]>>>>  >(vd,it2,bc,box,dom_ext,l_cnt,nl_cnt,n_out);

		// No particles should be out of domain + ghost
		BOOST_REQUIRE_EQUAL(n_out,0ul);

		// Ghost must populated because we synchronized them
		if (k > 524288)
		{
			BOOST_REQUIRE(nl_cnt != 0);
			BOOST_REQUIRE(l_cnt > nl_cnt);
		}

		// Sum all the particles inside the domain
		v_cl.sum(l_cnt);
		v_cl.execute();
		BOOST_REQUIRE_EQUAL(l_cnt,(size_t)k);

		l_cnt = 0;
		nl_cnt = 0;

		// Iterate only on the ghost particles
		auto itg = vd.getGhostIterator();
		count_local_n_local<3, vector_dist<3,float, aggregate<openfpm::vector<float>,grid_cpu<3,aggregate<double,double[3]>>>> >(vd,itg,bc,box,dom_ext,l_cnt,nl_cnt,n_out);

		// No particle on the ghost must be inside the domain
		BOOST_REQUIRE_EQUAL(l_cnt,0ul);

		// Ghost must be populated
		if (k > 524288)
		{
			BOOST_REQUIRE(nl_cnt != 0);
		}

		// check that the particles contain the correct information

		bool ret = true;
		auto it3 = vd.getDomainAndGhostIterator();

		while (it3.isNext())
		{
			auto key = it3.get();

			ret &= vd.getProp<0>(key).size() == 6;

			ret &= vd.getProp<0>(key).get(0) == vd.getPos(key)[0];
			ret &= vd.getProp<0>(key).get(1) == vd.getPos(key)[1];
			ret &= vd.getProp<0>(key).get(2) == vd.getPos(key)[2];

			ret &= vd.getProp<0>(key).get(3) == vd.getPos(key)[0]+vd.getPos(key)[1];
			ret &= vd.getProp<0>(key).get(4) == vd.getPos(key)[1]+vd.getPos(key)[2];
			ret &= vd.getProp<0>(key).get(5) == vd.getPos(key)[0]+vd.getPos(key)[2];

//			BOOST_REQUIRE_EQUAL(vd.getProp<0>(key).get(0),vd.getPos(key)[0]);

			ret &= vd.getProp<1>(key).size() == 3*3*3;
			ret &= vd.getProp<1>(key).get<0>({0,0,0}) == vd.getPos(key)[0];
			ret &= vd.getProp<1>(key).get<0>({0,0,1}) == vd.getPos(key)[1];
			ret &= vd.getProp<1>(key).get<0>({0,0,2}) == vd.getPos(key)[2];

			ret &= vd.getProp<1>(key).get<0>({0,1,0}) == 2.0*vd.getPos(key)[0];
			ret &= vd.getProp<1>(key).get<0>({0,1,1}) == 2.0*vd.getPos(key)[1];
			ret &= vd.getProp<1>(key).get<0>({0,1,2}) == 2.0*vd.getPos(key)[2];

			ret &= vd.getProp<1>(key).get<0>({0,2,0}) == 3.0*vd.getPos(key)[0];
			ret &= vd.getProp<1>(key).get<0>({0,2,1}) == 3.0*vd.getPos(key)[1];
			ret &= vd.getProp<1>(key).get<0>({0,2,2}) == 3.0*vd.getPos(key)[2];

			ret &= vd.getProp<1>(key).get<0>({1,0,0}) == 4.0*vd.getPos(key)[0];
			ret &= vd.getProp<1>(key).get<0>({1,0,1}) == 4.0*vd.getPos(key)[1];
			ret &= vd.getProp<1>(key).get<0>({1,0,2}) == 4.0*vd.getPos(key)[2];

			ret &= vd.getProp<1>(key).get<0>({1,1,0}) == 5.0*vd.getPos(key)[0];
			ret &= vd.getProp<1>(key).get<0>({1,1,1}) == 5.0*vd.getPos(key)[1];
			ret &= vd.getProp<1>(key).get<0>({1,1,2}) == 5.0*vd.getPos(key)[2];

			ret &= vd.getProp<1>(key).get<0>({1,2,0}) == 6.0*vd.getPos(key)[0];
			ret &= vd.getProp<1>(key).get<0>({1,2,1}) == 6.0*vd.getPos(key)[1];
			ret &= vd.getProp<1>(key).get<0>({1,2,2}) == 6.0*vd.getPos(key)[2];

			ret &= vd.getProp<1>(key).get<0>({2,0,0}) == 7.0*vd.getPos(key)[0];
			ret &= vd.getProp<1>(key).get<0>({2,0,1}) == 7.0*vd.getPos(key)[1];
			ret &= vd.getProp<1>(key).get<0>({2,0,2}) == 7.0*vd.getPos(key)[2];

			ret &= vd.getProp<1>(key).get<0>({2,1,0}) == 8.0*vd.getPos(key)[0];
			ret &= vd.getProp<1>(key).get<0>({2,1,1}) == 8.0*vd.getPos(key)[1];
			ret &= vd.getProp<1>(key).get<0>({2,1,2}) == 8.0*vd.getPos(key)[2];

			ret &= vd.getProp<1>(key).get<0>({2,2,0}) == 9.0*vd.getPos(key)[0];
			ret &= vd.getProp<1>(key).get<0>({2,2,1}) == 9.0*vd.getPos(key)[1];
			ret &= vd.getProp<1>(key).get<0>({2,2,2}) == 9.0*vd.getPos(key)[2];

			++it3;
		}

		BOOST_REQUIRE_EQUAL(ret,true);
	}
}


