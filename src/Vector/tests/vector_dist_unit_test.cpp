/*
 * vector_dist_unit_test.hpp
 *
 *  Created on: Mar 6, 2015
 *      Author: Pietro Incardona
 */

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "config.h"

#include <random>
#include "Vector/vector_dist.hpp"
#include "data_type/aggregate.hpp"
#include "vector_dist_util_unit_tests.hpp"
#include "Point_test.hpp"
#include "Vector/performance/vector_dist_performance_common.hpp"

/*! \brief Print a string about the test
 *
 * \param test string to print
 * \param sz size
 *
 */
void print_test_v(std::string test, size_t sz)
{
	if (create_vcluster().getProcessUnitID() == 0)
		std::cout << test << " " << sz << "\n";
}

/*! \brief Get next testing step decrementing the size
 *
 * \param k actual size
 * \param step
 *
 * \return the next step
 *
 */
long int decrement(long int k, long int step)
{
	if (k <= 32)
	{
		return 1;
	}
	else if (k - 2*step+1 <= 0)
	{
		return k - 32;
	}
	else
		return step;
}

/*! \brief Count the total number of particles
 *
 * \param vd distributed vector
 * \param bc boundary conditions
 *
 */
template<unsigned int dim, template <typename> class layout>
size_t total_n_part_lc(vector_dist<dim,float, Point_test<float>, CartDecomposition<dim,float>, HeapMemory, layout > & vd, size_t (& bc)[dim])
{
	Vcluster<> & v_cl = vd.getVC();
	auto it2 = vd.getDomainIterator();
	const CartDecomposition<3,float> & ct = vd.getDecomposition();

	bool noOut = true;

	size_t cnt = 0;
	while (it2.isNext())
	{
		auto key = it2.get();

		noOut &= ct.isLocal(vd.getPos(key));

		cnt++;

		++it2;
	}

	BOOST_REQUIRE_EQUAL(noOut,true);

	//
	v_cl.sum(cnt);
	v_cl.execute();

	return cnt;
}


BOOST_AUTO_TEST_SUITE( vector_dist_test )

void print_test(std::string test, size_t sz)
{
	if (create_vcluster().getProcessUnitID() == 0)
		std::cout << test << " " << sz << "\n";
}

template<typename vector>
void Test2D_ghost(Box<2,float> & box)
{
	// Communication object
	Vcluster<> & v_cl = create_vcluster();

	typedef Point_test<float> p;

	// Get the default minimum number of sub-sub-domain per processor (granularity of the decomposition)
	size_t n_sub = 64 * v_cl.getProcessingUnits();
	// Convert the request of having a minimum n_sub number of sub-sub domain into grid decompsition of the space
	size_t sz = CartDecomposition<2,float>::getDefaultGrid(n_sub);

	//! [Create a vector of elements distributed on a grid like way]

	size_t g_div[]= {sz,sz};

	// number of particles
	size_t np = sz * sz;

	// Calculate the number of elements this processor is going to obtain
	size_t p_np = np / v_cl.getProcessingUnits();

	// Get non divisible part
	size_t r = np % v_cl.getProcessingUnits();

	// Get the offset
	size_t offset = v_cl.getProcessUnitID() * p_np + std::min(v_cl.getProcessUnitID(),r);

	// Distribute the remain elements
	if (v_cl.getProcessUnitID() < r)
		p_np++;

	// Create a grid info
	grid_sm<2,void> g_info(g_div);

	// Calculate the grid spacing
	Point<2,float> spacing = box.getP2() - box.getP1();
	spacing = spacing / g_div;

	// middle spacing
	Point<2,float> m_spacing = spacing / 2.0;

	// set the ghost based on the radius cut off (make just a little bit smaller than the spacing)
	Ghost<2,float> g(spacing.get(0) - spacing .get(0) * 0.0001);

	// Boundary conditions
	size_t bc[2]={NON_PERIODIC,NON_PERIODIC};

	// Vector of particles
	vector vd(g_info.size(),box,bc,g);

	// size_t
	size_t cobj = 0;

	grid_key_dx_iterator_sp<2> it(g_info,offset,offset+p_np-1);
	auto v_it = vd.getIterator();

	while (v_it.isNext() && it.isNext())
	{
		auto key = it.get();
		auto key_v = v_it.get();

		// set the particle position

		vd.getPos(key_v)[0] = key.get(0) * spacing[0] + m_spacing[0] + box.getLow(0);
		vd.getPos(key_v)[1] = key.get(1) * spacing[1] + m_spacing[1] + box.getLow(1);

		cobj++;

		++v_it;
		++it;
	}

	//! [Create a vector of elements distributed on a grid like way]

	// Both iterators must signal the end, and the number of elements in the vector, must the equal to the
	// predicted one
	BOOST_REQUIRE_EQUAL(v_it.isNext(),false);
	BOOST_REQUIRE_EQUAL(it.isNext(),false);
	BOOST_REQUIRE_EQUAL(cobj,p_np);

	//! [Redistribute the particles and sync the ghost properties]

	// redistribute the particles according to the decomposition
	vd.map();

	auto v_it2 = vd.getIterator();

	while (v_it2.isNext())
	{
		auto key = v_it2.get();

		// fill with the processor ID where these particle live
		vd.template getProp<p::s>(key) = vd.getPos(key)[0] + vd.getPos(key)[1] * 16.0f;
		vd.template getProp<p::v>(key)[0] = v_cl.getProcessUnitID();
		vd.template getProp<p::v>(key)[1] = v_cl.getProcessUnitID();
		vd.template getProp<p::v>(key)[2] = v_cl.getProcessUnitID();

		++v_it2;
	}

	// do a ghost get
	vd.template ghost_get<p::s,p::v>();

	//! [Redistribute the particles and sync the ghost properties]

	// Get the decomposition
	const auto & dec = vd.getDecomposition();

	// Get the ghost external boxes
	openfpm::vector<size_t> vb(dec.getNEGhostBox());

	// Get the ghost iterator
	auto g_it = vd.getGhostIterator();

	size_t n_part = 0;

	// Check if the ghost particles contain the correct information
	while (g_it.isNext())
	{
		auto key = g_it.get();

		// Check the received data
		float prp = vd.template getProp<p::s>(key);
		float prp2 = vd.getPos(key)[0] + vd.getPos(key)[1] * 16.0f;
		BOOST_REQUIRE_EQUAL(prp2,prp);

		bool is_in = false;
		size_t b = 0;
		size_t lb = 0;

		// check if the received data are in one of the ghost boxes
		for ( ; b < dec.getNEGhostBox() ; b++)
		{
			Point<2,float> xp = vd.getPos(key);

			if (dec.getEGhostBox(b).isInside(xp) == true )
			{
				is_in = true;

				// Add
				vb.get(b)++;
				lb = b;
			}
		}
		BOOST_REQUIRE_EQUAL(is_in,true);

		// Check that the particle come from the correct processor
		BOOST_REQUIRE_EQUAL(vd.template getProp<p::v>(key)[0],dec.getEGhostBoxProcessor(lb));

		n_part++;
		++g_it;
	}

	if (v_cl.getProcessingUnits() > 1)
	{
		BOOST_REQUIRE(n_part != 0);
	}

    CellDecomposer_sm<2,float,shift<2,float>> cd(Box<2,float>(box),g_div,0);

	for (size_t i = 0 ; i < vb.size() ; i++)
	{
		// Calculate how many particle should be in the box
		size_t n_point = cd.getGridPoints(dec.getEGhostBox(i)).getVolumeKey();

		BOOST_REQUIRE_EQUAL(n_point,vb.get(i));
	}
}

BOOST_AUTO_TEST_CASE( vector_dist_ghost )
{
	typedef vector_dist<2,float, Point_test<float>> vector;

	Box<2,float> box({0.0,0.0},{1.0,1.0});
	Test2D_ghost<vector>(box);

	Box<2,float> box2({-1.0,-1.0},{2.5,2.5});
	Test2D_ghost<vector>(box2);
}

BOOST_AUTO_TEST_CASE( vector_dist_ghost_inte )
{
	typedef vector_dist_soa<2,float, Point_test<float>> vector;

	Box<2,float> box({0.0,0.0},{1.0,1.0});
	Test2D_ghost<vector>(box);

	Box<2,float> box2({-1.0,-1.0},{2.5,2.5});
	Test2D_ghost<vector>(box2);
}



BOOST_AUTO_TEST_CASE( vector_dist_iterator_test_use_2d )
{
	Vcluster<> & v_cl = create_vcluster();

    // set the seed
	// create the random generator engine
	std::srand(v_cl.getProcessUnitID());
    std::default_random_engine eg;
    std::uniform_real_distribution<float> ud(0.0f, 1.0f);

#ifdef TEST_COVERAGE_MODE
    long int k = 24288 * v_cl.getProcessingUnits();
#else
    long int k = 524288 * v_cl.getProcessingUnits();
#endif

	long int big_step = k / 4;
	big_step = (big_step == 0)?1:big_step;

	print_test_v( "Testing 2D vector k<=",k);

	// 2D test
	for ( ; k >= 2 ; k-= decrement(k,big_step) )
	{
		BOOST_TEST_CHECKPOINT( "Testing 2D vector k=" << k );

		//! [Create a vector of random elements on each processor 2D]

		Box<2,float> box({0.0,0.0},{1.0,1.0});

		// Boundary conditions
		size_t bc[2]={NON_PERIODIC,NON_PERIODIC};

		vector_dist<2,float, Point_test<float> > vd(k,box,bc,Ghost<2,float>(0.0));

		auto it = vd.getIterator();

		while (it.isNext())
		{
			auto key = it.get();

			vd.getPos(key)[0] = ud(eg);
			vd.getPos(key)[1] = ud(eg);

			++it;
		}

		vd.map();

		//! [Create a vector of random elements on each processor 2D]

		// Check if we have all the local particles
		size_t cnt = 0;
		const CartDecomposition<2,float> & ct = vd.getDecomposition();
		auto it2 = vd.getIterator();

		while (it2.isNext())
		{
			auto key = it2.get();

			// Check if local
			BOOST_REQUIRE_EQUAL(ct.isLocal(vd.getPos(key)),true);

			cnt++;

			++it2;
		}

		//
		v_cl.sum(cnt);
		v_cl.execute();
		BOOST_REQUIRE_EQUAL((long int)cnt,k);
	}
}

BOOST_AUTO_TEST_CASE( vector_dist_iterator_test_use_3d )
{
	Vcluster<> & v_cl = create_vcluster();

    // set the seed
	// create the random generator engine
	std::srand(v_cl.getProcessUnitID());
    std::default_random_engine eg;
    std::uniform_real_distribution<float> ud(0.0f, 1.0f);

#ifdef TEST_COVERAGE_MODE
    long int k = 24288 * v_cl.getProcessingUnits();
#else
    long int k = 524288 * v_cl.getProcessingUnits();
#endif

	long int big_step = k / 4;
	big_step = (big_step == 0)?1:big_step;

	print_test_v( "Testing 3D vector k<=",k);

	// 3D test
	for ( ; k >= 2 ; k-= decrement(k,big_step) )
	{
		BOOST_TEST_CHECKPOINT( "Testing 3D vector k=" << k );

		//! [Create a vector of random elements on each processor 3D]

		Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

		// Boundary conditions
		size_t bc[3]={NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

		vector_dist<3,float, Point_test<float> > vd(k,box,bc,Ghost<3,float>(0.0));

		auto it = vd.getIterator();

		while (it.isNext())
		{
			auto key = it.get();

			vd.getPos(key)[0] = ud(eg);
			vd.getPos(key)[1] = ud(eg);
			vd.getPos(key)[2] = ud(eg);

			++it;
		}

		vd.map();

		//! [Create a vector of random elements on each processor 3D]

		// Check if we have all the local particles
		size_t cnt = 0;
		const CartDecomposition<3,float> & ct = vd.getDecomposition();
		auto it2 = vd.getIterator();

		while (it2.isNext())
		{
			auto key = it2.get();

			// Check if local
			BOOST_REQUIRE_EQUAL(ct.isLocal(vd.getPos(key)),true);

			cnt++;

			++it2;
		}

		//
		v_cl.sum(cnt);
		v_cl.execute();
		BOOST_REQUIRE_EQUAL(cnt,(size_t)k);
	}
}


BOOST_AUTO_TEST_CASE( vector_dist_iterator_fixed_dec_3d )
{
	Vcluster<> & v_cl = create_vcluster();

    // set the seed
	// create the random generator engine
	std::srand(v_cl.getProcessUnitID());
    std::default_random_engine eg;
    std::uniform_real_distribution<float> ud(0.0f, 1.0f);

#ifdef TEST_COVERAGE_MODE
    long int k = 2428 * v_cl.getProcessingUnits();
#else
    long int k = 52428 * v_cl.getProcessingUnits();
#endif

	long int big_step = k / 4;
	big_step = (big_step == 0)?1:big_step;

	print_test_v( "Testing 3D vector copy decomposition k<=",k);

	// 3D test
	for ( ; k >= 2 ; k-= decrement(k,big_step) )
	{
		BOOST_TEST_CHECKPOINT( "Testing 3D vector copy decomposition k=" << k );

		Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

		// Boundary conditions
		size_t bc[3]={NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

		vector_dist<3,float, aggregate<double,double> > vd(k,box,bc,Ghost<3,float>(0.05));
		vector_dist<3,float, aggregate<double,double> > vd2(vd.getDecomposition(),k);

		auto it = vd.getIterator();

		while (it.isNext())
		{
			auto key = it.get();

			vd.getPos(key)[0] = ud(eg);
			vd.getPos(key)[1] = ud(eg);
			vd.getPos(key)[2] = ud(eg);

			vd2.getPos(key)[0] = vd.getPos(key)[0];
			vd2.getPos(key)[1] = vd.getPos(key)[1];
			vd2.getPos(key)[2] = vd.getPos(key)[2];

			++it;
		}

		vd.map();
		vd2.map();

		vd.ghost_get();
		vd2.ghost_get();

		auto NN = vd.getCellList(0.05);
		auto NN2 = vd2.getCellList(0.05);

		cross_calc<3,0>(NN,NN2,vd,vd2);
		cross_calc<3,1>(NN,NN,vd,vd);


		auto it3 = vd.getIterator();

		while (it3.isNext())
		{
			auto key = it3.get();

			BOOST_REQUIRE_EQUAL(vd.getProp<0>(key),vd.getProp<1>(key));

			++it3;
		}
	}
}

BOOST_AUTO_TEST_CASE( vector_dist_periodic_test_use_2d )
{
	Vcluster<> & v_cl = create_vcluster();

    // set the seed
	// create the random generator engine
	std::srand(v_cl.getProcessUnitID());
    std::default_random_engine eg;
    std::uniform_real_distribution<float> ud(0.0f, 1.0f);

#ifdef TEST_COVERAGE_MODE
    long int k = 24288 * v_cl.getProcessingUnits();
#else
    long int k = 524288 * v_cl.getProcessingUnits();
#endif

	long int big_step = k / 4;
	big_step = (big_step == 0)?1:big_step;

	print_test_v( "Testing 2D periodic vector k<=",k);

	// 2D test
	for ( ; k >= 2 ; k-= decrement(k,big_step) )
	{
		BOOST_TEST_CHECKPOINT( "Testing 2D periodic vector k=" << k );

		Box<2,float> box({0.0,0.0},{1.0,1.0});

		// Boundary conditions
		size_t bc[2]={PERIODIC,PERIODIC};

		// factor
		float factor = pow(create_vcluster().getProcessingUnits()/2.0f,1.0f/3.0f);

		// ghost
		Ghost<2,float> ghost(0.01 / factor);

		// ghost2 (a little bigger because of round off error)
		Ghost<2,float> ghost2(0.05001 / factor);

		// Distributed vector
		vector_dist<2,float, Point_test<float> > vd(k,box,bc,ghost);

		auto it = vd.getIterator();

		while (it.isNext())
		{
			auto key = it.get();

			vd.getPos(key)[0] = ud(eg);
			vd.getPos(key)[1] = ud(eg);

			++it;
		}

		vd.map();

		// sync the ghost, only the property zero
		vd.ghost_get<0>();

		// Domain + ghost box
		Box<2,float> dom_ext = box;
		dom_ext.enlarge(ghost2);

		// Iterate on all particles domain + ghost
		size_t l_cnt = 0;
		size_t nl_cnt = 0;
		size_t n_out = 0;


		auto it2 = vd.getIterator();
		count_local_n_local<2,vector_dist<2,float, Point_test<float> >>(vd,it2,bc,box,dom_ext,l_cnt,nl_cnt,n_out);

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

		// count that they are equal to the initial total number
		BOOST_REQUIRE_EQUAL((long int)l_cnt,k);

		l_cnt = 0;
		nl_cnt = 0;

		// Iterate only on the ghost particles
		auto itg = vd.getGhostIterator();
		count_local_n_local<2,vector_dist<2,float, Point_test<float> > >(vd,itg,bc,box,dom_ext,l_cnt,nl_cnt,n_out);

		// No particle on the ghost must be inside the domain
		BOOST_REQUIRE_EQUAL(l_cnt,0ul);

		// Ghost must be populated
		if (k > 524288)
		{
			BOOST_REQUIRE(nl_cnt != 0);
		}
	}
}

BOOST_AUTO_TEST_CASE( vector_dist_periodic_test_use_3d )
{
	Vcluster<> & v_cl = create_vcluster();

    // set the seed
	// create the random generator engine
	std::srand(v_cl.getProcessUnitID());
    std::default_random_engine eg;
    std::uniform_real_distribution<float> ud(0.0f, 1.0f);

#ifdef TEST_COVERAGE_MODE
    long int k = 24288 * v_cl.getProcessingUnits();
#else
    long int k = 524288 * v_cl.getProcessingUnits();
#endif

	long int big_step = k / 4;
	big_step = (big_step == 0)?1:big_step;

	print_test_v( "Testing 3D periodic vector k<=",k);

	// 3D test
	for ( ; k >= 2 ; k-= decrement(k,big_step) )
	{
		BOOST_TEST_CHECKPOINT( "Testing 3D periodic vector k=" << k );

		Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

		// Boundary conditions
		size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

		// factor
		float factor = pow(create_vcluster().getProcessingUnits()/2.0f,1.0f/3.0f);

		// ghost
		Ghost<3,float> ghost(0.05 / factor);

		// ghost2 (a little bigger because of round off error)
		Ghost<3,float> ghost2(0.05001 / factor);

		// Distributed vector
		vector_dist<3,float, Point_test<float> > vd(k,box,bc,ghost);

		auto it = vd.getIterator();

		while (it.isNext())
		{
			auto key = it.get();

			vd.getPos(key)[0] = ud(eg);
			vd.getPos(key)[1] = ud(eg);
			vd.getPos(key)[2] = ud(eg);

			++it;
		}

		vd.map();

		// sync the ghost
		vd.ghost_get<0>();

		// Domain + ghost
		Box<3,float> dom_ext = box;
		dom_ext.enlarge(ghost2);

		// Iterate on all particles domain + ghost
		size_t l_cnt = 0;
		size_t nl_cnt = 0;
		size_t n_out = 0;

		auto it2 = vd.getIterator();
		count_local_n_local<3,vector_dist<3,float, Point_test<float> >>(vd,it2,bc,box,dom_ext,l_cnt,nl_cnt,n_out);

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
		count_local_n_local<3,vector_dist<3,float, Point_test<float> > >(vd,itg,bc,box,dom_ext,l_cnt,nl_cnt,n_out);

		// No particle on the ghost must be inside the domain
		BOOST_REQUIRE_EQUAL(l_cnt,0ul);

		// Ghost must be populated
		if (k > 524288)
		{
			BOOST_REQUIRE(nl_cnt != 0);
		}
	}
}

void test_random_walk(size_t opt)
{
	Vcluster<> & v_cl = create_vcluster();

    // set the seed
	// create the random generator engine
	std::srand(v_cl.getProcessUnitID());
    std::default_random_engine eg;
    std::uniform_real_distribution<float> ud(0.0f, 1.0f);
	
	size_t nsz[] = {0,32,4};
	nsz[0] = 65536 * v_cl.getProcessingUnits();

	print_test_v( "Testing 3D random walk vector k<=",nsz[0]);

	// 3D test
	for (size_t i = 0 ; i < 3 ; i++ )
	{
		size_t k = nsz[i];

		BOOST_TEST_CHECKPOINT( "Testing 3D random walk vector k=" << k );

		Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

		// Boundary conditions
		size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

		// factor
		float factor = pow(create_vcluster().getProcessingUnits()/2.0f,1.0f/3.0f);

		// ghost
		Ghost<3,float> ghost(0.01 / factor);

		// Distributed vector
		vector_dist<3,float, Point_test<float> > vd(k,box,bc,ghost);

		auto it = vd.getIterator();

		while (it.isNext())
		{
			auto key = it.get();

			vd.getPos(key)[0] = ud(eg);
			vd.getPos(key)[1] = ud(eg);
			vd.getPos(key)[2] = ud(eg);

			++it;
		}

		vd.map();

		// 10 step random walk

		for (size_t j = 0 ; j < 4 ; j++)
		{
			auto it = vd.getDomainIterator();

			while (it.isNext())
			{
				auto key = it.get();

				vd.getPos(key)[0] += 0.02 * ud(eg);
				vd.getPos(key)[1] += 0.02 * ud(eg);
				vd.getPos(key)[2] += 0.02 * ud(eg);

				++it;
			}

			vd.map(opt);

			vd.ghost_get<0>();

			// Count the local particles and check that the total number is consistent
			size_t cnt = total_n_part_lc(vd,bc);

			BOOST_REQUIRE_EQUAL((size_t)k,cnt);
		}
	}
}

BOOST_AUTO_TEST_CASE( vector_dist_periodic_test_random_walk )
{
	test_random_walk(NONE);
}

BOOST_AUTO_TEST_CASE( vector_dist_periodic_test_random_walk_local_map )
{
	test_random_walk(MAP_LOCAL);
}

BOOST_AUTO_TEST_CASE( vector_dist_periodic_map )
{
	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// factor
	float factor = pow(create_vcluster().getProcessingUnits()/2.0f,1.0f/3.0f);

	// ghost
	Ghost<3,float> ghost(0.05 / factor);

	// Distributed vector
	vector_dist<3,float, Point_test<float> > vd(1,box,bc,ghost);

	// put particles al 1.0, check that they go to 0.0

	auto it = vd.getIterator();

	while (it.isNext())
	{
		auto key = it.get();

		vd.getPos(key)[0] = 1.0;
		vd.getPos(key)[1] = 1.0;
		vd.getPos(key)[2] = 1.0;

		++it;
	}

	vd.map();

	auto it2 = vd.getIterator();

	while (it2.isNext())
	{
		auto key = it2.get();

		float f = vd.getPos(key)[0];
		BOOST_REQUIRE_EQUAL(f, 0.0);
		f = vd.getPos(key)[1];
		BOOST_REQUIRE_EQUAL(f, 0.0);
		f = vd.getPos(key)[2];
		BOOST_REQUIRE_EQUAL(f, 0.0);

		++it2;
	}
}


BOOST_AUTO_TEST_CASE( vector_dist_not_periodic_map )
{
	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

	// factor
	float factor = pow(create_vcluster().getProcessingUnits()/2.0f,1.0f/3.0f);

	// ghost
	Ghost<3,float> ghost(0.05 / factor);

	// Distributed vector
	vector_dist<3,float, Point_test<float> > vd(1,box,bc,ghost);

	// put particles al 1.0, check that they go to 0.0

	auto it = vd.getIterator();

	while (it.isNext())
	{
		auto key = it.get();

		vd.getPos(key)[0] = 1.0;
		vd.getPos(key)[1] = 1.0;
		vd.getPos(key)[2] = 1.0;

		++it;
	}

	vd.map();

	auto it2 = vd.getIterator();

	while (it2.isNext())
	{
		auto key = it2.get();

		float f = vd.getPos(key)[0];
		BOOST_REQUIRE_EQUAL(f, 1.0);
		f = vd.getPos(key)[1];
		BOOST_REQUIRE_EQUAL(f, 1.0);
		f = vd.getPos(key)[2];
		BOOST_REQUIRE_EQUAL(f, 1.0);

		++it2;
	}
}

BOOST_AUTO_TEST_CASE( vector_dist_out_of_bound_policy )
{
	Vcluster<> & v_cl = create_vcluster();

	if (v_cl.getProcessingUnits() > 8)
		return;

	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

	// factor
	float factor = pow(create_vcluster().getProcessingUnits()/2.0f,1.0f/3.0f);

	// ghost
	Ghost<3,float> ghost(0.05 / factor);

	// Distributed vector
	vector_dist<3,float, Point_test<float> > vd(100,box,bc,ghost);

	// put particles at out of the boundary, they must be detected and and killed

	auto it = vd.getIterator();

	size_t cnt = 0;

	while (it.isNext())
	{
		auto key = it.get();

		if (cnt < 1)
		{
			vd.getPos(key)[0] = -0.06;
			vd.getPos(key)[1] = -0.06;
			vd.getPos(key)[2] = -0.06;
		}
		else
		{
			vd.getPos(key)[0] = 0.06;
			vd.getPos(key)[1] = 0.06;
			vd.getPos(key)[2] = 0.06;
		}

		cnt++;
		++it;
	}

	vd.map();

	// Particles out of the boundary are killed

	size_t cnt_l = vd.size_local();

	v_cl.sum(cnt_l);
	v_cl.execute();

	BOOST_REQUIRE_EQUAL(cnt_l,100-v_cl.getProcessingUnits());
}

void Test_interacting(Box<3,float> & box)
{
	Vcluster<> & v_cl = create_vcluster();

	if (v_cl.getProcessingUnits() > 8)
		return;

    // set the seed
	// create the random generator engine
	std::srand(v_cl.getProcessUnitID());
    std::default_random_engine eg;
    std::uniform_real_distribution<float> ud(-0.5f, 0.5f);

	size_t nsz[] = {0,32,4};

#ifdef TEST_COVERAGE_MODE
	nsz[0] = 5536 * v_cl.getProcessingUnits();
#else
	nsz[0] = 65536 * v_cl.getProcessingUnits();
#endif

	print_test_v("Testing 3D random walk interacting particles vector k=", nsz[0]);

	// 3D test
	for (size_t i = 0 ; i < 3 ; i++ )
	{
		size_t k = nsz[i];

		BOOST_TEST_CHECKPOINT( "Testing 3D random walk interacting particles vector k=" << k );

		// Boundary conditions
		size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

		// factor
		float factor = pow(create_vcluster().getProcessingUnits()/2.0f,1.0f/3.0f);

		// interaction radius
		float r_cut = 0.01 / factor;

		// ghost
		Ghost<3,float> ghost(r_cut);

		// Distributed vector
		vector_dist<3,float, Point_test<float> > vd(k,box,bc,ghost);

		auto it = vd.getIterator();

		while (it.isNext())
		{
			auto key = it.get();

			vd.getPos(key)[0] = ud(eg);
			vd.getPos(key)[1] = ud(eg);
			vd.getPos(key)[2] = ud(eg);

			++it;
		}

		vd.map();

		// 4 step random walk

		for (size_t j = 0 ; j < 4 ; j++)
		{
			auto it = vd.getDomainIterator();

			// Move the particles

			while (it.isNext())
			{
				auto key = it.get();

				vd.getPos(key)[0] += 0.02 * ud(eg);
				vd.getPos(key)[1] += 0.02 * ud(eg);
				vd.getPos(key)[2] += 0.02 * ud(eg);

				++it;
			}

			vd.map();
			vd.ghost_get<0>();

			// get the cell list with a cutoff radius

			bool error = false;

			auto NN = vd.getCellList(0.01 / factor);

			// iterate across the domain particle

			auto it2 = vd.getDomainIterator();

			while (it2.isNext())
			{
				auto p = it2.get();

				Point<3,float> xp = vd.getPos(p);

				auto Np = NN.getCellIterator(NN.getCell(xp));

				while (Np.isNext())
				{
					auto q = Np.get();

					// repulsive

					Point<3,float> xq = vd.getPos(q);
					Point<3,float> f = (xp - xq);

					float distance = f.norm();

					// Particle should be inside 2 * r_cut range

					if (distance > 2*r_cut*sqrt(2))
						error = true;

					++Np;
				}

				++it2;
			}

			// Error

			BOOST_REQUIRE_EQUAL(error,false);

			// Count the local particles and check that the total number is consistent
			size_t cnt = total_n_part_lc(vd,bc);

			BOOST_REQUIRE_EQUAL((size_t)k,cnt);
		}
	}
}

BOOST_AUTO_TEST_CASE( vector_dist_periodic_test_interacting_particles )
{
	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});
	Test_interacting(box);

	Box<3,float> box2({-0.5,-0.5,-0.5},{0.5,0.5,0.5});
	Test_interacting(box2);
}

BOOST_AUTO_TEST_CASE( vector_dist_grid_iterator )
{
#ifdef TEST_COVERAGE_MODE
	long int k = 32*32*32*create_vcluster().getProcessingUnits();
#else
	long int k = 64*64*64*create_vcluster().getProcessingUnits();
#endif
	k = std::pow(k, 1/3.);

	long int big_step = k / 30;
	big_step = (big_step == 0)?1:big_step;
	long int small_step = 21;

	print_test( "Testing vector grid iterator list k<=",k);

	// 3D test
	for ( ; k > 8*big_step ; k-= (k > 2*big_step)?big_step:small_step )
	{
		Vcluster<> & v_cl = create_vcluster();

		const size_t Ng = k;

		// we create a 128x128x128 Grid iterator
		size_t sz[3] = {Ng,Ng,Ng};

		Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

		// Boundary conditions
		size_t bc[3]={NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

		// ghost
		Ghost<3,float> ghost(1.0/(Ng-2));

		// Distributed vector
		vector_dist<3,float, Point_test<float> > vd(0,box,bc,ghost);

		// Put particles on a grid creating a Grid iterator
		auto it = vd.getGridIterator(sz);

		while (it.isNext())
		{
			vd.add();

			auto key = it.get();

			vd.getLastPos()[0] = key.get(0) * it.getSpacing(0);
			vd.getLastPos()[1] = key.get(1) * it.getSpacing(1);
			vd.getLastPos()[2] = key.get(2) * it.getSpacing(2);

			++it;
		}

		BOOST_REQUIRE_EQUAL(it.getSpacing(0),1.0f/(Ng-1));
		BOOST_REQUIRE_EQUAL(it.getSpacing(1),1.0f/(Ng-1));
		BOOST_REQUIRE_EQUAL(it.getSpacing(2),1.0f/(Ng-1));

		// distribute particles and sync ghost
		vd.map();


		// Check that the sum of all the particles is the grid size
		size_t total = vd.size_local();
		v_cl.sum(total);
		v_cl.execute();

		BOOST_REQUIRE_EQUAL(total,(Ng) * (Ng) * (Ng));
	}
}

BOOST_AUTO_TEST_CASE( vector_dist_cell_verlet_test )
{
#ifdef TEST_COVERAGE_MODE
	long int k = 16*16*16*create_vcluster().getProcessingUnits();
#else
	long int k = 64*64*64*create_vcluster().getProcessingUnits();
#endif
	k = std::pow(k, 1/3.);

	long int big_step = k / 30;
	big_step = (big_step == 0)?1:big_step;
	long int small_step = 21;

	print_test( "Testing cell and verlet list k<=",k);

	// 3D test
	for ( ; k > 8*big_step ; k-= (k > 2*big_step)?big_step:small_step )
	{
		Vcluster<> & v_cl = create_vcluster();

		const size_t Ng = k;

		// we create a 128x128x128 Grid iterator
		size_t sz[3] = {Ng,Ng,Ng};

		Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

		// Boundary conditions
		size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

		float spacing = 1.0/Ng;
		float first_dist = spacing;
		float second_dist = sqrt(2.0*spacing*spacing);
		float third_dist = sqrt(3.0 * spacing*spacing);

		// ghost
		Ghost<3,float> ghost(third_dist*1.1);

		// Distributed vector
		vector_dist<3,float, Point_test<float> > vd(0,box,bc,ghost);

		// Put particles on a grid creating a Grid iterator
		auto it = vd.getGridIterator(sz);

		while (it.isNext())
		{
			vd.add();

			auto key = it.get();

			vd.getLastPos()[0] = key.get(0) * it.getSpacing(0);
			vd.getLastPos()[1] = key.get(1) * it.getSpacing(1);
			vd.getLastPos()[2] = key.get(2) * it.getSpacing(2);

			++it;
		}

		BOOST_REQUIRE_EQUAL(it.getSpacing(0),1.0f/Ng);
		BOOST_REQUIRE_EQUAL(it.getSpacing(1),1.0f/Ng);
		BOOST_REQUIRE_EQUAL(it.getSpacing(2),1.0f/Ng);

		// distribute particles and sync ghost
		vd.map();

		// Check that the sum of all the particles is the grid size
		size_t total = vd.size_local();
		v_cl.sum(total);
		v_cl.execute();

		BOOST_REQUIRE_EQUAL(total,(Ng) * (Ng) * (Ng));

		vd.ghost_get<0>();

		// calculate the distance of the first, second and third neighborhood particle
		// Consider that they are on a regular grid

		// add a 5% to dist

		first_dist += first_dist * 0.05;
		second_dist += second_dist * 0.05;
		third_dist += third_dist * 0.05;

		// Create a verlet list for each particle

		VerletList<3,float,Mem_fast<>,shift<3,float>> verlet = vd.getVerlet(third_dist);

		bool correct = true;

		BOOST_REQUIRE_EQUAL(vd.size_local(),verlet.size());

		// for each particle
		for (size_t i = 0 ; i < verlet.size() ; i++)
		{
			// first NN
			size_t first_NN = 0;
			size_t second_NN = 0;
			size_t third_NN = 0;

			Point<3,float> p = vd.getPos(i);

			// for each neighborhood particle
			for (size_t j = 0 ; j < verlet.getNNPart(i) ; j++)
			{
				Point<3,float> q = vd.getPos(verlet.get(i,j));

				float dist = p.distance(q);

				if (dist <= first_dist)
					first_NN++;
				else if (dist <= second_dist)
					second_NN++;
				else
					third_NN++;
			}

			correct &= (first_NN == 7);
			correct &= (second_NN == 12);
			correct &= (third_NN == 8);
		}

		BOOST_REQUIRE_EQUAL(correct,true);
	}
}


BOOST_AUTO_TEST_CASE( vector_dist_periodic_map_list )
{
	Vcluster<> & v_cl = create_vcluster();

	if (v_cl.getProcessingUnits() > 3)
		return;

    // set the seed
	// create the random generator engine
	std::srand(v_cl.getProcessUnitID());
    std::default_random_engine eg;
    std::uniform_real_distribution<float> ud(0.0f, 1.0f);

#ifdef TEST_COVERAGE_MODE
    long int k = 24288 * v_cl.getProcessingUnits();
#else
    long int k = 524288 * v_cl.getProcessingUnits();
#endif

	long int big_step = k / 4;
	big_step = (big_step == 0)?1:big_step;

	print_test("Testing 3D periodic vector with map_list k=",k);
	BOOST_TEST_CHECKPOINT( "Testing 3D periodic vector with map_list k=" << k );

	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// factor
	float factor = pow(create_vcluster().getProcessingUnits()/2.0f,1.0f/3.0f);

	// ghost
	Ghost<3,float> ghost(0.05 / factor);

	// ghost2 (a little bigger because of round off error)
	Ghost<3,float> ghost2(0.05001 / factor);

	typedef  aggregate<float,float,std::list<int>,openfpm::vector<size_t>,openfpm::vector<Point_test<float>>> part_prop;

	// Distributed vector
	vector_dist<3,float, part_prop > vd(k,box,bc,ghost);

	auto it = vd.getIterator();

	while (it.isNext())
	{
		auto key = it.get();

		vd.getPos(key)[0] = ud(eg);
		vd.getPos(key)[1] = ud(eg);
		vd.getPos(key)[2] = ud(eg);

		// Fill some properties randomly

		vd.getProp<2>(key).push_back(1);
		vd.getProp<2>(key).push_back(2);
		vd.getProp<2>(key).push_back(3);
		vd.getProp<2>(key).push_back(4);

		vd.getProp<3>(key).add(1);
		vd.getProp<3>(key).add(2);
		vd.getProp<3>(key).add(3);
		vd.getProp<3>(key).add(4);

		vd.getProp<4>(key).add();
		vd.getProp<4>(key).add();
		vd.getProp<4>(key).add();
		vd.getProp<4>(key).add();

		++it;
	}

	vd.map_list<0,1>();

	// sync the ghost
	vd.ghost_get<0>();

	// Domain + ghost
	Box<3,float> dom_ext = box;
	dom_ext.enlarge(ghost2);

	// Iterate on all particles domain + ghost
	size_t l_cnt = 0;
	size_t nl_cnt = 0;
	size_t n_out = 0;

	auto it2 = vd.getIterator();
	count_local_n_local<3,vector_dist<3,float, part_prop>>(vd,it2,bc,box,dom_ext,l_cnt,nl_cnt,n_out);

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
	count_local_n_local<3, vector_dist<3,float,part_prop> >(vd,itg,bc,box,dom_ext,l_cnt,nl_cnt,n_out);

	// No particle on the ghost must be inside the domain
	BOOST_REQUIRE_EQUAL(l_cnt,0ul);

	// Ghost must be populated
	if (k > 524288)
	{
		BOOST_REQUIRE(nl_cnt != 0);
	}
}


BOOST_AUTO_TEST_CASE( vector_dist_ghost_with_ghost_buffering )
{
	Vcluster<> & v_cl = create_vcluster();

	if (v_cl.getProcessingUnits() > 3)
		return;

    // set the seed
	// create the random generator engine
	std::srand(v_cl.getProcessUnitID());
    std::default_random_engine eg;
    std::uniform_real_distribution<float> ud(0.0f, 1.0f);

#ifdef TEST_COVERAGE_MODE
    long int k = 24288 * v_cl.getProcessingUnits();
#else
    long int k = 524288 * v_cl.getProcessingUnits();
#endif

	long int big_step = k / 4;
	big_step = (big_step == 0)?1:big_step;

	print_test("Testing 3D periodic vector with ghost buffering k=",k);
	BOOST_TEST_CHECKPOINT( "Testing 3D periodic with ghost buffering k=" << k );

	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

	// ghost
	Ghost<3,float> ghost(0.1);

	typedef  aggregate<float,float,float> part_prop;

	// Distributed vector
	vector_dist<3,float, part_prop > vd(k,box,bc,ghost);

	auto it = vd.getIterator();

	while (it.isNext())
	{
		auto key = it.get();

		vd.getPos(key)[0] = ud(eg);
		vd.getPos(key)[1] = ud(eg);
		vd.getPos(key)[2] = ud(eg);

		// Fill some properties randomly

		vd.getProp<0>(key) = 0.0;
		vd.getProp<1>(key) = vd.getPos(key)[0];
		vd.getProp<2>(key) = vd.getPos(key)[0]*vd.getPos(key)[0];

		++it;
	}

	vd.map();

	// sync the ghost
	vd.ghost_get<0,1,2>();

	bool ret = true;
	auto it2 = vd.getGhostIterator();
	while (it2.isNext())
	{
		auto key = it2.get();

		ret &= vd.getProp<1>(key) == vd.getPos(key)[0];
		ret &= vd.getProp<2>(key) == vd.getPos(key)[0] * vd.getPos(key)[0];

		++it2;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	for (size_t i = 0 ; i < 10 ; i++)
	{
		auto it = vd.getDomainIterator();

		while (it.isNext())
		{
			auto key = it.get();

			vd.getPos(key)[1] = ud(eg);
			vd.getPos(key)[2] = ud(eg);

			// Fill some properties randomly

			vd.getProp<0>(key) = i;

			++it;
		}

		if (i % 2 == 0)
			vd.ghost_get<0>(SKIP_LABELLING);
		else
			vd.ghost_get<0>(SKIP_LABELLING | NO_CHANGE_ELEMENTS );

		auto it2 = vd.getGhostIterator();
		bool ret = true;

		while (it2.isNext())
		{
			auto key = it2.get();

			ret &= vd.getProp<0>(key) == i;
			ret &= vd.getProp<1>(key) == vd.getPos(key)[0];
			ret &= vd.getProp<2>(key) == vd.getPos(key)[0] * vd.getPos(key)[0];

			++it2;
		}

		BOOST_REQUIRE_EQUAL(ret,true);
	}

	vd.map();
	vd.ghost_get<0,1,2>();

	// shift the particle position by 1.0

	it = vd.getGhostIterator();
	while (it.isNext())
	{
		// Particle p
		auto p = it.get();

		// we shift down he particles
		vd.getPos(p)[0] = 10.0;

		// we shift
		vd.getPos(p)[1] = 17.0;

		// next particle
		++it;
	}

	for (size_t i = 0 ; i < 10 ; i++)
	{
		auto it = vd.getDomainIterator();

		while (it.isNext())
		{
			auto key = it.get();

			vd.getPos(key)[1] = ud(eg);
			vd.getPos(key)[2] = ud(eg);

			// Fill some properties randomly

			vd.getProp<0>(key) = i;
			vd.getProp<1>(key) = vd.getPos(key)[0];
			vd.getProp<2>(key) = vd.getPos(key)[0]*vd.getPos(key)[0];

			++it;
		}

		vd.ghost_get<0>(SKIP_LABELLING | NO_POSITION);

		auto it2 = vd.getGhostIterator();
		bool ret = true;

		while (it2.isNext())
		{
			// Particle p
			auto p = it.get();

			ret &= vd.getPos(p)[0] == 10.0;

			// we shift
			ret &= vd.getPos(p)[1] == 17.0;

			// next particle
			++it2;
		}

		BOOST_REQUIRE_EQUAL(ret,true);
	}
}



BOOST_AUTO_TEST_CASE( vector_dist_ghost_put )
{
	Vcluster<> & v_cl = create_vcluster();

	long int k = 25*25*25*create_vcluster().getProcessingUnits();
	k = std::pow(k, 1/3.);

	if (v_cl.getProcessingUnits() > 48)
		return;

	print_test("Testing 3D periodic ghost put k=",k);
	BOOST_TEST_CHECKPOINT( "Testing 3D periodic ghost put k=" << k );

	long int big_step = k / 30;
	big_step = (big_step == 0)?1:big_step;
	long int small_step = 21;

	// 3D test
	for ( ; k >= 2 ; k-= (k > 2*big_step)?big_step:small_step )
	{
		float r_cut = 1.3 / k;
		float r_g = 1.5 / k;

		Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

		// Boundary conditions
		size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

		// ghost
		Ghost<3,float> ghost(r_g);

		typedef  aggregate<float> part_prop;

		// Distributed vector
		vector_dist<3,float, part_prop > vd(0,box,bc,ghost);

		auto it = vd.getGridIterator({(size_t)k,(size_t)k,(size_t)k});

		while (it.isNext())
		{
			auto key = it.get();

			vd.add();

			vd.getLastPosWrite()[0] = key.get(0)*it.getSpacing(0);
			vd.getLastPosWrite()[1] = key.get(1)*it.getSpacing(1);
			vd.getLastPosWrite()[2] = key.get(2)*it.getSpacing(2);

			// Fill some properties randomly

			vd.getLastPropWrite<0>() = 0.0;

			++it;
		}

		vd.map();

		// sync the ghost
		vd.ghost_get<0>();

		{
			auto NN = vd.getCellList(r_cut);
			float a = 1.0f*k*k;

			// run trough all the particles + ghost

			auto it2 = vd.getDomainIterator();

			while (it2.isNext())
			{
				// particle p
				auto p = it2.get();
				Point<3,float> xp = vd.getPos(p);

				// Get an iterator over the neighborhood particles of p
				auto Np = NN.getNNIterator(NN.getCell(xp));

				// For each neighborhood particle ...
				while (Np.isNext())
				{
					auto q = Np.get();
					Point<3,float> xq = vd.getPosRead(q);

					float dist = xp.distance(xq);

					if (dist < r_cut)
						vd.getPropWrite<0>(q) += a*(-dist*dist+r_cut*r_cut);

					++Np;
				}

				++it2;
			}

			vd.ghost_put<add_,0>();

			bool ret = true;
			auto it3 = vd.getDomainIterator();

			float constant = vd.getProp<0>(it3.get());
			float eps = 0.001;

			while (it3.isNext())
			{
				float constant2 = vd.getProp<0>(it3.get());
				if (fabs(constant - constant2)/constant > eps)
				{
					Point<3,float> p = vd.getPosRead(it3.get());

					std::cout << p.toString() << "    " <<  constant2 << "/" << constant << "    " << v_cl.getProcessUnitID() << std::endl;
					ret = false;
					break;
				}

				++it3;
			}
			BOOST_REQUIRE_EQUAL(ret,true);
		}

		auto itp = vd.getDomainAndGhostIterator();
		while (itp.isNext())
		{
			auto key = itp.get();

			vd.getPropWrite<0>(key) = 0.0;

			++itp;
		}

		{
			auto NN = vd.getCellList(r_cut);
			float a = 1.0f*k*k;

			// run trough all the particles + ghost

			auto it2 = vd.getDomainIterator();

			while (it2.isNext())
			{
				// particle p
				auto p = it2.get();
				Point<3,float> xp = vd.getPosRead(p);

				// Get an iterator over the neighborhood particles of p
				auto Np = NN.getNNIterator(NN.getCell(xp));

				// For each neighborhood particle ...
				while (Np.isNext())
				{
					auto q = Np.get();
					Point<3,float> xq = vd.getPosRead(q);

					float dist = xp.distance(xq);

					if (dist < r_cut)
						vd.getPropWrite<0>(q) += a*(-dist*dist+r_cut*r_cut);

					++Np;
				}

				++it2;
			}

			vd.ghost_put<add_,0>();

			bool ret = true;
			auto it3 = vd.getDomainIterator();

			float constant = vd.getPropRead<0>(it3.get());
			float eps = 0.001;

			while (it3.isNext())
			{
				float constant2 = vd.getPropRead<0>(it3.get());
				if (fabs(constant - constant2)/constant > eps)
				{
					Point<3,float> p = vd.getPosRead(it3.get());

					std::cout << p.toString() << "    " <<  constant2 << "/" << constant << "    " << v_cl.getProcessUnitID() << std::endl;
					ret = false;
					break;
				}

				++it3;
			}
			BOOST_REQUIRE_EQUAL(ret,true);
		}
	}
}

BOOST_AUTO_TEST_CASE( vector_fixing_noposition_and_keep_prop )
{
	Vcluster<> & v_cl = create_vcluster();

	if (v_cl.getProcessingUnits() > 48)
		return;

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// Box
	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// ghost
	Ghost<3,float> ghost(0.1);

	vector_dist<3,float, aggregate<double,double>> vd(4096,box,bc,ghost);

	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		vd.getPos(key)[0] = ((double)rand())/RAND_MAX;
		vd.getPos(key)[1] = ((double)rand())/RAND_MAX;
		vd.getPos(key)[2] = ((double)rand())/RAND_MAX;

		++it;
	}

	vd.map();

	vd.ghost_get<>();
	size_t local = vd.getPosVector().size();

	vd.ghost_get<>(KEEP_PROPERTIES | NO_POSITION);

	size_t local2 = vd.getPosVector().size();

	BOOST_REQUIRE_EQUAL(local,local2);

	// Check now that map reset

	vd.map();

	local = vd.getPosVector().size();
	BOOST_REQUIRE_EQUAL(local,vd.size_local());
	vd.ghost_get<>(KEEP_PROPERTIES  | NO_POSITION);

	local2 = vd.getPosVector().size();

	BOOST_REQUIRE_EQUAL(local,local2);

	vd.ghost_get<>(KEEP_PROPERTIES);
	BOOST_REQUIRE_EQUAL(local,vd.getPosVector().size());
	BOOST_REQUIRE_EQUAL(vd.getPropVector().size(),local);

	vd.ghost_get<0>(KEEP_PROPERTIES);
	BOOST_REQUIRE_EQUAL(local,vd.getPosVector().size());
	BOOST_REQUIRE_EQUAL(vd.getPropVector().size(),local);
}


BOOST_AUTO_TEST_CASE( vector_of_vector_dist )
{
	Vcluster<> & v_cl = create_vcluster();

	if (v_cl.getProcessingUnits() > 48)
		return;

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// Box
	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// ghost
	Ghost<3,float> ghost(0.1);

	openfpm::vector< vector_dist<3,float, aggregate<double,double>> > phases;

	// first phase
	phases.add( vector_dist<3,float, aggregate<double,double>>(4096,box,bc,ghost) );

	// The other 3 phases
	phases.add( vector_dist<3,float, aggregate<double,double>>(phases.get(0).getDecomposition(),4096) );
	phases.add( vector_dist<3,float, aggregate<double,double>>(phases.get(0).getDecomposition(),4096) );
	phases.add( vector_dist<3,float, aggregate<double,double>>(phases.get(0).getDecomposition(),4096) );

	phases.get(0).map();
	phases.get(0).ghost_get<>();
	phases.get(1).map();
	phases.get(1).ghost_get<>();
	phases.get(2).map();
	phases.get(2).ghost_get<>();
	phases.get(3).map();
	phases.get(3).ghost_get<>();

	size_t cnt = 0;

	for (size_t i = 0 ; i < phases.size() ; i++)
		cnt += phases.get(i).size_local();

	v_cl.sum(cnt);
	v_cl.execute();

	BOOST_REQUIRE_EQUAL(cnt,4*4096ul);
}

BOOST_AUTO_TEST_CASE( vector_high_dimension )
{
	// Here we define our domain a 2D box with internals from 0 to 1.0 for x and y
	Box<10,double> domain;

	for (size_t i = 0 ; i < 10 ; i++)
	{
		domain.setLow(i,0.0);
		domain.setHigh(i,1.0);
	}

	// Here we define the boundary conditions of our problem
	size_t bc[10];
	for (size_t i = 0 ; i < 10 ; i++)
    {bc[i] = NON_PERIODIC;};

	// extended boundary around the domain, and the processor domain
	Ghost<10,double> g(0.0);

	// we check if the constructor does not stuck
	vector_dist<10,double, aggregate<double,double[10]> > vd(16,domain,bc,g);
}

BOOST_AUTO_TEST_CASE ( vector_of_cell_list_compile_test )
{
	auto & v_cl = create_vcluster();

    // set the seed
	// create the random generator engine
	std::srand(v_cl.getProcessUnitID());
    std::default_random_engine eg;
    std::uniform_real_distribution<float> ud(0.0f, 1.0f);

	Box<3,double> domain({0.0,0.0,0.0},{1.0,1.0,1.0});
	Ghost<3,double> g(0.1);
	size_t bc[3] = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

	vector_dist<3,double,aggregate<float,float[3]>> vd(100,domain,bc,g);

	auto it = vd.getIterator();

	while (it.isNext())
	{
		auto key = it.get();

		vd.getPos(key)[0] = ud(eg);
		vd.getPos(key)[1] = ud(eg);
		vd.getPos(key)[2] = ud(eg);

		++it;
	}

	vd.map();

	std::vector<decltype(vd.getCellList(0.1))> vector_of_celllist;

	typedef vector_dist<3,double,aggregate<float,float[3]>> my_particles;
	std::vector<decltype(std::declval<my_particles>().getCellList(0.0))> vector_of_celllist2;

	vector_of_celllist.push_back(vd.getCellList(0.1));

	vector_of_celllist2.push_back(vd.getCellList(0.1));
}


BOOST_AUTO_TEST_SUITE_END()

