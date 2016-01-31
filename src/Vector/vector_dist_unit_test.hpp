/*
 * vector_dist_unit_test.hpp
 *
 *  Created on: Mar 6, 2015
 *      Author: Pietro Incardona
 */

#ifndef VECTOR_DIST_UNIT_TEST_HPP_
#define VECTOR_DIST_UNIT_TEST_HPP_

#include <random>
#include "Vector/vector_dist.hpp"

/*! \brief Count the total number of particles
 *
 * \param vd distributed vector
 * \param bc boundary conditions
 *
 */
template<unsigned int dim> size_t total_n_part_lc(vector_dist<dim,float, Point_test<float>, CartDecomposition<dim,float> > & vd, size_t (& bc)[dim])
{
	typedef Point<dim,float> s;

	Vcluster & v_cl = vd.getVC();
	auto it2 = vd.getDomainIterator();
	const CartDecomposition<3,float> & ct = vd.getDecomposition();

	bool noOut = true;

	size_t cnt = 0;
	while (it2.isNext())
	{
		auto key = it2.get();

		noOut &= ct.isLocal(vd.template getPos<s::x>(key));

		cnt++;

		++it2;
	}

	BOOST_REQUIRE_EQUAL(noOut,true);

	//
	v_cl.sum(cnt);
	v_cl.execute();

	return cnt;
}

/*! \brief Count local and non local
 *
 * \param vd distributed vector
 * \param it iterator
 * \param bc boundary conditions
 * \param box domain box
 * \param dom_ext domain + ghost box
 * \param l_cnt local particles counter
 * \param nl_cnt non local particles counter
 * \param n_out out of domain + ghost particles counter
 *
 */
template<unsigned int dim> inline void count_local_n_local(vector_dist<dim,float, Point_test<float>, CartDecomposition<dim,float> > & vd, vector_dist_iterator & it, size_t (& bc)[dim] , Box<dim,float> & box, Box<dim,float> & dom_ext, size_t & l_cnt, size_t & nl_cnt, size_t & n_out)
{
	typedef Point<dim,float> s;
	const CartDecomposition<dim,float> & ct = vd.getDecomposition();

	while (it.isNext())
	{
		auto key = it.get();
		// Check if it is in the domain
		if (box.isInsideNP(vd.template getPos<s::x>(key)) == true)
		{
			// Check if local
			if (ct.isLocalBC(vd.template getPos<s::x>(key),bc) == true)
				l_cnt++;
			else
				nl_cnt++;
		}
		else
		{
			nl_cnt++;
		}

		// Check that all particles are inside the Domain + Ghost part
		if (dom_ext.isInside(vd.template getPos<s::x>(key)) == false)
				n_out++;

		++it;
	}
}

BOOST_AUTO_TEST_SUITE( vector_dist_test )

BOOST_AUTO_TEST_CASE( vector_dist_ghost )
{
	// Communication object
	Vcluster & v_cl = *global_v_cluster;

	typedef Point_test<float> p;
	typedef Point<2,float> s;

	// Get the default minimum number of sub-sub-domain per processor (granularity of the decomposition)
	size_t n_sub = vector_dist<2,float, Point_test<float>, CartDecomposition<2,float> >::getDefaultNsubsub() * v_cl.getProcessingUnits();
	// Convert the request of having a minimum n_sub number of sub-sub domain into grid decompsition of the space
	size_t sz = CartDecomposition<2,float>::getDefaultGrid(n_sub);

	//! [Create a vector of elements distributed on a grid like way]

	Box<2,float> box({0.0,0.0},{1.0,1.0});
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
	Point<2,float> spacing = box.getP2();
	spacing = spacing / g_div;

	// middle spacing
	Point<2,float> m_spacing = spacing / 2;

	// set the ghost based on the radius cut off (make just a little bit smaller than the spacing)
	Ghost<2,float> g(spacing.get(0) - spacing .get(0) * 0.0001);

	// Boundary conditions
	size_t bc[2]={NON_PERIODIC,NON_PERIODIC};

	// Vector of particles
	vector_dist<2,float, Point_test<float>, CartDecomposition<2,float> > vd(g_info.size(),box,bc,g);

	// size_t
	size_t cobj = 0;

	grid_key_dx_iterator_sp<2> it(g_info,offset,offset+p_np-1);
	auto v_it = vd.getIterator();

	while (v_it.isNext() && it.isNext())
	{
		auto key = it.get();
		auto key_v = v_it.get();

		// set the particle position

		vd.template getPos<s::x>(key_v)[0] = key.get(0) * spacing[0] + m_spacing[0];
		vd.template getPos<s::x>(key_v)[1] = key.get(1) * spacing[1] + m_spacing[1];

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
		vd.template getProp<p::s>(key) = vd.getPos<s::x>(key)[0] + vd.getPos<s::x>(key)[1] * 16;
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
		BOOST_REQUIRE_EQUAL(vd.getPos<s::x>(key)[0] + vd.getPos<s::x>(key)[1] * 16,vd.template getProp<p::s>(key));

		bool is_in = false;
		size_t b = 0;
		size_t lb = 0;

		// check if the received data are in one of the ghost boxes
		for ( ; b < dec.getNEGhostBox() ; b++)
		{
			if (dec.getEGhostBox(b).isInside(vd.getPos<s::x>(key)) == true )
			{
				is_in = true;

				// Add
				vb.get(b)++;
				lb = b;
			}
		}
		BOOST_REQUIRE_EQUAL(is_in,true);

		// Check that the particle come from the correct processor
		BOOST_REQUIRE_EQUAL(vd.getProp<p::v>(key)[0],dec.getEGhostBoxProcessor(lb));

		n_part++;
		++g_it;
	}

	BOOST_REQUIRE(n_part != 0);

    CellDecomposer_sm<2,float> cd(SpaceBox<2,float>(box),g_div,0);

	for (size_t i = 0 ; i < vb.size() ; i++)
	{
		// Calculate how many particle should be in the box
		size_t n_point = cd.getGridPoints(dec.getEGhostBox(i)).getVolumeKey();

		if (n_point != vb.get(i))
		{
			std::cout << n_point << "  " << dec.getEGhostBoxProcessor(i) << "  " << v_cl.getProcessUnitID() << dec.getEGhostBox(i).toString() << "\n";
		}
		//BOOST_REQUIRE_EQUAL(n_point,vb.get(i));
	}
}

void print_test_v(std::string test, size_t sz)
{
	if (global_v_cluster->getProcessUnitID() == 0)
		std::cout << test << " " << sz << "\n";
}

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

BOOST_AUTO_TEST_CASE( vector_dist_iterator_test_use_2d )
{
	typedef Point<2,float> s;

	Vcluster & v_cl = *global_v_cluster;

    // set the seed
	// create the random generator engine
	std::srand(v_cl.getProcessUnitID());
    std::default_random_engine eg;
    std::uniform_real_distribution<float> ud(0.0f, 1.0f);

    long int k = 524288 * v_cl.getProcessingUnits();

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

		vector_dist<2,float, Point_test<float>, CartDecomposition<2,float> > vd(k,box,bc,Ghost<2,float>(0.0));

		auto it = vd.getIterator();

		while (it.isNext())
		{
			auto key = it.get();

			vd.template getPos<s::x>(key)[0] = ud(eg);
			vd.template getPos<s::x>(key)[1] = ud(eg);

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
			BOOST_REQUIRE_EQUAL(ct.isLocal(vd.template getPos<s::x>(key)),true);

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
	typedef Point<3,float> s;

	Vcluster & v_cl = *global_v_cluster;

    // set the seed
	// create the random generator engine
	std::srand(v_cl.getProcessUnitID());
    std::default_random_engine eg;
    std::uniform_real_distribution<float> ud(0.0f, 1.0f);

    long int k = 524288 * v_cl.getProcessingUnits();

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

		vector_dist<3,float, Point_test<float>, CartDecomposition<3,float> > vd(k,box,bc,Ghost<3,float>(0.0));

		auto it = vd.getIterator();

		while (it.isNext())
		{
			auto key = it.get();

			vd.template getPos<s::x>(key)[0] = ud(eg);
			vd.template getPos<s::x>(key)[1] = ud(eg);
			vd.template getPos<s::x>(key)[2] = ud(eg);

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
			BOOST_REQUIRE_EQUAL(ct.isLocal(vd.template getPos<s::x>(key)),true);

			cnt++;

			++it2;
		}

		//
		v_cl.sum(cnt);
		v_cl.execute();
		BOOST_REQUIRE_EQUAL(cnt,(size_t)k);
	}
}

BOOST_AUTO_TEST_CASE( vector_dist_periodic_test_use_2d )
{
	typedef Point<2,float> s;

	Vcluster & v_cl = *global_v_cluster;

    // set the seed
	// create the random generator engine
	std::srand(v_cl.getProcessUnitID());
    std::default_random_engine eg;
    std::uniform_real_distribution<float> ud(0.0f, 1.0f);

    long int k = 524288 * v_cl.getProcessingUnits();

	long int big_step = k / 4;
	big_step = (big_step == 0)?1:big_step;

	print_test_v( "Testing 2D periodic vector k<=",k);

	// 2D test
	for ( ; k >= 2 ; k-= decrement(k,big_step) )
	{
		BOOST_TEST_CHECKPOINT( "Testing 2D periodic vector k=" << k );

		//! [Create a vector of random elements on each processor 2D]

		Box<2,float> box({0.0,0.0},{1.0,1.0});

		// Boundary conditions
		size_t bc[2]={PERIODIC,PERIODIC};

		// factor
		float factor = pow(global_v_cluster->getProcessingUnits()/2.0f,1.0f/3.0f);

		// ghost
		Ghost<2,float> ghost(0.01 / factor);

		// ghost2 (a little bigger because of round off error)
		Ghost<2,float> ghost2(0.05001 / factor);

		// Distributed vector
		vector_dist<2,float, Point_test<float>, CartDecomposition<2,float> > vd(k,box,bc,ghost);

		auto it = vd.getIterator();

		while (it.isNext())
		{
			auto key = it.get();

			vd.template getPos<s::x>(key)[0] = ud(eg);
			vd.template getPos<s::x>(key)[1] = ud(eg);

			++it;
		}

		vd.map();

		// sync the ghost, only the property zero
		vd.ghost_get<0>();

		//! [Create a vector of random elements on each processor 2D]

		// Domain + ghost box
		Box<2,float> dom_ext = box;
		dom_ext.enlarge(ghost2);

		// Iterate on all particles domain + ghost
		size_t l_cnt = 0;
		size_t nl_cnt = 0;
		size_t n_out = 0;


		auto it2 = vd.getIterator();
		count_local_n_local(vd,it2,bc,box,dom_ext,l_cnt,nl_cnt,n_out);

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
		count_local_n_local(vd,itg,bc,box,dom_ext,l_cnt,nl_cnt,n_out);

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
	typedef Point<3,float> s;

	Vcluster & v_cl = *global_v_cluster;

    // set the seed
	// create the random generator engine
	std::srand(v_cl.getProcessUnitID());
    std::default_random_engine eg;
    std::uniform_real_distribution<float> ud(0.0f, 1.0f);

    long int k = 524288 * v_cl.getProcessingUnits();

	long int big_step = k / 4;
	big_step = (big_step == 0)?1:big_step;

	print_test_v( "Testing 3D periodic vector k<=",k);

	// 3D test
	for ( ; k >= 2 ; k-= decrement(k,big_step) )
	{
		BOOST_TEST_CHECKPOINT( "Testing 3D periodic vector k=" << k );

		//! [Create a vector of random elements on each processor 3D]

		Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

		// Boundary conditions
		size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

		// factor
		float factor = pow(global_v_cluster->getProcessingUnits()/2.0f,1.0f/3.0f);

		// ghost
		Ghost<3,float> ghost(0.05 / factor);

		// ghost2 (a little bigger because of round off error)
		Ghost<3,float> ghost2(0.05001 / factor);

		// Distributed vector
		vector_dist<3,float, Point_test<float>, CartDecomposition<3,float> > vd(k,box,bc,ghost);

		auto it = vd.getIterator();

		while (it.isNext())
		{
			auto key = it.get();

			vd.template getPos<s::x>(key)[0] = ud(eg);
			vd.template getPos<s::x>(key)[1] = ud(eg);
			vd.template getPos<s::x>(key)[2] = ud(eg);

			++it;
		}

		vd.map();

		// sync the ghost
		vd.ghost_get<0>();

		//! [Create a vector of random elements on each processor 3D]

		// Domain + ghost
		Box<3,float> dom_ext = box;
		dom_ext.enlarge(ghost2);

		// Iterate on all particles domain + ghost
		size_t l_cnt = 0;
		size_t nl_cnt = 0;
		size_t n_out = 0;

		auto it2 = vd.getIterator();
		count_local_n_local(vd,it2,bc,box,dom_ext,l_cnt,nl_cnt,n_out);

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
		count_local_n_local(vd,itg,bc,box,dom_ext,l_cnt,nl_cnt,n_out);

		// No particle on the ghost must be inside the domain
		BOOST_REQUIRE_EQUAL(l_cnt,0ul);

		// Ghost must be populated
		if (k > 524288)
		{
			BOOST_REQUIRE(nl_cnt != 0);
		}
	}
}

BOOST_AUTO_TEST_CASE( vector_dist_periodic_test_random_walk )
{
	typedef Point<3,float> s;

	Vcluster & v_cl = *global_v_cluster;

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
		float factor = pow(global_v_cluster->getProcessingUnits()/2.0f,1.0f/3.0f);

		// ghost
		Ghost<3,float> ghost(0.01 / factor);

		// Distributed vector
		vector_dist<3,float, Point_test<float>, CartDecomposition<3,float> > vd(k,box,bc,ghost);

		auto it = vd.getIterator();

		while (it.isNext())
		{
			auto key = it.get();

			vd.template getPos<s::x>(key)[0] = ud(eg);
			vd.template getPos<s::x>(key)[1] = ud(eg);
			vd.template getPos<s::x>(key)[2] = ud(eg);

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

				vd.template getPos<s::x>(key)[0] += 0.02 * ud(eg);
				vd.template getPos<s::x>(key)[1] += 0.02 * ud(eg);
				vd.template getPos<s::x>(key)[2] += 0.02 * ud(eg);

				++it;
			}

			vd.map();

			vd.ghost_get<0>();

			// Count the local particles and check that the total number is consistent
			size_t cnt = total_n_part_lc(vd,bc);

			BOOST_REQUIRE_EQUAL((size_t)k,cnt);
		}
	}
}

BOOST_AUTO_TEST_CASE( vector_dist_periodic_map )
{
	typedef Point<3,float> s;

	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// factor
	float factor = pow(global_v_cluster->getProcessingUnits()/2.0f,1.0f/3.0f);

	// ghost
	Ghost<3,float> ghost(0.05 / factor);

	// Distributed vector
	vector_dist<3,float, Point_test<float>, CartDecomposition<3,float> > vd(1,box,bc,ghost);

	// put particles al 1.0, check that they go to 0.0

	auto it = vd.getIterator();

	while (it.isNext())
	{
		auto key = it.get();

		vd.template getPos<s::x>(key)[0] = 1.0;
		vd.template getPos<s::x>(key)[1] = 1.0;
		vd.template getPos<s::x>(key)[2] = 1.0;

		++it;
	}

	vd.map();

	auto it2 = vd.getIterator();

	while (it2.isNext())
	{
		auto key = it2.get();

		float f = vd.template getPos<s::x>(key)[0];
		BOOST_REQUIRE_EQUAL(f, 0.0);
		f = vd.template getPos<s::x>(key)[1];
		BOOST_REQUIRE_EQUAL(f, 0.0);
		f = vd.template getPos<s::x>(key)[2];
		BOOST_REQUIRE_EQUAL(f, 0.0);

		++it2;
	}
}

BOOST_AUTO_TEST_CASE( vector_dist_not_periodic_map )
{
	typedef Point<3,float> s;

	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

	// factor
	float factor = pow(global_v_cluster->getProcessingUnits()/2.0f,1.0f/3.0f);

	// ghost
	Ghost<3,float> ghost(0.05 / factor);

	// Distributed vector
	vector_dist<3,float, Point_test<float>, CartDecomposition<3,float> > vd(1,box,bc,ghost);

	// put particles al 1.0, check that they go to 0.0

	auto it = vd.getIterator();

	while (it.isNext())
	{
		auto key = it.get();

		vd.template getPos<s::x>(key)[0] = 1.0;
		vd.template getPos<s::x>(key)[1] = 1.0;
		vd.template getPos<s::x>(key)[2] = 1.0;

		++it;
	}

	vd.map();

	auto it2 = vd.getIterator();

	while (it2.isNext())
	{
		auto key = it2.get();

		float f = vd.template getPos<s::x>(key)[0];
		BOOST_REQUIRE_EQUAL(f, 1.0);
		f = vd.template getPos<s::x>(key)[1];
		BOOST_REQUIRE_EQUAL(f, 1.0);
		f = vd.template getPos<s::x>(key)[2];
		BOOST_REQUIRE_EQUAL(f, 1.0);

		++it2;
	}
}

BOOST_AUTO_TEST_CASE( vector_dist_out_of_bound_policy )
{
	Vcluster & v_cl = *global_v_cluster;

	if (v_cl.getProcessingUnits() > 8)
		return;

	typedef Point<3,float> s;

	size_t bc[3]={NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

	vector_dist<3,float, Point_test<float>, CartDecomposition<3,float> > vd(100,box,bc,ghost);

	// put particles at out of the boundary, they must be detected and and killed

	auto it = vd.getIterator();

	size_t cnt = 0;

	while (it.isNext())
	{
		auto key = it.get();

		if (cnt < 1)
		{
			vd.template getPos<s::x>(key)[0] = -0.06;
			vd.template getPos<s::x>(key)[1] = -0.06;
			vd.template getPos<s::x>(key)[2] = -0.06;
		}
		else
		{
			vd.template getPos<s::x>(key)[0] = 0.06;
			vd.template getPos<s::x>(key)[1] = 0.06;
			vd.template getPos<s::x>(key)[2] = 0.06;
		}

		cnt++;
		++it;
	}

	// Particles out of the boundary are killed

	size_t cnt_l = vd.size_local();

	v_cl.sum(cnt_l);
	v_cl.execute();

	BOOST_REQUIRE_EQUAL(cnt_l,100-v_cl.getProcessingUnits());
}


BOOST_AUTO_TEST_CASE( vector_dist_periodic_test_interacting_particles )
{
	typedef Point<3,float> s;

	Vcluster & v_cl = *global_v_cluster;

	if (v_cl.getProcessingUnits() > 8)
		return;

	// set the seed
	// create the random generator engine
	std::srand(v_cl.getProcessUnitID());
	std::default_random_engine eg;
	std::uniform_real_distribution<float> ud(0.0f, 1.0f);

	size_t nsz[] = {0,32,4};
	nsz[0] = 65536 * v_cl.getProcessingUnits();

	print_test_v("Testing 3D random walk interacting particles vector k=", nsz[0]);

	// 3D test
	for (size_t i = 0 ; i < 3 ; i++ )
	{
		size_t k = nsz[i];

		BOOST_TEST_CHECKPOINT( "Testing 3D random walk interacting particles vector k=" << k );

		Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

		// Boundary conditions
		size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

		// factor
		float factor = pow(global_v_cluster->getProcessingUnits()/2.0f,1.0f/3.0f);

		// interaction radius
		float r_cut = 0.01 / factor;

		// ghost
		Ghost<3,float> ghost(r_cut);

		// Distributed vector
		vector_dist<3,float, Point_test<float>, CartDecomposition<3,float> > vd(k,box,bc,ghost);

		auto it = vd.getIterator();

		while (it.isNext())
		{
			auto key = it.get();

			vd.template getPos<s::x>(key)[0] = ud(eg);
			vd.template getPos<s::x>(key)[1] = ud(eg);
			vd.template getPos<s::x>(key)[2] = ud(eg);

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

				vd.template getPos<s::x>(key)[0] += 0.02 * ud(eg);
				vd.template getPos<s::x>(key)[1] += 0.02 * ud(eg);
				vd.template getPos<s::x>(key)[2] += 0.02 * ud(eg);

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

				Point<3,float> xp = vd.getPos<0>(p);

				auto Np = NN.getIterator(NN.getCell(vd.getPos<0>(p)));

				while (Np.isNext())
				{
					auto q = Np.get();

					// repulsive

					Point<3,float> xq = vd.getPos<0>(q);
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

BOOST_AUTO_TEST_CASE( vector_dist_cell_verlet_test )
{
	typedef Point<3,float> s;

	// we create a 128x128x128 Grid iterator
	size_t sz[3] = {128,128,128};
	size_t total = sz[0]*sz[1]*sz[2];

	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// factor
	float factor = pow(global_v_cluster->getProcessingUnits()/2.0f,1.0f/3.0f);

	// ghost
	Ghost<3,float> ghost(0.05 / factor);

	// Distributed vector
	vector_dist<3,float, Point_test<float>, CartDecomposition<3,float> > vd(total,box,bc,ghost);

	// Put particles on a grid creating a Grid iterator
	auto it = vd.getGridIterator(sz);
	auto it_p = vd.getDomainIterator();

	while (it_p.isNext())
	{
		auto key_p = it_p.get();
		auto key = it.get();

		vd.template getPos<s::x>(key_p)[0] = key.get(0) * it.getSpacing(0);
		vd.template getPos<s::x>(key_p)[1] = key.get(1) * it.getSpacing(1);
		vd.template getPos<s::x>(key_p)[2] = key.get(2) * it.getSpacing(2);

		++it;
		++it_p;
	}

	vd.map();

	// calculate the distance of the first, second and third neighborhood particle
	// Consider that they are on a regular grid

	float spacing = it.getSpacing(0);
	float first_dist = spacing;
	float second_dist = sqrt(2.0*spacing*spacing);
	float third_dist = sqrt(3.0 * spacing*spacing);

	// add a 5% to dist

	first_dist += first_dist * 0.05;
	second_dist += second_dist * 0.05;
	third_dist += third_dist * 0.05;

	// Create a verlet list for each particle

	openfpm::vector<openfpm::vector<size_t>> verlet;
	vd.getVerlet(verlet);

	bool correct = true;

	// for each particle
	for (size_t i = 0 ; i < verlet.size() ; i++)
	{

		// first NN
		size_t first_NN = 0;
		size_t second_NN = 0;
		size_t third_NN = 0;

		Point<3,float> p = vd.getPos<0>(i);

		// for each neighborhood particle
		for (size_t j = 0 ; j < verlet.get(i).size() ; j++)
		{
			auto & NN = verlet.get(i);

			Point<3,float> q = vc.getPos<0>(NN.get(j));

			float dist = p.distance(q);

			if (dist <= first_dist)
				first_NN++;
			else if (dist <= second_dist)
				second_NN++;
			else
				third_NN++;
		}

		correct &= (first_NN == 6);
		correct &= (second_NN == 12);
		correct &= (third_NN = 8);
	}

	BOOST_REQUIRE_EQUAL(correct,true);
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* VECTOR_DIST_UNIT_TEST_HPP_ */
