#ifndef GRID_DIST_UNIT_TEST_HPP
#define GRID_DIST_UNIT_TEST_HPP

#include "grid_dist_id.hpp"
#include "data_type/scalar.hpp"

BOOST_AUTO_TEST_SUITE( grid_dist_id_test )

void print_test(std::string test, size_t sz)
{
	if (global_v_cluster->getProcessUnitID() == 0)
		std::cout << test << " " << sz << "\n";
}

BOOST_AUTO_TEST_CASE( grid_dist_id_domain_grid_unit_converter_test)
{
	// Domain
	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	// Initialize the global VCluster
	init_global_v_cluster(&boost::unit_test::framework::master_test_suite().argc,&boost::unit_test::framework::master_test_suite().argv);

	Vcluster & v_cl = *global_v_cluster;

	// Skip this test on big scale
	if (v_cl.getProcessingUnits() >= 32)
		return;

	// Test several grid dimensions

	for (size_t k = 1024 ; k >= 2 ; k--)
	{
		BOOST_TEST_CHECKPOINT( "Testing grid k=" << k );

		// grid size
		size_t sz[2];
		sz[0] = k;
		sz[1] = k;

		// Ghost
		Ghost<2,float> g(0.01);

		// Distributed grid with id decomposition
		grid_dist_id<2, float, scalar<float>, CartDecomposition<2,float>> g_dist(sz,domain,g);

		// get the decomposition
		auto & dec = g_dist.getDecomposition();

		// check the consistency of the decomposition
		bool val = dec.check_consistency();
		BOOST_REQUIRE_EQUAL(val,true);

		// for each local volume
		// Get the number of local grid needed
		size_t n_grid = dec.getNLocalHyperCube();

		size_t vol = 0;

		// Allocate the grids
		for (size_t i = 0 ; i < n_grid ; i++)
		{
			// Get the local hyper-cube
			SpaceBox<2,float> sub = dec.getLocalHyperCube(i);

			Box<2,size_t> g_box = g_dist.getCellDecomposer().convertDomainSpaceIntoGridUnits(sub);

			vol += g_box.getVolumeKey();
		}

		v_cl.sum(vol);
		v_cl.execute();

		BOOST_REQUIRE_EQUAL(vol,sz[0]*sz[1]);
	}
}

void Test2D_sub(const Box<2,float> & domain, long int k)
{
	long int big_step = k / 30;
	big_step = (big_step == 0)?1:big_step;
	long int small_step = 21;

	// this test is only performed when the number of processor is <= 32
	if (global_v_cluster->getProcessingUnits() > 32)
		return;

	print_test( "Testing 2D grid sub iterator k<=",k);

	// 2D test
	for ( ; k >= 2 ; k-= (k > 2*big_step)?big_step:small_step )
	{
		BOOST_TEST_CHECKPOINT( "Testing 2D grid k=" << k );

		//! [Create and access a distributed grid]

		// grid size
		size_t sz[2];
		sz[0] = k;
		sz[1] = k;

		float factor = pow(global_v_cluster->getProcessingUnits()/2.0f,1.0f/2.0f);

		// Ghost
		Ghost<2,float> g(0.01 / factor);

		// Distributed grid with id decomposition
		grid_dist_id<2, float, scalar<float>, CartDecomposition<2,float>> g_dist(sz,domain,g);

		// check the consistency of the decomposition
		bool val = g_dist.getDecomposition().check_consistency();
		BOOST_REQUIRE_EQUAL(val,true);

		size_t count;

		// Grid sm
		grid_sm<2,void> info(sz);

		{
		//! [Usage of a sub_grid iterator]

		grid_key_dx<2> one(1,1);
		grid_key_dx<2> one_end(k-2,k-2);

		bool check = true;
		count = 0;

		// get the sub-domain iterator
		auto dom = g_dist.getSubDomainIterator(one,one_end);

		while (dom.isNext())
		{
			auto key = dom.get();
			auto key_g = g_dist.getGKey(key);

			// key_g should never be 1 or k-1
			check &= (key_g.get(0) == 0 || key_g.get(0) == k-1)?false:true;
			check &= (key_g.get(1) == 0 || key_g.get(1) == k-1)?false:true;

			g_dist.template get<0>(key) = info.LinId(key_g);

			// Count the point
			count++;

			++dom;
		}

		BOOST_REQUIRE_EQUAL(check,true);

		//! [Usage of a sub_grid iterator]

		}

		// Get the virtual cluster machine
		Vcluster & vcl = g_dist.getVC();

		// reduce
		vcl.sum(count);
		vcl.execute();

		// Check
		BOOST_REQUIRE_EQUAL(count,(size_t)(k-2)*(k-2));

		// check with a 1x1 square

		{

		grid_key_dx<2> one(k/2,k/2);
		grid_key_dx<2> one_end(k/2,k/2);

		count = 0;

		// get the sub-domain iterator
		auto dom = g_dist.getSubDomainIterator(one,one_end);

		while (dom.isNext())
		{
			auto key = dom.get();
			auto key_g = g_dist.getGKey(key);

			// key_g
			BOOST_REQUIRE_EQUAL(key_g.get(0),k/2);
			BOOST_REQUIRE_EQUAL(key_g.get(1),k/2);

			auto key_s_it = dom.getGKey(key);

			BOOST_REQUIRE_EQUAL(key_g.get(0),key_s_it.get(0));
			BOOST_REQUIRE_EQUAL(key_g.get(1),key_s_it.get(1));

			// Count the point
			count++;

			++dom;
		}

		// reduce
		vcl.sum(count);
		vcl.execute();

		BOOST_REQUIRE_EQUAL(count,1ul);
		}
	}
}

void Test2D(const Box<2,float> & domain, long int k)
{
	long int big_step = k / 30;
	big_step = (big_step == 0)?1:big_step;
	long int small_step = 21;

	print_test( "Testing 2D grid k<=",k);

	// 2D test
	for ( ; k >= 2 ; k-= (k > 2*big_step)?big_step:small_step )
	{
		BOOST_TEST_CHECKPOINT( "Testing 2D grid k=" << k );

		//! [Create and access a distributed grid]

		// grid size
		size_t sz[2];
		sz[0] = k;
		sz[1] = k;

		float factor = pow(global_v_cluster->getProcessingUnits()/2.0f,1.0f/2.0f);

		// Ghost
		Ghost<2,float> g(0.01 / factor);

		// Distributed grid with id decomposition
		grid_dist_id<2, float, scalar<float>, CartDecomposition<2,float>> g_dist(sz,domain,g);

		// check the consistency of the decomposition
		bool val = g_dist.getDecomposition().check_consistency();
		BOOST_REQUIRE_EQUAL(val,true);

		// Grid sm
		grid_sm<2,void> info(sz);

		// get the domain iterator
		size_t count = 0;

		auto dom = g_dist.getDomainIterator();

		while (dom.isNext())
		{
			auto key = dom.get();
			auto key_g = g_dist.getGKey(key);

			g_dist.template get<0>(key) = info.LinId(key_g);

			// Count the point
			count++;

			++dom;
		}

		//! [Create and access a distributed grid]

		// Get the virtual cluster machine
		Vcluster & vcl = g_dist.getVC();

		// reduce
		vcl.sum(count);
		vcl.execute();

		// Check
		BOOST_REQUIRE_EQUAL(count,(size_t)k*k);

		auto dom2 = g_dist.getDomainIterator();

		bool match = true;

		// check that the grid store the correct information
		while (dom2.isNext())
		{
			auto key = dom2.get();
			auto key_g = g_dist.getGKey(key);

			match &= (g_dist.template get<0>(key) == info.LinId(key_g))?true:false;

			++dom2;
		}

		BOOST_REQUIRE_EQUAL(match,true);

		g_dist.template ghost_get<0>();

		// check that the communication is correctly completed

		auto domg = g_dist.getDomainGhostIterator();

		// check that the grid with the ghost past store the correct information
		while (domg.isNext())
		{
			auto key = domg.get();
			auto key_g = g_dist.getGKey(key);

			// In this case the boundary condition are non periodic
			if (g_dist.isInside(key_g))
			{
				match &= (g_dist.template get<0>(key),info.LinId(key_g));
			}

			++domg;
		}
	}
}

void Test3D_sub(const Box<3,float> & domain, long int k)
{
	long int big_step = k / 30;
	big_step = (big_step == 0)?1:big_step;
	long int small_step = 21;

	// this test is only performed when the number of processor is <= 32
	if (global_v_cluster->getProcessingUnits() > 32)
		return;

	print_test( "Testing 3D grid sub k<=",k);

	// 3D test
	for ( ; k >= 2 ; k-= (k > 2*big_step)?big_step:small_step )
	{
		BOOST_TEST_CHECKPOINT( "Testing 3D grid sub k=" << k );

		// grid size
		size_t sz[3];
		sz[0] = k;
		sz[1] = k;
		sz[2] = k;

		// factor
		float factor = pow(global_v_cluster->getProcessingUnits()/2.0f,1.0f/3.0f);

		// Ghost
		Ghost<3,float> g(0.01 / factor);

		// Distributed grid with id decomposition
		grid_dist_id<3, float, scalar<float>, CartDecomposition<3,float>> g_dist(sz,domain,g);

		// check the consistency of the decomposition
		bool val = g_dist.getDecomposition().check_consistency();
		BOOST_REQUIRE_EQUAL(val,true);

		// Grid sm
		grid_sm<3,void> info(sz);

		// get the domain iterator
		size_t count = 0;

		grid_key_dx<3> one(1,1,1);
		grid_key_dx<3> one_end(k-2,k-2,k-2);

		// Sub-domain iterator
		auto dom = g_dist.getSubDomainIterator(one,one_end);

		while (dom.isNext())
		{
			auto key = dom.get();
			auto key_g = g_dist.getGKey(key);

			g_dist.template get<0>(key) = info.LinId(key_g);

			// Count the point
			count++;

			++dom;
		}

		// Get the virtual cluster machine
		Vcluster & vcl = g_dist.getVC();

		// reduce
		vcl.sum(count);
		vcl.execute();

		// Check
		BOOST_REQUIRE_EQUAL(count,(size_t)(k-2)*(k-2)*(k-2));

		// check with a 1x1x1 square
		{

		grid_key_dx<3> one(k/2,k/2,k/2);
		grid_key_dx<3> one_end(k/2,k/2,k/2);

		count = 0;

		// get the sub-domain iterator
		auto dom = g_dist.getSubDomainIterator(one,one_end);

		while (dom.isNext())
		{
			auto key = dom.get();
			auto key_g = g_dist.getGKey(key);

			// key_g
			BOOST_REQUIRE_EQUAL(key_g.get(0),k/2);
			BOOST_REQUIRE_EQUAL(key_g.get(1),k/2);
			BOOST_REQUIRE_EQUAL(key_g.get(2),k/2);

			auto key_s_it = dom.getGKey(key);

			BOOST_REQUIRE_EQUAL(key_g.get(0),key_s_it.get(0));
			BOOST_REQUIRE_EQUAL(key_g.get(1),key_s_it.get(1));
			BOOST_REQUIRE_EQUAL(key_g.get(2),key_s_it.get(2));

			// Count the point
			count++;

			++dom;
		}

		// reduce
		vcl.sum(count);
		vcl.execute();

		BOOST_REQUIRE_EQUAL(count,1ul);
		}
	}
}

void Test3D(const Box<3,float> & domain, long int k)
{
	long int big_step = k / 30;
	big_step = (big_step == 0)?1:big_step;
	long int small_step = 21;

	print_test( "Testing 3D grid k<=",k);

	// 3D test
	for ( ; k >= 2 ; k-= (k > 2*big_step)?big_step:small_step )
	{
		BOOST_TEST_CHECKPOINT( "Testing 3D grid k=" << k );

		// grid size
		size_t sz[3];
		sz[0] = k;
		sz[1] = k;
		sz[2] = k;

		// factor
		float factor = pow(global_v_cluster->getProcessingUnits()/2.0f,1.0f/3.0f);

		// Ghost
		Ghost<3,float> g(0.01 / factor);

		// Distributed grid with id decomposition
		grid_dist_id<3, float, scalar<float>, CartDecomposition<3,float>> g_dist(sz,domain,g);

		// check the consistency of the decomposition
		bool val = g_dist.getDecomposition().check_consistency();
		BOOST_REQUIRE_EQUAL(val,true);

		// Grid sm
		grid_sm<3,void> info(sz);

		// get the domain iterator
		size_t count = 0;

		auto dom = g_dist.getDomainIterator();

		while (dom.isNext())
		{
			auto key = dom.get();
			auto key_g = g_dist.getGKey(key);

			g_dist.template get<0>(key) = info.LinId(key_g);

			// Count the point
			count++;

			++dom;
		}

		// Get the virtual cluster machine
		Vcluster & vcl = g_dist.getVC();

		// reduce
		vcl.sum(count);
		vcl.execute();

		// Check
		BOOST_REQUIRE_EQUAL(count,(size_t)k*k*k);

		bool match = true;

		auto dom2 = g_dist.getDomainIterator();

		// check that the grid store the correct information
		while (dom2.isNext())
		{
			auto key = dom2.get();
			auto key_g = g_dist.getGKey(key);

			match &= (g_dist.template get<0>(key) == info.LinId(key_g))?true:false;

			++dom2;
		}

		BOOST_REQUIRE_EQUAL(match,true);

		//! [Synchronize the ghost and check the information]

		g_dist.template ghost_get<0>();

		// check that the communication is correctly completed

		auto domg = g_dist.getDomainGhostIterator();

		// check that the grid with the ghost past store the correct information
		while (domg.isNext())
		{
			auto key = domg.get();
			auto key_g = g_dist.getGKey(key);

			// In this case the boundary condition are non periodic
			if (g_dist.isInside(key_g))
			{
				match &= (g_dist.template get<0>(key) == info.LinId(key_g))?true:false;
			}

			++domg;
		}

		BOOST_REQUIRE_EQUAL(match,true);

		//! [Synchronize the ghost and check the information]
	}
}


void Test3D_gg(const Box<3,float> & domain, long int k, long int gk)
{
	long int big_step = k / 30;
	big_step = (big_step == 0)?1:big_step;

	// this test is only performed when the number of processor is <= 32
	if (global_v_cluster->getProcessingUnits() > 32)
		return;

	print_test( "Testing 3D grid k<=",k);

	// 3D test
	for ( ; k > 64 ; k /= 2 )
	{
		BOOST_TEST_CHECKPOINT( "Testing 3D grid ghost integer k=" << k );

		// grid size
		size_t sz[3];
		sz[0] = k;
		sz[1] = k;
		sz[2] = k;

		// Ghost
		Ghost<3,long int> g(gk);

		// Distributed grid with id decomposition
		grid_dist_id<3, float, scalar<float>, CartDecomposition<3,float>> g_dist(sz,domain,g);

		// check the consistency of the decomposition
		bool val = g_dist.getDecomposition().check_consistency();
		BOOST_REQUIRE_EQUAL(val,true);

		auto lg = g_dist.getLocalGridsInfo();

		// for each local grid check that the border is 1 point
		// (Warning this property can only be ensured with k is a multiple of 2)
		// in the other case it will be mostly like that but cannot be ensured

		for (size_t i = 0 ; i < lg.size() ; i++)
		{
			for (size_t j = 0 ; j < 3 ; j++)
			{
				BOOST_REQUIRE(lg.get(i).Dbox.getLow(j) >= gk);
				BOOST_REQUIRE((lg.get(i).GDbox.getHigh(j) - lg.get(i).Dbox.getHigh(j)) >= gk);
			}
		}
	}
}

void Test2D_complex(const Box<2,float> & domain, long int k)
{
	typedef Point_test<float> p;

	long int big_step = k / 30;
	big_step = (big_step == 0)?1:big_step;
	long int small_step = 21;

	print_test( "Testing 2D complex grid k<=",k);

	// 2D test
	for ( ; k >= 2 ; k-= (k > 2*big_step)?big_step:small_step )
	{
		BOOST_TEST_CHECKPOINT( "Testing 2D complex grid k=" << k );

		//! [Create and access a distributed grid complex]

		// grid size
		size_t sz[2];
		sz[0] = k;
		sz[1] = k;

		float factor = pow(global_v_cluster->getProcessingUnits()/2.0f,1.0f/2.0f);

		// Ghost
		Ghost<2,float> g(0.01 / factor);

		// Distributed grid with id decomposition
		grid_dist_id<2, float, Point_test<float>, CartDecomposition<2,float>> g_dist(sz,domain,g);

		// check the consistency of the decomposition
		bool val = g_dist.getDecomposition().check_consistency();
		BOOST_REQUIRE_EQUAL(val,true);

		// Grid sm
		grid_sm<2,void> info(sz);

		// get the domain iterator
		size_t count = 0;

		auto dom = g_dist.getDomainIterator();

		while (dom.isNext())
		{
			auto key = dom.get();
			auto key_g = g_dist.getGKey(key);

			size_t k = info.LinId(key_g);

			g_dist.template get<p::x>(key) = 1 + k;
			g_dist.template get<p::y>(key) = 567 + k;
			g_dist.template get<p::z>(key) = 341 + k;
			g_dist.template get<p::s>(key) = 5670 + k;
			g_dist.template get<p::v>(key)[0] = 921 + k;
			g_dist.template get<p::v>(key)[1] = 5675 + k;
			g_dist.template get<p::v>(key)[2] = 117 + k;
			g_dist.template get<p::t>(key)[0][0] = 1921 + k;
			g_dist.template get<p::t>(key)[0][1] = 25675 + k;
			g_dist.template get<p::t>(key)[0][2] = 3117 + k;
			g_dist.template get<p::t>(key)[1][0] = 4921 + k;
			g_dist.template get<p::t>(key)[1][1] = 55675 + k;
			g_dist.template get<p::t>(key)[1][2] = 6117 + k;
			g_dist.template get<p::t>(key)[2][0] = 7921 + k;
			g_dist.template get<p::t>(key)[2][1] = 85675 + k;
			g_dist.template get<p::t>(key)[2][2] = 9117 + k;

			// Count the point
			count++;

			++dom;
		}

		//! [Create and access a distributed grid complex]

		// Get the virtual cluster machine
		Vcluster & vcl = g_dist.getVC();

		// reduce
		vcl.sum(count);
		vcl.execute();

		// Check
		BOOST_REQUIRE_EQUAL(count,(size_t)k*k);

		auto dom2 = g_dist.getDomainIterator();

		bool match = true;

		// check that the grid store the correct information
		while (dom2.isNext())
		{
			auto key = dom2.get();
			auto key_g = g_dist.getGKey(key);

			size_t k = info.LinId(key_g);

			match &= (g_dist.template get<p::x>(key) == 1 + k)?true:false;
			match &= (g_dist.template get<p::y>(key) == 567 + k)?true:false;
			match &= (g_dist.template get<p::z>(key) == 341 + k)?true:false;
			match &= (g_dist.template get<p::s>(key) == 5670 + k)?true:false;
			match &= (g_dist.template get<p::v>(key)[0] == 921 + k)?true:false;
			match &= (g_dist.template get<p::v>(key)[1] == 5675 + k)?true:false;
			match &= (g_dist.template get<p::v>(key)[2] == 117 + k)?true:false;
			match &= (g_dist.template get<p::t>(key)[0][0] == 1921 + k)?true:false;
			match &= (g_dist.template get<p::t>(key)[0][1] == 25675 + k)?true:false;
			match &= (g_dist.template get<p::t>(key)[0][2] == 3117 + k)?true:false;
			match &= (g_dist.template get<p::t>(key)[1][0] == 4921 + k)?true:false;
			match &= (g_dist.template get<p::t>(key)[1][1] == 55675 + k)?true:false;
			match &= (g_dist.template get<p::t>(key)[1][2] == 6117 + k)?true:false;
			match &= (g_dist.template get<p::t>(key)[2][0] == 7921 + k)?true:false;
			match &= (g_dist.template get<p::t>(key)[2][1] == 85675 + k)?true:false;
			match &= (g_dist.template get<p::t>(key)[2][2] == 9117 + k)?true:false;

			++dom2;
		}

		BOOST_REQUIRE_EQUAL(match,true);

		//! [Synchronized distributed grid complex]

		g_dist.template ghost_get<p::x,p::y,p::z,p::s,p::v,p::t>();

		// check that the communication is correctly completed

		auto domg = g_dist.getDomainGhostIterator();

		// check that the grid with the ghost past store the correct information
		while (domg.isNext())
		{
			auto key = domg.get();
			auto key_g = g_dist.getGKey(key);

			// In this case the boundary condition are non periodic
			if (g_dist.isInside(key_g))
			{
				size_t k = info.LinId(key_g);

				match &= (g_dist.template get<p::x>(key) == 1 + k)?true:false;
				match &= (g_dist.template get<p::y>(key) == 567 + k)?true:false;
				match &= (g_dist.template get<p::z>(key) == 341 + k)?true:false;
				match &= (g_dist.template get<p::s>(key) == 5670 + k)?true:false;

				match &= (g_dist.template get<p::v>(key)[0] == 921 + k)?true:false;
				match &= (g_dist.template get<p::v>(key)[1] == 5675 + k)?true:false;
				match &= (g_dist.template get<p::v>(key)[2] == 117 + k)?true:false;

				match &= (g_dist.template get<p::t>(key)[0][0] == 1921 + k)?true:false;
				match &= (g_dist.template get<p::t>(key)[0][1] == 25675 + k)?true:false;
				match &= (g_dist.template get<p::t>(key)[0][2] == 3117 + k)?true:false;
				match &= (g_dist.template get<p::t>(key)[1][0] == 4921 + k)?true:false;
				match &= (g_dist.template get<p::t>(key)[1][1] == 55675 + k)?true:false;
				match &= (g_dist.template get<p::t>(key)[1][2] == 6117 + k)?true:false;
				match &= (g_dist.template get<p::t>(key)[2][0] == 7921 + k)?true:false;
				match &= (g_dist.template get<p::t>(key)[2][1] == 85675 + k)?true:false;
				match &= (g_dist.template get<p::t>(key)[2][2] == 9117 + k)?true:false;
			}

			++domg;
		}

		//! [Synchronized distributed grid complex]
	}
}

void Test3D_complex(const Box<3,float> & domain, long int k)
{
	typedef Point_test<float> p;

	long int big_step = k / 30;
	big_step = (big_step == 0)?1:big_step;
	long int small_step = 21;

	print_test( "Testing 3D grid complex k<=",k);

	// 2D test
	for ( ; k >= 2 ; k-= (k > 2*big_step)?big_step:small_step )
	{
		BOOST_TEST_CHECKPOINT( "Testing 3D complex grid k=" << k );

		// grid size
		size_t sz[3];
		sz[0] = k;
		sz[1] = k;
		sz[2] = k;

		// factor
		float factor = pow(global_v_cluster->getProcessingUnits()/2.0f,1.0f/3.0f);

		// Ghost
		Ghost<3,float> g(0.01 / factor);

		// Distributed grid with id decomposition
		grid_dist_id<3, float, Point_test<float>, CartDecomposition<3,float>> g_dist(sz,domain,g);

		// check the consistency of the decomposition
		bool val = g_dist.getDecomposition().check_consistency();
		BOOST_REQUIRE_EQUAL(val,true);

		// Grid sm
		grid_sm<3,void> info(sz);

		// get the domain iterator
		size_t count = 0;

		auto dom = g_dist.getDomainIterator();

		while (dom.isNext())
		{
			auto key = dom.get();
			auto key_g = g_dist.getGKey(key);

			size_t k = info.LinId(key_g);

			g_dist.template get<p::x>(key) = 1 + k;
			g_dist.template get<p::y>(key) = 567 + k;
			g_dist.template get<p::z>(key) = 341 + k;
			g_dist.template get<p::s>(key) = 5670 + k;
			g_dist.template get<p::v>(key)[0] = 921 + k;
			g_dist.template get<p::v>(key)[1] = 5675 + k;
			g_dist.template get<p::v>(key)[2] = 117 + k;
			g_dist.template get<p::t>(key)[0][0] = 1921 + k;
			g_dist.template get<p::t>(key)[0][1] = 25675 + k;
			g_dist.template get<p::t>(key)[0][2] = 3117 + k;
			g_dist.template get<p::t>(key)[1][0] = 4921 + k;
			g_dist.template get<p::t>(key)[1][1] = 55675 + k;
			g_dist.template get<p::t>(key)[1][2] = 6117 + k;
			g_dist.template get<p::t>(key)[2][0] = 7921 + k;
			g_dist.template get<p::t>(key)[2][1] = 85675 + k;
			g_dist.template get<p::t>(key)[2][2] = 9117 + k;

			// Count the point
			count++;

			++dom;
		}

		// Get the virtual cluster machine
		Vcluster & vcl = g_dist.getVC();

		// reduce
		vcl.sum(count);
		vcl.execute();

		// Check
		BOOST_REQUIRE_EQUAL(count,(size_t)k*k*k);

		bool match = true;

		auto dom2 = g_dist.getDomainIterator();

		// check that the grid store the correct information
		while (dom2.isNext())
		{
			auto key = dom2.get();
			auto key_g = g_dist.getGKey(key);

			size_t k = info.LinId(key_g);

			match &= (g_dist.template get<p::x>(key) == 1 + k)?true:false;
			match &= (g_dist.template get<p::y>(key) == 567 + k)?true:false;
			match &= (g_dist.template get<p::z>(key) == 341 + k)?true:false;
			match &= (g_dist.template get<p::s>(key) == 5670 + k)?true:false;
			match &= (g_dist.template get<p::v>(key)[0] == 921 + k)?true:false;
			match &= (g_dist.template get<p::v>(key)[1] == 5675 + k)?true:false;
			match &= (g_dist.template get<p::v>(key)[2] == 117 + k)?true:false;
			match &= (g_dist.template get<p::t>(key)[0][0] == 1921 + k)?true:false;
			match &= (g_dist.template get<p::t>(key)[0][1] == 25675 + k)?true:false;
			match &= (g_dist.template get<p::t>(key)[0][2] == 3117 + k)?true:false;
			match &= (g_dist.template get<p::t>(key)[1][0] == 4921 + k)?true:false;
			match &= (g_dist.template get<p::t>(key)[1][1] == 55675 + k)?true:false;
			match &= (g_dist.template get<p::t>(key)[1][2] == 6117 + k)?true:false;
			match &= (g_dist.template get<p::t>(key)[2][0] == 7921 + k)?true:false;
			match &= (g_dist.template get<p::t>(key)[2][1] == 85675 + k)?true:false;
			match &= (g_dist.template get<p::t>(key)[2][2] == 9117 + k)?true:false;

			++dom2;
		}

		BOOST_REQUIRE_EQUAL(match,true);

		g_dist.template ghost_get<p::x,p::y,p::z,p::s,p::v,p::t>();

		// check that the communication is correctly completed

		auto domg = g_dist.getDomainGhostIterator();

		// check that the grid with the ghost past store the correct information
		while (domg.isNext())
		{
			auto key = domg.get();
			auto key_g = g_dist.getGKey(key);

			size_t k = info.LinId(key_g);

			// In this case the boundary condition are non periodic
			if (g_dist.isInside(key_g))
			{
				match &= (g_dist.template get<p::x>(key) == 1 + k)?true:false;
				match &= (g_dist.template get<p::y>(key) == 567 + k)?true:false;
				match &= (g_dist.template get<p::z>(key) == 341 + k)?true:false;
				match &= (g_dist.template get<p::s>(key) == 5670 + k)?true:false;

				match &= (g_dist.template get<p::v>(key)[0] == 921 + k)?true:false;
				match &= (g_dist.template get<p::v>(key)[1] == 5675 + k)?true:false;
				match &= (g_dist.template get<p::v>(key)[2] == 117 + k)?true:false;

				match &= (g_dist.template get<p::t>(key)[0][0] == 1921 + k)?true:false;
				match &= (g_dist.template get<p::t>(key)[0][1] == 25675 + k)?true:false;
				match &= (g_dist.template get<p::t>(key)[0][2] == 3117 + k)?true:false;
				match &= (g_dist.template get<p::t>(key)[1][0] == 4921 + k)?true:false;
				match &= (g_dist.template get<p::t>(key)[1][1] == 55675 + k)?true:false;
				match &= (g_dist.template get<p::t>(key)[1][2] == 6117 + k)?true:false;
				match &= (g_dist.template get<p::t>(key)[2][0] == 7921 + k)?true:false;
				match &= (g_dist.template get<p::t>(key)[2][1] == 85675 + k)?true:false;
				match &= (g_dist.template get<p::t>(key)[2][2] == 9117 + k)?true:false;
			}

			++domg;
		}

		BOOST_REQUIRE_EQUAL(match,true);
	}
}

// Test duplicated topology

void Test3D_dup(const Box<3,float> & domain, long int k)
{
	long int big_step = k / 30;
	big_step = (big_step == 0)?1:big_step;
	long int small_step = 21;
	long int k_old = k;

	Vcluster & v_cl = *global_v_cluster;

	if ( v_cl.getProcessingUnits() > 32 )
		return;

	print_test( "Testing 3D duplicate topology complex k<=",k);

	// 3D test
	for ( ; k >= 2 ; k-= (k > 2*big_step)?big_step:small_step )
	{
		BOOST_TEST_CHECKPOINT( "Testing 3D copy decomposition grid k=" << k );

		// grid size
		size_t sz[3];
		sz[0] = k;
		sz[1] = k;
		sz[2] = k;

		// factor
		float factor = pow(global_v_cluster->getProcessingUnits()/2.0f,1.0f/3.0f);

		// Ghost
		Ghost<3,float> g(0.01 / factor);

		//! [Construct two grid with the same decomposition]

		// Distributed grid with id decomposition
		grid_dist_id<3, float, Point_test<float>, CartDecomposition<3,float>> g_dist1(sz,domain,g);

		// another grid with the same decomposition
		grid_dist_id<3, float, Point_test<float>, CartDecomposition<3,float>> g_dist2(g_dist1.getDecomposition(),sz,domain,g);

		//! [Construct two grid with the same decomposition]

		BOOST_REQUIRE_EQUAL(g_dist2.getDecomposition().ref(),2);

		auto dom_g1 = g_dist1.getDomainIterator();
		auto dom_g2 = g_dist2.getDomainIterator();

		bool check = true;

		while (dom_g1.isNext())
		{
			auto key1 = dom_g1.get();
			auto key2 = dom_g2.get();

			check &= (key1 == key2)?true:false;

			++dom_g1;
			++dom_g2;
		}

		BOOST_REQUIRE_EQUAL(check,true);
	}

	k = k_old;

	// 3D test
	for ( ; k >= 2 ; k-= (k > 2*big_step)?big_step:small_step )
	{
		BOOST_TEST_CHECKPOINT( "Testing 3D copy decomposition grid k=" << k );

		// grid size
		size_t sz[3];
		sz[0] = k;
		sz[1] = k;
		sz[2] = k;

		// factor
		float factor = pow(global_v_cluster->getProcessingUnits()/2.0f,1.0f/3.0f);

		// Ghost
		Ghost<3,float> g(0.01 / factor);

		//! [Construct two grid with the same decomposition]

		// Distributed grid with id decomposition
		grid_dist_id<3, float, Point_test<float>, CartDecomposition<3,float>> * g_dist1 = new grid_dist_id<3, float, Point_test<float>, CartDecomposition<3,float>>(sz,domain,g);

		// another grid with the same decomposition
		grid_dist_id<3, float, Point_test<float>, CartDecomposition<3,float>> * g_dist2 = new grid_dist_id<3, float, Point_test<float>, CartDecomposition<3,float>>(g_dist1->getDecomposition(),sz,domain,g);

		//! [Construct two grid with the same decomposition]

		BOOST_REQUIRE_EQUAL(g_dist2->getDecomposition().ref(),2);

		delete g_dist1;

		BOOST_REQUIRE_EQUAL(g_dist2->getDecomposition().ref(),1);
		bool ret = g_dist2->getDecomposition().check_consistency();
		BOOST_REQUIRE_EQUAL(ret,true);

		delete g_dist2;
	}
}

BOOST_AUTO_TEST_CASE( grid_dist_id_iterator_test_use)
{
	// Domain
	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	// Initialize the global VCluster
	init_global_v_cluster(&boost::unit_test::framework::master_test_suite().argc,&boost::unit_test::framework::master_test_suite().argv);

	long int k = 1024*1024*global_v_cluster->getProcessingUnits();
	k = std::pow(k, 1/2.);

	Test2D(domain,k);
	Test2D_complex(domain,k);
	// Domain
	Box<3,float> domain3({0.0,0.0,0.0},{1.0,1.0,1.0});

	k = 128*128*128*global_v_cluster->getProcessingUnits();
	k = std::pow(k, 1/3.);
	Test3D(domain3,k);
	Test3D_complex(domain3,k);
}

BOOST_AUTO_TEST_CASE( grid_dist_id_dup)
{
	// Domain
	Box<3,float> domain3({0.0,0.0,0.0},{1.0,1.0,1.0});

	long int k = 128*128*128*global_v_cluster->getProcessingUnits();
	k = std::pow(k, 1/3.);
	Test3D_dup(domain3,k);
}

BOOST_AUTO_TEST_CASE( grid_dist_id_sub)
{
	// Domain
	Box<3,float> domain3({0.0,0.0,0.0},{1.0,1.0,1.0});

	long int k = 128*128*128*global_v_cluster->getProcessingUnits();
	k = std::pow(k, 1/3.);
	Test3D_sub(domain3,k);
}

BOOST_AUTO_TEST_CASE( grid_dist_id_sub_iterator_test_use)
{
	// Domain
	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	// Initialize the global VCluster
	init_global_v_cluster(&boost::unit_test::framework::master_test_suite().argc,&boost::unit_test::framework::master_test_suite().argv);

	long int k = 1024*1024*global_v_cluster->getProcessingUnits();
	k = std::pow(k, 1/2.);

	Test2D_sub(domain,k);
}

BOOST_AUTO_TEST_CASE( grid_dist_id_with_grid_unit_ghost )
{
	// Domain
	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	// Initialize the global VCluster
	init_global_v_cluster(&boost::unit_test::framework::master_test_suite().argc,&boost::unit_test::framework::master_test_suite().argv);

	long int k = 1024*1024*global_v_cluster->getProcessingUnits();
	k = std::pow(k, 1/2.);

//	Test2D_gg(domain,k);
	// Domain
	Box<3,float> domain3({0.0,0.0,0.0},{1.0,1.0,1.0});

	k = 128*128*128*global_v_cluster->getProcessingUnits();
	k = std::pow(k, 1/3.);
	Test3D_gg(domain3,k,1);
}

BOOST_AUTO_TEST_SUITE_END()

#endif
