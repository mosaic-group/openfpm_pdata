#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Point_test.hpp"
#include "Grid/grid_dist_id.hpp"
#include "data_type/scalar.hpp"
#include "data_type/aggregate.hpp"
#include "grid_dist_id_util_tests.hpp"

BOOST_AUTO_TEST_SUITE( grid_dist_id_test )


BOOST_AUTO_TEST_CASE( grid_dist_id_domain_grid_unit_converter3D_test)
{
	size_t bc[3] = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};

	// Domain
	Box<3,float> domain({-0.3,-0.3,-0.3},{1.0,1.0,1.0});

	Vcluster & v_cl = create_vcluster();

	// Skip this test on big scale
	if (v_cl.getProcessingUnits() >= 32)
		return;

	// Test several grid dimensions

	long int k = 293;
	long int big_step = k / 30;
	/* coverity[dead_error_line] */
	big_step = (big_step == 0)?1:big_step;
	long int small_step = 21;

	print_test( "Testing 3D grid converter k<=",k);

	// 3D test
	for ( ; k >= 2 ; k-= (k > 2*big_step)?big_step:small_step )
	{
		BOOST_TEST_CHECKPOINT( "Testing 3D grid converter k=" << k );

		// grid size
		size_t sz[3];
		sz[0] = k;
		sz[1] = k;
		sz[2] = k;

		// Ghost
		Ghost<3,float> g(0.01);

		// Distributed grid with id decomposition
		grid_dist_id<3, float, scalar<float>, CartDecomposition<3,float>> g_dist(sz,domain,g);

		// get the decomposition
		auto & dec = g_dist.getDecomposition();

		// check the consistency of the decomposition
		bool val = dec.check_consistency();
		BOOST_REQUIRE_EQUAL(val,true);

		// for each local volume
		// Get the number of local grid needed
		size_t n_grid = dec.getNSubDomain();

		size_t vol = 0;

		// vector of boxes
		openfpm::vector<Box<3,size_t>> vb;

		// Allocate the grids
		for (size_t i = 0 ; i < n_grid ; i++)
		{
			// Get the local hyper-cube
			SpaceBox<3,float> sub = dec.getSubDomain(i);
//			sub -= domain.getP1();

			Box<3,size_t> g_box = g_dist.getCellDecomposer().convertDomainSpaceIntoGridUnits(sub,bc);

			vb.add(g_box);

			vol += g_box.getVolumeKey();
		}

		// Create a writer and write
		VTKWriter<openfpm::vector<Box<3,size_t>>,VECTOR_BOX> vtk_box2;
		vtk_box2.add(vb);
		vtk_box2.write(std::to_string(v_cl.getProcessUnitID()) + "vtk_box_3D.vtk");

		v_cl.sum(vol);
		v_cl.execute();

		BOOST_REQUIRE_EQUAL(vol,sz[0]*sz[1]*sz[2]);
	}
}


BOOST_AUTO_TEST_CASE( grid_dist_id_domain_grid_unit_converter_test)
{
	size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};

	// Domain
	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	Vcluster & v_cl = create_vcluster();

	// Skip this test on big scale
	if (v_cl.getProcessingUnits() >= 32)
		return;

	for (size_t k = 1024 ; k >= 2 ; k--)
	{
		BOOST_TEST_CHECKPOINT( "Testing grid converter 3D k=" << k );

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
		size_t n_grid = dec.getNSubDomain();

		size_t vol = 0;

		// Allocate the grids
		for (size_t i = 0 ; i < n_grid ; i++)
		{
			// Get the local hyper-cube
			SpaceBox<2,float> sub = dec.getSubDomain(i);

			Box<2,size_t> g_box = g_dist.getCellDecomposer().convertDomainSpaceIntoGridUnits(sub,bc);

			vol += g_box.getVolumeKey();
		}

		v_cl.sum(vol);
		v_cl.execute();

		BOOST_REQUIRE_EQUAL(vol,sz[0]*sz[1]);
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

		float factor = pow(create_vcluster().getProcessingUnits()/2.0f,1.0f/2.0f);

		// Ghost
		Ghost<2,float> g(0.01 / factor);

		// Distributed grid with id decomposition
		grid_dist_id<2, float, scalar<float>> g_dist(sz,domain,g);

		Test2D_core(g_dist,sz,k);
	}
}


void Test1D(const Box<1,float> & domain, long int k)
{
	Vcluster & v_cl = create_vcluster();
	long int big_step = k / 30;
	big_step = (big_step == 0)?1:big_step;
	long int small_step = 21;

	if (v_cl.getProcessingUnits() > 48)
		return;

	print_test( "Testing 1D grid k<=",k);

	// 1D test
	for ( ; k >= 2 ; k-= (k > 2*big_step)?big_step:small_step )
	{
		BOOST_TEST_CHECKPOINT( "Testing 1D grid k=" << k );

		//! [Create and access a distributed grid]

		// grid size
		size_t sz[1];
		sz[0] = k;

		float factor = pow(create_vcluster().getProcessingUnits()/2.0f,1.0f);

		// Ghost
		Ghost<1,float> g(0.01 / factor);

		// Distributed grid with id decomposition
		grid_dist_id<1, float, scalar<float>> g_dist(sz,domain,g);

		// check the consistency of the decomposition
		bool val = g_dist.getDecomposition().check_consistency();
		BOOST_REQUIRE_EQUAL(val,true);

		// Grid sm
		grid_sm<1,void> info(sz);

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
		BOOST_REQUIRE_EQUAL(count,(size_t)k);

		auto dom2 = g_dist.getDomainIterator();

		grid_key_dx<1> start = dom2.getStart();
		grid_key_dx<1> stop = dom2.getStop();

		BOOST_REQUIRE_EQUAL((long int)stop.get(0),(long int)g_dist.size(0)-1);

		BOOST_REQUIRE_EQUAL(start.get(0),0);

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
				match &= (g_dist.template get<0>(key) == info.LinId(key_g))?true:false;
			}

			++domg;
		}

		BOOST_REQUIRE_EQUAL(match,true);
	}
}

void Test3D_sub(const Box<3,float> & domain, long int k)
{
	long int big_step = k / 30;
	big_step = (big_step == 0)?1:big_step;
	long int small_step = 21;

	// this test is only performed when the number of processor is <= 32
	if (create_vcluster().getProcessingUnits() > 32)
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
		float factor = pow(create_vcluster().getProcessingUnits()/2.0f,1.0f/3.0f);

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
		float factor = pow(create_vcluster().getProcessingUnits()/2.0f,1.0f/3.0f);

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

		// check that the grid with the ghost part store the correct information
		while (domg.isNext())
		{
			auto key = domg.get();
			auto key_g = g_dist.getGKey(key);

			// In this case the boundary condition are non periodic
			if (g_dist.isInside(key_g))
			{
				match &= (g_dist.template get<0>(key) == info.LinId(key_g))?true:false;
				if (match == false)
				{std::cout << "ERROR IN: " << key_g.to_string() << "   " << info.LinId(key_g) << " != " << g_dist.template get<0>(key) << std::endl; break;}
			}

			++domg;
		}

//		if (match == false)
//		{
			g_dist.write("Error_grid");

			g_dist.getDecomposition().write("Error_dec");
//		}

		BOOST_REQUIRE_EQUAL(match,true);

		//! [Synchronize the ghost and check the information]
	}
}


void Test3D_gg(const Box<3,float> & domain, long int k, long int gk)
{
	long int big_step = k / 30;
	big_step = (big_step == 0)?1:big_step;

	// this test is only performed when the number of processor is <= 32
	if (create_vcluster().getProcessingUnits() > 32)
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

/*! \brief Test when the domain is not from 0.0 to 1.0
 *
 *
 */

void Test3D_domain(const Box<3,float> & domain, long int k, const periodicity<3> & pr)
{
	long int big_step = k / 30;
	big_step = (big_step == 0)?1:big_step;
	long int small_step = 21;

	print_test( "Testing 3D grid shift domain k<=",k);

	// 3D test
	for ( ; k >= 2 ; k-= (k > 2*big_step)?big_step:small_step )
	{
		BOOST_TEST_CHECKPOINT( "Testing 3D grid shift domain k=" << k );

		// grid size
		size_t sz[3];
		sz[0] = k;
		sz[1] = k;
		sz[2] = k;

		// factor
		float factor = pow(create_vcluster().getProcessingUnits()/2.0f,1.0f/3.0f);

		// Ghost
		Ghost<3,float> g(0.01 / factor);

		// Distributed grid with id decomposition
		grid_dist_id<3, float, aggregate<long int,long int>, CartDecomposition<3,float>> g_dist(sz,domain,g,pr);

		auto & v_cl = create_vcluster();

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

			g_dist.template get<0>(key) = count;
			g_dist.template get<1>(key) = info.LinId(key_g);

			// Count the point
			count++;

			++dom;
		}

		size_t count2 = count;
		openfpm::vector<size_t> pnt;

		// Get the total size of the local grids on each processors
		// and the total size
		v_cl.sum(count2);
		v_cl.allGather(count,pnt);
		v_cl.execute();
		size_t s_pnt = 0;

		// calculate the starting point for this processor
		for (size_t i = 0 ; i < v_cl.getProcessUnitID() ; i++)
			s_pnt += pnt.get(i);

		// Check
		BOOST_REQUIRE_EQUAL(count2,(size_t)k*k*k);

		// sync the ghost
		g_dist.template ghost_get<0,1>();

		bool match = true;

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
				match &= (g_dist.template get<1>(key) == info.LinId(key_g))?true:false;
			}

			++domg;
		}

		BOOST_REQUIRE_EQUAL(match,true);
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

		float factor = pow(create_vcluster().getProcessingUnits()/2.0f,1.0f/2.0f);

		// Ghost
		Ghost<2,float> g(0.01 / factor);

		// Distributed grid with id decomposition
		grid_dist_id<2, float, Point_test<float>> g_dist(sz,domain,g);

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
		float factor = pow(create_vcluster().getProcessingUnits()/2.0f,1.0f/3.0f);

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

	Vcluster & v_cl = create_vcluster();

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
		float factor = pow(create_vcluster().getProcessingUnits()/2.0f,1.0f/3.0f);

		// Ghost
		Ghost<3,float> g(0.01 / factor);

		//! [Construct two grid with the same decomposition]

		// Distributed grid with id decomposition (It work also without the third template parameter)
		// Here is given to show that the 2 grid MUST have the same decomposition strategy
		grid_dist_id<3, float, Point_test<float>, CartDecomposition<3,float>> g_dist1(sz,domain,g);

		// another grid with the same decomposition
		grid_dist_id<3, float, Point_test<float>, CartDecomposition<3,float>> g_dist2(g_dist1.getDecomposition(),sz,g);

		//! [Construct two grid with the same decomposition]

		bool ret = g_dist2.getDecomposition().check_consistency();
		BOOST_REQUIRE_EQUAL(ret,true);
		ret = g_dist2.getDecomposition().is_equal(g_dist2.getDecomposition());
		BOOST_REQUIRE_EQUAL(ret,true);


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
		float factor = pow(create_vcluster().getProcessingUnits()/2.0f,1.0f/3.0f);

		// Ghost
		Ghost<3,float> g(0.01 / factor);

		// Distributed grid with id decomposition
		grid_dist_id<3, float, Point_test<float>, CartDecomposition<3,float>> * g_dist1 = new grid_dist_id<3, float, Point_test<float>, CartDecomposition<3,float>>(sz,domain,g);

		// another grid with the same decomposition
		grid_dist_id<3, float, Point_test<float>, CartDecomposition<3,float>> * g_dist2 = new grid_dist_id<3, float, Point_test<float>, CartDecomposition<3,float>>(g_dist1->getDecomposition(),sz,g);

		delete g_dist1;

		bool ret = g_dist2->getDecomposition().check_consistency();
		BOOST_REQUIRE_EQUAL(ret,true);

		delete g_dist2;
	}
}


// Test grid periodic

void Test3D_periodic(const Box<3,float> & domain, long int k)
{
	Vcluster & v_cl = create_vcluster();

	if ( v_cl.getProcessingUnits() > 32 )
		return;

	long int big_step = k / 30;
	big_step = (big_step == 0)?1:big_step;
	long int small_step = 21;

	print_test( "Testing grid periodic k<=",k);

	// 3D test
	for ( ; k >= 2 ; k-= (k > 2*big_step)?big_step:small_step )
	{
		BOOST_TEST_CHECKPOINT( "Testing grid periodick<=" << k );

		// grid size
		size_t sz[3];
		sz[0] = k;
		sz[1] = k;
		sz[2] = k;

		// factor
		float factor = pow(create_vcluster().getProcessingUnits()/2.0f,1.0f/3.0f);

		// Ghost
		Ghost<3,float> g(0.01 / factor);

		// periodicity
		periodicity<3> pr = {{PERIODIC,PERIODIC,PERIODIC}};

		// Distributed grid with id decomposition
		grid_dist_id<3, float, aggregate<long int>, CartDecomposition<3,float>> g_dist(sz,domain,g,pr);

		// check the consistency of the decomposition
		bool val = g_dist.getDecomposition().check_consistency();
		BOOST_REQUIRE_EQUAL(val,true);

		// Grid sm
		grid_sm<3,void> info(sz);

		size_t count = 0;

		// Set to zero the full grid

		auto dom1 = g_dist.getDomainGhostIterator();

		while (dom1.isNext())
		{
			auto key = dom1.get();

			g_dist.template get<0>(key) = -1;

			++dom1;
		}

		auto dom = g_dist.getDomainIterator();

		while (dom.isNext())
		{
			auto key = dom.get();
			auto key_g = g_dist.getGKey(key);

			g_dist.template get<0>(key) = info.LinId(key_g);

			// Count the points
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

		size_t tot = g_dist.getLocalDomainSize();
		// reduce
		vcl.sum(tot);
		vcl.execute();

		BOOST_REQUIRE_EQUAL(count,tot);

		// sync the ghosts
		g_dist.ghost_get<0>();

		bool match = true;

		// Domain + Ghost iterator
		auto dom_gi = g_dist.getDomainGhostIterator();

		size_t out_cnt = 0;

		while (dom_gi.isNext())
		{
			bool out_p = false;

			auto key = dom_gi.get();
			auto key_g = g_dist.getGKey(key);

			// Return the external boxes
			auto & gb = dom_gi.getGBoxes();

			// transform the key to be periodic
			for (size_t i = 0 ; i < 3 ; i++)
			{
				if (key_g.get(i) < 0)
				{key_g.set_d(i,key_g.get(i) + k);out_p = true;}
				else if (key_g.get(i) >= k)
				{key_g.set_d(i,key_g.get(i) - k);out_p = true;}
			}

			if (g_dist.template get<0>(key) != -1 && out_p == true)
				out_cnt++;

			// The last points can be invalid because of rounding off problems
			bool can_invalid = false;
			if (key.getKey().get(0) == 0 || key.getKey().get(1) == 0 || key.getKey().get(2) == 0)
				can_invalid = true;
			else if (key.getKey().get(0) == gb.get(key.getSub()).GDbox.getHigh(0) ||
					 key.getKey().get(1) == gb.get(key.getSub()).GDbox.getHigh(1) ||
					 key.getKey().get(2) == gb.get(key.getSub()).GDbox.getHigh(2))
				can_invalid = true;

			if (can_invalid == true)
			{
				if ( g_dist.template get<0>(key) != -1 && info.LinId(key_g) != g_dist.template get<0>(key) )
					match &= false;
			}
			else
			{
				if (info.LinId(key_g) != g_dist.template get<0>(key) )
					match &= false;
			}

			++dom_gi;
		}

		BOOST_REQUIRE_EQUAL(match, true);
		if (k > 83)
		{
			vcl.sum(out_cnt);
			vcl.execute();
			BOOST_REQUIRE(out_cnt != 0ul);
		}
	}
}

// Test grid periodic

void Test3D_periodic_put(const Box<3,float> & domain, long int k)
{
	Vcluster & v_cl = create_vcluster();

	if ( v_cl.getProcessingUnits() > 32 )
	{return;}

	long int big_step = k / 30;
	big_step = (big_step == 0)?1:big_step;
	long int small_step = 21;

	print_test( "Testing grid periodic put k<=",k);

	// 3D test
	for ( ; k >= 2 ; k-= (k > 2*big_step)?big_step:small_step )
	{
		BOOST_TEST_CHECKPOINT( "Testing grid periodick<=" << k );

		// grid size
		size_t sz[3];
		sz[0] = k;
		sz[1] = k;
		sz[2] = k;

		// Ghost
		Ghost<3,long int> g(1);

		// periodicity
		periodicity<3> pr = {{PERIODIC,PERIODIC,PERIODIC}};

		// Distributed grid with id decomposition
		grid_dist_id<3, float, aggregate<long int>, CartDecomposition<3,float>> g_dist(sz,domain,g,pr);

		// check the consistency of the decomposition
		bool val = g_dist.getDecomposition().check_consistency();
		BOOST_REQUIRE_EQUAL(val,true);

		// Grid sm
		grid_sm<3,void> info(sz);

		size_t count = 0;

		{
		auto dom = g_dist.getDomainIterator();

		while (dom.isNext())
		{
			auto key = dom.get();

			g_dist.template get<0>(key) = -6.0;

			// Count the points
			count++;

			++dom;
		}
		}

		// Set to zero the full grid

		{
		auto dom = g_dist.getDomainIterator();

		while (dom.isNext())
		{
			auto key = dom.get();

			g_dist.template get<0>(key.move(0,1)) += 1.0;
			g_dist.template get<0>(key.move(0,-1)) += 1.0;
			g_dist.template get<0>(key.move(1,1)) += 1.0;
			g_dist.template get<0>(key.move(1,-1)) += 1.0;
			g_dist.template get<0>(key.move(2,1)) += 1.0;
			g_dist.template get<0>(key.move(2,-1)) += 1.0;

			++dom;
		}
		}

		bool correct = true;

		// Domain + Ghost iterator
		auto dom_gi = g_dist.getDomainIterator();

		while (dom_gi.isNext())
		{
			auto key = dom_gi.get();

			correct &= (g_dist.template get<0>(key) == 0);

			++dom_gi;
		}

		g_dist.ghost_put<add_,0>();

		if (count != 0)
			BOOST_REQUIRE_EQUAL(correct, false);

		// sync the ghosts
		g_dist.ghost_get<0>();

		correct = true;

		// Domain + Ghost iterator
		auto dom_gi2 = g_dist.getDomainIterator();

		while (dom_gi2.isNext())
		{
			auto key = dom_gi2.get();

			correct &= (g_dist.template get<0>(key) == 0);

			++dom_gi2;
		}

		BOOST_REQUIRE_EQUAL(correct, true);
	}
}

void Test_grid_copy(const Box<3,float> & domain, long int k)
{
	typedef Point_test<float> p;

	Vcluster & v_cl = create_vcluster();

	if ( v_cl.getProcessingUnits() > 32 )
		return;

	long int big_step = k / 30;
	big_step = (big_step == 0)?1:big_step;
	long int small_step = 21;

	print_test( "Testing grid copy k<=",k);

	// 3D test
	for ( ; k >= 2 ; k-= (k > 2*big_step)?big_step:small_step )
	{
		BOOST_TEST_CHECKPOINT( "Testing grid periodick<=" << k );

		// grid size
		size_t sz[3];
		sz[0] = k;
		sz[1] = k;
		sz[2] = k;

		// factor
		float factor = pow(create_vcluster().getProcessingUnits()/2.0f,1.0f/3.0f);

		// Ghost
		Ghost<3,float> g(0.01 / factor);

		// periodicity
		periodicity<3> pr = {{PERIODIC,PERIODIC,PERIODIC}};

		// Distributed grid with id decomposition
		grid_dist_id<3,float,Point_test<float>> g_dist(sz,domain,g,pr);
		grid_dist_id<3,float,Point_test<float>> g_dist2(g_dist.getDecomposition(),sz,g);

		// Grid sm
		grid_sm<3,void> info(sz);

		// Set to zero the full grid
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

			++dom;
		}

		g_dist2.copy(g_dist);

		auto dom2 = g_dist2.getDomainIterator();

		bool match = true;

		// check that the grid store the correct information
		while (dom2.isNext())
		{
			auto key = dom2.get();
			auto key_g = g_dist.getGKey(key);

			size_t k = info.LinId(key_g);

			match &= (g_dist2.template get<p::x>(key) == 1 + k)?true:false;
			match &= (g_dist2.template get<p::y>(key) == 567 + k)?true:false;
			match &= (g_dist2.template get<p::z>(key) == 341 + k)?true:false;
			match &= (g_dist2.template get<p::s>(key) == 5670 + k)?true:false;
			match &= (g_dist2.template get<p::v>(key)[0] == 921 + k)?true:false;
			match &= (g_dist2.template get<p::v>(key)[1] == 5675 + k)?true:false;
			match &= (g_dist2.template get<p::v>(key)[2] == 117 + k)?true:false;
			match &= (g_dist2.template get<p::t>(key)[0][0] == 1921 + k)?true:false;
			match &= (g_dist2.template get<p::t>(key)[0][1] == 25675 + k)?true:false;
			match &= (g_dist2.template get<p::t>(key)[0][2] == 3117 + k)?true:false;
			match &= (g_dist2.template get<p::t>(key)[1][0] == 4921 + k)?true:false;
			match &= (g_dist2.template get<p::t>(key)[1][1] == 55675 + k)?true:false;
			match &= (g_dist2.template get<p::t>(key)[1][2] == 6117 + k)?true:false;
			match &= (g_dist2.template get<p::t>(key)[2][0] == 7921 + k)?true:false;
			match &= (g_dist2.template get<p::t>(key)[2][1] == 85675 + k)?true:false;
			match &= (g_dist2.template get<p::t>(key)[2][2] == 9117 + k)?true:false;

			++dom2;
		}

		BOOST_REQUIRE_EQUAL(match,true);
	}
}

void Test3D_copy(const Box<3,float> & domain, long int k)
{
	typedef Point_test<float> p;

	Vcluster & v_cl = create_vcluster();

	if ( v_cl.getProcessingUnits() > 32 )
		return;

	long int big_step = k / 30;
	big_step = (big_step == 0)?1:big_step;
	long int small_step = 21;

	print_test( "Testing grid copy k<=",k);

	// 3D test
	for ( ; k >= 2 ; k-= (k > 2*big_step)?big_step:small_step )
	{
		BOOST_TEST_CHECKPOINT( "Testing grid periodick<=" << k );

		// grid size
		size_t sz[3];
		sz[0] = k;
		sz[1] = k;
		sz[2] = k;

		// factor
		float factor = pow(create_vcluster().getProcessingUnits()/2.0f,1.0f/3.0f);

		// Ghost
		Ghost<3,float> g(0.01 / factor);

		// periodicity
		periodicity<3> pr = {{PERIODIC,PERIODIC,PERIODIC}};

		// Distributed grid with id decomposition
		grid_dist_id<3,float,Point_test<float>> g_dist(sz,domain,g,pr);

		// Grid sm
		grid_sm<3,void> info(sz);

		// Set to zero the full grid
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

			++dom;
		}

		grid_dist_id<3,float,Point_test<float>> g_dist2 = g_dist;
		g_dist2.template ghost_get<0>();

		auto dom2 = g_dist2.getDomainIterator();

		bool match = true;

		// check that the grid store the correct information
		while (dom2.isNext())
		{
			auto key = dom2.get();
			auto key_g = g_dist.getGKey(key);

			size_t k = info.LinId(key_g);

			match &= (g_dist2.template get<p::x>(key) == 1 + k)?true:false;
			match &= (g_dist2.template get<p::y>(key) == 567 + k)?true:false;
			match &= (g_dist2.template get<p::z>(key) == 341 + k)?true:false;
			match &= (g_dist2.template get<p::s>(key) == 5670 + k)?true:false;
			match &= (g_dist2.template get<p::v>(key)[0] == 921 + k)?true:false;
			match &= (g_dist2.template get<p::v>(key)[1] == 5675 + k)?true:false;
			match &= (g_dist2.template get<p::v>(key)[2] == 117 + k)?true:false;
			match &= (g_dist2.template get<p::t>(key)[0][0] == 1921 + k)?true:false;
			match &= (g_dist2.template get<p::t>(key)[0][1] == 25675 + k)?true:false;
			match &= (g_dist2.template get<p::t>(key)[0][2] == 3117 + k)?true:false;
			match &= (g_dist2.template get<p::t>(key)[1][0] == 4921 + k)?true:false;
			match &= (g_dist2.template get<p::t>(key)[1][1] == 55675 + k)?true:false;
			match &= (g_dist2.template get<p::t>(key)[1][2] == 6117 + k)?true:false;
			match &= (g_dist2.template get<p::t>(key)[2][0] == 7921 + k)?true:false;
			match &= (g_dist2.template get<p::t>(key)[2][1] == 85675 + k)?true:false;
			match &= (g_dist2.template get<p::t>(key)[2][2] == 9117 + k)?true:false;

			++dom2;
		}

		BOOST_REQUIRE_EQUAL(match,true);
	}
}

#include "grid_dist_id_unit_test_ext_dom.hpp"
#include "grid_dist_id_unit_test_unb_ghost.hpp"

BOOST_AUTO_TEST_CASE( grid_dist_id_iterator_test_use_2D)
{
	// Domain
	Box<2,float> domain({0.0,0.0},{1.0,1.0});

#ifdef TEST_COVERAGE_MODE
	long int k = 256*256*create_vcluster().getProcessingUnits();
#else
	long int k = 1024*1024*create_vcluster().getProcessingUnits();
#endif
	k = std::pow(k, 1/2.);

	Test2D(domain,k);
	Test2D_complex(domain,k);
}

BOOST_AUTO_TEST_CASE( grid_dist_id_iterator_test_use_3D)
{
	// Domain
	Box<3,float> domain3({0.0,0.0,0.0},{1.0,1.0,1.0});

	size_t k = 128*128*128*create_vcluster().getProcessingUnits();
	k = std::pow(k, 1/3.);
	Test3D(domain3,k);
	Test3D_complex(domain3,k);
}

BOOST_AUTO_TEST_CASE( grid_dist_id_dup)
{
	// Domain
	Box<3,float> domain3({0.0,0.0,0.0},{1.0,1.0,1.0});

	long int k = 128*128*128*create_vcluster().getProcessingUnits();
	k = std::pow(k, 1/3.);
	Test3D_dup(domain3,k);
}

BOOST_AUTO_TEST_CASE( grid_dist_id_sub)
{
	// Domain
	Box<3,float> domain3({0.0,0.0,0.0},{1.0,1.0,1.0});

	long int k = 128*128*128*create_vcluster().getProcessingUnits();
	k = std::pow(k, 1/3.);
	Test3D_sub(domain3,k);
}


BOOST_AUTO_TEST_CASE( grid_dist_id_with_grid_unit_ghost )
{
	// Domain
	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	long int k = 1024*1024*create_vcluster().getProcessingUnits();
	k = std::pow(k, 1/2.);

	// Domain
	Box<3,float> domain3({0.0,0.0,0.0},{1.0,1.0,1.0});

	k = 128*128*128*create_vcluster().getProcessingUnits();
	k = std::pow(k, 1/3.);
	Test3D_gg(domain3,k,1);
}


BOOST_AUTO_TEST_CASE( grid_dist_id_domain_test_use)
{
	// Domain
	Box<3,float> domain3({-0.3,-0.3,-0.3},{1.1,1.1,1.1});

	periodicity<3> np({{NON_PERIODIC,NON_PERIODIC,NON_PERIODIC}});
	periodicity<3> p({{PERIODIC,PERIODIC,PERIODIC}});

	long int k = 128*128*128*create_vcluster().getProcessingUnits();
	k = std::pow(k, 1/3.);
	Test3D_domain(domain3,k,np);

	auto & v_cl = create_vcluster();
	if (v_cl.getProcessingUnits() > 32)
		return;

	// We use a 128x128x128 and we move tha domain

	for (size_t i = 0 ; i < 10 ; i++)
	{
		Box<3,float> exp({0.0,0.0,0.0},{1.3,1.3,1.3});
		domain3.enlarge(exp);
		Test3D_domain(domain3,128,p);
	}
}

BOOST_AUTO_TEST_CASE( grid_dist_id_extended )
{
	// Domain
	Box<3,float> domain3({0.1,0.1,0.1},{1.1,1.1,1.1});

	long int k = 128*128*128*create_vcluster().getProcessingUnits();
	k = std::pow(k, 1/3.);

	Test3D_extended_grid(domain3,k);
}

BOOST_AUTO_TEST_CASE( grid_dist_id_periodic )
{
	// Domain
	Box<3,float> domain3({0.0,0.0,0.0},{1.0,1.0,1.0});

	long int k = 128*128*128*create_vcluster().getProcessingUnits();
	k = std::pow(k, 1/3.);

	Test3D_periodic(domain3,k);
}

BOOST_AUTO_TEST_CASE( grid_dist_id_unbound_ghost )
{
	// Domain
	Box<3,float> domain3({0.0,0.0,0.0},{1.0,1.0,1.0});

	long int k = 28*28*28*create_vcluster().getProcessingUnits();
	k = std::pow(k, 1/3.);

	Test3D_unb_ghost(domain3,k);
}

BOOST_AUTO_TEST_CASE( grid_dist_id_unbound_ghost_periodic )
{
	// Domain
	Box<3,float> domain3({0.0,0.0,0.0},{1.0,1.0,1.0});

	long int k = 25*25*25*create_vcluster().getProcessingUnits();
	k = std::pow(k, 1/3.);

	Test3D_unb_ghost_periodic(domain3,k);
}

BOOST_AUTO_TEST_CASE( grid_dist_id_copy )
{
	// Domain
	Box<3,float> domain3({0.0,0.0,0.0},{1.0,1.0,1.0});

	long int k = 32*32*32*create_vcluster().getProcessingUnits();
	k = std::pow(k, 1/3.);

	Test_grid_copy(domain3,k);
}

BOOST_AUTO_TEST_CASE( grid_1d_test )
{
	// Domain
	Box<1,float> domain1({-1.0},{1.0});

	long int k = 32*32*32*create_vcluster().getProcessingUnits();

	Test1D(domain1,k);
}

BOOST_AUTO_TEST_CASE( grid_dist_id_periodic_put_test )
{
	// Domain
	Box<3,float> domain3({0.0,0.0,0.0},{1.0,1.0,1.0});

	long int k = 128*128*128*create_vcluster().getProcessingUnits();
	k = std::pow(k, 1/3.);

	Test3D_periodic_put(domain3,k);
}

BOOST_AUTO_TEST_CASE( grid_dist_id_copy_test )
{
	// Domain
	Box<3,float> domain3({0.0,0.0,0.0},{1.0,1.0,1.0});

	long int k = 128*128*128*create_vcluster().getProcessingUnits();
	k = std::pow(k, 1/3.);

	Test3D_copy(domain3,k);
}


template<typename grid_amr>
void Test3D_ghost_put(grid_amr & g_dist_amr, long int k)
{
	// check the consistency of the decomposition
	bool val = g_dist_amr.getDecomposition().check_consistency();
	BOOST_REQUIRE_EQUAL(val,true);

	size_t sz[3] = {k,k,k};

	// Grid sm
	grid_sm<3,void> info(sz);

	size_t count = 0;

	auto dom = g_dist_amr.getGridIterator();

	while (dom.isNext())
	{
		auto key = dom.get_dist();

		g_dist_amr.template insert<0>(key) = -6.0;

		// Count the points
		count++;

		++dom;
	}

	// Set to zero the full grid

	{
	auto dom = g_dist_amr.getDomainIterator();

	while (dom.isNext())
	{
		auto key = dom.get();

		g_dist_amr.template insert<0>(key.move(0,1)) += 1.0;
		g_dist_amr.template insert<0>(key.move(0,-1)) += 1.0;
		g_dist_amr.template insert<0>(key.move(1,1)) += 1.0;
		g_dist_amr.template insert<0>(key.move(1,-1)) += 1.0;
		g_dist_amr.template insert<0>(key.move(2,1)) += 1.0;
		g_dist_amr.template insert<0>(key.move(2,-1)) += 1.0;

		++dom;
	}
	}

	bool correct = true;

	// Domain + Ghost iterator
	auto dom_gi = g_dist_amr.getDomainIterator();

	while (dom_gi.isNext())
	{
		auto key = dom_gi.get();

		correct &= (g_dist_amr.template get<0>(key) == 0);

		++dom_gi;
	}

	g_dist_amr.template ghost_put<add_,0>();

	if (count != 0)
	{BOOST_REQUIRE_EQUAL(correct, false);}

	// sync the ghosts
	g_dist_amr.template ghost_get<0>();

	correct = true;

	// Domain + Ghost iterator
	auto dom_gi2 = g_dist_amr.getDomainIterator();

	while (dom_gi2.isNext())
	{
		auto key = dom_gi2.get();

		correct &= (g_dist_amr.template get<0>(key) == 0);

		++dom_gi2;
	}

	BOOST_REQUIRE_EQUAL(correct, true);
}

BOOST_AUTO_TEST_CASE( grid_dist_domain_ghost_put_check )
{
	// Test grid periodic

	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	Vcluster & v_cl = create_vcluster();

	if ( v_cl.getProcessingUnits() > 32 )
	{return;}

	long int k = 13;

	BOOST_TEST_CHECKPOINT( "Testing grid periodic k<=" << k );

	// grid size
	size_t sz[3];
	sz[0] = k;
	sz[1] = k;
	sz[2] = k;

	// Ghost
	Ghost<3,long int> g(1);

	// periodicity
	periodicity<3> pr = {{PERIODIC,PERIODIC,PERIODIC}};

	// Distributed grid with id decomposition
	grid_dist_id<3, float, aggregate<long int>> g_dist(sz,domain,g,pr);

	Test3D_ghost_put(g_dist,k);

	// Distributed grid with id decomposition
	sgrid_dist_id<3, float, aggregate<long int>> sg_dist(sz,domain,g,pr);

	Test3D_ghost_put(sg_dist,k);
}


template<typename grid_amr>
void TestXD_ghost_put_create(grid_amr & g_dist_amr, long int k)
{
	// check the consistency of the decomposition
	bool val = g_dist_amr.getDecomposition().check_consistency();
	BOOST_REQUIRE_EQUAL(val,true);

	size_t count = 0;

	auto dom = g_dist_amr.getGridIterator();

	while (dom.isNext())
	{
		auto key = dom.get_dist();

		g_dist_amr.template insert<1>(key) = 1;

		// Count the points
		count++;

		++dom;
	}

	// Fill the ghost
	g_dist_amr.template ghost_get<1>();

	// Now we count the ghost point

	size_t g_point = 0;

	auto itg = g_dist_amr.getDomainGhostIterator();

	while (itg.isNext())
	{
		g_point++;

		++itg;
	}

	{
	auto it = g_dist_amr.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		g_dist_amr.remove(p);
		g_point--;

		++it;
	}
	}

	// A domain iterator should not produce points

	{
	auto it = g_dist_amr.getDomainIterator();

	size_t cnt = 0;
	while (it.isNext())
	{
		cnt++;

		++it;
	}

	BOOST_REQUIRE_EQUAL(cnt,0);
	}

	g_dist_amr.template ghost_put<add_,1>();

	{
	auto it = g_dist_amr.getDomainIterator();

	bool check = true;

	size_t cnt = 0;
	while (it.isNext())
	{
		auto p = it.get();

		cnt += g_dist_amr.template get<1>(p);

		check &= (g_dist_amr.template get<1>(p) >= 1);

		++it;
	}

	// Sum all the points
	auto & v_cl = create_vcluster();

	v_cl.sum(g_point);
	v_cl.sum(cnt);
	v_cl.execute();


	BOOST_REQUIRE_EQUAL(g_point,cnt);
	BOOST_REQUIRE_EQUAL(check,true);
	}
}

BOOST_AUTO_TEST_CASE( grid_dist_domain_ghost_2D_put_create_check )
{
	// Test grid periodic

	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	Vcluster & v_cl = create_vcluster();

	if ( v_cl.getProcessingUnits() > 32 )
	{return;}

	long int k = 13;

	BOOST_TEST_CHECKPOINT( "Testing grid periodic k<=" << k );

	// grid size
	size_t sz[2];
	sz[0] = k;
	sz[1] = k;

	// Ghost
	Ghost<2,long int> g(1);

	// periodicity
	periodicity<2> pr = {{PERIODIC,PERIODIC}};

	// Distributed grid with id decomposition
	sgrid_dist_id<2, float, aggregate<long int, int>> sg_dist(sz,domain,g,pr);

	TestXD_ghost_put_create(sg_dist,k);

	k = 7;
	sz[0] = k;
	sz[1] = k;

	// Distributed grid with id decomposition
	sgrid_dist_id<2, float, aggregate<long int, int>> sg_dist2(sz,domain,g,pr);

	TestXD_ghost_put_create(sg_dist2,k);

	k = 23;
	sz[0] = k;
	sz[1] = k;

	// Distributed grid with id decomposition
	sgrid_dist_id<2, float, aggregate<long int, int>> sg_dist3(sz,domain,g,pr);

	TestXD_ghost_put_create(sg_dist3,k);
}

BOOST_AUTO_TEST_CASE( grid_dist_domain_ghost_3D_put_create_check )
{
	// Test grid periodic

	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	Vcluster & v_cl = create_vcluster();

	if ( v_cl.getProcessingUnits() > 32 )
	{return;}

	long int k = 13;

	BOOST_TEST_CHECKPOINT( "Testing grid periodic k<=" << k );

	// grid size
	size_t sz[3];
	sz[0] = k;
	sz[1] = k;
	sz[2] = k;

	// Ghost
	Ghost<3,long int> g(1);

	// periodicity
	periodicity<3> pr = {{PERIODIC,PERIODIC,PERIODIC}};

	// Distributed grid with id decomposition
	sgrid_dist_id<3, float, aggregate<long int, int>> sg_dist(sz,domain,g,pr);

	TestXD_ghost_put_create(sg_dist,k);

	k = 7;
	sz[0] = k;
	sz[1] = k;
	sz[2] = k;

	// Distributed grid with id decomposition
	sgrid_dist_id<3, float, aggregate<long int, int>> sg_dist2(sz,domain,g,pr);

	TestXD_ghost_put_create(sg_dist2,k);

	k = 23;
	sz[0] = k;
	sz[1] = k;
	sz[2] = k;

	// Distributed grid with id decomposition
	sgrid_dist_id<3, float, aggregate<long int, int>> sg_dist3(sz,domain,g,pr);

	TestXD_ghost_put_create(sg_dist3,k);
}

BOOST_AUTO_TEST_SUITE_END()

