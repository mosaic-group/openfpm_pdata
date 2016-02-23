#ifndef CARTDECOMPOSITION_UNIT_TEST_HPP
#define CARTDECOMPOSITION_UNIT_TEST_HPP

#include "CartDecomposition.hpp"
#include "util/mathutil.hpp"

BOOST_AUTO_TEST_SUITE( CartDecomposition_test )

#define SUB_UNIT_FACTOR 64

BOOST_AUTO_TEST_CASE( CartDecomposition_non_periodic_test)
{
	// Vcluster
	Vcluster & vcl = *global_v_cluster;

	// Initialize the global VCluster
	init_global_v_cluster(&boost::unit_test::framework::master_test_suite().argc,&boost::unit_test::framework::master_test_suite().argv);

	//! [Create CartDecomposition]
	CartDecomposition<3,float> dec(vcl);

	// Physical domain
	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});
	size_t div[3];

	// Get the number of processor and calculate the number of sub-domain
	// for each processor (SUB_UNIT_FACTOR=64)
	size_t n_proc = vcl.getProcessingUnits();
	size_t n_sub = n_proc * SUB_UNIT_FACTOR;

	// Set the number of sub-domains on each dimension (in a scalable way)
	for (int i = 0 ; i < 3 ; i++)
	{div[i] = openfpm::math::round_big_2(pow(n_sub,1.0/3));}

	// Define ghost
	Ghost<3,float> g(0.01);

	// Boundary conditions
	size_t bc[] = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

	// Decompose
	dec.setParameters(div,box,bc,g);

	// create a ghost border
	dec.calculateGhostBoxes();

	//! [Create CartDecomposition]

	// For each calculated ghost box
	for (size_t i = 0 ; i < dec.getNIGhostBox() ; i++)
	{
		SpaceBox<3,float> b = dec.getIGhostBox(i);
		size_t proc = dec.getIGhostBoxProcessor(i);

		// sample one point inside the box
		Point<3,float> p = b.rnd();

		// Check that ghost_processorsID return that processor number
		const openfpm::vector<size_t> & pr = dec.template ghost_processorID<CartDecomposition<3,float>::processor_id>(p);

		bool found = false;

		for (size_t j = 0; j < pr.size() ; j++)
		{
			if (pr.get(j) == proc)
			{found = true; break;}
		}

		if (found == false)
		{
			const openfpm::vector<size_t> pr2 = dec.template ghost_processorID<CartDecomposition<3,float>::processor_id>(p);
		}

		BOOST_REQUIRE_EQUAL(found,true);
	}

	// Check the consistency

	bool val = dec.check_consistency();
	BOOST_REQUIRE_EQUAL(val,true);

	// We duplicate the decomposition
	CartDecomposition<3,float> dec2 = dec.duplicate();
	dec2.check_consistency();

	// check that dec and dec2 contain the same information
	bool ret = dec.is_equal(dec2);

	// We check if the two decomposition are equal
	BOOST_REQUIRE_EQUAL(ret,true);

	// We duplicate the decomposition redefining the ghost

	// Define ghost
	Ghost<3,float> g3(0.005);

	// We duplicate the decomposition redefining the ghost
	CartDecomposition<3,float> dec3 = dec.duplicate(g3);

	ret = dec3.check_consistency();
	BOOST_REQUIRE_EQUAL(ret,true);

	// Check that dec3 is equal to dec2 with the exception of the ghost part
	ret = dec3.is_equal_ng(dec2);
	BOOST_REQUIRE_EQUAL(ret,true);
}


BOOST_AUTO_TEST_CASE( CartDecomposition_periodic_test)
{
	// Vcluster
	Vcluster & vcl = *global_v_cluster;

	// Initialize the global VCluster
	init_global_v_cluster(&boost::unit_test::framework::master_test_suite().argc,&boost::unit_test::framework::master_test_suite().argv);

	//! [Create CartDecomposition]
	CartDecomposition<3,float> dec(vcl);

	// Physical domain
	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});
	size_t div[3];

	// Get the number of processor and calculate the number of sub-domain
	// for each processor (SUB_UNIT_FACTOR=64)
	size_t n_proc = vcl.getProcessingUnits();
	size_t n_sub = n_proc * SUB_UNIT_FACTOR;

	// Set the number of sub-domains on each dimension (in a scalable way)
	for (int i = 0 ; i < 3 ; i++)
	{div[i] = openfpm::math::round_big_2(pow(n_sub,1.0/3));}

	// Define ghost
	Ghost<3,float> g(0.01);

	// Boundary conditions
	size_t bc[] = {PERIODIC,PERIODIC,PERIODIC};

	// Decompose
	dec.setParameters(div,box,bc,g);

	// create a ghost border
	dec.calculateGhostBoxes();

	//! [Create CartDecomposition]

	// For each calculated ghost box
	for (size_t i = 0 ; i < dec.getNIGhostBox() ; i++)
	{
		SpaceBox<3,float> b = dec.getIGhostBox(i);
		size_t proc = dec.getIGhostBoxProcessor(i);

		// sample one point inside the box
		Point<3,float> p = b.rnd();

		// Check that ghost_processorsID return that processor number
		const openfpm::vector<size_t> & pr = dec.template ghost_processorID<CartDecomposition<3,float>::processor_id>(p);

		bool found = false;

		for (size_t j = 0; j < pr.size() ; j++)
		{
			if (pr.get(j) == proc)
			{found = true; break;}
		}

		if (found == false)
		{
			const openfpm::vector<size_t> pr2 = dec.template ghost_processorID<CartDecomposition<3,float>::processor_id>(p);
		}

		BOOST_REQUIRE_EQUAL(found,true);
	}

	// Check the consistency

	bool val = dec.check_consistency();
	BOOST_REQUIRE_EQUAL(val,true);

	// We duplicate the decomposition
	CartDecomposition<3,float> dec2 = dec.duplicate();
	dec2.check_consistency();

	bool ret = dec.is_equal(dec2);

	// We check if the two decomposition are equal
	BOOST_REQUIRE_EQUAL(ret,true);

	// check that dec and dec2 contain the same information

	// We duplicate the decomposition redefining the ghost

	// Define ghost
	Ghost<3,float> g3(0.005);

	// We duplicate the decomposition refefining the ghost
	CartDecomposition<3,float> dec3 = dec.duplicate(g3);

	ret = dec3.check_consistency();
	BOOST_REQUIRE_EQUAL(ret,true);

	// Check that g3 is equal to dec2 with the exception of the ghost part
	ret = dec3.is_equal_ng(dec2);
	BOOST_REQUIRE_EQUAL(ret,true);
}

BOOST_AUTO_TEST_CASE( CartDecomposition_extend_test)
{
	// Initialize the global VCluster
	init_global_v_cluster(&boost::unit_test::framework::master_test_suite().argc,&boost::unit_test::framework::master_test_suite().argv);

	// Vcluster
	Vcluster & vcl = *global_v_cluster;

	//! [Create CartDecomposition]
	CartDecomposition<3,float> dec(vcl);

	// Physical domain
	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});
	Box<3,float> bulk({0.05,0.05,0.05},{0.95,0.95,0.95});
	Box<3,float> bulk_g({0.01,0.01,0.01},{0.95,0.95,0.95});
	size_t div[3];

	// Get the number of processor and calculate the number of sub-domain
	// for each processor (SUB_UNIT_FACTOR=64)
	size_t n_proc = vcl.getProcessingUnits();
	size_t n_sub = n_proc * SUB_UNIT_FACTOR;

	// Set the number of sub-sub-domains on each dimension (in a scalable way)
	for (int i = 0 ; i < 3 ; i++)
	{div[i] = openfpm::math::round_big_2(pow(n_sub,1.0/3));}

	// Define ghost size
	Ghost<3,float> g(0.01);

	// Boundary conditions
	size_t bc[] = {PERIODIC,PERIODIC,PERIODIC};

	// Decompose
	dec.setParameters(div,box,bc,g);

	// Internally create the ghosts
	dec.calculateGhostBoxes();

	// Duplicate the decomposition extending the domain
	Box<3,float> box_ext({-0.1,-0.1,-0.1},{1.1,1.1,1.1});
	CartDecomposition<3,float> dec2 = dec.duplicate(g,box_ext);

	// than we create a grid iterator that iterate across the

	// Create a grid based iterator on the smaller decomposition
	size_t sz[3] = {128,128,128};
	grid_dist_id_iterator_dec<CartDecomposition<3,float>> git(dec,sz);

	// iterate across the grid points and we check that the results are consistent
	while (git.isNext())
	{
		auto key = git.get();

		Point<3,float> p = key.toPoint();

		if (bulk_g.isInside(p) == false)
		{
			++git;
			continue;
		}

		for (size_t i = 0 ; i < 3 ; i++)
			p.get(i) = key.get(i) * git.getSpacing(i);

		const openfpm::vector<size_t> & pr_pid = dec.template ghost_processorID<CartDecomposition<3,float>::processor_id>(p);
		const openfpm::vector<size_t> & pr2_pid = dec2.template ghost_processorID<CartDecomposition<3,float>::processor_id>(p);

		BOOST_REQUIRE(pr_pid == pr2_pid);

		const openfpm::vector<size_t> & pr_bid = dec.template ghost_processorID<CartDecomposition<3,float>::box_id>(p);
		const openfpm::vector<size_t> & pr2_bid = dec2.template ghost_processorID<CartDecomposition<3,float>::box_id>(p);

		BOOST_REQUIRE(pr_bid == pr2_bid);

		const openfpm::vector<size_t> & pr_lid = dec.template ghost_processorID<CartDecomposition<3,float>::lc_processor_id>(p);
		const openfpm::vector<size_t> & pr2_lid = dec2.template ghost_processorID<CartDecomposition<3,float>::lc_processor_id>(p);

		BOOST_REQUIRE(pr_lid == pr2_lid);

		size_t n_test = 0;

		// Given a point p, if the point is near the boundary
		if (bulk.isInside(p) == false)
		{
			// we can move p in one direction and check that dec2 return the same result
			for (size_t i = 0 ; i < 3 ; i++)
			{
				Point<3,float> q = p;

				q.get(i) = p.get(i) + 0.05;

				if (box.isInside(q) == false)
				{
					const openfpm::vector<size_t> & pr3_pid = dec2.template ghost_processorID<CartDecomposition<3,float>::processor_id>(q);
					const openfpm::vector<size_t> & pr3_bid = dec2.template ghost_processorID<CartDecomposition<3,float>::box_id>(q);
					const openfpm::vector<size_t> & pr3_lid = dec2.template ghost_processorID<CartDecomposition<3,float>::lc_processor_id>(q);
					BOOST_REQUIRE(pr2_pid == pr3_pid);
					BOOST_REQUIRE(pr2_bid == pr3_bid);
					BOOST_REQUIRE(pr2_lid == pr3_lid);

					n_test++;
				}

				q = p;
				q.get(i) = p.get(i) - 0.05;

				if (box.isInside(q) == false)
				{
					const openfpm::vector<size_t> & pr3_pid = dec2.template ghost_processorID<CartDecomposition<3,float>::processor_id>(q);
					const openfpm::vector<size_t> & pr3_bid = dec2.template ghost_processorID<CartDecomposition<3,float>::box_id>(q);
					const openfpm::vector<size_t> & pr3_lid = dec2.template ghost_processorID<CartDecomposition<3,float>::lc_processor_id>(q);

					if (pr2_pid != pr3_pid)
					{
						const openfpm::vector<size_t> & pr3_pid = dec2.template ghost_processorID<CartDecomposition<3,float>::processor_id>(q);

						int debug = 0;
						debug++;
					}

					BOOST_REQUIRE(pr2_pid == pr3_pid);
					BOOST_REQUIRE(pr2_bid == pr3_bid);
					BOOST_REQUIRE(pr2_lid == pr3_lid);

					n_test++;
				}
			}
		}

		BOOST_REQUIRE(n_test != 0);

		++git;
	}
}

BOOST_AUTO_TEST_SUITE_END()


#endif
