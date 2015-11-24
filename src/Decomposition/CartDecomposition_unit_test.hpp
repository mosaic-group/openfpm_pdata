#ifndef CARTDECOMPOSITION_UNIT_TEST_HPP
#define CARTDECOMPOSITION_UNIT_TEST_HPP

#include "CartDecomposition.hpp"
#include "util/mathutil.hpp"
#include "DLB.hpp"
#include <boost/algorithm/string.hpp>

BOOST_AUTO_TEST_SUITE (CartDecomposition_test)

#define SUB_UNIT_FACTOR 128
#define DIM 2

BOOST_AUTO_TEST_CASE( CartDecomposition_test_use)
{
	// Vcluster
	Vcluster & vcl = *global_v_cluster;

	// Initialize the global VCluster
	init_global_v_cluster(&boost::unit_test::framework::master_test_suite().argc,&boost::unit_test::framework::master_test_suite().argv);

	//! [Create CartDecomposition]
	CartDecomposition<DIM, float> dec(vcl);

	// Init DLB tool
	DLB dlb(vcl);

	// Physical domain
	Box<DIM, float> box( { 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
	size_t div[DIM];

	// Get the number of processor and calculate the number of sub-domain
	// for each processor (SUB_UNIT_FACTOR=64)
	size_t n_proc = vcl.getProcessingUnits();
	size_t n_sub = n_proc * SUB_UNIT_FACTOR;

	// Set the number of sub-domains on each dimension (in a scalable way)
	for (int i = 0; i < DIM; i++)
	{
		div[i] = openfpm::math::round_big_2(pow(n_sub,1.0/DIM));
	}

	// Define ghost
	Ghost<DIM, float> g(0.01);

	// Decompose
	dec.setParameters(div, box, g);

	//Add weights
	float c_x = 0, c_y = 0, c_z = 0, radius2, eq;
	float x = 0, y = 0, z = 0;
	float thickness = 2;
	int weight = 5;

	size_t n_v = pow(div[0], DIM);

	// Set center and radius of the distribution shape (Sphere)
	radius2 = pow(2,2);
	c_x = 2;
	c_y = 2;

	if (DIM == 3)
	c_z = 2;

	openfpm::vector<real_t> pos(DIM);
	for (int i = 0; i < n_v; i++)
	{
		dec.getSubSubDomainPosition(i, pos);
		x = pos.get(0) * 10;
		y = pos.get(1) * 10;

		if (DIM == 3)
		z = pos.get(2) * 10;

		eq = pow((x - c_x), 2) + pow((y - c_y), 2) + pow((z - c_z), 2);

		if (eq <= radius2)
		{
			dec.setSubSubDomainComputationCost(i, weight);
		}
		else
		{
			dec.setSubSubDomainComputationCost(i, 1);
		}
	}

	dec.decompose();

	dec.printCurrentDecomposition(0);

	float stime = 0.1, etime = 5, tstep = 0.1;

	dlb.setSimulationStartTime(0);
	dlb.setSimulationEndTime(5);

	for(real_t t = stime, i = 1, t_sim = 1;
			t < etime;
			t = t + tstep, i++, t_sim++)
	{
		dlb.setIterationStartTime(clock());

		if(t < etime/2)
		{
			c_x += tstep;
			c_y += tstep;
		}
		else
		{
			c_x -= tstep;
			c_y -= tstep;
		}

		if (DIM == 3)
		c_z += tstep;

		openfpm::vector<real_t> pos(DIM);
		for (int i = 0; i < n_v; i++)
		{
			dec.getSubSubDomainPosition(i, pos);
			x = pos.get(0) * 10;
			y = pos.get(1) * 10;

			if (DIM == 3)
			z = pos.get(2) * 10;

			eq = pow((x - c_x), 2) + pow((y - c_y), 2) + pow((z - c_z), 2);

			if (eq <= radius2)
			{
				dec.setSubSubDomainComputationCost(i, weight);
			}
			else
			{
				dec.setSubSubDomainComputationCost(i, 1);
			}
		}

		usleep(1000*t_sim);

		dlb.setIterationEndTime(clock());

		dec.rebalance(dlb);

		if(dlb.rebalanceNeeded())
		{
			t_sim = 1;
		}

		dec.printCurrentDecomposition(i);
	}

	//print statistics
	/*
	 if(vcl.getProcessUnitID() == 0)
	 {
	 float avg = dec.getTotalMovedV()/((etime-stime)/tstep);

	 std::cout << "Moved vertices average: " << avg << "\n";
	 std::cout << "Max number of moved vertices: " << dec.getMaxMovedV() << "\n";
	 }
	 */

	// create a ghost border
	dec.calculateGhostBoxes();

	// For each calculated ghost box
	for (size_t i = 0; i < dec.getNIGhostBox(); i++)
	{
		SpaceBox<DIM,float> b = dec.getIGhostBox(i);
		size_t proc = dec.getIGhostBoxProcessor(i);

		// sample one point inside the box
		Point<DIM,float> p = b.rnd();

		// Check that ghost_processorsID return that processor number
		const openfpm::vector<size_t> & pr = dec.template ghost_processorID<CartDecomposition<DIM,float>::processor_id>(p);

		bool found = false;

		for (size_t j = 0; j < pr.size(); j++)
		{
			if (pr.get(j) == proc)
			{	found = true; break;}
		}

		if (found == false)
		{
			const openfpm::vector<size_t> pr2 = dec.template ghost_processorID<CartDecomposition<DIM,float>::processor_id>(p);
		}

		BOOST_REQUIRE_EQUAL(found,true);
	}

	// Check the consistency

	bool val = dec.check_consistency();
	BOOST_REQUIRE_EQUAL(val,true);
}

BOOST_AUTO_TEST_SUITE_END()

#endif
