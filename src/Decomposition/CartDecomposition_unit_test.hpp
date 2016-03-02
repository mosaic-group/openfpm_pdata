#ifndef CARTDECOMPOSITION_UNIT_TEST_HPP
#define CARTDECOMPOSITION_UNIT_TEST_HPP

#include "CartDecomposition.hpp"
#include "util/mathutil.hpp"
#include "DLB.hpp"
#include <boost/algorithm/string.hpp>

BOOST_AUTO_TEST_SUITE (CartDecomposition_test)

#define SUB_UNIT_FACTOR 64
#define DIM 2

void setComputationCosts(CartDecomposition<DIM, float> &dec, size_t n_v, Point<DIM, float> center, float radius, size_t weight_h, size_t weight_l)
{
	float radius2 = pow(radius, 2);
	float eq;

	// Position structure for the single vertex
	float pos[DIM];

	for (int i = 0; i < n_v; i++)
	{
		dec.getSubSubDomainPosition(i, pos);

		eq = pow((pos[0] - center.get(0)), 2) + pow((pos[1] - center.get(1)), 2);

		if (eq <= radius2)
		{
			dec.setSubSubDomainComputationCost(i, weight_h);
		}
		else
		{
			dec.setSubSubDomainComputationCost(i, weight_l);
		}
	}
}

void setComputationCosts3D(CartDecomposition<3, float> &dec, size_t n_v, Point<3, float> center, float radius, size_t weight_h, size_t weight_l)
{
	float radius2 = pow(radius, 2);
	float eq;

	// Position structure for the single vertex
	float pos[3];

	for (int i = 0; i < n_v; i++)
	{
		dec.getSubSubDomainPosition(i, pos);

		eq = pow((pos[0] - center.get(0)), 2) + pow((pos[1] - center.get(1)), 2) + pow((pos[2] - center.get(2)), 2);

		if (eq <= radius2)
		{
			dec.setSubSubDomainComputationCost(i, weight_h);
		}
		else
		{
			dec.setSubSubDomainComputationCost(i, weight_l);
		}
	}
}

BOOST_AUTO_TEST_CASE( CartDecomposition_test_2D )
{
	//size_t balance_values_4p_64[] = {2.86,2.86,2.86,6.7,7.43,4.9,8.07,1.82,1.82,4.47,5.3};

	// Vcluster
	Vcluster & vcl = *global_v_cluster;

	// Initialize the global VCluster
	init_global_v_cluster(&boost::unit_test::framework::master_test_suite().argc,&boost::unit_test::framework::master_test_suite().argv);

	//! [Create CartDecomposition]
	CartDecomposition<DIM, float> dec(vcl);

	// Init DLB tool
	DLB dlb(vcl);

	// Physical domain
	Box<DIM, float> box( { 0.0, 0.0 }, { 10.0, 10.0 });
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

	// Set unbalance threshold
	dlb.setHeurisitc(DLB::Heuristic::UNBALANCE_THRLD);
	dlb.setThresholdLevel(DLB::ThresholdLevel::THRLD_MEDIUM);

	// Add weights to points

	// First create the center of the weights distribution, check it is coherent to the size of the domain
	Point<DIM, float> center( { 2.0, 2.0 });

	// Radius of the weights distribution
	float radius = 2.0;

	// Weight if the distribution (high)
	size_t weight_h = 5, weight_l = 1;

	setComputationCosts(dec, dec.getNSubSubDomains(), center, radius, weight_h, weight_l);

	dec.printCurrentDecomposition(0);

	dec.decompose();

	dec.printCurrentDecomposition(1);

	float stime = 0.0, etime = 10.0, tstep = 0.1;

	for(float t = stime, i = 1; t < etime; t = t + tstep, i++)
	{

		if(t < etime/2)
		{
			center.get(0) += tstep;
			center.get(1) += tstep;
		}
		else
		{
			center.get(0) -= tstep;
			center.get(1) -= tstep;
		}

		setComputationCosts(dec, dec.getNSubSubDomains(), center, radius, weight_h, weight_l);

		dlb.endIteration();

		dec.rebalance(dlb);

		dec.printCurrentDecomposition(i+1);

	}
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

BOOST_AUTO_TEST_CASE( CartDecomposition_test_2D_sar)
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
	Box<DIM, float> box( { 0.0, 0.0 }, { 10.0, 10.0 });
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

	// Set type of heuristic
	dlb.setHeurisitc(DLB::Heuristic::SAR_HEURISTIC);

	// Add weights to points

	// First create the center of the weights distribution, check it is coherent to the size of the domain
	Point<DIM, float> center( { 2.0, 2.0 });

	// Radius of the weights distribution
	float radius = 2.0;

	// Weight if the distribution (high)
	size_t weight_h = 5, weight_l = 1;

	size_t n_v = pow(div[0], DIM);

	setComputationCosts(dec, n_v, center, radius, weight_h, weight_l);

	dec.decompose();

	dec.printCurrentDecomposition(0);

	float stime = 0.0, etime = 10.0, tstep = 0.1;

	dlb.setSimulationStartTime(0);
	dlb.setSimulationEndTime(10);

	for(float t = stime, i = 1; t < etime; t = t + tstep, i++)
	{
		dlb.startIteration();

		if(t < etime/2)
		{
			center.get(0) += tstep;
			center.get(1) += tstep;
		}
		else
		{
			center.get(0) -= tstep;
			center.get(1) -= tstep;
		}

		setComputationCosts(dec, n_v, center, radius, weight_h, weight_l);

		sleep((n_v/dec.getProcessorLoad())/vcl.getProcessingUnits());

		dlb.endIteration();

		dec.rebalance(dlb);

		dec.printCurrentDecomposition(i);
	}

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

BOOST_AUTO_TEST_CASE( CartDecomposition_test_3D)
{
	// Vcluster
	Vcluster & vcl = *global_v_cluster;

	// Initialize the global VCluster
	init_global_v_cluster(&boost::unit_test::framework::master_test_suite().argc,&boost::unit_test::framework::master_test_suite().argv);

	//! [Create CartDecomposition]
	CartDecomposition<3, float> dec(vcl);

	// Init DLB tool
	DLB dlb(vcl);

	// Physical domain
	Box<3, float> box( { 0.0, 0.0, 0.0 }, { 10.0, 10.0, 10.0 });
	size_t div[3];

	// Get the number of processor and calculate the number of sub-domain
	// for each processor (SUB_UNIT_FACTOR=64)
	size_t n_proc = vcl.getProcessingUnits();
	size_t n_sub = n_proc * SUB_UNIT_FACTOR;

	// Set the number of sub-domains on each dimension (in a scalable way)
	for (int i = 0; i < 3; i++)
	{
		div[i] = openfpm::math::round_big_2(pow(n_sub,1.0/3));
	}

	// Define ghost
	Ghost<3, float> g(0.01);

	// Decompose
	dec.setParameters(div, box, g);

	// Set unbalance threshold
	dlb.setHeurisitc(DLB::Heuristic::UNBALANCE_THRLD);
	dlb.setThresholdLevel(DLB::ThresholdLevel::THRLD_MEDIUM);

	// Add weights to points

	// First create the center of the weights distribution, check it is coherent to the size of the domain
	Point<3, float> center( { 2.0, 2.0, 2.0 });

	// Radius of the weights distribution
	float radius = 2.0;

	// Weight if the distribution (high)
	size_t weight_h = 5, weight_l = 1;

	size_t n_v = pow(div[0], 3);

	setComputationCosts3D(dec, n_v, center, radius, weight_h, weight_l);

	dec.decompose();

	dec.printCurrentDecomposition(0);

	float stime = 0.0, etime = 10.0, tstep = 0.1;

	for(float t = stime, i = 1; t < etime; t = t + tstep, i++)
	{

		if(t < etime/2)
		{
			center.get(0) += tstep;
			center.get(1) += tstep;
			center.get(2) += tstep;
		}
		else
		{
			center.get(0) -= tstep;
			center.get(1) -= tstep;
			center.get(2) -= tstep;
		}

		setComputationCosts3D(dec, n_v, center, radius, weight_h, weight_l);

		dlb.endIteration();

		dec.rebalance(dlb);

		dec.printCurrentDecomposition(i);
	}

	// create a ghost border
	dec.calculateGhostBoxes();

	// For each calculated ghost box
	for (size_t i = 0; i < dec.getNIGhostBox(); i++)
	{
		SpaceBox<3,float> b = dec.getIGhostBox(i);
		size_t proc = dec.getIGhostBoxProcessor(i);

		// sample one point inside the box
		Point<3,float> p = b.rnd();

		// Check that ghost_processorsID return that processor number
		const openfpm::vector<size_t> & pr = dec.template ghost_processorID<CartDecomposition<3,float>::processor_id>(p);

		bool found = false;

		for (size_t j = 0; j < pr.size(); j++)
		{
			if (pr.get(j) == proc)
			{	found = true; break;}
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
}

BOOST_AUTO_TEST_SUITE_END()

#endif
