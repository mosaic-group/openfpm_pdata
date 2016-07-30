/*
 * test_isolation.hpp
 *
 *  Created on: Mar 11, 2016
 *      Author: i-bird
 */

#ifndef SRC_TEST_ISOLATION_HPP_
#define SRC_TEST_ISOLATION_HPP_

/*
 * Test that are not test but more are isolated here
 *
 *
 */

BOOST_AUTO_TEST_SUITE( antoniol_test_isolation )


BOOST_AUTO_TEST_CASE( CartDecomposition_test_2D )
{
	//size_t balance_values_4p_64[] = {2.86,2.86,2.86,6.7,7.43,4.9,8.07,1.82,1.82,4.47,5.3};

	// Vcluster
	Vcluster & vcl = *global_v_cluster;

	// non-periodic boundary condition
	size_t bc[2] = { NON_PERIODIC, NON_PERIODIC };

	// Initialize the global VCluster
	init_global_v_cluster(&boost::unit_test::framework::master_test_suite().argc,&boost::unit_test::framework::master_test_suite().argv);

	//! [Create CartDecomposition]
	CartDecomposition<2, float> dec(vcl);

	// Init DLB tool
	DLB dlb(vcl);

	// Physical domain
	Box<2, float> box( { 0.0, 0.0 }, { 10.0, 10.0 });
	size_t div[2];

	// Get the number of processor and calculate the number of sub-domain
	// for each processor (SUB_UNIT_FACTOR=64)
	size_t n_proc = vcl.getProcessingUnits();
	size_t n_sub = n_proc * SUB_UNIT_FACTOR;

	// Set the number of sub-domains on each dimension (in a scalable way)
	for (int i = 0; i < 2; i++)
	{
		div[i] = openfpm::math::round_big_2(pow(n_sub,1.0/2));
	}

	// Define ghost
	Ghost<2, float> g(0.01);

	// Decompose
	dec.setParameters(div, box, bc, g);

	// Set unbalance threshold
	dlb.setHeurisitc(DLB::Heuristic::UNBALANCE_THRLD);
	dlb.setThresholdLevel(DLB::ThresholdLevel::THRLD_MEDIUM);

	// Add weights to points

	// First create the center of the weights distribution, check it is coherent to the size of the domain
	Point<2, float> center( { 2.0, 2.0 });

	// Radius of the weights distribution
	float radius = 2.0;

	// Weight if the distribution (high)
	size_t weight_h = 5, weight_l = 1;

	setComputationCosts(dec, dec.getNSubSubDomains(), center, radius, weight_h, weight_l);

	dec.getDistribution().write("DLB_test_graph_0.vtk");

	dec.decompose();

	dec.getDistribution().write("DLB_test_graph_1.vtk");

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

		if(dec.rebalance(dlb))
			dec.write("DLB_test_graph_cart_" + std::to_string(i+1) + "_");

		std::stringstream str;
		str << "DLB_test_graph_" << i + 1 << ".vtk";
		dec.getDistribution().write(str.str());
	}

	// For each calculated ghost box
	for (size_t i = 0; i < dec.getNIGhostBox(); i++)
	{
		SpaceBox<2,float> b = dec.getIGhostBox(i);
		size_t proc = dec.getIGhostBoxProcessor(i);

		// sample one point inside the box
		Point<2,float> p = b.rnd();

		// Check that ghost_processorsID return that processor number
		const openfpm::vector<size_t> & pr = dec.template ghost_processorID<CartDecomposition<2,float>::processor_id>(p);

		bool found = false;

		for (size_t j = 0; j < pr.size(); j++)
		{
			if (pr.get(j) == proc)
			{	found = true; break;}
		}

		if (found == false)
		{
			const openfpm::vector<size_t> pr2 = dec.template ghost_processorID<CartDecomposition<2,float>::processor_id>(p);
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

	// non-periodic boundary condition
	size_t bc[2] = { NON_PERIODIC, NON_PERIODIC };

	// Initialize the global VCluster
	init_global_v_cluster(&boost::unit_test::framework::master_test_suite().argc,&boost::unit_test::framework::master_test_suite().argv);

	//! [Create CartDecomposition]
	CartDecomposition<2, float> dec(vcl);

	// Init DLB tool
	DLB dlb(vcl);

	// Physical domain
	Box<2, float> box( { 0.0, 0.0 }, { 10.0, 10.0 });
	size_t div[2];

	// Get the number of processor and calculate the number of sub-domain
	// for each processor (SUB_UNIT_FACTOR=64)
	size_t n_proc = vcl.getProcessingUnits();
	size_t n_sub = n_proc * SUB_UNIT_FACTOR;

	// Set the number of sub-domains on each dimension (in a scalable way)
	for (int i = 0; i < 2; i++)
	{
		div[i] = openfpm::math::round_big_2(pow(n_sub,1.0/2));
	}

	// Define ghost
	Ghost<2, float> g(0.01);

	// Decompose
	dec.setParameters(div, box, bc, g);

	// Set type of heuristic
	dlb.setHeurisitc(DLB::Heuristic::SAR_HEURISTIC);

	// Add weights to points

	// First create the center of the weights distribution, check it is coherent to the size of the domain
	Point<2, float> center( { 2.0, 2.0 });

	// Radius of the weights distribution
	float radius = 2.0;

	// Weight if the distribution (high)
	size_t weight_h = 5, weight_l = 1;

	size_t n_v = pow(div[0], 2);

	setComputationCosts(dec, n_v, center, radius, weight_h, weight_l);

	dec.decompose();

	dec.getDistribution().write("DLB_test_graph_0.vtk");

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

		if(dec.rebalance(dlb))
		{
			dec.write("DLB_test_graph_cart_" + std::to_string(i) + "_");
		}

		std::stringstream str;
		str << "DLB_test_graph_" << i << ".vtk";
		dec.getDistribution().write(str.str());
	}

	// For each calculated ghost box
	for (size_t i = 0; i < dec.getNIGhostBox(); i++)
	{
		SpaceBox<2,float> b = dec.getIGhostBox(i);
		size_t proc = dec.getIGhostBoxProcessor(i);

		// sample one point inside the box
		Point<2,float> p = b.rnd();

		// Check that ghost_processorsID return that processor number
		const openfpm::vector<size_t> & pr = dec.template ghost_processorID<CartDecomposition<2,float>::processor_id>(p);

		bool found = false;

		for (size_t j = 0; j < pr.size(); j++)
		{
			if (pr.get(j) == proc)
			{	found = true; break;}
		}

		if (found == false)
		{
			const openfpm::vector<size_t> pr2 = dec.template ghost_processorID<CartDecomposition<2,float>::processor_id>(p);
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

	// non-periodic boundary condition
	size_t bc[3] = { NON_PERIODIC, NON_PERIODIC, NON_PERIODIC };

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
	dec.setParameters(div, box, bc, g);

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

	dec.getDistribution().write("DLB_test_graph_0.vtk");

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

		std::stringstream str;
		str << "DLB_test_graph_" << i << ".vtk";
		dec.getDistribution().write(str.str());
	}

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


BOOST_AUTO_TEST_CASE( dist_map_graph_use_cartesian)
{
	//! Vcluster
	Vcluster & vcl = *global_v_cluster;

//	if(vcl.getProcessingUnits() != 3)
//	return;

	// non-periodic boundary condition
	size_t bc[3] = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

	//! [Create CartDecomposition]
	CartDecomposition<3, float> dec(vcl);

	// Physical domain
	Box<3, float> box( { 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
	size_t div[3] = {8,8,8};

	// Grid size and info
	size_t gsz[3] = {8,8,8};
	grid_sm<3,void> g_sm(gsz);

	// Define ghost
	Ghost<3, float> g(0.01);

	// Decompose
	dec.setParameters(div, box, bc, g);
	dec.decompose();

	grid_dist_id_iterator_dec<CartDecomposition<3,float>> it_dec(dec,gsz);

	size_t cnt = 0;
	while (it_dec.isNext())
	{
		cnt++;
		++it_dec;
	}

	openfpm::vector<size_t> v_cnt(vcl.getProcessingUnits());

	// Sent and receive the size of each subgraph
	vcl.allGather(cnt, v_cnt);
	vcl.execute();

	cnt = 0;
	for (long int i = 0; i <= ((long int)vcl.getProcessUnitID()) - 1 ; ++i)
		cnt += v_cnt.get(i);

	// count the points

	//! Distributed graph
	DistGraph_CSR<aggregate<size_t[3]>, aggregate<size_t>> dg;

	grid_dist_id_iterator_dec<CartDecomposition<3,float>> it_dec2(dec,gsz);
	while (it_dec2.isNext())
	{
		auto key = it_dec2.get();

		aggregate<size_t[3]> v;

		v.template get<0>()[0] = key.get(0);
		v.template get<0>()[1] = key.get(1);
		v.template get<0>()[2] = key.get(2);

		size_t gid = g_sm.LinId(key);
		dg.add_vertex(v, gid, cnt);

		cnt++;
		++it_dec2;
	}

	dg.init();

	// we ask for some random vertex

	std::default_random_engine rg;
	std::uniform_int_distribution<size_t> d(0,g_sm.size()-1);

	openfpm::vector<size_t> v_req;

/*	for (size_t i = 0 ; i < 16 ; i++)
	{
		size_t v = d(rg);*/

		if (vcl.getProcessUnitID() == 0)
			dg.reqVertex(450);

/*		dg.reqVertex(v);
	}*/

	dg.sync();

	if (vcl.getProcessUnitID() == 0)
	{
		grid_key_dx<3> key;
		// get the position information
		key.set_d(0,dg.getVertex(450).template get<0>()[0]);
		key.set_d(1,dg.getVertex(450).template get<0>()[1]);
		key.set_d(2,dg.getVertex(450).template get<0>()[2]);

		size_t lin_id = g_sm.LinId(key);

	//	BOOST_REQUIRE_EQUAL(lin_id,v_req.get(i));

		std::cout << "Error: " << "   " << lin_id << "    " << key.to_string() << "\n";
	}

/*	for (size_t i = 0 ; i < 16 ; i++)
	{
		grid_key_dx<3> key;
		// get the position information
		key.set_d(0,dg.getVertex(v_req.get(i)).template get<0>()[0]);
		key.set_d(1,dg.getVertex(v_req.get(i)).template get<0>()[1]);
		key.set_d(2,dg.getVertex(v_req.get(i)).template get<0>()[2]);

		size_t lin_id = g_sm.LinId(key);

	//	BOOST_REQUIRE_EQUAL(lin_id,v_req.get(i));

		std::cout << "Error: " << i << "   " << lin_id << "  " << v_req.get(i) << "\n";
	}*/

/*	if (vcl.getProcessUnitID() == 0)
		std::cout << "Error: " << i << "   " << lin_id << "  " << v_req.get(i) << "\n";*/
}

BOOST_AUTO_TEST_CASE( Parmetis_distribution_test_random_walk )
{
	typedef Point<3,float> s;

	Vcluster & v_cl = *global_v_cluster;

	// set the seed
	// create the random generator engine
	std::srand(v_cl.getProcessUnitID());
	std::default_random_engine eg;
	std::uniform_real_distribution<float> ud(0.0f, 1.0f);

	size_t nsz[] = { 0, 32, 4 };
	nsz[0] = 65536 * v_cl.getProcessingUnits();

	print_test_v( "Testing 3D random walk vector k<=",nsz[0]);

	// 3D test
	for (size_t i = 0; i < 3; i++ )
	{
		size_t k = nsz[i];

		BOOST_TEST_CHECKPOINT( "Testing 3D random walk k=" << k );

		Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

		// Grid info
		grid_sm<3, void> info( { GS_SIZE, GS_SIZE, GS_SIZE });

		// Boundary conditions
		size_t bc[3] = { NON_PERIODIC,NON_PERIODIC,NON_PERIODIC };

		// factor
		float factor = pow(global_v_cluster->getProcessingUnits()/2.0f,1.0f/3.0f);

		// ghost
		Ghost<3,float> ghost(0.01 / factor);

		// Distributed vector
		vector_dist<3,float, Point_test<float>, CartDecomposition<3, float, HeapMemory, ParMetisDistribution<3, float>>> vd(k,box,bc,ghost);

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

		vd.addComputationCosts();

		vd.getDecomposition().getDistribution().write("parmetis_random_walk_" + std::to_string(0) + ".vtk");

		// 10 step random walk

		for (size_t j = 0; j < 10; j++)
		{
			auto it = vd.getDomainIterator();

			while (it.isNext())
			{
				auto key = it.get();

				vd.template getPos<s::x>(key)[0] += 0.01 * ud(eg);
				vd.template getPos<s::x>(key)[1] += 0.01 * ud(eg);
				vd.template getPos<s::x>(key)[2] += 0.01 * ud(eg);

				++it;
			}

			vd.map();

			/////// Interactions ///

			//vd.ghost_get<>();
			//vd.getDomainIterator;

			////////////////////////

			vd.addComputationCosts();

			vd.getDecomposition().rebalance(10);

			vd.map();

			vd.getDecomposition().getDistribution().write("parmetis_random_walk_" + std::to_string(j+1) + ".vtk");

			size_t l = vd.size_local();
			v_cl.sum(l);
			v_cl.execute();

			// Count the local particles and check that the total number is consistent
			size_t cnt = total_n_part_lc(vd,bc);

			//BOOST_REQUIRE_EQUAL((size_t)k,cnt);
		}
	}
}

BOOST_AUTO_TEST_CASE( Parmetis_distribution_test_random_walk_2D )
{

	//Particle: position, type of poistion, type of animal (0 rabbit, 1 fox), dead or alive (0 or 1), time the fox stays alive without eating
	typedef Point<2,float> s;

	Vcluster & v_cl = *global_v_cluster;

	size_t JUST_EATEN = 5;

	// set the seed
	// create the random generator engine
	std::srand(v_cl.getProcessUnitID());
	std::default_random_engine eg;
	std::uniform_real_distribution<float> ud(0.0f, 0.5f);

	size_t k = 100000;

	print_test_v( "Testing 2D random walk vector k<=",k);

	BOOST_TEST_CHECKPOINT( "Testing 2D random walk k=" << k );

	Box<2,float> box({0.0,0.0},{1.0,1.0});

	// Grid info
	grid_sm<2, void> info( { GS_SIZE, GS_SIZE });

	// Boundary conditions
	size_t bc[2] = { PERIODIC, PERIODIC };

	// factor
	float factor = pow(global_v_cluster->getProcessingUnits()/2.0f,1.0f/3.0f);

	// ghost
	Ghost<2,float> ghost(0.01 / factor);

	// Distributed vector
	vector_dist<2,float, Point_test<float>, CartDecomposition<2, float, HeapMemory, ParMetisDistribution<2, float>>> vd(k,box,bc,ghost);

	// Init DLB tool
	DLB dlb(v_cl);

	// Set unbalance threshold
	dlb.setHeurisitc(DLB::Heuristic::UNBALANCE_THRLD);
	dlb.setThresholdLevel(DLB::ThresholdLevel::THRLD_MEDIUM);

	auto it = vd.getIterator();

	size_t c = 0;
	while (it.isNext())
	{
		auto key = it.get();
		if(c % 5)
		{
			vd.template getPos<s::x>(key)[0] = ud(eg);
			vd.template getPos<s::x>(key)[1] = ud(eg);
		}else{
			vd.template getPos<s::x>(key)[0] = ud(eg)*2;
			vd.template getPos<s::x>(key)[1] = ud(eg)*2;
		}
		++it;
		++c;
	}

	vd.map();

	vd.addComputationCosts();

	vd.getDecomposition().rebalance(dlb);

	vd.map();

	vd.getDecomposition().write("dec_init");
	vd.getDecomposition().getDistribution().write("parmetis_random_walk_" + std::to_string(0) + ".vtk");
	vd.write("particles_", 0, NO_GHOST);

	// 10 step random walk
	for (size_t j = 0; j < 50; j++)
	{
		std::cout << "Iteration " << (j+1) << "\n";

		auto it = vd.getDomainIterator();

		while (it.isNext())
		{
			auto key = it.get();

			vd.template getPos<s::x>(key)[0] += 0.01 * ud(eg);
			vd.template getPos<s::x>(key)[1] += 0.01 * ud(eg);

			++it;
		}

		vd.map();

		vd.addComputationCosts();

		vd.getDecomposition().rebalance(dlb);

		vd.map();

		vd.getDecomposition().getDistribution().write("parmetis_random_walk_" + std::to_string(j+1) + ".vtk");
		vd.write("particles_", j+1, NO_GHOST);
		vd.getDecomposition().write("dec_");

		size_t l = vd.size_local();
		v_cl.sum(l);
		v_cl.execute();

		// Count the local particles and check that the total number is consistent
		//size_t cnt = total_n_part_lc(vd,bc);

		//BOOST_REQUIRE_EQUAL((size_t)k,cnt);
	}

}


BOOST_AUTO_TEST_CASE( Parmetis_distribution_test_prey_and_predators )
{
	Vcluster & v_cl = *global_v_cluster;

	//time the animal stays alive without eating or reproducing
	size_t TIME_A = 5;

	size_t PREDATOR = 1, PREY = 0;
	size_t ALIVE = 1, DEAD = 0;

	// Predators reproducing probability
	float PRED_REPR = 0.1;

	// Predators eating probability
	float PRED_EAT = 0.2;

	// Prey reproducing probability
	float PREY_REPR = 0.1;

	// set the seed
	// create the random generator engine
	std::srand(v_cl.getProcessUnitID());
	std::default_random_engine eg;
	std::uniform_real_distribution<float> ud(0.0f, 1.0f);

	size_t k = 50000;

	print_test_v( "Testing 2D random walk vector k<=",k);

	BOOST_TEST_CHECKPOINT( "Testing 2D random walk k=" << k );

	Box<2,float> box({0.0,0.0},{1.0,1.0});

	// Grid info
	grid_sm<2, void> info( { GS_SIZE, GS_SIZE });

	// Boundary conditions
	size_t bc[2] = { PERIODIC, PERIODIC };

	// factor
	float factor = pow(global_v_cluster->getProcessingUnits()/2.0f,1.0f/3.0f);

	// interaction radius
	float r_cut = 0.01 / factor;

	// ghost
	Ghost<2,float> ghost(0.01 / factor);

	// Distributed vector
	vector_dist<2,float, animal, CartDecomposition<2, float, HeapMemory, ParMetisDistribution<2, float>>> vd(k,box,bc,ghost);

	// Init DLB tool
	DLB dlb(v_cl);

	// Set unbalance threshold
	dlb.setHeurisitc(DLB::Heuristic::UNBALANCE_THRLD);
	dlb.setThresholdLevel(DLB::ThresholdLevel::THRLD_MEDIUM);

	auto it = vd.getIterator();

	size_t c = 0;
	while (it.isNext())
	{
		auto key = it.get();
		if(c % 3)
		{
			vd.template getPos<animal::pos>(key)[0] = ud(eg);
			vd.template getPos<animal::pos>(key)[1] = ud(eg);
			vd.template getProp<animal::genre>(key) = 0; //prey
			vd.template getProp<animal::status>(key) = 1; //alive
			vd.template getProp<animal::time_a>(key) = TIME_A; //alive
		}else{
			vd.template getPos<animal::pos>(key)[0] = ud(eg);
			vd.template getPos<animal::pos>(key)[1] = ud(eg);
			vd.template getProp<animal::genre>(key) = 1; //predator
			vd.template getProp<animal::status>(key) = 1; //alive
			vd.template getProp<animal::time_a>(key) = TIME_A; //alive
		}
		++it;
		++c;
	}

	vd.map();

	vd.addComputationCosts();

	vd.getDecomposition().rebalance(dlb);

	vd.map();

	vd.getDecomposition().write("dec_init");
	vd.getDecomposition().getDistribution().write("parmetis_random_walk_" + std::to_string(0) + ".vtk");
	vd.write("particles_", 0, NO_GHOST);

	// 10 step random walk
	for (size_t j = 0; j < 50; j++)
	{
		std::cout << "Iteration " << (j+1) << "\n";

		auto it = vd.getDomainIterator();

		while (it.isNext())
		{
			auto key = it.get();

			vd.template getPos<animal::pos>(key)[0] += 0.01 * ud(eg);
			vd.template getPos<animal::pos>(key)[1] += 0.01 * ud(eg);

			++it;
		}

		vd.map();

		/////// Interactions ///

		vd.ghost_get<0>();

		openfpm::vector<size_t> deads;

		// get the cell list with a cutoff radius

		bool error = false;

		auto NN = vd.getCellList(0.01 / factor);

		// iterate across the domain particle

		auto it2 = vd.getDomainIterator();

		while (it2.isNext())
		{
			auto p = it2.get();

			Point<2,float> xp = vd.getPos<0>(p);

			size_t gp = vd.getProp<animal::genre>(p);
			size_t sp = vd.getProp<animal::status>(p);

			auto Np = NN.getIterator(NN.getCell(vd.getPos<0>(p)));

			while (Np.isNext())
			{
				auto q = Np.get();

				size_t gq = vd.getProp<animal::genre>(q);
				size_t sq = vd.getProp<animal::status>(q);

				// repulsive

				Point<2,float> xq = vd.getPos<0>(q);
				Point<2,float> f = (xp - xq);

				float distance = f.norm();

				//if p is a fox and q a rabit and they are both alive then the fox eats the rabbit
				if (distance < 2*r_cut*sqrt(2) && sp == ALIVE)
				{
					if(gp == PREDATOR && gq == PREY && sq == ALIVE)
					{
						vd.getProp<animal::status>(q) = DEAD;
						vd.getProp<animal::time_a>(q) = TIME_A;
					}
					else if (gp == PREY && gq == PREY && sq != DEAD)
					{
						vd.add();
						//vd.getLastPos<animal::pos>()[0] = ud(eg);
						//vd.getLastPos<animal::pos>()[1] = ud(eg);
						vd.getLastProp<animal::genre>() = 0;
					}
				}

				++Np;
			}

			if(vd.getProp<animal::status>(p) == DEAD)
			{
				deads.add(p.getKey());
			}

			++it2;
		}

		vd.remove(deads, 0);
		deads.resize(0);

		////////////////////////

		vd.addComputationCosts();

		vd.getDecomposition().rebalance(dlb);

		vd.map();

		vd.getDecomposition().getDistribution().write("parmetis_random_walk_" + std::to_string(j+1) + ".vtk");
		vd.write("particles_", j+1, NO_GHOST);
		vd.getDecomposition().write("dec_");

		size_t l = vd.size_local();
		v_cl.sum(l);
		v_cl.execute();

		// Count the local particles and check that the total number is consistent
		//size_t cnt = total_n_part_lc(vd,bc);

		//BOOST_REQUIRE_EQUAL((size_t)k,cnt);
	}

}

BOOST_AUTO_TEST_SUITE_END()

#endif /* SRC_TEST_ISOLATION_HPP_ */
