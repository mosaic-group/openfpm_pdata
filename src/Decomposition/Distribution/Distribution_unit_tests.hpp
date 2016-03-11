/*
 * Distribution_unit_tests.hpp
 *
 *  Created on: Feb 27, 2016
 *      Author: i-bird
 */

#ifndef SRC_DECOMPOSITION_DISTRIBUTION_DISTRIBUTION_UNIT_TESTS_HPP_
#define SRC_DECOMPOSITION_DISTRIBUTION_DISTRIBUTION_UNIT_TESTS_HPP_

/*! \brief Set a sphere as high computation cost
 *
 * \param dist Distribution structure
 * \param gr grid info
 * \param center of the sphere
 * \param radius radius of the sphere
 * \param max_l maximum load of the processor
 * \param min_l minimum load of the processor
 *
 */
template<unsigned int dim, typename Distribution> void setSphereComputationCosts(Distribution & dist, grid_sm<dim, void> & gr, Point<3, float> center, float radius, size_t max_l, size_t min_l)
{
	float radius2 = radius * radius;
	float eq;

	// Position structure for the single vertex
	float pos[dim];

	for (size_t i = 0; i < dist.getNSubSubDomains(); i++)
	{
		dist.getSubSubDomainPosition(i, pos);

		eq = 0;
		for (size_t j = 0; j < dim; j++)
			eq += (pos[j] - center.get(j)) * (pos[j] - center.get(j));

		if (eq <= radius2)
		{
			dist.setComputationCost(i, max_l);
			dist.setMigrationCost(i, max_l * 2);
		}
		else
		{
			dist.setComputationCost(i, min_l);
			dist.setMigrationCost(i, min_l * 2);
		}

		// set Migration cost and communication cost
		for (size_t j = 0; j < dist.getNSubSubDomainNeighbors(i); j++)
			dist.setCommunicationCost(i, j, 1);
	}
}

struct animal
{
	typedef boost::fusion::vector<float[2], size_t, size_t, size_t> type;

	//! Attributes name
	struct attributes
	{
		static const std::string name[];
	};

	//! type of the positional field
	typedef float s_type;

	//! The data
	type data;

	//! position property id in boost::fusion::vector
	static const unsigned int pos = 0;
	//! genre of animal property id in boost::fusion::vector
	static const unsigned int genre = 1;
	//! state property id in boost::fusion::vector
	static const unsigned int status = 2;
	//! alive time property id in boost::fusion::vector
	static const unsigned int time_a = 3;

	//! total number of properties boost::fusion::vector
	static const unsigned int max_prop = 4;

	animal()
	{
	}

	inline animal(const animal & p)
	{
		boost::fusion::at_c<0>(data)[0] = boost::fusion::at_c<0>(p.data)[0];
		boost::fusion::at_c<0>(data)[1] = boost::fusion::at_c<0>(p.data)[1];
		//boost::fusion::at_c<0>(data)[2] = boost::fusion::at_c<0>(p.data)[2];
		boost::fusion::at_c<1>(data) = boost::fusion::at_c<1>(p.data);
		boost::fusion::at_c<2>(data) = boost::fusion::at_c<2>(p.data);
		boost::fusion::at_c<3>(data) = boost::fusion::at_c<3>(p.data);
	}

	template<unsigned int id> inline auto get() -> decltype(boost::fusion::at_c < id > (data))
	{
		return boost::fusion::at_c<id>(data);
	}

	template<unsigned int id> inline auto get() const -> const decltype(boost::fusion::at_c < id > (data))
	{
		return boost::fusion::at_c<id>(data);
	}

	template<unsigned int dim, typename Mem> inline animal(const encapc<dim, animal, Mem> & p)
	{
		this->operator=(p);
	}

	template<unsigned int dim, typename Mem> inline animal & operator=(const encapc<dim, animal, Mem> & p)
	{
		boost::fusion::at_c<0>(data)[0] = p.template get<0>()[0];
		boost::fusion::at_c<0>(data)[1] = p.template get<0>()[1];
		//boost::fusion::at_c<0>(data)[2] = p.template get<0>()[2];
		boost::fusion::at_c<1>(data) = p.template get<1>();
		boost::fusion::at_c<2>(data) = p.template get<2>();
		boost::fusion::at_c<3>(data) = p.template get<3>();

		return *this;
	}

	static bool noPointers()
	{
		return true;
	}
};

const std::string animal::attributes::name[] = { "pos", "genre", "status", "time_a" };

BOOST_AUTO_TEST_SUITE (Distribution_test)

BOOST_AUTO_TEST_CASE( Metis_distribution_test)
{
	Vcluster & v_cl = *global_v_cluster;

	if (v_cl.getProcessingUnits() != 3)
	return;

	if (v_cl.getProcessUnitID() != 0)
	return;

	//! [Initialize a Metis Cartesian graph and decompose]

	MetisDistribution<3, float> met_dist(v_cl);

	// Cartesian grid
	size_t sz[3] = { GS_SIZE, GS_SIZE, GS_SIZE };

	// Box
	Box<3, float> box( { 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });

	// Grid info
	grid_sm<3, void> info(sz);

	// Set metis on test, It fix the seed (not required if we are not testing)
	met_dist.onTest();

	// Initialize Cart graph and decompose

	met_dist.createCartGraph(info,box);
	met_dist.decompose();

	//! [Initialize a Metis Cartesian graph and decompose]

	BOOST_REQUIRE(met_dist.getUnbalance() < 0.03);

	met_dist.write("vtk_metis_distribution.vtk");

	size_t b = GS_SIZE * GS_SIZE * GS_SIZE / 5;

	//! [Decomposition Metis with weights]

	// Initialize the weights to 1.0
	// not required, if we set ALL Computation,Migration,Communication cost
	met_dist.initWeights();

	// Change set some weight on the graph and re-decompose

	for (size_t i = 0; i < met_dist.getNSubSubDomains(); i++)
	{
		if (i == 0 || i == b || i == 2*b || i == 3*b || i == 4*b)
		met_dist.setComputationCost(i,10);
		else
		met_dist.setComputationCost(i,1);

		// We also show how to set some Communication and Migration cost

		met_dist.setMigrationCost(i,1);

		for (size_t j = 0; j < met_dist.getNSubSubDomainNeighbors(i); j++)
		met_dist.setCommunicationCost(i,j,1);
	}

	met_dist.decompose();

	//! [Decomposition Metis with weights]

	BOOST_REQUIRE(met_dist.getUnbalance() < 0.03);

	met_dist.write("vtk_metis_distribution_red.vtk");

	// check that match

	bool test = compare("vtk_metis_distribution.vtk", "src/Decomposition/Distribution/test_data/vtk_metis_distribution_test.vtk");
	BOOST_REQUIRE_EQUAL(true,test);

	test = compare("vtk_metis_distribution_red.vtk","src/Decomposition/Distribution/test_data/vtk_metis_distribution_red_test.vtk");
	BOOST_REQUIRE_EQUAL(true,test);

	// Copy the Metis distribution

	MetisDistribution<3, float> met_dist2(v_cl);

	met_dist2 = met_dist;

	test = (met_dist2 == met_dist);

	BOOST_REQUIRE_EQUAL(test,true);

	// We fix the size of MetisDistribution if you are gointg to change this number
	// please check the following
	// duplicate functions
	// swap functions
	// Copy constructors
	// operator= functions
	// operator== functions

	BOOST_REQUIRE_EQUAL(sizeof(MetisDistribution<3,float>),568ul);
}

BOOST_AUTO_TEST_CASE( Parmetis_distribution_test)
{
	Vcluster & v_cl = *global_v_cluster;

	if (v_cl.getProcessingUnits() != 3)
	return;

	//! [Initialize a ParMetis Cartesian graph and decompose]

	ParMetisDistribution<3, float> pmet_dist(v_cl);

	// Physical domain
	Box<3, float> box( { 0.0, 0.0, 0.0 }, { 10.0, 10.0, 10.0 });

	// Grid info
	grid_sm<3, void> info( { GS_SIZE, GS_SIZE, GS_SIZE });

	// Initialize Cart graph and decompose
	pmet_dist.createCartGraph(info,box);

	// First create the center of the weights distribution, check it is coherent to the size of the domain
	Point<3, float> center( { 2.0, 2.0, 2.0 });

	// It produces a sphere of radius 2.0
	// with high computation cost (5) inside the sphere and (1) outside
	setSphereComputationCosts(pmet_dist, info, center, 2.0f, 5ul, 1ul);

	// first decomposition
	pmet_dist.decompose();

	//! [Initialize a ParMetis Cartesian graph and decompose]

	if (v_cl.getProcessingUnits() == 0)
	{
		// write the first decomposition
		pmet_dist.write("vtk_parmetis_distribution_0.vtk");

		bool test = compare("vtk_parmetis_distribution_0.vtk","src/Decomposition/Distribution/test_data/vtk_parmetis_distribution_0_test.vtk");
		BOOST_REQUIRE_EQUAL(true,test);
	}

	//! [refine with parmetis the decomposition]

	float stime = 0.0, etime = 10.0, tstep = 0.1;

	// Shift of the sphere at each iteration
	Point<3, float> shift( { tstep, tstep, tstep });

	size_t iter = 1;

	for(float t = stime; t < etime; t = t + tstep, iter++)
	{
		if(t < etime/2)
		center += shift;
		else
		center -= shift;

		setSphereComputationCosts(pmet_dist, info, center, 2.0f, 5, 1);

		// With some regularity refine and write the parmetis distribution
		if ((size_t)iter % 10 == 0)
		{
			pmet_dist.refine();

			if (v_cl.getProcessUnitID() == 0)
			{
				std::stringstream str;
				str << "vtk_parmetis_distribution_" << iter;
				pmet_dist.write(str.str() + ".vtk");

				// Check
				bool test = compare(str.str() + ".vtk",std::string("src/Decomposition/Distribution/test_data/") + str.str() + "_test.vtk");
				BOOST_REQUIRE_EQUAL(true,test);
			}
		}
	}

	//! [refine with parmetis the decomposition]

	BOOST_REQUIRE_EQUAL(sizeof(MetisDistribution<3,float>),568ul);
}

BOOST_AUTO_TEST_CASE( DistParmetis_distribution_test)
{
	Vcluster & v_cl = *global_v_cluster;

	if (v_cl.getProcessingUnits() != 3)
	return;

	//! [Initialize a ParMetis Cartesian graph and decompose]

	DistParMetisDistribution<3, float> pmet_dist(v_cl);

	// Physical domain
	Box<3, float> box( { 0.0, 0.0, 0.0 }, { 10.0, 10.0, 10.0 });

	// Grid info
	grid_sm<3, void> info( { GS_SIZE, GS_SIZE, GS_SIZE });

	// Initialize Cart graph and decompose
	pmet_dist.createCartGraph(info,box);

	// First create the center of the weights distribution, check it is coherent to the size of the domain
	Point<3, float> center( { 2.0, 2.0, 2.0 });

	// It produces a sphere of radius 2.0
	// with high computation cost (5) inside the sphere and (1) outside
	setSphereComputationCosts(pmet_dist, info, center, 2.0f, 5ul, 1ul);

	// first decomposition
	pmet_dist.decompose();

	//! [Initialize a ParMetis Cartesian graph and decompose]

	// write the first decomposition
	pmet_dist.write("vtk_dist_parmetis_distribution_0.vtk");

	if (v_cl.getProcessingUnits() == 0)
	{
		bool test = compare("vtk_dist_parmetis_distribution_0.vtk","src/Decomposition/Distribution/test_data/vtk_dist_parmetis_distribution_0_test.vtk");
		BOOST_REQUIRE_EQUAL(true,test);
	}

	//! [refine with dist_parmetis the decomposition]

	float stime = 0.0, etime = 10.0, tstep = 0.1;

	// Shift of the sphere at each iteration
	Point<3, float> shift( { tstep, tstep, tstep });

	size_t iter = 1;

	for(float t = stime; t < etime; t = t + tstep, iter++)
	{
		if(t < etime/2)
		center += shift;
		else
		center -= shift;

		setSphereComputationCosts(pmet_dist, info, center, 2.0f, 5, 1);

		// With some regularity refine and write the parmetis distribution
		if ((size_t)iter % 10 == 0)
		{
			pmet_dist.refine();

			std::stringstream str;
			str << "vtk_dist_parmetis_distribution_" << iter;
			pmet_dist.write(str.str() + ".vtk");

			// Check
			if (v_cl.getProcessUnitID() == 0)
			{
				bool test = compare(str.str() + ".vtk",std::string("src/Decomposition/Distribution/test_data/") + str.str() + "_test.vtk");
				BOOST_REQUIRE_EQUAL(true,test);
			}
		}
	}

	//! [refine with dist_parmetis the decomposition]

	BOOST_REQUIRE_EQUAL(sizeof(DistParMetisDistribution<3,float>),1440ul);
}

void print_test_v(std::string test, size_t sz)
{
	if (global_v_cluster->getProcessUnitID() == 0)
		std::cout << test << " " << sz << "\n";
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
	Ghost<2,float> ghost(r_cut);

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

		auto NN = vd.getCellList(0.01  / factor);

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

#endif /* SRC_DECOMPOSITION_DISTRIBUTION_DISTRIBUTION_UNIT_TESTS_HPP_ */
