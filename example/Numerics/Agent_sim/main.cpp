/*
 * ### WIKI 1 ###
 *
 * ## Simple example
 *
 * In this example we show 1D PSE derivative function approximation
 *
 * ### WIKI END ###
 *
 */

#define SE_CLASS1
#define THROW_ON_ERROR

#include "Vector/vector_dist.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "PSE/Kernels.hpp"
#include "Plot/util.hpp"
#include "Plot/GoogleChart.hpp"
#include "data_type/aggregate.hpp"
#include <cmath>

struct animal
{
	typedef boost::fusion::vector<float[2], size_t, size_t, long int> type;

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

const std::string animal::attributes::name[] = { "pos", "genre", "status", "time_a", "j_repr" };

int main(int argc, char* argv[])
{
	init_global_v_cluster(&argc,&argv);

	Vcluster & v_cl = *global_v_cluster;

	//time the animal stays alive without eating
	size_t PRED_TIME_A = 14;

	size_t PREY_TIME_A = 7;

	size_t PREDATOR = 1, PREY = 0;
	size_t ALIVE = 1, DEAD = 0;

	// Predators reproducing probability
	float PRED_REPR = 0.2;

	// Predators eating probability
	float PRED_EAT = 0.6;

	// Prey reproducing probability
	float PREY_REPR = 0.5;

	// set the seed
	// create the random generator engine
	std::srand(v_cl.getProcessUnitID());
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine eg(seed);
	std::uniform_real_distribution<float> ud(0.0f, 1.0f);
	std::uniform_real_distribution<float> md(-1.0f, 1.0f);
	std::uniform_real_distribution<float> uc(0.0f, 0.7f);
	std::uniform_real_distribution<float> lc(0.3f, 1.0f);

	size_t k = 100000;

	Box<2, float> box( { 0.0, 0.0 }, { 1.0, 1.0 });

	// Grid info
	grid_sm<2, void> info( { 8, 8 });

	// Boundary conditions
	size_t bc[2] = { PERIODIC, PERIODIC };

	// factor
	float factor = pow(global_v_cluster->getProcessingUnits() / 2.0f, 1.0f / 3.0f);

	// interaction radius
	float r_cut = 0.01 / factor;

	// ghost
	Ghost<2, float> ghost(r_cut);

	// Distributed vector
	vector_dist<2, float, animal, CartDecomposition<2, float, HeapMemory, ParMetisDistribution<2, float>>> vd(k,box,bc,ghost);

	// Init DLB tool
	DLB dlb(v_cl);

	// Set unbalance threshold
	dlb.setHeurisitc(DLB::Heuristic::UNBALANCE_THRLD);
	dlb.setThresholdLevel(DLB::ThresholdLevel::THRLD_MEDIUM);

	auto it = vd.getIterator();

	while (it.isNext())
	{
		auto key = it.get();
		if(ud(eg) < 0.7 )
		{
			vd.template getPos<animal::pos>(key)[0] = lc(eg);
			vd.template getPos<animal::pos>(key)[1] = lc(eg);
			vd.template getProp<animal::genre>(key) = PREY;
			vd.template getProp<animal::status>(key) = ALIVE;
			vd.template getProp<animal::time_a>(key) = PREY_TIME_A;
		}
		else
		{
			vd.template getPos<animal::pos>(key)[0] = uc(eg);
			vd.template getPos<animal::pos>(key)[1] = uc(eg);
			vd.template getProp<animal::genre>(key) = PREDATOR;
			vd.template getProp<animal::status>(key) = ALIVE;
			vd.template getProp<animal::time_a>(key) = PRED_TIME_A;
		}
		++it;
	}

	vd.map();

	vd.addComputationCosts();

	vd.getDecomposition().rebalance(dlb);

	vd.map();

	//vd.getDecomposition().getDistribution().write("parmetis_prey_predators_" + std::to_string(0) + ".vtk");
	//vd.write("particles_", 0, NO_GHOST);

	// 100 step random walk
	for (size_t j = 0; j < 100; j++)
	{
		std::cout << "iteration: " << j << std::endl;

		size_t prey = 0, predators = 0;

		auto it = vd.getDomainIterator();

		while (it.isNext())
		{
			auto key = it.get();

			vd.template getPos<animal::pos>(key)[0] += 0.01 * md(eg);
			vd.template getPos<animal::pos>(key)[1] += 0.01 * md(eg);

			if(vd.template getProp<animal::genre>(key) == PREY)
			prey++;
			else
			predators++;

			++it;
		}

		vd.map();

		float tot_p = vd.size_local();

		v_cl.sum(tot_p);
		v_cl.sum(prey);
		v_cl.sum(predators);
		v_cl.execute();


	//	 if(v_cl.getProcessUnitID()==0)
	//	 {
	//		 std::cout << j << " -> tot: " << tot_p << " prey: " << prey << " predators: " << predators << "\n";
	//	 }


		/////// Interactions ///
		// get ghosts
		try
		{
			vd.ghost_get<0>();
		}
		catch (size_t e)
		{
			vd.write("err_particles",j);
			vd.getDecomposition().write("err_decomposition");
			MPI_Abort(MPI_COMM_WORLD,-1);
			return -1;
		}

		// vector of dead animals
		openfpm::vector<size_t> deads;
		openfpm::vector<vect_dist_key_dx> reps_prey;
		openfpm::vector<vect_dist_key_dx> reps_pred;

		// get the cell list with a cutoff radius

		bool error = false;

		CellList<2, float, FAST, shift<2, float> > NN;


			NN = vd.getCellList(0.01/factor);


		// iterate across the domain particle

		auto it2 = vd.getDomainIterator();

		while (it2.isNext())
		{
			auto p = it2.get();

			Point<2,float> xp = vd.getPos<0>(p);

			size_t gp = vd.getProp<animal::genre>(p);
			size_t sp = vd.getProp<animal::status>(p);

			if(sp == ALIVE)
			{
				if(gp == PREY)
				{
					if( prey < k/1.5 && ud(eg) < PREY_REPR )
						reps_prey.add(p);

					vd.getProp<animal::time_a>(p)--;

					if(vd.getProp<animal::time_a>(p) <= 0)
					{
						vd.getProp<animal::status>(p) = DEAD;
						prey--;
					}
				}
				else if(gp == PREDATOR)
				{
					vd.getProp<animal::time_a>(p)--;

					if(vd.getProp<animal::time_a>(p) <= 0)
					{
						vd.getProp<animal::status>(p) = DEAD;
					}
					else
					{
						auto Np = NN.getIterator(NN.getCell(xp));

						while (Np.isNext())
						{
							auto q = Np.get();

							size_t gq = vd.getProp<animal::genre>(q);
							size_t sq = vd.getProp<animal::status>(q);

							Point<2,float> xq = vd.getPos<0>(q);
							Point<2,float> f = (xp - xq);

							float distance = f.norm();

							if (distance < 2*r_cut*sqrt(2) && gq == PREY && sq == ALIVE)
							{
								if( ud(eg) < PRED_EAT )
								{
									vd.getProp<animal::status>(q) = DEAD;
									vd.getProp<animal::time_a>(p) = PRED_TIME_A;

									if( ud(eg) < PRED_REPR )
										reps_pred.add(p);
								}
							}
							++Np;
						}
					}
				}

			}

			++it2;
		}

		vd.deleteGhost();

		// Replicate

		for (size_t i = 0 ; i < reps_prey.size() ; i++)
		{
			vd.add();
			vd.getLastPos<animal::pos>()[0] = vd.getPos<0>(reps_prey.get(i))[0];
			vd.getLastPos<animal::pos>()[1] = vd.getPos<0>(reps_prey.get(i))[1];
			vd.getLastProp<animal::genre>() = PREY;
			vd.getLastProp<animal::status>() = ALIVE;
			vd.getLastProp<animal::time_a>() = PREY_TIME_A;
		}

		for (size_t i = 0 ; i < reps_pred.size() ; i++)
		{
			vd.add();
			vd.getLastPos<animal::pos>()[0] = vd.getPos<0>(reps_pred.get(i))[0];
			vd.getLastPos<animal::pos>()[1] = vd.getPos<0>(reps_pred.get(i))[1];
			vd.getLastProp<animal::genre>() = PREDATOR;
			vd.getLastProp<animal::status>() = ALIVE;
			vd.getLastProp<animal::time_a>() = PRED_TIME_A;
		}

		auto it3 = vd.getDomainIterator();
		while (it3.isNext())
		{
			auto key = it3.get();
			if(vd.getProp<animal::status>(key.getKey()) == DEAD)
			{
				deads.add(key.getKey());
			}
			++it3;
		}

		deads.sort();
		vd.remove(deads, 0);

/*		for (size_t i = 0 ; i < deads.size() ; i++)
		{
			size_t d = deads.get(i);

			vd.remove(d);
		}*/

		deads.resize(0);

		vd.deleteGhost();

		////////////////////////

		vd.map();

		vd.addComputationCosts();

//		vd.write("before_ref_particles" + std::to_string(j),j);
//		vd.getDecomposition().write("before_ref_decomposition" + std::to_string(j));

		vd.getDecomposition().rebalance(1);

		vd.map();

//		vd.write("after_ref_particles" + std::to_string(j),j);
//		vd.getDecomposition().write("after_ref_decomposition" + std::to_string(j));

		//vd.getDecomposition().getDistribution().write("parmetis_prey_predators_" + std::to_string(j+1) + ".vtk");
		//vd.write("particles_", j+1, NO_GHOST);
		//vd.getDecomposition().write("dec_");

		size_t l = vd.size_local();
		v_cl.sum(l);
		v_cl.execute();

		if (v_cl.getProcessUnitID() == 0)
		{
			std::cout << "Particles: " << l  << "\n";
		}

		// Count the local particles and check that the total number is consistent
		//size_t cnt = total_n_part_lc(vd,bc);

		//BOOST_REQUIRE_EQUAL((size_t)k,cnt);
	}

	//
	// ### WIKI 10 ###
	//
	// Deinitialize the library
	//
	delete_global_v_cluster();
}
