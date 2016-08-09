
#include "Vector/vector_dist.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "data_type/aggregate.hpp"
#include "NN/CellList/CellListM.hpp"

/*!
 * \page Vector_4_mp_cl Vector 4 Multi Phase cell-list
 *
 * [TOC]
 *
 *
 * # Vector Multi Phase cell-list # {#e4_ph_cl}
 *
 * This example show multi-phase cell lists for the distributed vector
 *
 * \warning BETA version
 *
 */

int main(int argc, char* argv[])
{
	/*!
	 * \page Vector_4_mp_cl Vector 4 Multi Phase cell-list
	 *
	 * ## Initialization ##
	 *
	 * Here we Initialize the library, and we create a set of distributed vectors all forced to have the same
	 * decomposition. Each vector identify one phase
	 *
	 * \snippet Vector/1_celllist/main.cpp Initialization and parameters
	 *
	 */

	//! \cond [Initialization and parameters] \endcond

	openfpm_init(&argc,&argv);
	Vcluster & v_cl = create_vcluster();

	// we will place the particles on a grid like way with 128 particles on each direction
	size_t sz[3] = {128,128,128};

	// The domain
	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// ghost, big enough to contain the interaction radius
	Ghost<3,float> ghost(1.0/(128-2));

	//! \cond [Initialization and parameters] \endcond

	/*!
	 * \page Vector_1_celllist Vector 1 Cell-list
	 *
	 * ## %Vector create ##
	 *
	 * Here we define a distributed vector in 3D, containing 3 properties, a
	 * scalar double, a vector double[3], and a tensor or rank 2 double[3][3].
	 * In this case the vector contain 0 particles initially
	 *
	 * \see \ref vector_inst
	 *
	 * \snippet Vector/1_celllist/main.cpp vector inst
	 *
	 */

	//! \cond [vector inst] \endcond

	openfpm::vector< vector_dist<3,float, aggregate<double,double>> > phases;

	// first phase
	phases.add( vector_dist<3,float, aggregate<double,double>>(4096,box,bc,ghost) );

	// The other 3 phases
	phases.add( vector_dist<3,float, aggregate<double,double>>(phases.get(1).getDecomposition(),4096) );
	phases.add( vector_dist<3,float, aggregate<double,double>>(phases.get(2).getDecomposition(),4096) );
	phases.add( vector_dist<3,float, aggregate<double,double>>(phases.get(3).getDecomposition(),4096) );

	//! \cond [grid like part] \endcond

	auto it = phases.get(0).getDomainIterator();

	// For all the particles
	while (it.isNext())
	{
		// for all phases
		for (size_t i = 0 ; i < phases.size() ; i++)
		{
			auto key = it.get();

			phases.get(i).getPos(key)[0] = (float)rand() / RAND_MAX;
			phases.get(i).getPos(key)[1] = (float)rand() / RAND_MAX;
			phases.get(i).getPos(key)[2] = (float)rand() / RAND_MAX;
		}
		// next point
		++it;
	}

	// Redistribute all the phases in the mean while also take the iterators of all the phases

	typedef decltype(phases.get(0).getIterator()) iterator;
	openfpm::vector<iterator> phase_it;

	for (size_t i = 0 ; i < phases.size() ; i++)
	{
		phases.get(i).map();
		phase_it.add(phases.get(i).getDomainIterator());
	}

	// Construct one single Multi-phase cell list to use in the computation

	CellListM<3,float,2> NN;

	while (it.isNext())
	{
		for (size_t i = 0; i < phases.size() ; i++)
		{
			auto key = it.get();

			NN.add(phases.get(i).getPos(key), key.getKey(), i);

			++it;
		}
	}

	vector_dist<3,float, aggregate<double,double> > & current_phase = phases.get(0);

	// Get the iterator of the particles of phase 0
	auto it2 = current_phase.getIterator();

	// For each particle ...
	while (it2.isNext())
	{
		// ... p
		auto p = it2.get();

		// Get the position of the particle p
		Point<3,float> xp = current_phase.getPos(p);

		// Get an iterator of all the particles neighborhood of p
		auto Np = NN.getNNIterator(NN.getCell(current_phase.getPos(p)));

		// For each particle near p
		while (Np.isNext())
		{
			// Get the particle q near to p
			auto q = Np.getP();

			// Get from which phase it come from
			auto ph_q = Np.getV();

			Point<3,float> xq = phases.get(ph_q).getPos(q);

			// we accumulate all the distances
			current_phase.getProp<0>(p) = norm(xp - xq);

			++Np;
		}

		// Next particle p
		++it2;
	}

	//! \cond [verletlist] \endcond

	/*!
	 * \page Vector_1_celllist Vector 1 Cell-list
	 *
	 * ## Finalize ## {#finalize}
	 *
	 *  At the very end of the program we have always to de-initialize the library
	 *
	 * \snippet Vector/1_celllist/main.cpp finalize
	 *
	 */

	//! \cond [finalize] \endcond

	openfpm_finalize();

	//! \cond [finalize] \endcond

	/*!
	 * \page Vector_1_celllist Vector 1 Cell-list
	 *
	 * # Full code # {#code}
	 *
	 * \include Vector/1_celllist/main.cpp
	 *
	 */
}




