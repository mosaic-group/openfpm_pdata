
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
	 * \snippet Vector/4_multiphase_celllist/main.cpp Initialization and parameters
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

	float r_cut = 0.05;

	// ghost, big enough to contain the interaction radius
	Ghost<3,float> ghost(r_cut);

	openfpm::vector< vector_dist<3,float, aggregate<double,double>> > phases;

	// first phase
	phases.add( vector_dist<3,float, aggregate<double,double>>(4096,box,bc,ghost) );

	// The other 3 phases
	phases.add( vector_dist<3,float, aggregate<double,double>>(phases.get(0).getDecomposition(),4096) );
	phases.add( vector_dist<3,float, aggregate<double,double>>(phases.get(0).getDecomposition(),4096) );
	phases.add( vector_dist<3,float, aggregate<double,double>>(phases.get(0).getDecomposition(),4096) );

	//! \cond [Initialization and parameters] \endcond


	/*!
	 * \page Vector_4_mp_cl Vector 4 Multi Phase cell-list
	 *
	 * ## Initialization ##
	 *
	 * We initialize all the phases with particle randomly positioned in the space
	 *
	 * \snippet Vector/4_multiphase_celllist/main.cpp rand dist
	 *
	 */

	//! \cond [rand dist] \endcond

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
		phases.get(i).map();

	//! \cond [rand dist] \endcond

	/*!
	 * \page Vector_4_mp_cl Vector 4 Multi Phase cell-list
	 *
	 * ## Multi-phase cell-list construction ##
	 *
	 * In this part we construct the Multi-phase cell list. The multiphase cell list has 3 parameters
	 * * one is the dimensionality (3)
	 * * The precision of the coordinates (float),
	 * * How many bit to use for the phase.
	 *   Multi-phase cell-list try to pack into a 64bit number information about the particle id and the
	 *   the phase id.
	 *
	 * \snippet Vector/4_multiphase_celllist/main.cpp cl construction
	 *
	 */

	//! \cond [cl construction] \endcond

	// Construct one single Multi-phase cell list to use in the computation
	// in 3d, precision float, 2 bit dedicated to the phase for a maximum of 2^2 = 4 (Maximum number of phase)
	//
	//
	
	size_t div[3];
	Box<3,float> box_cl;
	phases.get(0).getCellListParams(r_cut,div,box_cl);

	CellListM<3,float,2> NN;
	NN.Initialize(box_cl,div);

	// for all the phases i
	for (size_t i = 0; i < phases.size() ; i++)
	{
		// iterate across all the particle of the phase i
		auto it = phases.get(i).getDomainIterator();
		while (it.isNext())
		{
			auto key = it.get();

			// Add the particle of the phase i to the cell list
			NN.add(phases.get(i).getPos(key), key.getKey(), i);

			++it;
		}
	}

	//! \cond [cl construction] \endcond

	/*!
	 * \page Vector_4_mp_cl Vector 4 Multi Phase cell-list
	 *
	 * ## Multi-phase cell-list usage ##
	 *
	 * After construction we show how to use the Cell-list. In this case we accumulate on the property
	 * 0 of the phase 0 the distance of the near particles from all the phases
	 *
	 * \snippet Vector/4_multiphase_celllist/main.cpp cl usage
	 *
	 */

	//! \cond [cl usage] \endcond

/*	vector_dist<3,float, aggregate<double,double> > & current_phase = phases.get(0);

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
	}*/

	//! \cond [cl usage] \endcond

	/*!
	 * \page Vector_4_mp_cl Vector 4 Multi Phase cell-list
	 *
	 * ## Finalize ## {#finalize}
	 *
	 *  At the very end of the program we have always de-initialize the library
	 *
	 * \snippet Vector/4_multiphase_celllist/main.cpp finalize
	 *
	 */

	//! \cond [finalize] \endcond

	openfpm_finalize();

	//! \cond [finalize] \endcond

	/*!
	 * \page Vector_4_mp_cl Vector 4 Multi Phase cell-list
	 *
	 * # Full code # {#code}
	 *
	 * \include Vector/4_multiphase_celllist/main.cpp
	 *
	 */
}




