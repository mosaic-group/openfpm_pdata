#include "Vector/vector_dist.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "data_type/aggregate.hpp"
#include "Plot/GoogleChart.hpp"
#include "Plot/util.hpp"
#include "timer.hpp"

/*!
 * \page Vector_3_md_vl Vector 3 molecular dynamic with Verlet list
 *
 *
 * ## Variations for symmetric interaction ##
 *
 * In this example we show the variation in case of symmetric interations.
 *
 */

//! \cond [constants] \endcond

constexpr int velocity = 0;
constexpr int force = 1;

//! \cond [constants] \endcond

/*!
 *
 * \page Vector_3_md_vl Vector 3 molecular dynamic with Verlet list symmetric
 *
 * ## Calculate forces ##
 *
 * In this function we calculate the forces between particles. It require the vector of particles
 * Cell list and sigma factor for the Lennard-Jhones potential. The function is exactly the same
 * as the original with the following changes
 *
 *  \see \ref e3_md_vl_cf
 *
 * * We have to reset for force counter for each particles
 *  \snippet Vector/3_molecular_dynamic/main_vl_sym.cpp reset
 *
 * * The iterator go through real AND ghost.
 *  \snippet Vector/3_molecular_dynamic/main_vl_sym.cpp real and ghost
 *
 * * If we calculate the force for p-q we are also adding this force to q-p
 *  \snippet Vector/3_molecular_dynamic/main_vl_sym.cpp add to q
 *
 * ** Explanation**
 *
 * The reason why symmetric Verlet exist is to reduce the computation. In a symmetric interaction
 * we have that the interaction p-q is is the same as q-p. This mean that we are calculating 2 times
 * the same quantity one when we are on particle p and one when we are on particle q. Symmetric verlet
 *  avoid this. Every time we calculate a force p-q we also increment the force q-p.
 *
 * The reason instead why we have to go through ghost is the following
 *
 *
	\verbatim

	|      |
	| Ghost|
	|      |
	|      |
	|      |
	|      |
	|      |
	|   +----+----+----+
	|   |  | |    |    |
	|   | 1| | 2  | 3  |
	|   |  | |    |    |
	|   +--------------+
	|      | |  * | 4  |
	|      +--------------------
	|    X   | 5  |    |
	|        +---------+
	+---------------------------

	\endverbatim
 *
 * Where * is a particle and X is another particle on ghost
 *
 * The symmetric verlet list is constructed from the cell-List 1,2,3,4,5.
 * As you can see the interaction X-* is not present in the neighborhood of *
 * (cells 1,2,3,4,5). But is present if we draw the neighborhood of X. This
 * mean we have to iterate also across X if we want the interation X-*
 *
 *
 * \snippet Vector/3_molecular_dynamic/main_vl_sym.cpp calc forces vl
 *
 * The symmetric force calculation is faster than the normal one. but is not 2
 * time faster. The reason can be founded in 2 factors.
 *
 * * We have to reset the force
 *
 * * Every particle in the ghost
 * that interact with one in the domain is still calculated twice by the same
 * processor or some other.
 *
 * * The access pattern is worst we alternate read and write stressing more the Memory
 *
 * * We have to iterate across domain + ghost
 *
 * Overall the function calc_force is faster
 *
 * \snippet Vector/3_molecular_dynamic/main_vl_sym.cpp calc forces vl
 *
 */

//! \cond [calc forces vl] \endcond

void calc_forces(vector_dist<3,double, aggregate<double[3],double[3]> > & vd, VerletList<3, double, FAST, shift<3, double> > & NN, double sigma12, double sigma6, double r_cut)
{
	//! \cond [reset] \endcond

	// Reset the force for all particles
	auto it3 = vd.getDomainIterator();

	// For each particle p ...
	while (it3.isNext())
	{
		// ... get the particle p
		auto p = it3.get();

		// Reset the forice counter
		vd.template getProp<force>(p)[0] = 0.0;
		vd.template getProp<force>(p)[1] = 0.0;
		vd.template getProp<force>(p)[2] = 0.0;

		// Next particle p
		++it3;
	}

	//! \cond [reset] \endcond

	//! \cond [real and ghost] \endcond

	// Get an iterator over particles
	auto it2 = vd.getDomainAndGhostIterator();

	//! \cond [real and ghost] \endcond

	// For each particle p ...
	while (it2.isNext())
	{
		// ... get the particle p
		auto p = it2.get();

		// Get the position xp of the particle
		Point<3,double> xp = vd.getPos(p);

		// Get an iterator over the neighborhood particles of p
		// Note that in case of symmetric
		auto Np = NN.template getNNIterator<NO_CHECK>(p.getKey());

		// For each neighborhood particle ...
		while (Np.isNext())
		{
			// ... q
			auto q = Np.get();

			// if (p == q) skip this particle
			if (q == p.getKey() || (p.getKey() >= vd.size_local() && q >= vd.size_local()))	{++Np; continue;};

			// Get the position of p
			Point<3,double> xq = vd.getPos(q);

			// Get the distance between p and q
			Point<3,double> r = xp - xq;

			// take the norm of this vector
			double rn = norm2(r);

			if (rn > r_cut * r_cut) {++Np;continue;}

			// Calculate the force, using pow is slower
			Point<3,double> f = 24.0*(2.0 *sigma12 / (rn*rn*rn*rn*rn*rn*rn) -  sigma6 / (rn*rn*rn*rn)) * r;

			// we sum the force produced by q on p
			vd.template getProp<force>(p)[0] += f.get(0);
			vd.template getProp<force>(p)[1] += f.get(1);
			vd.template getProp<force>(p)[2] += f.get(2);

			//! \cond [add to q] \endcond

			// we sum the force produced by p on q
			vd.template getProp<force>(q)[0] -= f.get(0);
			vd.template getProp<force>(q)[1] -= f.get(1);
			vd.template getProp<force>(q)[2] -= f.get(2);

			//! \cond [add to q] \endcond

			// Next neighborhood
			++Np;
		}

		// Next particle
		++it2;
	}
}

//! \cond [calc forces vl] \endcond


/*!
 * \page Vector_3_md_vl Vector 3 molecular dynamic with Verlet list symmetric
 *
 * ## Calculate energy ##
 *
 * No changes that been made to this code compared to the molecular dynamic example
 * with Verlet-list
 *
 * \see \ref e3_md_vl
 *
 *
 * \snippet Vector/3_molecular_dynamic/main_vl_sym.cpp calc energy vl
 *
 */

//! \cond [calc energy vl] \endcond

double calc_energy(vector_dist<3,double, aggregate<double[3],double[3]> > & vd, VerletList<3, double, FAST, shift<3, double> > & NN, double sigma12, double sigma6, double r_cut)
{
	double E = 0.0;

	double rc = r_cut*r_cut;
	double shift = 2.0 * ( sigma12 / (rc*rc*rc*rc*rc*rc) - sigma6 / ( rc*rc*rc) );

	// Get the iterator
	auto it2 = vd.getDomainIterator();

	// For each particle ...
	while (it2.isNext())
	{
		// ... p
		auto p = it2.get();

		// Get the position of the particle p
		Point<3,double> xp = vd.getPos(p);

		// Get an iterator over the neighborhood of the particle p
		auto Np = NN.template getNNIterator<NO_CHECK>(p.getKey());

		double Ep = E;

		// For each neighborhood of the particle p
		while (Np.isNext())
		{
			// Neighborhood particle q
			auto q = Np.get();

			// if p == q skip this particle
			if (q == p.getKey() && q < vd.size_local())	{++Np; continue;};

			// Get position of the particle q
			Point<3,double> xq = vd.getPos(q);

			// take the normalized direction
			double rn = norm2(xp - xq);

			if (rn >= r_cut*r_cut)
			{++Np;continue;}

			// potential energy (using pow is slower)
			E += 2.0 * ( sigma12 / (rn*rn*rn*rn*rn*rn) - sigma6 / ( rn*rn*rn) ) - shift;

			// Next neighborhood
			++Np;
		}

		// Kinetic energy of the particle given by its actual speed
		E +=   (vd.template getProp<velocity>(p)[0]*vd.template getProp<velocity>(p)[0] +
				vd.template getProp<velocity>(p)[1]*vd.template getProp<velocity>(p)[1] +
				vd.template getProp<velocity>(p)[2]*vd.template getProp<velocity>(p)[2]) / 2;

		// Next Particle
		++it2;
	}

	// Calculated energy
	return E;
}

//! \cond [calc energy vl] \endcond

int main(int argc, char* argv[])
{
	/*!
	 * \page Vector_3_md_vl Vector 3 molecular dynamic with Verlet list symmetric
	 *
	 * ## Simulation ##
	 *
	 * The simulation is equal to the simulation explained in the example molecular dynamic
	 *
	 * \see \ref e3_md_vl
	 *
	 * The differences are that:
	 *
	 * * Create symmetric Verlet-list
	 *   \snippet Vector/3_molecular_dynamic/main_vl_sym.cpp sim verlet
	 *
	 * * For Energy calculation we have to create another verlet in this case a normal Verlet-list
	 *  \snippet Vector/3_molecular_dynamic/main_vl_sym.cpp norm verlet
	 *
	 * We need 2 verlet-list with one normal because we cannot use
	 * symmetric Verlet to calculate total reductions like total energy
	 *
	 * \snippet Vector/3_molecular_dynamic/main_vl_sym.cpp simulation
	 *
	 */

	//! \cond [simulation] \endcond

	double dt = 0.00025;
	double sigma = 0.1;
	double r_cut = 3.0*sigma;
	double r_gskin = 1.3*r_cut;
	double sigma12 = pow(sigma,12);
	double sigma6 = pow(sigma,6);

	openfpm::vector<double> x;
	openfpm::vector<openfpm::vector<double>> y;

	openfpm_init(&argc,&argv);
	Vcluster & v_cl = create_vcluster();

	// we will use it do place particles on a 10x10x10 Grid like
	size_t sz[3] = {10,10,10};

	// domain
	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// ghost, big enough to contain the interaction radius
	Ghost<3,float> ghost(r_gskin);

	vector_dist<3,double, aggregate<double[3],double[3]> > vd(0,box,bc,ghost);

	auto it = vd.getGridIterator(sz);

	while (it.isNext())
	{
		vd.add();

		auto key = it.get();

		vd.getLastPos()[0] = key.get(0) * it.getSpacing(0);
		vd.getLastPos()[1] = key.get(1) * it.getSpacing(1);
		vd.getLastPos()[2] = key.get(2) * it.getSpacing(2);

		vd.template getLastProp<velocity>()[0] = 0.0;
		vd.template getLastProp<velocity>()[1] = 0.0;
		vd.template getLastProp<velocity>()[2] = 0.0;

		vd.template getLastProp<force>()[0] = 0.0;
		vd.template getLastProp<force>()[1] = 0.0;
		vd.template getLastProp<force>()[2] = 0.0;

		++it;
	}

	timer tsim;
	tsim.start();

	//! \cond [sim verlet] \endcond

	// Get the Cell list structure
	auto NN = vd.getVerlet(r_gskin,VL_SYMMETRIC);

	//! \cond [sim verlet] \endcond

	//! \cond [norm verlet] \endcond

	auto NNE = vd.getVerlet(r_gskin);

	//! \cond [norm verlet] \endcond

	// calculate forces
	calc_forces(vd,NN,sigma12,sigma6,r_cut);
	unsigned long int f = 0;

	int cnt = 0;
	double max_disp = 0.0;

	// MD time stepping
	for (size_t i = 0; i < 10000 ; i++)
	{
		// Get the iterator
		auto it3 = vd.getDomainIterator();

		double max_displ = 0.0;

		// integrate velicity and space based on the calculated forces (Step1)
		while (it3.isNext())
		{
			auto p = it3.get();

			// here we calculate v(tn + 0.5)
			vd.template getProp<velocity>(p)[0] += 0.5*dt*vd.template getProp<force>(p)[0];
			vd.template getProp<velocity>(p)[1] += 0.5*dt*vd.template getProp<force>(p)[1];
			vd.template getProp<velocity>(p)[2] += 0.5*dt*vd.template getProp<force>(p)[2];

			Point<3,double> disp({vd.template getProp<velocity>(p)[0]*dt,vd.template getProp<velocity>(p)[1]*dt,vd.template getProp<velocity>(p)[2]*dt});

			// here we calculate x(tn + 1)
			vd.getPos(p)[0] += disp.get(0);
			vd.getPos(p)[1] += disp.get(1);
			vd.getPos(p)[2] += disp.get(2);

			if (disp.norm() > max_displ)
				max_displ = disp.norm();

			++it3;
		}

		if (max_disp < max_displ)
			max_disp = max_displ;

		// Because we moved the particles in space we have to map them and re-sync the ghost
		if (cnt % 10 == 0)
		{
			vd.map();
			vd.template ghost_get<>();
			// Get the Cell list structure
			vd.updateVerlet(NN,r_gskin,VL_SYMMETRIC);
		}
		else
		{
			vd.template ghost_get<>(SKIP_LABELLING);
		}

		cnt++;

		// calculate forces or a(tn + 1) Step 2
		calc_forces(vd,NN,sigma12,sigma6,r_cut);

		// Integrate the velocity Step 3
		auto it4 = vd.getDomainIterator();

		while (it4.isNext())
		{
			auto p = it4.get();

			// here we calculate v(tn + 1)
			vd.template getProp<velocity>(p)[0] += 0.5*dt*vd.template getProp<force>(p)[0];
			vd.template getProp<velocity>(p)[1] += 0.5*dt*vd.template getProp<force>(p)[1];
			vd.template getProp<velocity>(p)[2] += 0.5*dt*vd.template getProp<force>(p)[2];

			++it4;
		}

		// After every iteration collect some statistic about the confoguration
		if (i % 100 == 0)
		{
			// We write the particle position for visualization (Without ghost)
			vd.deleteGhost();
			vd.write("particles_",f);

			// we resync the ghost
			vd.ghost_get<>(SKIP_LABELLING);

			// update the verlet for energy calculation
			vd.updateVerlet(NNE,r_gskin);

			// We calculate the energy
			double energy = calc_energy(vd,NNE,sigma12,sigma6,r_cut);
			auto & vcl = create_vcluster();
			vcl.sum(energy);
			vcl.max(max_disp);
			vcl.execute();

			// we save the energy calculated at time step i c contain the time-step y contain the energy
			x.add(i);
			y.add({energy});

			// We also print on terminal the value of the energy
			// only one processor (master) write on terminal
			if (vcl.getProcessUnitID() == 0)
				std::cout << "Energy: " << energy << "   " << max_disp << "   " << std::endl;

			max_disp = 0.0;

			f++;
		}
	}

	tsim.stop();
	std::cout << "Time: " << tsim.getwct()  << std::endl;

	//! \cond [simulation] \endcond

	// Google charts options, it store the options to draw the X Y graph
	GCoptions options;

	// Title of the graph
	options.title = std::string("Energy with time");

	// Y axis name
	options.yAxis = std::string("Energy");

	// X axis name
	options.xAxis = std::string("iteration");

	// width of the line
	options.lineWidth = 1.0;

	// Object that draw the X Y graph
	GoogleChart cg;

	// Add the graph
	// The graph that it produce is in svg format that can be opened on browser
	cg.AddLinesGraph(x,y,options);

	// Write into html format
	cg.write("gc_plot2_out.html");

	//! \cond [google chart] \endcond

	/*!
	 * \page Vector_3_md_vl Vector 3 molecular dynamic with Verlet list symmetric
	 *
	 * ## Finalize ## {#finalize_v_e3_md_sym}
	 *
	 *  At the very end of the program we have always to de-initialize the library
	 *
	 * \snippet Vector/3_molecular_dynamic/main_vl_sym.cpp finalize
	 *
	 */

	//! \cond [finalize] \endcond

	openfpm_finalize();

	//! \cond [finalize] \endcond

	/*!
	 * \page Vector_3_md_vl Vector 3 molecular dynamic
	 *
	 * ## Full code ## {#code_v_e3_md_vl}
	 *
	 * \include Vector/3_molecular_dynamic/main_vl.cpp
	 *
	 */

	/*!
	 * \page Vector_3_md_vl Vector 3 molecular dynamic with Verlet list symmetric
	 *
	 * ## Full code symmetric ## {#code_v_e3_md_sym}
	 *
	 * \include Vector/3_molecular_dynamic/main_vl_sym.cpp
	 *
	 */
}
