#include "Vector/vector_dist.hpp"
#include "data_type/aggregate.hpp"
#include "Plot/GoogleChart.hpp"
#include "Plot/util.hpp"
#include "timer.hpp"
#include "energy_force.hpp"

/*!
 * \page Vector_4_reo Vector 4 data-reordering and cache friendliness
 *
 * [TOC]
 *
 * # Data reordering, computation-reordering and cache friendly computation # {#e4_reo}
 *
 * In this example we show how reordering the data can significantly improve the computation speed.
 * In order to do this we will re-work the molecular dynamic example.
 *
 *
 */


//! \cond [calc energy2] \endcond

int main(int argc, char* argv[])
{
	/*!
	 * \page Vector_4_reo Vector 4 data-reordering and cache friendliness
	 *
	 * ## Initialization ##
	 *
	 * The initialization is the same as the molecular dynamic example. The difference are in the
	 * parameters. We will use a bigger system, with more particles. The delta time for integration
	 * is chosen in order to keep the system stable.
	 *
	 * \see \ref e3_md_init
	 *
	 * \snippet Vector/4_reorder/main_data_ord.cpp vect create
	 *
	 */

	//! \cond [vect create] \endcond

	double dt = 0.0001;
	float r_cut = 0.03;
	double sigma = r_cut/3.0;
	double sigma12 = pow(sigma,12);
	double sigma6 = pow(sigma,6);

	openfpm::vector<double> x;
	openfpm::vector<openfpm::vector<double>> y;

	openfpm_init(&argc,&argv);
	Vcluster<> & v_cl = create_vcluster();

	// we will use it do place particles on a 10x10x10 Grid like
	size_t sz[3] = {40,40,40};

	// domain
	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// ghost, big enough to contain the interaction radius
	Ghost<3,float> ghost(r_cut);

	// Create vector
	vector_dist<3,double, aggregate<double[3],double[3]> > vd(0,box,bc,ghost);

	//! \cond [vect create] \endcond

	/*!
	 * \page Vector_4_reo Vector 4 data-reordering and cache friendliness
	 *
	 * ## Particles on a grid like position ##
	 *
	 * Here we place the particles on a grid like manner
	 *
	 * \see \ref e3_md_gl
	 *
	 * \snippet Vector/4_reorder/main_data_ord.cpp vect grid
	 *
	 */

	//! \cond [vect grid] \endcond

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

	//! \cond [vect grid] \endcond

	/*!
	 * \page Vector_4_reo Vector 4 data-reordering and cache friendliness
	 *
	 * ## Molecular dynamic steps ##
	 *
	 * Here we do 30000 MD steps using verlet integrator the cycle is the same as the
	 * molecular dynamic example. with the following changes.
	 *
	 * ### Reordering ###
	 *
	 * Every 200 iterations we reorder the data.
	 *
	 * \snippet Vector/4_reorder/main_data_ord.cpp vector reorder
	 *
	 * This function reorder the particles internally using an hilbert curve of order m
	 * (in this case m=5).
	 * More in detail, a cell list is created with \f$ 2^m \times 2^m \f$ Cells in 2D or
	 * \f$ 2^m \times 2^m \times 2^m \f$ Cells in 3D. The particles are reordered in the
	 * vector following an Hilbert curve passing through the cells, in the way described
	 * in the figure below
	 *
	 * \verbatim

	+------+------+------+------+     Example of Hilbert curve for m = 2
	|      |      |      |      |
	|  6+---->7   |  10+--->11  |
	|  ^   |  +   |  ^   |   +  |
	+--|------|------|-------|--+
	|  +   |  v   |  +   |   v  |
	|  5   |  8+---->9   |  12  |
	|  ^   |      |      |   +  |
	+--|---------------------|--+
	|  +   |      |      |   v  |
	|  4<----+3   |  14<---+13  |
	|      |  ^   |   +  |      |
	+---------|-------|---------+
	|      |  +   |   v  |      |
	|  1+---->2   |  15+--->16  |
	|      |      |      |      |
	+------+------+------+------+


	 Suppose now that the particles are ordered like the situation (BEFORE).
	 After the call to reorder they will be reorder like (AFTER)


                   BEFORE                                AFTER

    Particles   id      Cell                         id        Cell
                 0         1                          0           1
                 1         7                          3           1
                 2         8                         16           2
				 3		   1                         10           3
				 4		   9                         20           3
				 5		   9                         18           4
				 6		   6                          6           6
				 7		   7                          1           7
				 8		  12                          7           7
				 9		  10                          2           8
				10		   3                          4           9
				11		  13                          5           9
				12		  13                          9          10
				13		  15                         15          10
				14		  14                          8          12
				15		  10                         11          13
				16		   2                         12          13
				17		  16                         14          14
				18		   4                         19          14
				19		  14                         13          15
				20		   3                          2          16



	 * \endverbatim
	 *
	 * We cannot explain here what is a cache, but in practice is a fast memory in the CPU able
	 * to store chunks of memory. The cache in general is much smaller than RAM, but the big advantage
	 * is its speed. Retrieve data from the cache is much faster than RAM. Unfortunately the factors
	 *  that determine what is on cache and what is not are multiples: Type of cache, algorithm ... .
	 * Qualitatively all caches will tend to load chunks of data that you read multiple-time, or chunks
	 *  of data that probably you will read based on pattern analysis. A small example is a linear memory copy where
	 * you read consecutively memory and you write on consecutive memory.
	 * Modern CPU recognize such pattern and decide to load on cache the consecutive memory before
	 *  you actually require it.
	 *
	 *
	 * Reordering the vector in the way described above has the advantage that when we do computation on particles
	 * and its neighborhood consecutively (from the first to the end) we will have the tendency to:
	 *
	 * * read always the same particles
	 * * read the memory more consecutively or in predictable way
	 *
	 * That are the 2 factor required to take advantage of the cache
	 *
	 *
	 * \snippet Vector/4_reorder/main_data_ord.cpp vector reorder
	 *
	 * In order to show in practice what happen we first show the graph when we do not reorder
	 *
	 * \htmlinclude Vector/4_reorder/no_reorder.html
	 *
	 * The measure has oscillation but we see an asymptotic behavior from 0.04 in the initial condition to
	 * 0.124 . Below we show what happen when we reorder every 10000 iterations
	 *
	 * \htmlinclude Vector/4_reorder/reorder_10000.html
	 *
	 * As we can see the reorder at iteration 0
	 * does not produce any effect or performance improve in force calculation. This is because the particles
	 *  on a grid are already ordered enough. As soon as
	 * the particle move the situation degrade and the calculation on the force increase. At 10000 we try
	 * to reorder the particles and as we can see this reduce drastically the time to compute the forces.
	 *
	 *  \note There is still a gap between the initial condition and the reordered situation. This is because
	 *        the particles from the initial conditions are getting near. This increase
	 *        the average number of neighborhood per particle and so the overall computational
	 *        cost of the force calculation
	 *
	 * This show how the reordering of the data can significantly improve the performance. As soon as the
	 * the particles move we see a progressive degrade of the performance. The particles
	 * will have the tendency to disorder again. At 20000 we do the same things, and as we can see after
	 * reordering we get again a dramatic increase in performance. From this graph we clearly see that wait
	 * 10000 iteration is too much for reordering. We have to increase the frequency of reordering. Finding
	 * a good frequency depend from
	 *
	 *  * How fast the performance degrade (Or how fast particles move in a disordered way)
	 *  * How reorder is heavy compared to the full iteration
	 *
	 *  In this configuration
	 *
	 *  * Every 600 time step we get around 10% performance degrade
	 *  * Reordering is in general one order of magnitude than the calculation of the force
	 *
	 *  Even if the optimal is an higher frequency, we decide to reorder every 200 iteration. Getting
	 *  the following graph where the force calculation improve of a factor a little smaller than 2.
	 *
	 * \htmlinclude Vector/4_reorder/reorder_200.html
	 *
	 * In cases where particles move a lot, and so when degradation of performance is fast, consider to
	 * use computation reordering
	 *
	 * \see \ref e4_comp_reo
	 *
	 *  ## Timers ##
	 *
	 *  In order to collect the time of the force calculation we insert two timers around the function
	 *  calc_force. The overall performance is instead calculated with another timer around the time stepping
	 *
	 * \snippet Vector/4_reorder/main_data_ord.cpp timer start
	 * \snippet Vector/4_reorder/main_data_ord.cpp timer stop
	 *
	 * \see \ref e3_md_vi
	 *
	 * \snippet Vector/4_reorder/main_data_ord.cpp md steps
	 *
	 */

	//! \cond [md steps] \endcond

	// Get the Cell list structure
	auto NN = vd.getCellList(r_cut);

	// calculate forces
	calc_forces(vd,NN,sigma12,sigma6);
	unsigned long int f = 0;

	timer time2;
	time2.start();

	#ifndef TEST_RUN
	size_t Nstep = 30000;
	#else
	size_t Nstep = 300;
	#endif

	// MD time stepping
	for (size_t i = 0; i < Nstep ; i++)
	{
		// Get the iterator
		auto it3 = vd.getDomainIterator();

		// integrate velicity and space based on the calculated forces (Step1)
		while (it3.isNext())
		{
			auto p = it3.get();

			// here we calculate v(tn + 0.5)
			vd.template getProp<velocity>(p)[0] += 0.5*dt*vd.template getProp<force>(p)[0];
			vd.template getProp<velocity>(p)[1] += 0.5*dt*vd.template getProp<force>(p)[1];
			vd.template getProp<velocity>(p)[2] += 0.5*dt*vd.template getProp<force>(p)[2];

			// here we calculate x(tn + 1)
			vd.getPos(p)[0] += vd.template getProp<velocity>(p)[0]*dt;
			vd.getPos(p)[1] += vd.template getProp<velocity>(p)[1]*dt;
			vd.getPos(p)[2] += vd.template getProp<velocity>(p)[2]*dt;

			++it3;
		}

		// Because we mooved the particles in space we have to map them and re-sync the ghost
		vd.map();
		vd.template ghost_get<>();

		//! \cond [vector reorder] \endcond

		if (i % 200 == 0)
			vd.reorder(5);

		//! \cond [vector reorder] \endcond

		//! \cond [timer start] \endcond

		timer time;
		if (i % 10 == 0)
			time.start();

		//! \cond [timer start] \endcond

		// calculate forces or a(tn + 1) Step 2
		calc_forces(vd,NN,sigma12,sigma6);

		//! \cond [timer stop] \endcond

		if (i % 10 == 0)
		{
			time.stop();
			x.add(i);
			y.add({time.getwct()});
		}

		//! \cond [timer stop] \endcond

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
			// We calculate the energy
			double energy = calc_energy(vd,NN,sigma12,sigma6);
			auto & vcl = create_vcluster();
			vcl.sum(energy);
			vcl.execute();

			// We also print on terminal the value of the energy
			// only one processor (master) write on terminal
			if (vcl.getProcessUnitID() == 0)
				std::cout << "Energy: " << energy << std::endl;

			f++;
		}
	}

	time2.stop();
	std::cout << "Performance: " << time2.getwct() << std::endl;

	//! \cond [md steps] \endcond

	/*!
	 * \page Vector_4_reo Vector 4 data-reordering and cache friendliness
	 *
	 * ## Plotting graphs ##
	 *
	 * This code follow the same as the one in molecular dynamic code. The difference is that in this case
	 * the output is the computation time of the force for each iteration
	 *
	 * \see \ref e3_md_pg
	 *
	 * \snippet Vector/4_reorder/main_data_ord.cpp google chart
	 *
	 */

	//! \cond [google chart] \endcond

	// Google charts options, it store the options to draw the X Y graph
	GCoptions options;

	// Title of the graph
	options.title = std::string("Force calculation time");

	// Y axis name
	options.yAxis = std::string("Time");

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
	 * \page Vector_4_reo Vector 4 data-reordering and cache friendliness
	 *
	 * ## Finalize ## {#finalize}
	 *
	 *  At the very end of the program we have always de-initialize the library
	 *
	 * \snippet Vector/4_reorder/main_data_ord.cpp finalize
	 *
	 */

	//! \cond [finalize] \endcond

	openfpm_finalize();

	//! \cond [finalize] \endcond

	/*!
	 * \page Vector_4_reo Vector 4 data-reordering and cache friendliness
	 *
	 * # Full code # {#code}
	 *
	 * \include Vector/4_reorder/main_data_ord.cpp
	 *
	 */
}




