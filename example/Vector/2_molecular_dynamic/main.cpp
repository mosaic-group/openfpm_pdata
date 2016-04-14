
#include "Vector/vector_dist.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "data_type/aggregate.hpp"
#include "Plot/GoogleChart.hpp"
#include "Plot/util.hpp"

/*
 * ### WIKI 1 ###
 *
 * ## Molecular Dynamic with Lennard-Jones potential
 *
 * This example show a simple Lennard-Jones molecular dynamic simulation in a stable regime
 *
 * ### WIKI END ###
 *
 */

constexpr int velocity = 0;
constexpr int force = 1;

/* ### WIKI 11 ###
 *
 * The function to calculate the forces between particles. It require the vector of particles
 * Cell list and scaling factor.
 *
 */
void calc_forces(vector_dist<3,double, aggregate<double[3],double[3]> > & vd, CellList<3, double, FAST, shift<3, double> > & NN, double L)
{

	// ### WIKI 12 ###
	//
	// Update the cell list from the actual particle configuration
	//
	//
	vd.updateCellList(NN);

	// ### WIKI 13 ###
	//
	// Calculate the forces
	//
	auto it2 = vd.getDomainIterator();
	while (it2.isNext())
	{
		auto p = it2.get();

		Point<3,double> xp = vd.getPos<0>(p);

		vd.template getProp<force>(p)[0] = 0.0;
		vd.template getProp<force>(p)[1] = 0.0;
		vd.template getProp<force>(p)[2] = 0.0;

		// For each neighborhood particle
		auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos<0>(p)));

		while (Np.isNext())
		{
			// Neighborhood particle q
			auto q = Np.get();

			if (q == p.getKey())	{++Np; continue;};

			// repulsive
			Point<3,double> xq = vd.getPos<0>(q);
			Point<3,double> r = xp - xq;

			// take the norm, normalize
			float rn = r.norm();
			r /= rn;
			rn *= L;

			// Calculate the force, using pow is slower
			Point<3,double> f = 24.0*(2.0 / (rn*rn*rn*rn*rn*rn*rn*rn*rn*rn*rn*rn*rn) -  1.0 / (rn*rn*rn*rn*rn*rn*rn)) * r;

			// we sum the force produced by q on p
			vd.template getProp<force>(p)[0] += f.get(0);
			vd.template getProp<force>(p)[1] += f.get(1);
			vd.template getProp<force>(p)[2] += f.get(2);

			++Np;
		}

		++it2;
	}
}

/* ### WIKI 14 ###
 *
 * The function to calculate the total energy. It require the same parameter as calculate forces
 *
 */
double calc_energy(vector_dist<3,double, aggregate<double[3],double[3]> > & vd, CellList<3, double, FAST, shift<3, double> > & NN, double L)
{
	// ### WIKI 15 ###
	//
	// Reset the counter for the energy and
	// update the cell list from the actual particle configuration
	//
	//
	double E = 0.0;
	vd.updateCellList(NN);

	// ### WIKI 16 ###
	//
	// Calculate the forces
	//
	auto it2 = vd.getDomainIterator();
	while (it2.isNext())
	{
		auto p = it2.get();

		Point<3,double> xp = vd.getPos<0>(p);

		vd.template getProp<force>(p)[0] = 0.0;
		vd.template getProp<force>(p)[1] = 0.0;
		vd.template getProp<force>(p)[2] = 0.0;

		// For each neighborhood particle
		auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos<0>(p)));

		while (Np.isNext())
		{
			// Neighborhood particle q
			auto q = Np.get();

			if (q == p.getKey())	{++Np; continue;};

			Point<3,double> xq = vd.getPos<0>(q);
			Point<3,double> r = xp - xq;

			// take the normalized direction
			float rn = r.norm();
			r /= rn;

			rn *= L;

			// potential energy (using pow is slower)
			E += 4.0 * ( 1.0 / (rn*rn*rn*rn*rn*rn*rn*rn*rn*rn*rn*rn) + 1.0 / ( rn*rn*rn*rn*rn*rn) );

			// Kinetic energy
			E += sqrt(vd.template getProp<force>(p)[0]*vd.template getProp<force>(p)[0] +
					vd.template getProp<force>(p)[1]*vd.template getProp<force>(p)[1] +
					vd.template getProp<force>(p)[2]*vd.template getProp<force>(p)[2]);

			++Np;
		}

		++it2;
	}

	return E;
}

int main(int argc, char* argv[])
{
	//
	// ### WIKI 2 ###
	//
	// Here we define important parameters or the simulation, time step integration,
	// size of the box, and cut-off radius of the interaction
	//
	double dt = 0.005;
	double L = 10.0;
	float r_cut = 0.3;

	openfpm::vector<double> x;
	openfpm::vector<openfpm::vector<double>> y;

	//
	// ### WIKI 2 ###
	//
	// Here we Initialize the library, we create a Box that define our domain, boundary conditions, ghost
	// and the grid size
	//
	openfpm_init(&argc,&argv);
	Vcluster & v_cl = create_vcluster();

	// we create a 10x10x10 Grid iterator
	size_t sz[3] = {10,10,10};

	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// ghost, big enough to contain the interaction radius
	Ghost<3,float> ghost(r_cut);

	//
	// ### WIKI 3 ###
	//
	// Here we define a distributed vector in 3D, containing 3 properties, a
	// scalar double, a vector double[3], and a tensor or rank 2 double[3][3].
	// In this case the vector contain 0 particles in total
	//
	vector_dist<3,double, aggregate<double[3],double[3]> > vd(0,box,bc,ghost);

	//
	// ### WIKI 4 ###
	//
	// We define a grid iterator, to create particles on a grid like way.
	// An important note is that the grid iterator, iterate only on the
	// local nodes for each processor for example suppose to have a domain like
	// the one in figure
	//
	//   +---------+
	//   |* * *|* *|
	//   |  2  |   |
	//   |* * *|* *|
	//   |   ---   |
	//   |* *|* * *|
	//   |   |     |
	//   |* *|* * *|
	//   |   |  1  |
	//   |* *|* * *|
	//   +---------+
	//
	// divided in 2 processors, the processor 1 will iterate only on the points
	// inside the portion of space marked with one. Also the grid iterator follow the
	// boundary condition specified in vector. For a perdiodic 2D 5x5 grid we have
	//
	//   +---------+
	//   * * * * * |
	//   |         |
	//   * * * * * |
	//   |         |
	//   * * * * * |
	//   |         |
	//   * * * * * |
	//   |         |
	//   *-*-*-*-*-+
	//
	// Because the right border is equivalent to the left border, while for a non periodic we have the
	// following distribution of points
	//
	//   *-*-*-*-*
	//   |       |
	//   * * * * *
	//   |       |
	//   * * * * *
	//   |       |
	//   * * * * *
	//   |       |
	//   *-*-*-*-*
	//
	// The loop will place particles on a grid like on each processor
	//
	auto it = vd.getGridIterator(sz);

	while (it.isNext())
	{
		vd.add();

		auto key = it.get();

		vd.template getLastPos<0>()[0] = key.get(0) * it.getSpacing(0);
		vd.template getLastPos<0>()[1] = key.get(1) * it.getSpacing(1);
		vd.template getLastPos<0>()[2] = key.get(2) * it.getSpacing(2);

		vd.template getLastProp<velocity>()[0] = 0.0;
		vd.template getLastProp<velocity>()[1] = 0.0;
		vd.template getLastProp<velocity>()[2] = 0.0;

		vd.template getLastProp<force>()[0] = 0.0;
		vd.template getLastProp<force>()[1] = 0.0;
		vd.template getLastProp<force>()[2] = 0.0;

		++it;
	}

	vd.map();

	//
	// ### WIKI 5 ###
	//
	// we get the cell list to compute the neighborhood of the particles
	auto NN = vd.getCellList(r_cut);

	// calculate forces a(tn)
	calc_forces(vd,NN,L);
	unsigned long int f = 0;

	//
	// ### WIKI 6 ###
	//
	// Here we do 100 MD steps using verlet integrator
	//
	// $$ \vec{v}(t_{n+1/2}) = \vec{v}_p(t_n) + \frac{1}{2} \delta t \vec{a}(t_n) $$
	// $$ \vec{x}(t_{n}) = \vec{x}_p(t_n) + \delta t \vec{v}(t_n+1/2) $$
	//
	// calculate the forces $$ \vec{a} (t_{n}) $$ from $$ \vec{x} (t_{n}) $$
	//
	// $$ \vec{v}(t_{n+1}) = \vec{v}_p(t_n+1/2) + \frac{1}{2} \delta t \vec{a}(t_n+1) $$
	//
	//
	for (size_t i = 0; i < 10000 ; i++)
	{
		auto it3 = vd.getDomainIterator();

		while (it3.isNext())
		{
			auto p = it3.get();

			// here we calculate v(tn + 0.5)
			vd.template getProp<velocity>(p)[0] += 0.5*dt*vd.template getProp<force>(p)[0];
			vd.template getProp<velocity>(p)[1] += 0.5*dt*vd.template getProp<force>(p)[1];
			vd.template getProp<velocity>(p)[2] += 0.5*dt*vd.template getProp<force>(p)[2];

			// here we calculate x(tn + 1)
			vd.template getPos<0>(p)[0] += vd.template getProp<velocity>(p)[0]*dt;
			vd.template getPos<0>(p)[1] += vd.template getProp<velocity>(p)[1]*dt;
			vd.template getPos<0>(p)[2] += vd.template getProp<velocity>(p)[2]*dt;

			++it3;
		}

		vd.map();
		vd.template ghost_get<>();

		// calculate forces or a(tn + 1)
		calc_forces(vd,NN,L);

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

		if (i % 100 == 0)
		{
			vd.write("particles",f);
			double energy = calc_energy(vd,NN,L);
			auto & vcl = create_vcluster();
			vcl.sum(energy);
			vcl.execute();

			x.add(i);
			y.add({energy});
			if (vcl.getProcessUnitID() == 0)
				std::cout << "Energy: " << energy << std::endl;

			f++;
		}
	}

	// Google charts options
	GCoptions options;

	options.title = std::string("Energy with time");
	options.yAxis = std::string("Energy");
	options.xAxis = std::string("iteration");
	options.lineWidth = 1.0;

	GoogleChart cg;
	cg.AddLinesGraph(x,y,options);
	cg.write("gc_plot2_out.html");

	//
	// ### WIKI 10 ###
	//
	// Deinitialize the library
	//
	openfpm_finalize();
}




