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

#include "Vector/vector_dist.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "PSE/Kernels.hpp"
#include "Plot/util.hpp"
#include "Plot/GoogleChart.hpp"
#include "data_type/aggregate.hpp"
#include <cmath>


/*
 * ### WIKI 2 ###
 *
 * Here we define the function x*e^(-x*x) and its
 *
 */

double f_xex2(double x)
{
	return x*exp(-x*x);
}

double f_xex2(Point<1,double> & x)
{
	return x.get(0)*exp(-x.get(0)*x.get(0));
}

/*
 *
 * ### WIKI 3 ###
 *
 * It calculate the Laplacian of the function in one point using the PSE
 * Approximation
 *
 * \$ \frac{1}{\epsilon^{2}} h (u_{q} - u_{p}) \eta_{\epsilon}(x_q - x_p) \$
 *
 * \param p point
 * \param key vector
 * \param vd Distributed vector
 * \param eps epsilon of the kernel
 * \param spacing between particles
 * \param cl CellList
 *
 */
template<typename CellL> double calcLap(Point<1,double> p, vect_dist_key_dx key, vector_dist<1,double, aggregate<double>, CartDecomposition<1,double> > & vd, double eps, double spacing, CellL & cl)
{
	// Laplacian PSE kernel 1 dimension, on double, second order
	Lap_PSE<1,double,2> lker(eps);

	// double PSE integration accumulator
	double pse = 0;

	// Get f(x) at the position of the particle
	double prp_x = vd.template getProp<0>(key);

	// Get the neighborhood of the particles
	auto NN = cl.template getNNIterator<NO_CHECK>(cl.getCell(p));
	while(NN.isNext())
	{
		auto nnp = NN.get();

		// Calculate contribution given by the kernel value at position p,
		// given by the Near particle
		if (nnp != key.getKey())
		{
			// W(x-y)
			double ker = lker.value(p,vd.template getPos<0>(nnp));

			// f(y)
			double prp_y = vd.template getProp<0>(nnp);

			// 1.0/(eps)^2 [f(y)-f(x)]*W(x,y)*V_q
			double prp = 1.0/eps/eps * (prp_y - prp_x) * spacing;
			pse += prp * ker;
		}

    	// Next particle
    	++NN;
	}

	return pse;
}

/*
 *
 * ### WIKI END ###
 *
 */

int main(int argc, char* argv[])
{
	//
	// ### WIKI 3 ###
	//
	// Some useful parameters. Like
	//
	// * Number of particles
	// * Minimum number of padding particles
	// * The computational domain
	// * The spacing
	// * The mollification length
	// * Second order Laplacian kernel in 1D
	//

	// Number of particles
	const size_t Npart = 1000;

	// Number of step
	const size_t Nstep = 1000;

	// Delta t
	const double dt = 10.0 / Nstep;

	// Number of padding particles (At least)
	const size_t Npad = 20;

	// The domain
	Box<1,double> box({0.0},{4.0});

	// Calculated spacing
	double spacing = box.getHigh(0) / Npart;

	// The mollification length
	const double eps = 2*spacing;

	//
	// ### WIKI 2 ###
	//
	// Here we Initialize the library and we define Ghost size
	// and non-periodic boundary conditions
	//
	init_global_v_cluster(&argc,&argv);
	Vcluster & v_cl = *global_v_cluster;

    size_t bc[1]={NON_PERIODIC};
	Ghost<1,double> g(12*eps);

	//
	// ### WIKI 3 ###
	//
	// Here we are creating a distributed vector defined by the following parameters
	//
	// we create a set of N+1 particles to have a fully covered domain of particles between 0.0 and 4.0
	// Suppose we have a spacing given by 1.0 you need 4 +1 particles to cover your domain
	//
	vector_dist<1,double, aggregate<double>, CartDecomposition<1,double> > vd(Npart+1,box,bc,g);

	//
	// ### WIKI 4 ###
	//
	// We assign the position to the particles, the scalar property is set to
	// the function x*e^(-x*x) value.
	// Each processor has parts of the particles and fill part of the space, the
	// position is assigned independently from the decomposition.
	//
	// In this case, if there are 1001 particles and 3 processors the in the
	// domain from 0.0 to 4.0
	//
	// * processor 0 place particles from 0.0 to 1.332 (334 particles)
	// * processor 1 place particles from 1.336 to 2.668 (334 particles)
	// * processor 2 place particles from 2.672 to 4.0 (333 particles)
	//

	// It return how many particles the processors with id < rank has in total
	size_t base = vd.init_size_accum(Npart+1);
	auto it2 = vd.getIterator();

	while (it2.isNext())
	{
		auto key = it2.get();

		// set the position of the particles
		vd.template getPos<0>(key)[0] = (key.getKey() + base) * spacing;
		//set the property of the particles
		vd.template getProp<0>(key) = f_xex2((key.getKey() + base) * spacing);

		++it2;
	}

	//
	// ### WIKI 5 ###
	//
	// Once defined the position, we distribute them across the processors
	// following the decomposition, finally we get the ghost part
	//
	vd.map();
	vd.ghost_get<0>();

	//
	// ### WIKI 6 ###
	//
	// near the boundary we have two options, or we use one sided kernels,
	// or we add additional particles, it is required that such particles
	// produces a 2 time differentiable function. In order to obtain such
	// result we extend for x < 0.0 and x > 4.0 with the test function xe^(-x*x).
	//
	// Note that for x < 0.0 such extension is equivalent to mirror the
	// particles changing the sign of the strength
	//
	// \verbatim
	//
	//-0.6  -0.5      0.5  0.6   Strength
	//  +----+----*----*----*-
	//          0.0              Position
	//
	//  with * = real particle
	//       + = mirror particle
	//
	// \endverbatim
	//
	//
	Box<1,double> m_pad({0.0},{0.1});
	Box<1,double> m_pad2({3.9},{4.0});
	double enlarge = 0.1;

	// Create a box
	if (Npad * spacing > 0.1)
	{
		m_pad.setHigh(0,Npad * spacing);
		m_pad2.setLow(0,4.0 - Npad*spacing);
		enlarge = Npad * spacing;
	}

	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		// set the position of the particles
		if (m_pad.isInsideNB(vd.template getPos<0>(key)) == true)
		{
			vd.add();
			vd.template getLastPos<0>()[0] = - vd.template getPos<0>(key)[0];
			vd.template getLastProp<0>() = - vd.template getProp<0>(key);
		}

		// set the position of the particles
		if (m_pad2.isInsideNB(vd.template getPos<0>(key)) == true)
		{
			vd.add();
			vd.template getLastPos<0>()[0] = 2.0 * box.getHigh(0) - vd.template getPos<0>(key)[0];
			vd.template getLastProp<0>() = f_xex2(vd.template getLastPos<0>()[0]);
		}

		++it;
	}

	//
	// ### WIKI 6 ###
	//
	// We create a CellList with cell spacing 12 sigma
	//

    // get and construct the Cell list
	Ghost<1,double> gp(enlarge);
    auto cl = vd.getCellList(8*eps,gp);

    // Maximum infinity norm
    double linf = 0.0;

	//
	// ### WIKI 7 ###
	//
    // Euler time integration
    //
    //

    for (double t = 0; t <= 10 ; t += dt)
    {
    	auto it_p = vd.getDomainIterator();
    	while (it_p.isNext())
    	{
    		// key
    		vect_dist_key_dx key = it_p.get();

    		Point<1,double> p = vd.template getPos<0>(key);

    		// We are not interested in calculating out the domain
    		// note added padding particle are considered domain particles
    		if (p.get(0) < 0.0 || p.get(0) >= 4.0)
    		{
    			++it_p;
    			continue;
    		}

    		double pse = calcLap(p,key,vd,eps,spacing,cl);

    		// Euler update step
    		vd.template getProp<0>(key) = vd.template getProp<0>(key) + pse * dt;

    		++it_p;
    	}
    }

    // Once we have the value of the Laplacian we update with the eulerian step


	//
	// ### WIKI 8 ###
	//
    // Collect the solution in one processor
    //

    // Resize a vector with the local size of the vector
    openfpm::vector<double> v;
    v.resize(vd.size_local());

    // y contain 2 vectors the first is the calculated solution
    // the second contain the analytical solution
    openfpm::vector<openfpm::vector<double>> y;

    // Copy the property we want to Plot
    for (size_t i = 0 ; i < vd.size_local() ; i++)
    	y.get(i) = vd.template getProp<0>(i);

    // Gather the solution on master
    v_cl.SGather(v,y.get(0),0);

    //
    // ### WIKI 9 ###
    //
    // Here we plot
    //

    // If there are too many points terminate
    if (y.get(0).size() > 16000)
    {
    	std::cerr << "Too many points for plotting " << std::endl;
    	return 1;
    }

    openfpm::vector<double> x;

    // It fill the vector x with vg.size() points with interval values between 0.0 and 4.0
    // and fill the vector y;

    // Total points to plot
    size_t tot_p = y.get(0).size();

    Fill1D(0.0,4.0,tot_p,x);

    // It fill y.get(1) with the Analytical solution
    Fill1D(0.0,4.0,tot_p,y.get(1),f_xex2);

    // Plot with GooglePlot
	// Google charts options
	GCoptions options;

	options.yAxis = std::string("Y Axis");
	options.xAxis = std::string("X Axis");
	options.lineWidth = 1.0;

	GoogleChart cg;
	cg.AddPointsGraph(x,y,options);
	cg.write("PSE_plot.html");

	//
	// ### WIKI 8 ###
	//
	// Deinitialize the library
	//
	delete_global_v_cluster();
}
