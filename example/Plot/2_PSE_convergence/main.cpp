/*
 * ### WIKI 1 ###
 *
 * ## Simple example
 *
 * In this example we show 1D PSE derivative function approximation
 * convergence plot
 *
 * ### WIKI END ###
 *
 */

#include "config/config.h"

// Some compiler like clang does not have libquadmath this example
// require libquadmath. So is active only if in openfpm installation
// such library has been detected

#ifdef HAVE_LIBQUADMATH

#include "VCluster.hpp"
#include "PSE/Kernels_test_util.hpp"
#include "Plot/GoogleChart.hpp"
#include <boost/multiprecision/float128.hpp>

/*
 *
 * ### WIKI 2 ###
 *
 * In this header we encapsulate the example 0_Derivative_approx_1D.
 * In particular instead of printing the l2 and l_infinity norm, it
 * expose a function PSE_test that given a number of particles it
 * return the calculated l2 and l_inf error for the second order Laplacian
 * approximation of the function x*e^{-x^2} in the interval 0.0 and 4.0
 *
 */

/*
 *
 * ### WIKI END ###
 *
 */

int main(int argc, char* argv[])
{
	//
	// ### WIKI 2 ###
	//
	// Here we Initialize the library and we define Ghost size
	// and non-periodic boundary conditions
	//
	openfpm_init(&argc,&argv);
	Vcluster & v_cl = create_vcluster();

	// This test is not made to run in parallel
	if (v_cl.getProcessingUnits() > 1)
	{
		std::cerr << "Error: only one processor is allowed" << "\n";
		return 1;
	}

	openfpm::vector<size_t> x;
	openfpm::vector<openfpm::vector<double>> y;

	openfpm::vector<std::string> yn;
	yn.add("(order 2) float128 overlap=2");
	yn.add("(order 2) float128 overlap=4");
	yn.add("(order 2) overlap 2 double ");
	yn.add("(order 2) overlap 4 double (order 2)");
	yn.add("(order 2) overlap 2 float (order 2)");
	yn.add("(order 2) overlap 4 float (order 2)");

	yn.add("(order 4) float128 overlap=2");
	yn.add("(order 4) float128 overlap=4");

	yn.add("(order 6) float128 overlap=2");
	yn.add("(order 6) float128 overlap=4");

	yn.add("(order 8) float128 overlap=2");
	yn.add("(order 8) float128 overlap=4");

	// Every time increase the number of particles by 2
	for (size_t i = 125 ; i <= 2097152000 ; i*=2)
	{
		x.add(i);
		y.add();

		PSEError err;

		/////// Order 2 //////////////

		PSE_test<boost::multiprecision::float128,Lap_PSE<1,boost::multiprecision::float128,2>>(i,2,err);
		std::cout << "Particles: " << i << " error: " << err.linf_error << std::endl;
		y.last().add(err.linf_error);

		PSE_test<boost::multiprecision::float128,Lap_PSE<1,boost::multiprecision::float128,2>>(i,4,err);
		std::cout << "Particles: " << i << " error: " << err.linf_error << std::endl;
		y.last().add(err.linf_error);

		PSE_test<double,Lap_PSE<1,double,2>>(i,2,err);
		std::cout << "Particles: " << i << " error: " << err.linf_error << std::endl;
		y.last().add(err.linf_error);

		PSE_test<double,Lap_PSE<1,double,2>>(i,4,err);
		std::cout << "Particles: " << i << " error: " << err.linf_error << std::endl;
		y.last().add(err.linf_error);

		PSE_test<float,Lap_PSE<1,float,2>>(i,2,err);
		std::cout << "Particles: " << i << " error: " << err.linf_error << std::endl;
		y.last().add(err.linf_error);

		PSE_test<float,Lap_PSE<1,float,2>>(i,4,err);
		std::cout << "Particles: " << i << " error: " << err.linf_error << std::endl;
		y.last().add(err.linf_error);

		//////// Order 4 /////////////

		PSE_test<boost::multiprecision::float128,Lap_PSE<1,boost::multiprecision::float128,4>>(i,2,err);
		std::cout << "Particles: " << i << " error: " << err.linf_error << std::endl;
		y.last().add(err.linf_error);

		PSE_test<boost::multiprecision::float128,Lap_PSE<1,boost::multiprecision::float128,4>>(i,4,err);
		std::cout << "Particles: " << i << " error: " << err.linf_error << std::endl;
		y.last().add(err.linf_error);


		//////// Order 6 /////////////

		PSE_test<boost::multiprecision::float128,Lap_PSE<1,boost::multiprecision::float128,6>>(i,2,err);
		std::cout << "Particles: " << i << " error: " << err.linf_error << std::endl;
		y.last().add(err.linf_error);

		PSE_test<boost::multiprecision::float128,Lap_PSE<1,boost::multiprecision::float128,6>>(i,4,err);
		std::cout << "Particles: " << i << " error: " << err.linf_error << std::endl;
		y.last().add(err.linf_error);

		//////// Order 8 /////////////

		PSE_test<boost::multiprecision::float128,Lap_PSE<1,boost::multiprecision::float128,8>>(i,8,err);
		std::cout << "Particles: " << i << " error: " << err.linf_error << std::endl;
		y.last().add(err.linf_error);

		PSE_test<boost::multiprecision::float128,Lap_PSE<1,boost::multiprecision::float128,8>>(i,16,err);
		std::cout << "Particles: " << i << " error: " << err.linf_error << std::endl;
		y.last().add(err.linf_error);
	}

	// Google charts options
	GCoptions options;

	options.title = std::string("PSE Laplacian convergence");
	options.yAxis = std::string("Y Axis");
	options.xAxis = std::string("X Axis");
	options.more = std::string("hAxis:{logScale: true},vAxis:{logScale: true,format: 'scientific'}");
	options.lineWidth = 1.0;

	GoogleChart cg;
	cg.AddLinesGraph(x,y,yn,options);
	cg.write("PSE_plot.html");

	//
	// ### WIKI 8 ###
	//
	// Deinitialize the library
	//
	openfpm_finalize();
}

#else

int main()
{

}

#endif
