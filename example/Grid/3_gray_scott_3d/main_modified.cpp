#include "Grid/grid_dist_id.hpp"
#include "data_type/aggregate.hpp"
#include "timer.hpp"
#include "SubdomainGraphNodes.hpp"
#include "Decomposition/Distribution/BoxDistribution.hpp"

/*!
 *
 * \page Grid_3_gs_3D Gray Scott in 3D
 *
 * [TOC]
 *
 * # Solving a gray scott-system in 3D # {#e3_gs_gray_scott}
 *
 * This example is just an extension of the 2D Gray scott example.
 * Here we show how to solve a non-linear reaction diffusion system in 3D
 *
 * In figure is the final solution of the problem
 *
 * \htmlonly
 * <img src="http://ppmcore.mpi-cbg.de/web/images/examples/gray_scott_3d/gs_alpha.png"/>
 * \endhtmlonly
 *
 * More or less this example is the adaptation of the previous example to 3D
 * with the improvement of using stencil iterator.
 *
 * ## Stencil iterator {#e3_gs_grat_scott_si}
 *
 * Stencil iterator require that you define a stencil,
 *
 * \snippet Grid/3_gray_scott_3d/main.cpp stencil def
 *
 * once is defined it is
 * possible get and use a stencil iterator
 *
 * \snippet Grid/3_gray_scott_3d/main.cpp stencil get and use
 *
 * The rest of the example remain the same with the exception
 * that the code has been extended in 3D.
 *
 * \see \ref Grid_2_solve_eq
 *
 * 
 */

//! \cond [constants] \endcond

constexpr int U = 0;
constexpr int V = 1;

constexpr int x = 0;
constexpr int y = 1;
constexpr int z = 2;

//! \cond [constants] \endcond

typedef grid_dist_id<3, double, aggregate<double,double>, CartDecomposition<3,double,HeapMemory,memory_traits_lin,BoxDistribution<3,double>>> grid_type;


void draw_box(double x_s, double y_s, double z_s, double x_e, double y_e, double z_e, grid_type & Old, Box<3, double> & domain)
//void draw_box(double x_s, double y_s, double z_s, double x_e, double y_e, double z_e, grid_dist_id<3,double,aggregate<double,double>> & Old, Box<3, double> & domain)
{
    long int x_start = Old.size(0)*x_s/domain.getHigh(0);
    long int y_start = Old.size(1)*y_s/domain.getHigh(1);
    long int z_start = Old.size(1)*z_s/domain.getHigh(2);

    long int x_stop = Old.size(0)*x_e/domain.getHigh(0);
    long int y_stop = Old.size(1)*y_e/domain.getHigh(1);
    long int z_stop = Old.size(1)*z_e/domain.getHigh(2);

    grid_key_dx<3> start({x_start,y_start,z_start});
    grid_key_dx<3> stop ({x_stop,y_stop,z_stop});
    auto it_init = Old.getSubDomainIterator(start,stop);

    while (it_init.isNext())
    {
        auto key = it_init.get();

        Old.template get<U>(key) = 0.5 + (((double)std::rand())/RAND_MAX -0.5)/10.0;
        Old.template get<V>(key) = 0.25 + (((double)std::rand())/RAND_MAX -0.5)/20.0;

        ++it_init;
    }
}

void init(grid_type & Old, grid_type & New, Box<3,double> & domain)
//void init(grid_dist_id<3,double,aggregate<double,double> > & Old, grid_dist_id<3,double,aggregate<double,double> > & New, Box<3,double> & domain)
{
	auto it = Old.getDomainIterator();

	while (it.isNext())
	{
		// Get the local grid key
		auto key = it.get();

		// Old values U and V
		Old.template get<U>(key) = 1.0;
		Old.template get<V>(key) = 0.0;

		// Old values U and V
		New.template get<U>(key) = 0.0;
		New.template get<V>(key) = 0.0;

		++it;
	}
//    draw_box(0.05, 0.05, 0.05, 1.95, 1.95, 1.95, Old, domain);

	draw_box(0.75, 0.75, 0.75, 1.95, 1.95, 1.95, Old, domain);
    draw_box(0.15, 0.15, 0.15, 1.15, 1.15, 1.15, Old, domain);
    draw_box(0.75, 0.75, 0.15, 1.95, 1.95, 0.75, Old, domain);
    draw_box(0.75, 0.15, 0.75, 1.95, 1.15, 1.95, Old, domain);
    draw_box(0.15, 0.75, 0.75, 0.75, 1.95, 1.95, Old, domain);
    draw_box(0.75, 0.75, 0.15, 1.95, 1.95, 1.15, Old, domain);

//	long int x_start = Old.size(0)*1.55f/domain.getHigh(0);
//	long int y_start = Old.size(1)*1.55f/domain.getHigh(1);
//	long int z_start = Old.size(1)*1.55f/domain.getHigh(2);
//
//	long int x_stop = Old.size(0)*1.85f/domain.getHigh(0);
//	long int y_stop = Old.size(1)*1.85f/domain.getHigh(1);
//	long int z_stop = Old.size(1)*1.85f/domain.getHigh(2);
//
//	grid_key_dx<3> start({x_start,y_start,z_start});
//	grid_key_dx<3> stop ({x_stop,y_stop,z_stop});
//	auto it_init = Old.getSubDomainIterator(start,stop);
//
//	while (it_init.isNext())
//	{
//		auto key = it_init.get();
//
//                Old.template get<U>(key) = 0.5 + (((double)std::rand())/RAND_MAX -0.5)/10.0;
//                Old.template get<V>(key) = 0.25 + (((double)std::rand())/RAND_MAX -0.5)/20.0;
//
//		++it_init;
//	}
}


int main(int argc, char* argv[])
{
	openfpm_init(&argc,&argv, init_options::in_situ_visualization);

	// domain
	Box<3,double> domain({0.0,0.0},{2.5,2.5,2.5});
	
	// grid size
        size_t sz[3] = {1500,1500,1500};

	// Define periodicity of the grid
	periodicity<3> bc = {PERIODIC,PERIODIC,PERIODIC};
	
	// Ghost in grid unit
	Ghost<3,long int> g(1);
    Ghost<3,long int> g_zero(0);

    // deltaT
	double deltaT = 1.0/50.0;

	// Diffusion constant for specie U
	double du = 2*1e-5;

	// Diffusion constant for specie V
	double dv = 1*1e-5;

	// Number of timesteps
        size_t timeSteps = 200;

	// K and F (Physical constant in the equation)
        double K = 0.053;
        double F = 0.014;

    grid_type Old(sz,domain,g,bc);

	// New grid with the decomposition of the old grid
    grid_type New(Old.getDecomposition(),sz,g);

//	grid_dist_id<3, double, aggregate<double,double>> Old(sz,domain,g,bc);

	// New grid with the decomposition of the old grid
//	grid_dist_id<3, double, aggregate<double,double>> New(Old.getDecomposition(),sz,g);
	
//	grid_dist_id<3, double, aggregate<unsigned short>> Vis_new(Old.getDecomposition(),sz,g);

    // spacing of the grid on x and y
	double spacing[3] = {Old.spacing(0),Old.spacing(1),Old.spacing(2)};

	init(Old,New,domain);
//	Old.write("InitGrayScott");

    auto &v_cl = create_vcluster();

    size_t i = 0;

    if (argc > 1)
    {
        // take that argument
        std::string restart(argv[1]);

        // convert it into number
        i = std::stoi(restart);

        // load the file
        Old.load("checkpoint" + std::to_string(i));

        // Print to inform that we are restarting from a
        // particular position
        if (v_cl.getProcessUnitID() == 0)
        {std::cout << "Restarting from " << i << std::endl;}
    }

	// sync the ghost
	size_t count = 0;
	Old.template ghost_get<U,V>();

	// because we assume that spacing[x] == spacing[y] we use formula 2
	// and we calculate the prefactor of Eq 2
	double uFactor = deltaT * du/(spacing[x]*spacing[x]);
	double vFactor = deltaT * dv/(spacing[x]*spacing[x]);

	timer tot_sim;
	tot_sim.start();

	//! \cond [stencil def] \endcond

	static grid_key_dx<3> star_stencil_3D[7] = {{0,0,0},
                                         	    {0,0,-1},
												{0,0,1},
												{0,-1,0},
												{0,1,0},
												{-1,0,0},
												{1,0,0}};

	//! \cond [stencil def] \endcond

    timer timercomp;
    timer timerghost;
    timer timervisualize;

    double visualizetime = 0.0;
    double comptime = 0.0;
    double ghosttime = 0.0;

    for (; i < timeSteps; ++i)
	{
        if (i % 30 == 0)
		{std::cout << "STEP: " << i << std::endl;}

	timercomp.start();
        //timer tfake;
        //tfake.start();

	//if (v_cl.rank() == 0)
        //{
        //        std::cout << "before" << std::endl;
        //}

	//double fake_comp = 0.0;
        //while(tfake.getwct() < 2.0)
	//{
	//	fake_comp += 1.0; 
	//}

        //if (v_cl.rank() == 0)
	//{
        //	std::cout << fake_comp << std::endl;
	//}

        //! \cond [stencil get and use] \endcond

        auto it = Old.getDomainIteratorStencil(star_stencil_3D);

        while (it.isNext())
		{
            // center point
            auto Cp = it.getStencil<0>();

            // plus,minus X,Y,Z
            auto mx = it.getStencil<1>();
            auto px = it.getStencil<2>();
            auto my = it.getStencil<3>();
            auto py = it.getStencil<4>();
            auto mz = it.getStencil<5>();
            auto pz = it.getStencil<6>();

            // update based on Eq 2
            New.get<U>(Cp) = Old.get<U>(Cp) + uFactor * (
										Old.get<U>(mz) +
										Old.get<U>(pz) +
										Old.get<U>(my) +
										Old.get<U>(py) +
										Old.get<U>(mx) +
										Old.get<U>(px) -
										6.0*Old.get<U>(Cp)) +
										- deltaT * Old.get<U>(Cp) * Old.get<V>(Cp) * Old.get<V>(Cp) +
										- deltaT * F * (Old.get<U>(Cp) - 1.0);


            // update based on Eq 2
            New.get<V>(Cp) = Old.get<V>(Cp) + vFactor * (
										Old.get<V>(mz) +
										Old.get<V>(pz) +
										Old.get<V>(my) +
										Old.get<V>(py) +
										Old.get<V>(mx) +
                                        Old.get<V>(px) -
										6*Old.get<V>(Cp)) +
										deltaT * Old.get<U>(Cp) * Old.get<V>(Cp) * Old.get<V>(Cp) +
										- deltaT * (F+K) * Old.get<V>(Cp);

            // Next point in the grid
            ++it;
        }

        //! \cond [stencil get and use] \endcond

        // Here we copy New into the old grid in preparation of the new step
        // It would be better to alternate, but using this we can show the usage
        // of the function copy. To note that copy work only on two grid of the same
        // decomposition. If you want to copy also the decomposition, or force to be
        // exactly the same, use Old = New
        Old.copy(New);

	comptime += timercomp.getwct();

	timerghost.start();
        // After copy we synchronize again the ghost part U and V
        Old.ghost_get<U,V>();

	ghosttime += timerghost.getwct();

        //if ( i % 1 == 0) std::cout<<"The maximum concentration is " << 20 << " and the minimum is " << 0 <<std::endl;

	//if(i==0)
	timervisualize.start();
	New.visualize<V>(scale_vis::fixed, 0.0, 1.0);
	visualizetime+=timervisualize.getwct();
//        if(i == (int)2000/deltaT || i == (int)3000/deltaT || i == (int)4000/deltaT )
//        {
//            Old.save("checkpoint" + std::to_string(i));
//        }
	}

	std::cout << "Compute time " << comptime << std::endl;	
	std::cout << "Ghost-get time: " << ghosttime << std::endl; 
	std::cout << "Time in visualize function " << visualizetime << std::endl;      
	double my_time = tot_sim.getwct();
        v_cl.max(my_time);
	v_cl.max(comptime);
	v_cl.sum(ghosttime);
	v_cl.sum(visualizetime);
        v_cl.execute();
        std::cout << "Total simulation: " << tot_sim.getwct() << std::endl;
        std::cout<< "Max time is: " << my_time <<std::endl;


	std::cout << "Max compute time " << comptime << std::endl;	
	std::cout << "Total ghost-get time: " << ghosttime << " and average: " << ghosttime/((double)v_cl.size()) << std::endl;       
	
	std::cout << "Total time in visualize: " << visualizetime << " and average: " << visualizetime/((double)v_cl.size()) << std::endl;
//! \cond [time stepping] \endcond

	/*!
	 * \page Grid_3_gs_3D Gray Scott in 3D
	 *
	 * ## Finalize ##
	 *
	 * Deinitialize the library
	 *
	 * \snippet Grid/3_gray_scott/main.cpp finalize
	 *
	 */

	//! \cond [finalize] \endcond

	openfpm_finalize();

	//! \cond [finalize] \endcond

	/*!
	 * \page Grid_3_gs_3D Gray Scott in 3D
	 *
	 * # Full code # {#code}
	 *
	 * \include Grid/3_gray_scott_3d/main.cpp
	 *
	 */
}
