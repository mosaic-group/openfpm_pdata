#include "Grid/grid_dist_id.hpp"
#include "data_type/aggregate.hpp"
#include "timer.hpp"

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

void draw_box(double x_s, double y_s, double z_s, double x_e, double y_e, double z_e, grid_dist_id<3,double,aggregate<double,double>> & Old, Box<3, double> & domain)
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

void init(grid_dist_id<3,double,aggregate<double,double> > & Old, grid_dist_id<3,double,aggregate<double,double> > & New, Box<3,double> & domain)
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

	draw_box(1.35, 1.35, 1.35, 1.95, 1.95, 1.95, Old, domain);
    draw_box(0.15, 0.15, 0.15, 0.75, 0.75, 0.75, Old, domain);
    draw_box(1.35, 1.35, 0.15, 1.95, 1.95, 0.75, Old, domain);
    draw_box(1.35, 0.15, 1.35, 1.95, 0.75, 1.95, Old, domain);
    draw_box(0.15, 1.35, 1.35, 0.35, 1.95, 1.95, Old, domain);
    draw_box(1.35, 1.35, 0.15, 1.95, 1.95, 0.75, Old, domain);

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
        size_t sz[3] = {128,128,128};

	// Define periodicity of the grid
	periodicity<3> bc = {PERIODIC,PERIODIC,PERIODIC};
	
	// Ghost in grid unit
	Ghost<3,long int> g(1);
    Ghost<3,long int> g_zero(0);

    // deltaT
	double deltaT = 1;

	// Diffusion constant for specie U
	double du = 2*1e-5;

	// Diffusion constant for specie V
	double dv = 1*1e-5;

	// Number of timesteps
        size_t timeSteps = 5000;

	// K and F (Physical constant in the equation)
        double K = 0.053;
        double F = 0.014;

	grid_dist_id<3, double, aggregate<double,double>> Old(sz,domain,g,bc);

	// New grid with the decomposition of the old grid
    grid_dist_id<3, double, aggregate<double,double>> New(Old.getDecomposition(),sz,g);

    grid_dist_id<3, double, aggregate<unsigned short>> Vis_new(Old.getDecomposition(),sz,g);

    Vis_new.visualize();

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


    for (; i < timeSteps; ++i)
	{
        if (i % 300 == 0)
		{std::cout << "STEP: " << i << std::endl;}

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

        // After copy we synchronize again the ghost part U and V
        Old.ghost_get<U,V>();

        // Every 500 time step we output the configuration for
        // visualization
//		if (i % 500 == 0)
//		{
//			Old.save("output_" + std::to_string(count));
//			Vis_new.write_frame("vis_output", count);
//			Old.write_frame("old_output", count);
//			count++;
//		}

        //find max and min velocity
        float maxConc = -1.0f;
        float minConc = 100000.0f;
        auto it1 = New.getDomainIterator();


        while(it1.isNext())
        {
            auto key = it1.get();

            float curConc = (float) New.template get<V>(key);

            if(curConc > maxConc)
            {
                maxConc = curConc;
            }

            if(curConc < minConc)
            {
                minConc = curConc;
            }
            ++it1;
        }

        if ( i % 100 == 0) std::cout<<"The maximum concentration is "<<maxConc << " and the minimum is " << minConc <<std::endl;

        // calculate the magnitude of velocity
        auto it2 = New.getDomainIterator();
        auto it2_vis = Vis_new.getDomainIterator();
        while (it2.isNext())
        {
            auto key = it2.get();
            auto key_vis = it2_vis.get();

            float curConc = (float) New.template get<V>(key);

            float scaled = (curConc / (maxConc - minConc)) * 65535;
            // copy
            Vis_new.get<0>(key) = (unsigned short)(scaled);

//                std::cout << key.to_string() << "lin " << loc_grid.get(i).get<0>(key) << "  " <<  lin.LinId(key);
//            auto &lin = g_vis.getGrid();

//			std::cout<<"Value at "<<key.to_string()<<" linearized as "<<g_vis.get<0>(key) << " " << lin.LinId(key);
//            std::cout<<"Value at "<<key.to_string()<<" is "<< (unsigned short)(scaled) <<std::endl;

            ++it2;
            ++it2_vis;
        }

//        if(i == (int)2000/deltaT || i == (int)3000/deltaT || i == (int)4000/deltaT )
//        {
//            Old.save("checkpoint" + std::to_string(i));
//        }
	}
	
	tot_sim.stop();
	std::cout << "Total simulation: " << tot_sim.getwct() << std::endl;

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
