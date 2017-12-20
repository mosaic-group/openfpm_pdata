#include "Grid/grid_dist_id.hpp"
#include "data_type/aggregate.hpp"
#include "timer.hpp"
#include "Vc/Vc"

/*!
 *
 * \page Grid_3_gs_3D_vector Gray Scott in 3D fast implementation with vectorization
 *
 * # Solving a gray scott-system in 3D # {#e3_gs_gray_scott_vector}
 *
 * This example is just an improved version of the previous 3D Gray scott example.
 * In particular we do the following improvements we separate U and V in two grids
 * in order to vectorize. Every loop now handle 4 double in case of AVX-256 and 2 double
 * in case of SSE. We also avoid to use the function copy and we alternate the use of the
 * fields New and Old. If at the first iteration we read from Old and we write on New in
 * the second iteration we read from New and we write on Old. The last improvement is write
 * on hdf5 rather that VTK. VTK writers are convenient but are slow for performances. HDF5
 * files can be saved with **save()** reload with **load()** and after loading can be written
 * on VTK with **write** this mean that HDF5 files can be easily converted into VTK in a second moment.
 * Not only but because HDF5 files can be saved on multiple processors and reloaded on a different
 * number of processors, you can use this method to stitch VTK files together.
 *
 *
 * In figure is the final solution of the problem
 *
 * \htmlonly
 * <img src="http://ppmcore.mpi-cbg.de/web/images/examples/gray_scott_3d/gs_alpha.png"/>
 * \endhtmlonly
 *
 * \see \ref Grid_2_solve_eq
 *
 * \snippet Grid/3_gray_scott_3d_vectorization/main.cpp constants
 * 
 */

//! \cond [constants] \endcond

//#define FORTRAN_UPDATE

constexpr int x = 0;
constexpr int y = 1;
constexpr int z = 2;

extern "C" void update_new(const int* lo, const int* hi,
                  double* u, const int* ulo, const int* uhi,
                  double* v, const int* vlo, const int* vhi,
                  double* flu, const int* fulo, const int* fuhi,
                  double* flv, const int* fvlo, const int* fvhi,
                  const double * dt, const double * uFactor, const double * vFactor, const double * F,
                  const double * K);


//! \cond [constants] \endcond

void init(grid_dist_id<3,double,aggregate<double> > & OldU,
		  grid_dist_id<3,double,aggregate<double> > & OldV,
		  grid_dist_id<3,double,aggregate<double> > & NewU,
		  grid_dist_id<3,double,aggregate<double> > & NewV,
		  Box<3,double> & domain)
{
	auto it = OldU.getDomainIterator();

	while (it.isNext())
	{
		// Get the local grid key
		auto key = it.get();

		// Old values U and V
		OldU.get(key) = 1.0;
		OldV.get(key) = 0.0;

		// Old values U and V
		NewU.get(key) = 0.0;
		NewV.get(key) = 0.0;

		++it;
	}

	long int x_start = OldU.size(0)*1.55f/domain.getHigh(0);
	long int y_start = OldU.size(1)*1.55f/domain.getHigh(1);
	long int z_start = OldU.size(1)*1.55f/domain.getHigh(2);

	long int x_stop = OldU.size(0)*1.85f/domain.getHigh(0);
	long int y_stop = OldU.size(1)*1.85f/domain.getHigh(1);
	long int z_stop = OldU.size(1)*1.85f/domain.getHigh(2);

	grid_key_dx<3> start({x_start,y_start,z_start});
	grid_key_dx<3> stop ({x_stop,y_stop,z_stop});
	auto it_init = OldU.getSubDomainIterator(start,stop);

	while (it_init.isNext())
	{
		auto key = it_init.get();

        OldU.get(key) = 0.5 + (((double)std::rand())/RAND_MAX -0.5)/10.0;
        OldV.get(key) = 0.25 + (((double)std::rand())/RAND_MAX -0.5)/20.0;

		++it_init;
	}
}


//! \cond [vectorization] \endcond

void step(grid_dist_id<3, double, aggregate<double>> & OldU,
		  grid_dist_id<3, double, aggregate<double>> & OldV,
		  grid_dist_id<3, double, aggregate<double>> & NewU,
		  grid_dist_id<3, double, aggregate<double>> & NewV,
		  grid_key_dx<3> (& star_stencil_3D)[7],
		  Vc::double_v uFactor, Vc::double_v vFactor, double deltaT, double F, double K)
{
#ifndef FORTRAN_UPDATE

	//! \cond [cpp_update] \endcond

	WHILE_M(OldU,star_stencil_3D)
			auto & U_old = GET_GRID_M(OldU);
			auto & V_old = GET_GRID_M(OldV);

			auto & U_new = GET_GRID_M(NewU);
			auto & V_new = GET_GRID_M(NewV);
	ITERATE_3D_M

			// center point
			auto Cp = it.getStencil<0>();

			// plus,minus X,Y,Z
			auto mx = it.getStencil<1>();
			auto px = it.getStencil<2>();
			auto my = it.getStencil<3>();
			auto py = it.getStencil<4>();
			auto mz = it.getStencil<5>();
			auto pz = it.getStencil<6>();

			//
			Vc::double_v u_c(&U_old.get<0>(Cp),Vc::Unaligned);
			Vc::double_v u_mz(&U_old.get<0>(mz),Vc::Unaligned);
			Vc::double_v u_pz(&U_old.get<0>(pz),Vc::Unaligned);
			Vc::double_v u_my(&U_old.get<0>(my),Vc::Unaligned);
			Vc::double_v u_py(&U_old.get<0>(py),Vc::Unaligned);
			Vc::double_v u_mx(&U_old.get<0>(mx),Vc::Unaligned);
			Vc::double_v u_px(&U_old.get<0>(px),Vc::Unaligned);


			Vc::double_v v_c(&V_old.get<0>(Cp),Vc::Unaligned);
			Vc::double_v v_mz(&V_old.get<0>(mz),Vc::Unaligned);
			Vc::double_v v_pz(&V_old.get<0>(pz),Vc::Unaligned);
			Vc::double_v v_my(&V_old.get<0>(my),Vc::Unaligned);
			Vc::double_v v_py(&V_old.get<0>(py),Vc::Unaligned);
			Vc::double_v v_mx(&V_old.get<0>(mx),Vc::Unaligned);
			Vc::double_v v_px(&V_old.get<0>(px),Vc::Unaligned);

			Vc::double_v out1 = u_c + uFactor * (u_mz + u_pz +
					                            u_my + u_py +
												u_mx + u_px +
											    - 6.0 * u_c) +
												- deltaT * u_c * v_c * v_c
												- deltaT * F * (u_c - 1.0);

			Vc::double_v out2 = v_c + vFactor * (v_mz + v_pz +
					               v_my + v_py +
								   v_mx + v_px +
								   - 6.0 * v_c ) +
								   deltaT * u_c * v_c * v_c +
								   - deltaT * (F+K) * v_c;

			out1.store(&U_new.get<0>(Cp),Vc::Unaligned);
			out2.store(&V_new.get<0>(Cp),Vc::Unaligned);
	END_LOOP_M

	//! \cond [cpp_update] \endcond

#else

	//! \cond [fort_update] \endcond

	double uFactor_s = uFactor[0];
	double vFactor_s = vFactor[0];

	auto & ginfo = OldU.getLocalGridsInfo();

	for (size_t i = 0 ; i < OldU.getN_loc_grid() ; i++)
	{
		auto & U_old = OldU.get_loc_grid(i);
		auto & V_old = OldV.get_loc_grid(i);

		auto & U_new = NewU.get_loc_grid(i);
		auto & V_new = NewV.get_loc_grid(i);

		int lo[3] = {(int)ginfo.get(i).Dbox.getLow(0),(int)ginfo.get(i).Dbox.getLow(1),(int)ginfo.get(i).Dbox.getLow(2)};
		int hi[3] = {(int)ginfo.get(i).Dbox.getHigh(0),(int)ginfo.get(i).Dbox.getHigh(1),(int)ginfo.get(i).Dbox.getHigh(2)};

		int ulo[3] = {0,0,0};
		int uhi[3] = {(int)ginfo.get(i).GDbox.getHigh(0),(int)ginfo.get(i).GDbox.getHigh(1),(int)ginfo.get(i).GDbox.getHigh(2)};
		int nulo[3] = {0,0,0};
		int nuhi[3] = {(int)ginfo.get(i).GDbox.getHigh(0),(int)ginfo.get(i).GDbox.getHigh(1),(int)ginfo.get(i).GDbox.getHigh(2)};

		int vlo[3] = {0,0,0};
		int vhi[3] = {(int)ginfo.get(i).GDbox.getHigh(0),(int)ginfo.get(i).GDbox.getHigh(1),(int)ginfo.get(i).GDbox.getHigh(2)};
		int nvlo[3] = {0,0,0};
		int nvhi[3] = {(int)ginfo.get(i).GDbox.getHigh(0),(int)ginfo.get(i).GDbox.getHigh(1),(int)ginfo.get(i).GDbox.getHigh(2)};

		update_new(lo,hi,
				   (double *)U_old.getPointer(),ulo,uhi,
				   (double *)V_old.getPointer(),vlo,vhi,
				   (double *)U_new.getPointer(),nulo,nuhi,
				   (double *)V_new.getPointer(),nulo,nvhi,
				   &deltaT, &uFactor_s, &vFactor_s,&F,&K);
	}

	//! \cond [fort_update] \endcond

#endif
}

//! \cond [vectorization] \endcond

int main(int argc, char* argv[])
{
	openfpm_init(&argc,&argv);

	// domain
	Box<3,double> domain({0.0,0.0},{2.5,2.5,2.5});
	
	// grid size
        size_t sz[3] = {128,128,128};

	// Define periodicity of the grid
	periodicity<3> bc = {PERIODIC,PERIODIC,PERIODIC};
	
	// Ghost in grid unit
	Ghost<3,long int> g(1);
	
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

	//! \cond [init lib] \endcond

	/*!
	 * \page Grid_3_gs_3D_vector Gray Scott in 3D fast implementation with vectorization
	 *
	 * Here we create 2 distributed grid in 3D Old and New splitting U and V in two different fields.
	 *  In particular because we want that all the grids are distributed across processors in the same
	 *   way we pass the decomposition of the first grid.
	 *
	 * \snippet Grid/3_gray_scott_3d_vectorization/main.cpp init grid
	 *
	 */

	//! \cond [init grid] \endcond

	grid_dist_id<3, double, aggregate<double>> OldU(sz,domain,g,bc);
	grid_dist_id<3, double, aggregate<double>> OldV(OldU.getDecomposition(),sz,g);

	// New grid with the decomposition of the old grid
    grid_dist_id<3, double, aggregate<double>> NewU(OldU.getDecomposition(),sz,g);
    grid_dist_id<3, double, aggregate<double>> NewV(OldV.getDecomposition(),sz,g);

	// spacing of the grid on x and y

	double spacing[3] = {OldU.spacing(0),OldU.spacing(1),OldU.spacing(2)};

	init(OldU,OldV,NewU,NewV,domain);

	//! \cond [init grid] \endcond

	// sync the ghost
	size_t count = 0;
	OldU.template ghost_get<0>();
	OldV.template ghost_get<0>();

	// because we assume that spacing[x] == spacing[y] we use formula 2
	// and we calculate the prefactor of Eq 2
	Vc::double_v uFactor = deltaT * du/(spacing[x]*spacing[x]);
	Vc::double_v vFactor = deltaT * dv/(spacing[x]*spacing[x]);

	timer tot_sim;
	tot_sim.start();

	static grid_key_dx<3> star_stencil_3D[7] = {{0,0,0},
                                         	    {0,0,-1},
						    {0,0,1},
						    {0,-1,0},
						    {0,1,0},
						    {-1,0,0},
						    {1,0,0}};

	for (size_t i = 0; i < timeSteps; ++i)
	{
		if (i % 300 == 0)
			std::cout << "STEP: " << i << std::endl;

		/*!
		 * \page Grid_3_gs_3D_vector Gray Scott in 3D fast implementation with vectorization
		 *
		 * Alternate New and Old field to run one step, switch between old and new if the iteration
		 * is even or odd. The function step is nothing else than the implementation of Gray-Scott
		 * 3D in the previous example but in a more optimized way.
		 *
		 * \snippet Grid/3_gray_scott_3d_vectorization/main.cpp alternate
		 *
		 * In this function we show two methods to optimize this function.
		 *
		 * * We can use the macro **WHILE_M** passing the stencil definition, **ITERATE_3D** to define the loop,
		 *  **END_LOOP** to close the loop, and use the function
		 * function **getStencil<0>()** to retrieve the stencil points. Additionaly we can use Vc::double_v instead
		 *  of double to vectorize the code. This method give the advantage to keep all the
		 * code in C++.
		 *
		 * \snippet Grid/3_gray_scott_3d_vectorization/main.cpp cpp_update
		 *
		 * * Another possibility is to use FORTRAN. Because FORTRAN has better
		 *  support for multi dimensional array another possibility is to process each local grid using
		 *  FORTRAN, this also give us the opportunity to show hybrid code. We can switch between
		 *   one and the other method commenting
		 *  and uncommeting the line #define FORTRAN_UPDATE in the code.
		 *
		 * \snippet Grid/3_gray_scott_3d_vectorization/main.cpp fort_update
		 *
		 * \include Grid/3_gray_scott_3d_vectorization/update_new.f90
		 *
		 */

		//! \cond [alternate] \endcond

		if (i % 2 == 0)
		{
			step(OldU,OldV,NewU,NewV,star_stencil_3D,uFactor,vFactor,deltaT,F,K);

			NewU.ghost_get<0>();
			NewV.ghost_get<0>();
		}
		else
		{
			step(NewU,NewV,OldU,OldV,star_stencil_3D,uFactor,vFactor,deltaT,F,K);

			OldU.ghost_get<0>();
			OldV.ghost_get<0>();
		}

		//! \cond [alternate] \endcond

		/*!
		 * \page Grid_3_gs_3D_vector Gray Scott in 3D fast implementation with vectorization
		 *
		 * Instead of using the function **write** we use the function **save** to save on HDF5
		 *
		 * \snippet Grid/3_gray_scott_3d_vectorization/main.cpp save hdf5
		 *
		 */

		//! \cond [save hdf5] \endcond

		// Every 2000 time step we output the configuration on hdf5
		if (i % 2000 == 0)
		{
			OldU.save("output_u_" + std::to_string(count));
			OldV.save("output_v_" + std::to_string(count));
			count++;
		}

		//! \cond [save hdf5] \endcond
	}
	
	tot_sim.stop();

	if (create_vcluster().rank() == 0)
	{std::cout << "Total simulation: " << tot_sim.getwct() << std::endl;}

	// We frite the final configuration
	OldV.write("final");

	//! \cond [time stepping] \endcond

	/*!
	 * \page Grid_3_gs_3D_vector Gray Scott in 3D fast implementation with vectorization
	 *
	 * ## Finalize ##
	 *
	 * Deinitialize the library
	 *
	 * \snippet Grid/3_gray_scott_3d_vectorization/main.cpp finalize
	 *
	 */

	//! \cond [finalize] \endcond

	openfpm_finalize();

	//! \cond [finalize] \endcond

	/*!
	 * \page Grid_3_gs_3D_vector Gray Scott in 3D fast implementation with vectorization
	 *
	 * # Full code # {#code}
	 *
	 * \include Grid/3_gray_scott_3d_vectorization/main.cpp
	 *
	 */
}


