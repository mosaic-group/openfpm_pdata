#include "Grid/grid_dist_id.hpp"
#include "data_type/aggregate.hpp"
#include "timer.hpp"

/*!
 *
 * \page Grid_3_gs_3D_sparse_cs Gray Scott in 3D using sparse grids with complex geometry
 *
 * [TOC]
 *
 * # Solving a gray scott-system in 3D using Sparse grids# {#e3_gs_gray_scott_cs}
 *
 * This example show how to solve a Gray-Scott system in 3D using sparse grids. In this example we will construct
 * a complex geometry
 *
 * In figure is the final solution of the problem
 *
 *
 * \htmlonly
<table border="1" bgcolor="black">
  <tr>
    <td>
      <img src="http://ppmcore.mpi-cbg.de/web/images/examples/1_gray_scott_3d_sparse_cs/gs_3d_sparse_cs_section.png" style="width: 500px;" />
    </td>
    <td>
      <img src="http://ppmcore.mpi-cbg.de/web/images/examples/1_gray_scott_3d_sparse_cs/gs_3d_sparse_cs.png" style="width: 500px;" />
    </td>
  </tr>
</table>
\endhtmlonly
 *
 *
 * This example follow the gray scott 3d sparse with the addition that the initialization is completely different to create a
 * complex geometry initialization
 *
 *
 * We recall here the main differences between sparse and dense.
 *
 * * Get function return now constant values, so cannot be used to get values, a get in write is an insert
 *   a get on a point position that has not been inserted yet return the background value
 *
 * * Insert function create/overwrite the points value
 *
 * * getDomainIterator return an iterator on the existing points
 *
 * * getGridIterator return an iterator on the dense version of the grid
 *
 * 
 * # Initialization
 *
 * The initialization involve the creation of 3 sphere and one cylinder channel connecting them in order to do it we
 * create an iterator over the grid (inserted and not inserted) point with **getGridIterator**
 *
 * \snippet SparseGrid/3_gray_scott_3d_sparse_cs/main.cpp init sphere channel
 *
 * After creating the domain we make a perturbation in the up sphere
 *
 * \snippet SparseGrid/3_gray_scott_3d_sparse_cs/main.cpp perturbation
 *
 * # Boundary conditions
 *
 * For this example we use mirror on direction X Y Z If the point is missing. If the point is missing in both direction than
 * the second derivative is considered zero
 *
 * \snippet SparseGrid/3_gray_scott_3d_sparse_cs/main.cpp boundary condition
 *
 */

constexpr int U = 0;
constexpr int V = 1;

constexpr int x = 0;
constexpr int y = 1;
constexpr int z = 2;

typedef sgrid_dist_id<3,double,aggregate<double,double,int> > sgrid_type;

void init(sgrid_type & Old, sgrid_type & New, Box<3,double> & domain)
{
	//! \cond [init sphere channel] \endcond

	auto it = Old.getGridIterator();
	Point<3,double> p1;
	Point<3,double> p2;
	Point<3,double> p3;

	// Shere1
	for (int i = 0 ; i < 3 ; i++)
	{p1.get(i) = 2.0;}

        // Shere2
        for (int i = 0 ; i < 3 ; i++)
        {p2.get(i) = 1.0;}

        // Shere3
        for (int i = 0 ; i < 3 ; i++)
        {p3.get(i) = 0.5;}

	Sphere<3,double> sph1(p1,0.3);
	Sphere<3,double> sph2(p2,0.3);
	Sphere<3,double> sph3(p3,0.3);

	Point<3,double> u({1.0,1.0,1.0});
	Box<3,double> channel_box(p3,p1);

	while (it.isNext())
	{
		// Get the local grid key
		auto key = it.get_dist();
		auto keyg = it.get();

		Point<3,double> pc;
		Point<3,double> vp;

		for (int i = 0 ; i < 3 ; i++)
        	{pc.get(i) = keyg.get(i) * it.getSpacing(i);}

                // calculate the distance from the diagonal
                vp.get(0) = pc.get(1)*u.get(2) - pc.get(2)*u.get(1);
                vp.get(1) = pc.get(2)*u.get(0) - pc.get(0)*u.get(2);
                vp.get(2) = pc.get(0)*u.get(1) - pc.get(1)*u.get(0);

                double distance = vp.norm() / sqrt(3);

		// Check if the point is in the first sphere
		if (sph1.isInside(pc) || sph2.isInside(pc) || sph3.isInside(pc) || (distance < 0.1 && channel_box.isInside(pc)) )
		{
			// Old values U and V
			Old.template insert<U>(key) = 1.0;
			Old.template insert<V>(key) = 0.0;

			// Old values U and V
			New.template insert<U>(key) = 0.0;
			New.template insert<V>(key) = 0.0;
		}

		++it;
	}

	//! \cond [init sphere channel] \endcond

	//! \cond [perturbation] \endcond

	long int x_start = Old.size(0)*1.95f/domain.getHigh(0);
	long int y_start = Old.size(1)*1.95f/domain.getHigh(1);
	long int z_start = Old.size(1)*1.95f/domain.getHigh(2);

	long int x_stop = Old.size(0)*2.05f/domain.getHigh(0);
	long int y_stop = Old.size(1)*2.05f/domain.getHigh(1);
	long int z_stop = Old.size(1)*2.05f/domain.getHigh(2);

	grid_key_dx<3> start({x_start,y_start,z_start});
	grid_key_dx<3> stop ({x_stop,y_stop,z_stop});
	auto it_init = Old.getGridIterator(start,stop);

	while (it_init.isNext())
	{
		auto key = it_init.get_dist();

                Old.template insert<U>(key) = 0.5 + (((double)std::rand())/RAND_MAX -0.5)/10.0;
                Old.template insert<V>(key) = 0.25 + (((double)std::rand())/RAND_MAX -0.5)/20.0;

		++it_init;
	}

	//! \cond [perturbation] \endcond
}


int main(int argc, char* argv[])
{
	openfpm_init(&argc,&argv);

	// domain
	Box<3,double> domain({0.0,0.0,0.0},{2.5,2.5,2.5});
	
	// grid size
        size_t sz[3] = {256,256,256};

	// Define periodicity of the grid
	periodicity<3> bc = {PERIODIC,PERIODIC,PERIODIC};
	
	// Ghost in grid unit
	Ghost<3,long int> g(1);
	
	// deltaT
	double deltaT = 0.1;

	// Diffusion constant for specie U
	double du = 2*1e-5;

	// Diffusion constant for specie V
	double dv = 1*1e-5;

#ifdef TEST_RUN
        // Number of timesteps
        size_t timeSteps = 200;
#else
	// Number of timesteps
        size_t timeSteps = 150000;
#endif

	// K and F (Physical constant in the equation)
        double K = 0.053;
        double F = 0.014;

	sgrid_type Old(sz,domain,g,bc);

	// New grid with the decomposition of the old grid
        sgrid_type New(Old.getDecomposition(),sz,g);

	
	// spacing of the grid on x and y
	double spacing[3] = {Old.spacing(0),Old.spacing(1),Old.spacing(2)};

	init(Old,New,domain);

	// sync the ghost
	size_t count = 0;
	Old.template ghost_get<U,V>();

	// because we assume that spacing[x] == spacing[y] we use formula 2
	// and we calculate the prefactor of Eq 2
	double uFactor = deltaT * du/(spacing[x]*spacing[x]);
	double vFactor = deltaT * dv/(spacing[x]*spacing[x]);

	auto & v_cl = create_vcluster();

	Old.write("Init_condition");

	timer tot_sim;
	tot_sim.start();

	for (size_t i = 0; i < timeSteps; ++i)
	{
		auto it = Old.getDomainIterator();

		while (it.isNext())
		{
			// center point
			auto Cp = it.get();

			// plus,minus X,Y,Z
			auto mx = Cp.move(0,-1);
			auto px = Cp.move(0,+1);
			auto my = Cp.move(1,-1);
			auto py = Cp.move(1,1);
			auto mz = Cp.move(2,-1);
			auto pz = Cp.move(2,1);

			//! \cond [boundary condition] \endcond

			double LapU = 0.0;
			double LapV = 0.0;

			// Mirror z

			if (Old.existPoint(mz) == true)
			{
				LapU += Old.get<U>(mz);
				LapV += Old.get<V>(mz);
			}
			else if (Old.existPoint(pz) == true)
			{
				LapU += Old.get<U>(pz);
				LapV += Old.get<V>(pz);
			}
			else
			{
				LapU += Old.get<U>(Cp);
				LapV += Old.get<V>(Cp);
			}

            if (Old.existPoint(pz) == true)
            {
				LapU += Old.get<U>(pz);
				LapV += Old.get<V>(pz);
			}
            else if (Old.existPoint(mz) == true)
            {
				LapU+= Old.get<U>(mz);
				LapV += Old.get<V>(mz);
			}
            else
            {
				LapU+= Old.get<U>(Cp);
				LapV += Old.get<V>(Cp);
            }


			// Mirror y

			if (Old.existPoint(my) == true)
			{
					LapU += Old.get<U>(my);
					LapV += Old.get<V>(my);
			}
			else if (Old.existPoint(py) == true)
			{
					LapU += Old.get<U>(py);
					LapV += Old.get<V>(py);
			}
            else
            {
				LapU+= Old.get<U>(Cp);
				LapV += Old.get<V>(Cp);
            }

			if (Old.existPoint(py) == true)
			{
					LapU += Old.get<U>(py);
					LapV += Old.get<V>(py);
			}
			else if (Old.existPoint(my) == true)
			{
					LapU+= Old.get<U>(my);
					LapV += Old.get<V>(my);
			}
			else
			{
				LapU+= Old.get<U>(Cp);
				LapV += Old.get<V>(Cp);
			}

			// Mirror x

			if (Old.existPoint(mx) == true)
			{
					LapU += Old.get<U>(mx);
					LapV += Old.get<V>(mx);
			}
			else if (Old.existPoint(px) == true)
			{
					LapU += Old.get<U>(px);
					LapV += Old.get<V>(px);
			}
			else
			{
				LapU+= Old.get<U>(Cp);
				LapV += Old.get<V>(Cp);
			}

			if (Old.existPoint(px) == true)
			{
					LapU += Old.get<U>(px);
					LapV += Old.get<V>(px);
			}
			else if (Old.existPoint(mx) == true)
			{
					LapU+= Old.get<U>(mx);
					LapV += Old.get<V>(mx);
			}
			else
			{
				LapU+= Old.get<U>(Cp);
				LapV += Old.get<V>(Cp);
			}

			LapU -= 6.0*Old.get<U>(Cp);
			LapV -= 6.0*Old.get<V>(Cp);

			//! \cond [boundary condition] \endcond

			// update based on Eq 2
			New.insert<U>(Cp) = Old.get<U>(Cp) + uFactor * LapU  +
										- deltaT * Old.get<U>(Cp) * Old.get<V>(Cp) * Old.get<V>(Cp) +
										- deltaT * F * (Old.get<U>(Cp) - 1.0);


			// update based on Eq 2
			New.insert<V>(Cp) = Old.get<V>(Cp) + vFactor * LapV  +
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
		Old.copy_sparse(New);

		// After copy we synchronize again the ghost part U and V
		Old.ghost_get<U,V>();

		// Every 500 time step we output the configuration for
		// visualization
		if (i % 500 == 0)
		{
			Old.save("output_" + std::to_string(count));
			count++;
		}

		if (v_cl.rank() == 0)
		{std::cout << "STEP: " << i  << "   " << std::endl;}
                if (i % 300 == 0)
                {
                        Old.write_frame("out",i);
                }
	}
	
	tot_sim.stop();
	std::cout << "Total simulation: " << tot_sim.getwct() << std::endl;

	//! \cond [time stepping] \endcond

	/*!
	 * \page Grid_3_gs_3D_sparse Gray Scott in 3D
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
	 * \page Grid_3_gs_3D_sparse Gray Scott in 3D
	 *
	 * # Full code # {#code}
	 *
	 * \include Grid/3_gray_scott_3d/main.cpp
	 *
	 */
}
