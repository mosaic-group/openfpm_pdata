#include "Grid/grid_dist_id.hpp"
#include "data_type/aggregate.hpp"
#include "timer.hpp"

/*!
 *
 *
 */

constexpr int U = 0;
constexpr int V = 1;
constexpr int phi = 2;
constexpr int normal = 3;
constexpr int tgrad_u = 4;
constexpr int tgrad_v = 5;
constexpr int U_next = 6;
constexpr int V_next = 7;

constexpr int x = 0;
constexpr int y = 1;
constexpr int z = 2;

typedef sgrid_dist_id<3,double,aggregate<double,double,double,double[3],double[3],double[3],double,double> > sgrid_type;

void init(sgrid_type & grid, Box<3,double> & domain)
{
	//! \cond [init sphere channel] \endcond

	auto it = grid.getGridIterator();
	Point<3,double> p1({0.5,0.5,0.5});

	double sx = grid.spacing(0);

	Sphere<3,double> sph1(p1,0.3);
	Sphere<3,double> sph2(p1,0.3 - sx*10);
	Sphere<3,double> sph_zero(p1,0.3 - sx*5);

	while (it.isNext())
	{
		// Get the local grid key
		auto key = it.get_dist();
		auto keyg = it.get();

		Point<3,double> pc;
		Point<3,double> vp;

		for (int i = 0 ; i < 3 ; i++)
        {pc.get(i) = keyg.get(i) * it.getSpacing(i);}

		// Check if the point is in the first sphere
		if (sph1.isInside(pc) == true && sph2.isInside(pc) == false)
		{
			Point<3,double> pn = pc - p1;
			pn /= pn.norm();
			double theta = acos(pn * Point<3,double>({0.0,0.0,1.0}));
			Point<3,double> pn_ = pn;
			pn_[2] = 0.0;
			pn_ /= pn_.norm();
			double aphi = acos(pn_ * Point<3,double>({1.0,0.0,0.0}));

			// Create a perturbation in the solid angle
			if (theta > 0.6 && theta < 0.8 && aphi > 0.0 && aphi < 0.2)
			{
				grid.template insert<U>(key) = 0.5;
				grid.template insert<V>(key) = 0.25;
			}
			else
			{
				grid.template insert<U>(key) = 1.0;
				grid.template insert<V>(key) = 0.0;
			}
			grid.template insert<phi>(key) = sph_zero.distance(pc);
			grid.template insert<normal>(key)[0] = pn[0];
			grid.template insert<normal>(key)[1] = pn[1];
			grid.template insert<normal>(key)[2] = pn[2];

			// Old values U and V
			grid.template insert<U_next>(key) = 0.0;
			grid.template insert<V_next>(key) = 0.0;
		}

		++it;
	}

	//! \cond [init sphere channel] \endcond
}

template<unsigned int U_src,unsigned int V_src,unsigned int U_dst, unsigned int V_dst>
void extend(sgrid_type & grid)
{
	double delta = 1e-10;
	double max = 0.0;
	auto it = grid.getDomainIterator();

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

		double s = grid.get<phi>(Cp) / sqrt(fabs(grid.get<phi>(Cp)) + delta);

		double Uext = 0.0;
		double Vext = 0.0;

		double dir = s*grid.get<normal>(Cp)[x];

		if (dir > 0)
		{
			Uext += dir * (grid.get<U_src>(Cp) - grid.get<U_src>(mx));
			Vext += dir * (grid.get<V_src>(Cp) - grid.get<V_src>(mx));
		}
		else if (dir < 0)
		{
			Uext += dir * (grid.get<U_src>(px) - grid.get<U_src>(Cp));
			Vext += dir * (grid.get<V_src>(px) - grid.get<V_src>(Cp));
		}


		dir = s*grid.get<normal>(Cp)[y];
		if (dir > 0)
		{
			Uext += dir * (grid.get<U_src>(Cp) - grid.get<U_src>(my));
			Vext += dir * (grid.get<V_src>(Cp) - grid.get<V_src>(my));
		}
		else if (dir < 0)
		{
			Uext += dir * (grid.get<U_src>(py) - grid.get<U_src>(Cp));
			Vext += dir * (grid.get<V_src>(py) - grid.get<V_src>(Cp));
		}

		dir = s*grid.get<normal>(Cp)[z];
		if (dir > 0)
		{
			Uext += dir * (grid.get<U_src>(Cp) - grid.get<U_src>(mz));
			Vext += dir * (grid.get<V_src>(Cp) - grid.get<V_src>(mz));
		}
		else if (dir < 0)
		{
			Uext += dir * (grid.get<U_src>(pz) - grid.get<U_src>(Cp));
			Vext += dir * (grid.get<V_src>(pz) - grid.get<V_src>(Cp));
		}

		if (Uext >= max)
		{
			max = Uext;
		}

		grid.insert<U_dst>(Cp) = grid.get<U_src>(Cp) - 1.0*Uext;
		grid.insert<V_dst>(Cp) = grid.get<V_src>(Cp) - 1.0*Vext;

		// Next point in the grid
		++it;
	}

	std::cout << "UEX max: " << max << std::endl;
}

int main(int argc, char* argv[])
{
	openfpm_init(&argc,&argv);

	// domain
	Box<3,double> domain({0.0,0.0,0.0},{2.5,2.5,2.5});
	
	// grid size
    size_t sz[3] = {512,512,512};

	// Define periodicity of the grid
	periodicity<3> bc = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};
	
	// Ghost in grid unit
	Ghost<3,long int> g(1);
	
	// deltaT
	double deltaT = 0.3;

	// Diffusion constant for specie U
	double du = 1*1e-5;

	// Diffusion constant for specie V
	double dv = 0.5*1e-5;

//#ifdef TEST_RUN
        // Number of timesteps
//        size_t timeSteps = 200;
//#else
	// Number of timesteps
        size_t timeSteps = 150000;
//#endif

	// K and F (Physical constant in the equation)
        double K = 0.053;
        double F = 0.014;

	sgrid_type grid(sz,domain,g,bc);

	
	// spacing of the grid on x and y
	double spacing[3] = {grid.spacing(0),grid.spacing(1),grid.spacing(2)};

	init(grid,domain);

	// sync the ghost
	size_t count = 0;
	grid.template ghost_get<U,V>();

	// because we assume that spacing[x] == spacing[y] we use formula 2
	// and we calculate the prefactor of Eq 2
	double uFactor = deltaT * du;
	double vFactor = deltaT * dv;

	auto & v_cl = create_vcluster();

	timer tot_sim;
	tot_sim.start();

	for (size_t i = 0; i < timeSteps ; ++i)
	{
		{
		auto it = grid.getDomainIterator();

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

			grid.insert<tgrad_u>(Cp)[0] = 0.0;
			grid.insert<tgrad_u>(Cp)[1] = 0.0;
			grid.insert<tgrad_u>(Cp)[2] = 0.0;
			grid.insert<tgrad_v>(Cp)[0] = 0.0;
			grid.insert<tgrad_v>(Cp)[1] = 0.0;
			grid.insert<tgrad_v>(Cp)[2] = 0.0;

			//! \cond [boundary condition] \endcond

			if (grid.existPoint(mz) == true && grid.existPoint(pz) == true &&
				grid.existPoint(my) == true && grid.existPoint(py) == true &&
				grid.existPoint(mx) == true && grid.existPoint(px) == true )
			{
				Point<3,double> gradU;
				gradU[x] = (grid.get<U>(Cp) - grid.get<U>(mx)) / grid.spacing(0);
				gradU[y] = (grid.get<U>(Cp) - grid.get<U>(my)) / grid.spacing(1);
				gradU[z] = (grid.get<U>(Cp) - grid.get<U>(mz)) / grid.spacing(2);

				Point<3,double> gradV;
				gradV[x] = (grid.get<V>(Cp) - grid.get<V>(mx)) / grid.spacing(0);
				gradV[y] = (grid.get<V>(Cp) - grid.get<V>(my)) / grid.spacing(1);
				gradV[z] = (grid.get<V>(Cp) - grid.get<V>(mz)) / grid.spacing(2);

				Point<3,double> PgradU;
				Point<3,double> PgradV;

				PgradU.zero();
				PgradV.zero();

				for (int i = 0 ; i < 3 ; i++)
				{
					for (int j = 0 ; j < 3 ; j++)
					{
						grid.insert<tgrad_u>(Cp)[i] += (((i == j)?1.0:0.0) - grid.get<normal>(Cp)[i]*grid.get<normal>(Cp)[j])*gradU[j];
						grid.insert<tgrad_v>(Cp)[i] += (((i == j)?1.0:0.0) - grid.get<normal>(Cp)[i]*grid.get<normal>(Cp)[j])*gradV[j];
					}
				}
			}
			++it;
		}
		}

//		Old.write_frame("Init_condition",i);

		{
		auto it = grid.getDomainIterator();

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

			// Mirror z

			if (grid.existPoint(mz) == true && grid.existPoint(pz) == true &&
				grid.existPoint(my) == true && grid.existPoint(py) == true &&
				grid.existPoint(mx) == true && grid.existPoint(px) == true )
			{
				double lapU = 0;
				double lapV = 0;

				//Div
				lapU += (grid.get<tgrad_u>(px)[0] - grid.get<tgrad_u>(Cp)[0]) / grid.spacing(0);
				lapV += (grid.get<tgrad_v>(px)[0] - grid.get<tgrad_v>(Cp)[0]) / grid.spacing(0);
				lapU += (grid.get<tgrad_u>(py)[1] - grid.get<tgrad_u>(Cp)[1]) / grid.spacing(1);
				lapV += (grid.get<tgrad_v>(py)[1] - grid.get<tgrad_v>(Cp)[1]) / grid.spacing(1);
				lapU += (grid.get<tgrad_u>(pz)[2] - grid.get<tgrad_u>(Cp)[2]) / grid.spacing(2);
				lapV += (grid.get<tgrad_v>(pz)[2] - grid.get<tgrad_v>(Cp)[2]) / grid.spacing(2);

				// update based on Eq 2
				grid.insert<U_next>(Cp) = grid.get<U>(Cp) + uFactor * lapU +
											- deltaT * grid.get<U>(Cp) * grid.get<V>(Cp) * grid.get<V>(Cp) +
											- deltaT * F * (grid.get<U>(Cp) - 1.0);


				// update based on Eq 2
				grid.insert<V_next>(Cp) = grid.get<V>(Cp) + vFactor * lapV +
											deltaT * grid.get<U>(Cp) * grid.get<V>(Cp) * grid.get<V>(Cp) +
											- deltaT * (F+K) * grid.get<V>(Cp);
			}

			// Next point in the grid
			++it;
		}
		}

//		New.write_frame("update",i);

		// Extend

		if (i % 5 == 0)
		{
		for (int j = 0 ; j < 2 ; j++)
		{
			if (j % 2 == 0)
			{extend<U_next,V_next,U,V>(grid);}
			else
			{extend<U,V,U_next,V_next>(grid);}

			// Here we copy New into the old grid in preparation of the new step
			// It would be better to alternate, but using this we can show the usage
			// of the function copy. To note that copy work only on two grid of the same
			// decomposition. If you want to copy also the decomposition, or force to be
			// exactly the same, use Old = New
			//New.copy_sparse(Old);
		}
		}

/*		auto it = grid.getDomainIterator();

		while (it.isNext())
		{
			// center point
			auto Cp = it.get();

			// update based on Eq 2
			grid.insert<U>(Cp) = grid.get<U_next>(Cp);
			grid.insert<V>(Cp) = grid.get<V_next>(Cp);

			++it;
		}*/

		//! \cond [stencil get and use] \endcond

		// After copy we synchronize again the ghost part U and V
		grid.ghost_get<U,V>();

		// Every 500 time step we output the configuration for
		// visualization
		if (i % 500 == 0)
		{
			grid.save("output_" + std::to_string(count));
			count++;
		}

		if (v_cl.rank() == 0)
		{std::cout << "STEP: " << i  << "   " << std::endl;}
		if (i % 100 == 0)
		{
			grid.write_frame("out",i);
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
