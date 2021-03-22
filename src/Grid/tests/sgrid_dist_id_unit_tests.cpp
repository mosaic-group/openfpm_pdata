/*
 * sgrid_dist_id_unit_tests.cpp
 *
 *  Created on: Nov 18, 2017
 *      Author: i-bird
 */


#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "Grid/grid_dist_id.hpp"
#include "Point_test.hpp"


const int x = 0;
const int y = 1;
const int z = 2;

BOOST_AUTO_TEST_SUITE( sgrid_dist_id_test )



BOOST_AUTO_TEST_CASE (sgrid_dist_id_soa )
{
	periodicity<3> bc = {PERIODIC, PERIODIC, PERIODIC};

	// Domain
	Box<3,double> domain({-0.3,-0.3,-0.3},{1.0,1.0,1.0});

	// grid size
	size_t sz[3];
	sz[0] = 1024;
	sz[1] = 1024;
	sz[2] = 1024;

	// Ghost
	Ghost<3,double> g(0.01);

	sgrid_dist_soa<3,double,Point_test<float>> sg(sz,domain,g,bc);

	// create a grid iterator over a bilion point

	auto it = sg.getGridIterator();

	while(it.isNext())
	{
		auto gkey = it.get();
		auto key = it.get_dist();

		size_t sx = gkey.get(0) - 512;
		size_t sy = gkey.get(1) - 512;
		size_t sz = gkey.get(2) - 512;

		if (sx*sx + sy*sy + sz*sz < 128*128)
		{
			sg.template insert<0>(key) = 1.0;
		}

		++it;
	}

	bool match = true;
	auto it2 = sg.getGridIterator();

	while(it2.isNext())
	{
		auto gkey = it2.get();
		auto key = it2.get_dist();

		size_t sx = gkey.get(0) - 512;
		size_t sy = gkey.get(1) - 512;
		size_t sz = gkey.get(2) - 512;

		if (sx*sx + sy*sy + sz*sz < 128*128)
		{
			match &= (sg.template get<0>(key) == 1.0);
		}

		++it2;
	}

	auto & gr = sg.getGridInfo();

	auto it3 = sg.getDomainIterator();

	while (it3.isNext())
	{
		auto key = it3.get();
		auto gkey = it3.getGKey(key);

		sg.template insert<0>(key) = gkey.get(0)*gkey.get(0) + gkey.get(1)*gkey.get(1) + gkey.get(2)*gkey.get(2);

		++it3;
	}

	sg.ghost_get<0>();
	// now we check the stencil

	bool good = true;
	auto it4 = sg.getDomainIterator();

	while (it4.isNext())
	{
		auto key = it4.get();
		auto gkey = it4.getGKey(key);

		size_t sx = gkey.get(0) - 512;
		size_t sy = gkey.get(1) - 512;
		size_t sz = gkey.get(2) - 512;

		double lap;

		lap = sg.template get<0>(key.move(x,1)) + sg.template get<0>(key.move(x,-1)) +
			  sg.template get<0>(key.move(y,1)) + sg.template get<0>(key.move(y,-1)) +
			  sg.template get<0>(key.move(z,1)) + sg.template get<0>(key.move(z,-1)) -
			  6.0*sg.template get<0>(key);

		good &= (lap == 6.0);

		++it4;
	}

	BOOST_REQUIRE_EQUAL(match,true);
}

BOOST_AUTO_TEST_CASE( sgrid_dist_id_basic_test_2D)
{
	periodicity<2> bc = {NON_PERIODIC, NON_PERIODIC};

	// Domain
	Box<2,double> domain({-0.3,-0.3},{1.0,1.0});

	// grid size
	size_t sz[2];
	sz[0] = 1024;
	sz[1] = 1024;

	// Ghost
	Ghost<2,double> g(0.01);

	sgrid_dist_id<2,double,Point_test<double>> sg(sz,domain,g,bc);

	// create a grid iterator

	auto it = sg.getGridIterator();

	while(it.isNext())
	{
		auto gkey = it.get();
		auto key = it.get_dist();


		long int sx = gkey.get(0) - 512;
		long int sy = gkey.get(1) - 512;

		if (sx*sx + sy*sy < 128*128)
		{
			sg.template insert<0>(key) = 1.0;
		}

		++it;
	}

	bool match = true;
	auto it2 = sg.getGridIterator();

	while(it2.isNext())
	{
		auto gkey = it2.get();
		auto key = it2.get_dist();

		long int sx = gkey.get(0) - 512;
		long int sy = gkey.get(1) - 512;

		if (sx*sx + sy*sy < 128*128)
		{
			match &= (sg.template get<0>(key) == 1.0);
		}

		++it2;
	}

	BOOST_REQUIRE_EQUAL(match,true);

	auto & gr = sg.getGridInfo();

	auto it3 = sg.getDomainIterator();

	while (it3.isNext())
	{
		auto key = it3.get();
		auto gkey = it3.getGKey(key);

		sg.template insert<0>(key) = gkey.get(0)*gkey.get(0) + gkey.get(1)*gkey.get(1);

		++it3;
	}

	sg.ghost_get<0>();

	// now we check the stencil

	bool good = true;
	auto it4 = sg.getDomainIterator();

	while (it4.isNext())
	{
		auto key = it4.get();
		auto gkey = it4.getGKey(key);

		double lap;

		// Here we check that all point of the stencil are inside*/

		long int sx = gkey.get(0) - 512;
		long int sy = gkey.get(1) - 512;

		if (sx*sx + sy*sy < 126*126)
		{
			lap = sg.template get<0>(key.move(x,1)) + sg.template get<0>(key.move(x,-1)) +
				  sg.template get<0>(key.move(y,1)) + sg.template get<0>(key.move(y,-1)) -
				  4.0*sg.template get<0>(key);

			good &= (lap == 4.0);
		}

		++it4;
	}

	BOOST_REQUIRE_EQUAL(good,true);
}

BOOST_AUTO_TEST_CASE( sgrid_dist_id_basic_test)
{
	periodicity<3> bc = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};

	// Domain
	Box<3,double> domain({-0.3,-0.3,-0.3},{1.0,1.0,1.0});

	// grid size
	size_t sz[3];
	sz[0] = 1024;
	sz[1] = 1024;
	sz[2] = 1024;

	// Ghost
	Ghost<3,double> g(0.01);

	sgrid_dist_id<3,double,Point_test<float>> sg(sz,domain,g,bc);

	// create a grid iterator over a bilion point

	auto it = sg.getGridIterator();

	while(it.isNext())
	{
		auto gkey = it.get();
		auto key = it.get_dist();

		size_t sx = gkey.get(0) - 512;
		size_t sy = gkey.get(1) - 512;
		size_t sz = gkey.get(2) - 512;

		if (sx*sx + sy*sy + sz*sz < 128*128)
		{
			sg.template insert<0>(key) = 1.0;
		}

		++it;
	}

	bool match = true;
	auto it2 = sg.getGridIterator();

	while(it2.isNext())
	{
		auto gkey = it2.get();
		auto key = it2.get_dist();

		size_t sx = gkey.get(0) - 512;
		size_t sy = gkey.get(1) - 512;
		size_t sz = gkey.get(2) - 512;

		if (sx*sx + sy*sy + sz*sz < 128*128)
		{
			match &= (sg.template get<0>(key) == 1.0);
		}

		++it2;
	}

	auto & gr = sg.getGridInfo();

	auto it3 = sg.getDomainIterator();

	while (it3.isNext())
	{
		auto key = it3.get();
		auto gkey = it3.getGKey(key);

		sg.template insert<0>(key) = gkey.get(0)*gkey.get(0) + gkey.get(1)*gkey.get(1) + gkey.get(2)*gkey.get(2);

		++it3;
	}

	sg.ghost_get<0>();
	// now we check the stencil

	bool good = true;
	auto it4 = sg.getDomainIterator();

	while (it4.isNext())
	{
		auto key = it4.get();
		auto gkey = it4.getGKey(key);

		size_t sx = gkey.get(0) - 512;
		size_t sy = gkey.get(1) - 512;
		size_t sz = gkey.get(2) - 512;

		if (sx*sx + sy*sy + sz*sz < 125*125)
		{
			double lap;

			lap = sg.template get<0>(key.move(x,1)) + sg.template get<0>(key.move(x,-1)) +
				  sg.template get<0>(key.move(y,1)) + sg.template get<0>(key.move(y,-1)) +
				  sg.template get<0>(key.move(z,1)) + sg.template get<0>(key.move(z,-1)) -
				  6.0*sg.template get<0>(key);

			good &= (lap == 6.0);
		}

		++it4;
	}

	BOOST_REQUIRE_EQUAL(match,true);
}


BOOST_AUTO_TEST_CASE( sparse_grid_fast_stencil_vectorized_simplified_conv2)
{
	constexpr int U = 0;
	constexpr int V = 1;

	constexpr int U_next = 2;
	constexpr int V_next = 3;

	constexpr int x = 0;
	constexpr int y = 1;
	constexpr int z = 2;

    Box<3,double> domain({0.0,0.0,0.0},{2.5,2.5,2.5});

    // grid size
    size_t sz[3] = {32,32,32};

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

    sgrid_dist_soa<3, double, aggregate<double,double,double,double>> grid(sz,domain,g,bc);

    auto it = grid.getGridIterator();

    while (it.isNext())
    {
            // Get the local grid key
            auto key = it.get_dist();

            // Old values U and V
            grid.template insert<U>(key) = 1.0;
            grid.template insert<V>(key) = 0.0;

            // Old values U and V
            grid.template insert<U_next>(key) = 0.0;
            grid.template insert<V_next>(key) = 0.0;

            ++it;
    }

    long int x_start = grid.size(0)*1.55f/domain.getHigh(0);
    long int y_start = grid.size(1)*1.55f/domain.getHigh(1);
    long int z_start = grid.size(1)*1.55f/domain.getHigh(2);

    long int x_stop = grid.size(0)*1.85f/domain.getHigh(0);
    long int y_stop = grid.size(1)*1.85f/domain.getHigh(1);
    long int z_stop = grid.size(1)*1.85f/domain.getHigh(2);

    grid_key_dx<3> start({x_start,y_start,z_start});
    grid_key_dx<3> stop ({x_stop,y_stop,z_stop});
    auto it_init = grid.getGridIterator(start,stop);

    while (it_init.isNext())
    {
            auto key = it_init.get_dist();

            grid.template insert<U>(key) = 0.5 + (((double)std::rand())/RAND_MAX -0.5)/10.0;
            grid.template insert<V>(key) = 0.25 + (((double)std::rand())/RAND_MAX -0.5)/20.0;

            ++it_init;
    }

    // spacing of the grid on x and y
    double spacing[3] = {grid.spacing(0),grid.spacing(1),grid.spacing(2)};
    // sync the ghost
    size_t count = 0;
    grid.template ghost_get<U,V>();

    // because we assume that spacing[x] == spacing[y] we use formula 2
    // and we calculate the prefactor of Eq 2
    double uFactor = deltaT * du/(spacing[x]*spacing[x]);
    double vFactor = deltaT * dv/(spacing[x]*spacing[x]);

    int stencil[6][3] = {{1,0,0},{-1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};


     //! \cond [stencil get and use] \endcond


    auto func = [uFactor,vFactor,deltaT,F,K](Vc::double_v & u_out,Vc::double_v & v_out,
                                                                Vc::double_v (& u)[7],Vc::double_v (& v)[7],
                                                                unsigned char * mask){

																													 u_out = u[0] + uFactor *(u[1] + u[2] +
																																								  u[3] + u[4] +
																																								  u[5] + u[6] - 6.0*u[0]) - deltaT * u[0]*v[0]*v[0]
																																								- deltaT * F * (u[0] - 1.0);

																													 v_out = v[0] + vFactor *(v[1] + v[2] +
																																								  v[3] + v[4] +
																																								  v[5] + v[6] - 6.0*v[0]) + deltaT * u[0]*v[0]*v[0]
																																								- deltaT * (F+K) * v[0];
                                                                                     };

    grid.conv2<U,V,U_next,V_next,1>(stencil,{0,0,0},{(long int)sz[0]-1,(long int)sz[1]-1,(long int)sz[2]-1},func);
    grid.conv2<U,V,U_next,V_next,1>(stencil,{0,0,0},{(long int)sz[0]-1,(long int)sz[1]-1,(long int)sz[2]-1},func);

    bool match = true;

    {
		auto it = grid.getDomainIterator();

		double max_U = 0.0;
		double max_V = 0.0;
		grid_dist_key_dx<3> k_max;
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

			// update based on Eq 2
			if ( fabs(grid.get<U>(Cp) + uFactor * (
																	grid.get<U>(mz) +
																	grid.get<U>(pz) +
																	grid.get<U>(my) +
																	grid.get<U>(py) +
																	grid.get<U>(mx) +
																	grid.get<U>(px) -
																	6.0*grid.get<U>(Cp)) +
																	- deltaT * grid.get<U>(Cp) * grid.get<V>(Cp) * grid.get<V>(Cp) +
																	- deltaT * F * (grid.get<U>(Cp) - 1.0) - grid.get<U_next>(Cp)) > 0.000000001 )
			{
				match = false;
				break;
			}

			// update based on Eq 2
			if ( fabs(grid.get<V>(Cp) + vFactor * (
																	grid.get<V>(mz) +
																	grid.get<V>(pz) +
																	grid.get<V>(my) +
																	grid.get<V>(py) +
																	grid.get<V>(mx) +
																	grid.get<V>(px) -
																	6*grid.get<V>(Cp)) +
																	deltaT * grid.get<U>(Cp) * grid.get<V>(Cp) * grid.get<V>(Cp) +
																	- deltaT * (F+K) * grid.get<V>(Cp) - grid.get<V_next>(Cp)) > 0.000000001 )
			{
				match = false;
				break;
			}

			++it;
		}
    }

    BOOST_REQUIRE_EQUAL(match,true);
}


BOOST_AUTO_TEST_CASE( sparse_grid_fast_stencil_vectorized_simplified_conv2_crossing)
{
	constexpr int U = 0;
	constexpr int V = 1;

	constexpr int U_next = 2;
	constexpr int V_next = 3;

	constexpr int x = 0;
	constexpr int y = 1;
	constexpr int z = 2;

    Box<3,double> domain({0.0,0.0,0.0},{2.5,2.5,2.5});

    // grid size
    size_t sz[3] = {32,32,32};

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

    sgrid_dist_soa<3, double, aggregate<double,double,double,double>> grid(sz,domain,g,bc);

    auto it = grid.getGridIterator();

    while (it.isNext())
    {
            // Get the local grid key
            auto key = it.get_dist();

            // Old values U and V
            grid.template insert<U>(key) = 1.0;
            grid.template insert<V>(key) = 0.0;

            // Old values U and V
            grid.template insert<U_next>(key) = 0.0;
            grid.template insert<V_next>(key) = 0.0;

            ++it;
    }

    long int x_start = grid.size(0)*1.55f/domain.getHigh(0);
    long int y_start = grid.size(1)*1.55f/domain.getHigh(1);
    long int z_start = grid.size(1)*1.55f/domain.getHigh(2);

    long int x_stop = grid.size(0)*1.85f/domain.getHigh(0);
    long int y_stop = grid.size(1)*1.85f/domain.getHigh(1);
    long int z_stop = grid.size(1)*1.85f/domain.getHigh(2);

    grid_key_dx<3> start({x_start,y_start,z_start});
    grid_key_dx<3> stop ({x_stop,y_stop,z_stop});
    auto it_init = grid.getGridIterator(start,stop);

    while (it_init.isNext())
    {
            auto key = it_init.get_dist();

            grid.template insert<U>(key) = 0.5 + (((double)std::rand())/RAND_MAX -0.5)/10.0;
            grid.template insert<V>(key) = 0.25 + (((double)std::rand())/RAND_MAX -0.5)/20.0;

            ++it_init;
    }

    // spacing of the grid on x and y
    double spacing[3] = {grid.spacing(0),grid.spacing(1),grid.spacing(2)};
    // sync the ghost
    size_t count = 0;
    grid.template ghost_get<U,V>();

    // because we assume that spacing[x] == spacing[y] we use formula 2
    // and we calculate the prefactor of Eq 2
    double uFactor = deltaT * du/(spacing[x]*spacing[x]);
    double vFactor = deltaT * dv/(spacing[x]*spacing[x]);


     //! \cond [stencil get and use] \endcond


    auto func = [uFactor,vFactor,deltaT,F,K](Vc::double_v & u_out,Vc::double_v & v_out,
    															Vc::double_v & u,Vc::double_v & v,
                                                                cross_stencil_v & us,cross_stencil_v & vs,
                                                                unsigned char * mask){

																														 u_out = u + uFactor *(us.xm + us.xp +
																																 	           us.ym + us.yp +
																																			   us.zm + us.zp - 6.0*u) - deltaT * u*v*v
																																									- deltaT * F * (u - 1.0);

																														 v_out = v + vFactor *(vs.xm + vs.xp +
																																	  	  	   vs.ym + vs.yp +
																																			   vs.zm + vs.zp - 6.0*v) + deltaT * u*v*v
																																									- deltaT * (F+K) * v;
                                                                                     };

    grid.conv_cross2<U,V,U_next,V_next,1>({0,0,0},{(long int)sz[0]-1,(long int)sz[1]-1,(long int)sz[2]-1},func);
    grid.conv_cross2<U,V,U_next,V_next,1>({0,0,0},{(long int)sz[0]-1,(long int)sz[1]-1,(long int)sz[2]-1},func);

    bool match = true;

    {
		auto it = grid.getDomainIterator();

		double max_U = 0.0;
		double max_V = 0.0;
		grid_dist_key_dx<3> k_max;
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

			// update based on Eq 2
			if ( fabs(grid.get<U>(Cp) + uFactor * (
																	grid.get<U>(mz) +
																	grid.get<U>(pz) +
																	grid.get<U>(my) +
																	grid.get<U>(py) +
																	grid.get<U>(mx) +
																	grid.get<U>(px) -
																	6.0*grid.get<U>(Cp)) +
																	- deltaT * grid.get<U>(Cp) * grid.get<V>(Cp) * grid.get<V>(Cp) +
																	- deltaT * F * (grid.get<U>(Cp) - 1.0) - grid.get<U_next>(Cp)) > 0.000000001 )
			{
				match = false;
				break;
			}

			// update based on Eq 2
			if ( fabs(grid.get<V>(Cp) + vFactor * (
																	grid.get<V>(mz) +
																	grid.get<V>(pz) +
																	grid.get<V>(my) +
																	grid.get<V>(py) +
																	grid.get<V>(mx) +
																	grid.get<V>(px) -
																	6*grid.get<V>(Cp)) +
																	deltaT * grid.get<U>(Cp) * grid.get<V>(Cp) * grid.get<V>(Cp) +
																	- deltaT * (F+K) * grid.get<V>(Cp) - grid.get<V_next>(Cp)) > 0.000000001 )
			{
				match = false;
				break;
			}

			++it;
		}
    }

    BOOST_REQUIRE_EQUAL(match,true);
}

BOOST_AUTO_TEST_CASE (sgrid_dist_id_soa_write )
{
	periodicity<3> bc = {PERIODIC, PERIODIC, PERIODIC};

	auto & v_cl = create_vcluster<>();

	if (v_cl.size() > 16)
	{return;}

	// Domain
	Box<3,double> domain({-0.3,-0.3,-0.3},{1.0,1.0,1.0});

	// grid size
	size_t sz[3];
	sz[0] = 256;
	sz[1] = 256;
	sz[2] = 256;

	// Ghost
	Ghost<3,long int> g(1);

	sgrid_dist_soa<3,double,aggregate<double,double[3]>> sg1(sz,domain,g,bc);
	sgrid_dist_id<3,double,aggregate<double,double[3]>> sg2(sg1.getDecomposition(),sz,g);

	// create a grid iterator over a bilion point

	auto it = sg1.getGridIterator();

	while(it.isNext())
	{
		auto gkey = it.get();
		auto key = it.get_dist();

		size_t sx = gkey.get(0) - 128;
		size_t sy = gkey.get(1) - 128;
		size_t sz = gkey.get(2) - 128;

		if (sx*sx + sy*sy + sz*sz < 32*32)
		{
			sg1.template insert<0>(key) = 1.0;
			sg1.template insert<1>(key)[0] = gkey.get(0);
			sg1.template insert<1>(key)[1] = gkey.get(1);
			sg1.template insert<1>(key)[2] = gkey.get(2);

			sg2.template insert<0>(key) = 1.0;
			sg2.template insert<1>(key)[0] = gkey.get(0);
			sg2.template insert<1>(key)[1] = gkey.get(1);
			sg2.template insert<1>(key)[2] = gkey.get(2);
		}

		++it;
	}

	sg1.write("sg1_test");
	sg2.write("sg2_test");

	bool test = compare("sg1_test_" + std::to_string(v_cl.rank()) + ".vtk","sg2_test_" + std::to_string(v_cl.rank()) + ".vtk");
	BOOST_REQUIRE_EQUAL(true,test);

	sg1.save("hdf5_w1_test");
	sg2.save("hdf5_w2_test");

	// To uncomment and check
//	sgrid_dist_soa<3,double,aggregate<double,double[3]>> sg1_(sz,domain,g,bc);
//	sgrid_dist_id<3,double,aggregate<double,double[3]>> sg2_(sg1.getDecomposition(),sz,g);

//	sg1.load("hdf5_w1_test");
//	sg2.load("hdf5_w2_test");
}

BOOST_AUTO_TEST_SUITE_END()
