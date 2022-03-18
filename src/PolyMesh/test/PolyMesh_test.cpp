/*
 * PolyMesh_tests.cpp
 *
 *  Created on: Aug 16, 2016
 *      Author: i-bird
 */

#include "config.h"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "PolyMesh/PolyMesh.hpp"
#include <stdlib.h>


BOOST_AUTO_TEST_SUITE( polymesh_test_suite )

static void poly_meshtest_case_2_cells(size_t (& bc)[3])
{
    constexpr int x = 0;
    constexpr int y = 1;
    constexpr int z = 2;

    Box<3,double> SimBox;

    constexpr int Volume = 0;
    constexpr int CellMinVolume = 1;
    constexpr int Pressure = 2;
    constexpr int OsmoticPressure = 3;
    constexpr int NumberOfIons = 4;

    constexpr int Area = 1;
    constexpr int SurfaceTension = 0;
    constexpr int SurfaceWaterPermeability = 2;
    constexpr int Flux = 3;

    Box<3,double> simBox({0,0,0},{1.0,1.0,1.0});

    Ghost<3,double> g(0.95);

    PolyMesh<3,double,aggregate<double,double,double,double,double>, // Physical quantities for Cells
                      aggregate<double,double,double,double>,        // Physical quantities for surface
                      aggregate<>,                                   // Physical quantities for edges
                      aggregate<> >                                  // Physical quantities for vertices
    voroModel(simBox,bc,g);

    voroModel.addCell({0.40,0.1,0.5});
    voroModel.addCell({0.90,0.9,0.5});

    voroModel.ghost_get_volumes<>();

    voroModel.createVoronoi();

    bool check = voroModel.check_consistent();
    BOOST_REQUIRE_EQUAL(check,true);

    voroModel.write("aaaaaaaaaa");
}


static void poly_meshtest_case_N_cells(size_t (& bc)[3], int N)
{
    constexpr int x = 0;
    constexpr int y = 1;
    constexpr int z = 2;

    Box<3,double> SimBox;

    constexpr int Volume = 0;
    constexpr int CellMinVolume = 1;
    constexpr int Pressure = 2;
    constexpr int OsmoticPressure = 3;
    constexpr int NumberOfIons = 4;

    constexpr int Area = 1;
    constexpr int SurfaceTension = 0;
    constexpr int SurfaceWaterPermeability = 2;
    constexpr int Flux = 3;

    Box<3,double> simBox({0,0,0},{1.0,1.0,1.0});

    Ghost<3,double> g(0.95);

    PolyMesh<3,double,aggregate<double,double,double,double,double>, // Physical quantities for Cells
                      aggregate<double,double,double,double>,        // Physical quantities for surface
                      aggregate<>,                                   // Physical quantities for edges
                      aggregate<> >                                  // Physical quantities for vertices
    voroModel(simBox,bc,g);

    for (int i = 0 ; i < N ; i++)
    {
        Point<3,double> p = {(double)rand()/RAND_MAX,(double)rand()/RAND_MAX,(double)rand()/RAND_MAX};

        voroModel.addCell(p);
    }

    voroModel.ghost_get_volumes<>();

    voroModel.createVoronoi();

    bool check = voroModel.check_consistent();
    BOOST_REQUIRE_EQUAL(check,true);

    std::cout << voroModel.toString() << std::endl;

    voroModel.write("aaaaaaaaaa");
}

BOOST_AUTO_TEST_CASE( polymesh_test_case )
{
    size_t bc[3] = {PERIODIC,NON_PERIODIC,NON_PERIODIC};

    poly_meshtest_case_2_cells(bc);

    bc[0] = PERIODIC;
    bc[1] = PERIODIC;
    bc[2] = PERIODIC;

    poly_meshtest_case_2_cells(bc);
}

BOOST_AUTO_TEST_CASE( polymesh_test_case_2 )
{
    constexpr int x = 0;
    constexpr int y = 1;
    constexpr int z = 2;

    Box<3,double> SimBox;

    constexpr int Volume = 0;
    constexpr int CellMinVolume = 1;
    constexpr int Pressure = 2;
    constexpr int OsmoticPressure = 3;
    constexpr int NumberOfIons = 4;

    constexpr int Area = 1;
    constexpr int SurfaceTension = 0;
    constexpr int SurfaceWaterPermeability = 2;
    constexpr int Flux = 3;

    Box<3,double> simBox({0,0,0},{1.0,1.0,1.0});
    size_t bc[3] = {PERIODIC,NON_PERIODIC,NON_PERIODIC};

    Ghost<3,double> g(0.99);

    PolyMesh<3,double,aggregate<double,double,double,double,double>, // Physical quantities for Cells
                      aggregate<double,double,double,double>,        // Physical quantities for surface
                      aggregate<>,                                   // Physical quantities for edges
                      aggregate<> >                                  // Physical quantities for vertices
    voroModel(simBox,bc,g);

    voroModel.addCell({0.84018771715471,0.394382926819093,0.783099223758606});
    voroModel.addCell({0.798440033476073,0.911647357936784,0.197551369293384});

    voroModel.ghost_get_volumes<>();

    voroModel.createVoronoi();

    bool check = voroModel.check_consistent();
    BOOST_REQUIRE_EQUAL(check,true);

    voroModel.write("output");

}

BOOST_AUTO_TEST_CASE( polymesh_test_case_100 )
{
    size_t bc[3] = {PERIODIC,NON_PERIODIC,NON_PERIODIC};
    poly_meshtest_case_N_cells(bc,100);

    bc[0] = PERIODIC;
    bc[1] = PERIODIC;
    bc[2] = PERIODIC;

    poly_meshtest_case_N_cells(bc,100);

/*    constexpr int x = 0;
    constexpr int y = 1;
    constexpr int z = 2;

    Box<3,double> SimBox;

    constexpr int Volume = 0;
    constexpr int CellMinVolume = 1;
    constexpr int Pressure = 2;
    constexpr int OsmoticPressure = 3;
    constexpr int NumberOfIons = 4;

    constexpr int Area = 1;
    constexpr int SurfaceTension = 0;
    constexpr int SurfaceWaterPermeability = 2;
    constexpr int Flux = 3;

    Box<3,double> simBox({0,0,0},{1.0,1.0,1.0});
    size_t bc[3] = {PERIODIC,NON_PERIODIC,NON_PERIODIC};

    Ghost<3,double> g(0.95);

    PolyMesh<3,double,aggregate<double,double,double,double,double>, // Physical quantities for Cells
                      aggregate<double,double,double,double>,        // Physical quantities for surface
                      aggregate<>,                                   // Physical quantities for edges
                      aggregate<> >                                  // Physical quantities for vertices
    voroModel(simBox,bc,g);

    // 
    for (int i = 0 ; i < 100 ; i++)
    {
        Point<3,double> p = {(double)rand()/RAND_MAX,(double)rand()/RAND_MAX,(double)rand()/RAND_MAX};

        voroModel.addCell(p);
    }

    voroModel.ghost_get_volumes<>();

    voroModel.createVoronoi();

    bool check = voroModel.check_consistent();
    BOOST_REQUIRE_EQUAL(check,true);

    voroModel.write("output");*/

}

BOOST_AUTO_TEST_CASE( polymesh_test_case_1000 )
{
//    size_t bc[3] = {PERIODIC,NON_PERIODIC,NON_PERIODIC};
//    poly_meshtest_case_N_cells(bc,1000);

    size_t bc[3];

    bc[0] = PERIODIC;
    bc[1] = PERIODIC;
    bc[2] = PERIODIC;

    poly_meshtest_case_N_cells(bc,1000);

/*    constexpr int x = 0;
    constexpr int y = 1;
    constexpr int z = 2;

    Box<3,double> SimBox;

    constexpr int Volume = 0;
    constexpr int CellMinVolume = 1;
    constexpr int Pressure = 2;
    constexpr int OsmoticPressure = 3;
    constexpr int NumberOfIons = 4;

    constexpr int Area = 1;
    constexpr int SurfaceTension = 0;
    constexpr int SurfaceWaterPermeability = 2;
    constexpr int Flux = 3;

    Box<3,double> simBox({0,0,0},{1.0,1.0,1.0});
    size_t bc[3] = {PERIODIC,NON_PERIODIC,NON_PERIODIC};

    Ghost<3,double> g(0.95);

    PolyMesh<3,double,aggregate<double,double,double,double,double>, // Physical quantities for Cells
                      aggregate<double,double,double,double>,        // Physical quantities for surface
                      aggregate<>,                                   // Physical quantities for edges
                      aggregate<> >                                  // Physical quantities for vertices
    voroModel(simBox,bc,g);

    // 
    for (int i = 0 ; i < 1000 ; i++)
    {
        Point<3,double> p = {(double)rand()/RAND_MAX,(double)rand()/RAND_MAX,(double)rand()/RAND_MAX};

        voroModel.addCell(p);
    }

    voroModel.ghost_get_volumes<>();

    voroModel.createVoronoi();

    bool check = voroModel.check_consistent();
    BOOST_REQUIRE_EQUAL(check,true);

    voroModel.write("output");*/

}

BOOST_AUTO_TEST_SUITE_END()

