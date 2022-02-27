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

BOOST_AUTO_TEST_CASE( polymesh_test_case )
{
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
    size_t bc[3] = {PERIODIC,PERIODIC,PERIODIC};

    Ghost<3,double> g(0.9);

    PolyMesh<3,double,aggregate<double,double,double,double,double>, // Physical quantities for Cells
                      aggregate<double,double,double,double>,        // Physical quantities for surface
                      aggregate<>,                                   // Physical quantities for edges
                      aggregate<> >                                  // Physical quantities for vertices
    voroModel(simBox,bc,g);

    // 
/*    for (int i = 0 ; i < 100 ; i++)
    {
        Point<3,double> p = {(double)rand()/RAND_MAX,(double)rand()/RAND_MAX,(double)rand()/RAND_MAX};

        voroModel.addCell(p);
    }*/

    voroModel.addCell({0.25,0.5,0.5});
    voroModel.addCell({0.75,0.5,0.5});

    voroModel.createVoronoi();

    

    voroModel.write("output");

}

BOOST_AUTO_TEST_SUITE_END()

