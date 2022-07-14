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

    voroModel.addVolume({0.40,0.1,0.5});
    voroModel.addVolume({0.90,0.9,0.5});

    voroModel.ghost_get_volumes<>();

    voroModel.createVoronoi();

    bool check = voroModel.check_consistent();
    BOOST_REQUIRE_EQUAL(check,true);
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

    Ghost<3,double> g(0.99);

    PolyMesh<3,double,aggregate<double>, // Physical quantities for Cells
                      aggregate<double>,        // Physical quantities for surface
                      aggregate<double>,                                   // Physical quantities for edges
                      aggregate<double> >                                  // Physical quantities for vertices
    voroModel(simBox,bc,g);

    for (int i = 0 ; i < N ; i++)
    {
        Point<3,double> p = {(double)rand()/RAND_MAX,(double)rand()/RAND_MAX,(double)rand()/RAND_MAX};

        voroModel.addVolume(p);
    }

    voroModel.ghost_get_volumes<>();

    voroModel.createVoronoi();

    voroModel.write("bbbbbbb");
    voroModel.getVolumesDist().write("bbbbbbb_vols");

    bool check = voroModel.check_consistent();
    BOOST_REQUIRE_EQUAL(check,true);
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

    voroModel.addVolume({0.84018771715471,0.394382926819093,0.783099223758606});
    voroModel.addVolume({0.798440033476073,0.911647357936784,0.197551369293384});

    voroModel.ghost_get_volumes<>();

    voroModel.createVoronoi();

    bool check = voroModel.check_consistent();
    BOOST_REQUIRE_EQUAL(check,true);

    voroModel.write("output");

}

BOOST_AUTO_TEST_CASE( polymesh_test_case_10 )
{
    size_t bc[3] = {PERIODIC,NON_PERIODIC,NON_PERIODIC};
//    poly_meshtest_case_N_cells(bc,100);

    bc[0] = PERIODIC;
    bc[1] = PERIODIC;
    bc[2] = PERIODIC;

    poly_meshtest_case_N_cells(bc,10);
}

BOOST_AUTO_TEST_CASE( polymesh_test_case_100 )
{
    size_t bc[3] = {PERIODIC,NON_PERIODIC,NON_PERIODIC};
//    poly_meshtest_case_N_cells(bc,100);

    bc[0] = PERIODIC;
    bc[1] = PERIODIC;
    bc[2] = PERIODIC;

    poly_meshtest_case_N_cells(bc,100);
}

BOOST_AUTO_TEST_CASE( polymesh_test_case_1000 )
{
    size_t bc[3];

    bc[0] = PERIODIC;
    bc[1] = PERIODIC;
    bc[2] = PERIODIC;

    poly_meshtest_case_N_cells(bc,1000);
}

BOOST_AUTO_TEST_CASE( polymesh_test_grad_volume )
{
    size_t bc[3];

    bc[0] = PERIODIC;
    bc[1] = PERIODIC;
    bc[2] = PERIODIC;

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

    voroModel.addVolume({0.40,0.1,0.5});
    voroModel.addVolume({0.90,0.9,0.1});

    voroModel.ghost_get_volumes<>();

    voroModel.createVoronoi();

    bool check = voroModel.check_consistent();
    BOOST_REQUIRE_EQUAL(check,true);

    // First we take 2 cells we create voronoi

    // we move cell 0 by dx 1e-6 ... 0.5e-7 ... 0.25e-7 ... 0.125e-7
    //                   dy 1e-6 ... 0.5e-7 ... 0.25e-7 ... 0.125e-7
    //                   dz 1e-6 ... 0.5e-7 ... 0.25e-7 ... 0.125e-7

    openfpm::vector<aggregate<double[4][3]>> Gv;
    openfpm::vector<aggregate<double[4][3]>> Gb;

    bool Mask[4][3];

    grad_voronoi_cell(voroModel,0,Gv);
    grad_barycenter(voroModel,0,Gv,Gb);
    Point<3,double> deriv = derivate_V_cell(voroModel,0,0,Gv,Gb,Mask);

    BOOST_REQUIRE(fabs(deriv[0]) < 1e-7);
    BOOST_REQUIRE(fabs(deriv[1]) < 1e-7);
    BOOST_REQUIRE(fabs(deriv[2]) < 1e-7);

    std::cout << deriv.toString() << std::endl;
}

template<typename points_type, unsigned int N>
void test_grad_center(points_type & points, int (& num)[N])
{
    size_t bc[3];

    bc[0] = PERIODIC;
    bc[1] = PERIODIC;
    bc[2] = PERIODIC;

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

    for (int i = 0 ; i < points.size() ; i++)
    {
        voroModel.addVolume(points.get(i));
    }

    voroModel.ghost_get_volumes<>();

    voroModel.createVoronoi();

    bool check = voroModel.check_consistent();
    BOOST_REQUIRE_EQUAL(check,true);

    // First we take 2 cells we create voronoi

    // we move cell 0 by dx 1e-6 ... 0.5e-7 ... 0.25e-7 ... 0.125e-7
    //                   dy 1e-6 ... 0.5e-7 ... 0.25e-7 ... 0.125e-7
    //                   dz 1e-6 ... 0.5e-7 ... 0.25e-7 ... 0.125e-7

    openfpm::vector<aggregate<double[4][3]>> Gv;
    openfpm::vector<aggregate<double[4][3]>> Gb;

    bool Mask[4][3];

    for (int j = 0 ; j < N; j++)
    {
        int i = num[j];
        grad_voronoi_cell(voroModel,i,Gv);
        grad_barycenter(voroModel,i,Gv,Gb);
        Point<3,double> deriv = derivate_V_cell(voroModel,i,i,Gv,Gb,Mask);

        std::cout << deriv.toString() << std::endl;

        for (int c = 0 ; c < 3 ; c++)
        {
            // Check convergence (Point i) Derivative X
            double Volume_now = voroModel.getVolume(i);

            voroModel.getVolumePos(i)[c] = voroModel.getVolumePos(i)[c] + 1e-6;

            voroModel.ghost_get_volumes<>();
            voroModel.createVoronoi();

            double Volume_after = voroModel.getVolume(i);

            double Derivative = (Volume_after - Volume_now) / 1.e-6;

            BOOST_REQUIRE_CLOSE(Derivative,deriv[c],0.01);

            voroModel.getVolumePos(i)[c] = voroModel.getVolumePos(i)[c] - 1e-6;

            voroModel.ghost_get_volumes<>();
            voroModel.createVoronoi();
        }

        // NN change dVolume


    }
}

template<typename points_type, unsigned int N>
void test_grad_nn(points_type & points, int (& num)[N])
{
    size_t bc[3];

    bc[0] = PERIODIC;
    bc[1] = PERIODIC;
    bc[2] = PERIODIC;

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

    for (int i = 0 ; i < points.size() ; i++)
    {
        voroModel.addVolume(points.get(i));
    }

    voroModel.ghost_get_volumes<>();

    voroModel.createVoronoi();

    bool check = voroModel.check_consistent();
    BOOST_REQUIRE_EQUAL(check,true);

    // First we take 2 cells we create voronoi

    // we move cell 0 by dx 1e-6 ... 0.5e-7 ... 0.25e-7 ... 0.125e-7
    //                   dy 1e-6 ... 0.5e-7 ... 0.25e-7 ... 0.125e-7
    //                   dz 1e-6 ... 0.5e-7 ... 0.25e-7 ... 0.125e-7

    openfpm::vector<aggregate<double[4][3]>> Gv;
    openfpm::vector<aggregate<double[4][3]>> Gb;
    openfpm::vector<aggregate<double[4][3]>> Gv2;

    bool Mask[4][3];

    for (int j = 0 ; j < N; j++)
    {
        int i = num[j];
        grad_voronoi_cell(voroModel,i,Gv);
        grad_barycenter(voroModel,i,Gv,Gb);
        Point<3,double> deriv = derivate_V_cell(voroModel,i,i,Gv,Gb,Mask);

        // get neighborhood cells

        voroModel.ForAllVolumeSurfaces(i,[&](int f_ind, int conn){

//            voroModel.get
            int k = voroModel.getFaceNNVolume(f_ind,i);
            grad_vor_bar_neighbour(voroModel,i,k,Gv,Gb,Gv2);
            Point<3,double> dV=derivate_V_cell(voroModel, k, i,Gv2,Gb, Mask);

            std::cout << "Grad: " << dV.toString() << std::endl;
        });

/*        for (int c = 0 ; c < 3 ; c++)
        {
            // Check convergence (Point i) Derivative X
            double Volume_now = voroModel.getVolume(i);

            voroModel.getVolumePos(i)[c] = voroModel.getVolumePos(i)[c] + 1e-6;

            voroModel.ghost_get_volumes<>();
            voroModel.createVoronoi();

            double Volume_after = voroModel.getVolume(i);

            double Derivative = (Volume_after - Volume_now) / 1.e-6;

            BOOST_REQUIRE_CLOSE(Derivative,deriv[c],0.01);

            voroModel.getVolumePos(i)[c] = voroModel.getVolumePos(i)[c] - 1e-6;

            voroModel.ghost_get_volumes<>();
            voroModel.createVoronoi();
        }

        // NN change dVolume

        grad_vor_bar_neighbour(voroModel,0,1,Gv,Gv2,Gb);*/
    }
}

BOOST_AUTO_TEST_CASE( polymesh_test_grad_volume_3_point )
{
    int num[3] = {0,1,2};

    openfpm::vector<Point<3,double>> points;

    points.add({0.40,0.1,0.5});
    points.add({0.90,0.9,0.1});
    points.add({0.5,0.5,0.3});

    test_grad_center(points,num);
    test_grad_nn(points,num);
}

BOOST_AUTO_TEST_CASE( polymesh_test_grad_volume_100_point )
{
    int num[3] = {11,27,82};

    openfpm::vector<Point<3,double>> points;

    for (int i = 0 ; i < 100 ; i++)
    {
        Point<3,double> p = {(double)rand()/RAND_MAX,(double)rand()/RAND_MAX,(double)rand()/RAND_MAX};

        points.add(p);
    }

    test_grad_center(points,num);
    test_grad_nn(points,num);
}

BOOST_AUTO_TEST_SUITE_END()

