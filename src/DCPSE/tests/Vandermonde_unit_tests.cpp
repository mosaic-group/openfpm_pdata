//
// Created by tommaso on 22/03/19.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <DCPSE/MonomialBasis.hpp>
#include <DCPSE/VandermondeRowBuilder.hpp>
#include <DCPSE/Vandermonde.hpp>
#include "../../openfpm_numerics/src/DMatrix/EMatrix.hpp"

BOOST_AUTO_TEST_SUITE(Vandermonde_tests)

// If EIGEN is not present, EMatrix is not available and we don't need to build this test
#ifdef HAVE_EIGEN

    BOOST_AUTO_TEST_CASE(VandermondeRowBuilder_AllOnes_test)
    {
        MonomialBasis<2> mb({1, 0}, 2);
        EMatrix<double, Eigen::Dynamic, Eigen::Dynamic> row(1, mb.size());
        Point<2, double> x({2, 2});
        double eps = 2;
        VandermondeRowBuilder<2, double> vrb(mb);
        vrb.buildRow(row, 0, x, eps);
        // For the way the row has been constructed, it should be composed of only 1s
        bool isRowAllOnes = true;
        for (int i = 0; i < mb.size(); ++i)
        {
            isRowAllOnes = isRowAllOnes && (row(0, i) == 1);
        }
        BOOST_REQUIRE(isRowAllOnes);
    }

    BOOST_AUTO_TEST_CASE(VandermondeRowBuilder_OneZero_test)
    {
        MonomialBasis<2> mb({1, 0}, 2);
        EMatrix<double, Eigen::Dynamic, Eigen::Dynamic> row(1, mb.size());
        Point<2, double> x({1, 0});
        double eps = 1;
        VandermondeRowBuilder<2, double> vrb(mb);
        vrb.buildRow(row, 0, x, eps);
        // For the way the row has been constructed, it should be composed of only 1s
        bool areValuesOk = true;
        for (int i = 0; i < mb.size(); ++i)
        {
            bool isThereY = mb.getElement(i).getExponent(1) > 0;
            bool curCheck = (row(0, i) == !isThereY);
            areValuesOk = areValuesOk && curCheck;
        }
        BOOST_REQUIRE(areValuesOk);
    }

    BOOST_AUTO_TEST_CASE(VandermondeRowBuilder_ZeroOne_test)
    {
        MonomialBasis<2> mb({1, 0}, 2);
        EMatrix<double, Eigen::Dynamic, Eigen::Dynamic> row(1, mb.size());
        Point<2, double> x({0, 1});
        double eps = 1;
        VandermondeRowBuilder<2, double> vrb(mb);
        vrb.buildRow(row, 0, x, eps);
        // For the way the row has been constructed, it should be composed of only 1s
        bool areValuesOk = true;
        for (int i = 0; i < mb.size(); ++i)
        {
            bool isThereX = mb.getElement(i).getExponent(0) > 0;
            bool curCheck = (row(0, i) == !isThereX);
            areValuesOk = areValuesOk && curCheck;
        }
        BOOST_REQUIRE(areValuesOk);
    }

    BOOST_AUTO_TEST_CASE(Vandermonde_KnownValues_test)
    {
        MonomialBasis<2> mb({1, 0}, 2);
//        std::cout << mb << std::endl; // Debug
        Point<2, double> x({0, 0});
        std::vector<Point<2, double>> neighbours;
        EMatrix<double, Eigen::Dynamic, Eigen::Dynamic> V(mb.size(), mb.size());
        EMatrix<double, Eigen::Dynamic, Eigen::Dynamic> ExpectedV(mb.size(), mb.size());
        // Now push 4 points on diagonal as neighbours
        neighbours.push_back({1, 1});
        neighbours.push_back({-1, -1});
        neighbours.push_back({1, -1});
        neighbours.push_back({-1, 1});
        neighbours.push_back({0, 1});
        neighbours.push_back({0, -1});

        // ...and get the matrix V
        Vandermonde<2, double, EMatrix<double, Eigen::Dynamic, Eigen::Dynamic>> vandermonde(x, neighbours, mb);
        vandermonde.getMatrix(V);

        // Now build the matrix of expected values
        ExpectedV(0, 0) = 1;
        ExpectedV(0, 1) = -0.3;
        ExpectedV(0, 2) = 0.09;
        ExpectedV(0, 3) = -0.3;
        ExpectedV(0, 4) = +0.09;
        ExpectedV(0, 5) = 0.09;
        ExpectedV(1, 0) = 1;
        ExpectedV(1, 1) = +0.3;
        ExpectedV(1, 2) = 0.09;
        ExpectedV(1, 3) = +0.3;
        ExpectedV(1, 4) = +0.09;
        ExpectedV(1, 5) = 0.09;
        ExpectedV(2, 0) = 1;
        ExpectedV(2, 1) = -0.3;
        ExpectedV(2, 2) = 0.09;
        ExpectedV(2, 3) = +0.3;
        ExpectedV(2, 4) = -0.09;
        ExpectedV(2, 5) = 0.09;
        ExpectedV(3, 0) = 1;
        ExpectedV(3, 1) = +0.3;
        ExpectedV(3, 2) = 0.09;
        ExpectedV(3, 3) = -0.3;
        ExpectedV(3, 4) = -0.09;
        ExpectedV(3, 5) = 0.09;
        ExpectedV(4, 0) = 1;
        ExpectedV(4, 1) = 0;
        ExpectedV(4, 2) = 0;
        ExpectedV(4, 3) = -0.3;
        ExpectedV(4, 4) = 0;
        ExpectedV(4, 5) = 0.09;
        ExpectedV(5, 0) = 1;
        ExpectedV(5, 1) = 0;
        ExpectedV(5, 2) = 0;
        ExpectedV(5, 3) = +0.3;
        ExpectedV(5, 4) = 0;
        ExpectedV(5, 5) = 0.09;

        // Now check that the values in V are all the expected ones
        for (int i = 0; i < mb.size(); ++i)
        {
//            std::cout << "N[" << i << "] = " << neighbours[i].toString() << std::endl;
            for (int j = 0; j < mb.size(); ++j)
            {
//                std::cout << ">> V[" << i << "," << j << "] = " << V(i,j) << std::endl;
                BOOST_REQUIRE_CLOSE(V(i, j), ExpectedV(i, j), 1e-6);
            }
        }
    }

    BOOST_AUTO_TEST_CASE(Vandermonde_TranslatedSetup_test)
    {
        MonomialBasis<2> mb({1, 0}, 2);
//        std::cout << mb << std::endl; // Debug
        Point<2, double> x({1, 1});
        std::vector<Point<2, double>> neighbours;
        EMatrix<double, Eigen::Dynamic, Eigen::Dynamic> V(mb.size(), mb.size());
        EMatrix<double, Eigen::Dynamic, Eigen::Dynamic> ExpectedV(mb.size(), mb.size());
        // Now push 4 points on diagonal as neighbours
        neighbours.push_back({2, 2});
        neighbours.push_back({0, 0});
        neighbours.push_back({2, 0});
        neighbours.push_back({0, 2});
        neighbours.push_back({1, 2});
        neighbours.push_back({1, 0});
        // ...and get the matrix V
        Vandermonde<2, double, EMatrix<double, Eigen::Dynamic, Eigen::Dynamic>> vandermonde(x, neighbours, mb);
        vandermonde.getMatrix(V);

        // Now build the matrix of expected values
        ExpectedV(0, 0) = 1;
        ExpectedV(0, 1) = -0.3;
        ExpectedV(0, 2) = 0.09;
        ExpectedV(0, 3) = -0.3;
        ExpectedV(0, 4) = +0.09;
        ExpectedV(0, 5) = 0.09;
        ExpectedV(1, 0) = 1;
        ExpectedV(1, 1) = +0.3;
        ExpectedV(1, 2) = 0.09;
        ExpectedV(1, 3) = +0.3;
        ExpectedV(1, 4) = +0.09;
        ExpectedV(1, 5) = 0.09;
        ExpectedV(2, 0) = 1;
        ExpectedV(2, 1) = -0.3;
        ExpectedV(2, 2) = 0.09;
        ExpectedV(2, 3) = +0.3;
        ExpectedV(2, 4) = -0.09;
        ExpectedV(2, 5) = 0.09;
        ExpectedV(3, 0) = 1;
        ExpectedV(3, 1) = +0.3;
        ExpectedV(3, 2) = 0.09;
        ExpectedV(3, 3) = -0.3;
        ExpectedV(3, 4) = -0.09;
        ExpectedV(3, 5) = 0.09;
        ExpectedV(4, 0) = 1;
        ExpectedV(4, 1) = 0;
        ExpectedV(4, 2) = 0;
        ExpectedV(4, 3) = -0.3;
        ExpectedV(4, 4) = 0;
        ExpectedV(4, 5) = 0.09;
        ExpectedV(5, 0) = 1;
        ExpectedV(5, 1) = 0;
        ExpectedV(5, 2) = 0;
        ExpectedV(5, 3) = +0.3;
        ExpectedV(5, 4) = 0;
        ExpectedV(5, 5) = 0.09;

        // Now check that the values in V are all the expected ones
        for (int i = 0; i < mb.size(); ++i)
        {
//            std::cout << "N[" << i << "] = " << neighbours[i].toString() << std::endl;
            for (int j = 0; j < mb.size(); ++j)
            {
//                std::cout << ">> V[" << i << "," << j << "] = " << V(i,j) << std::endl;
                BOOST_REQUIRE_CLOSE(V(i, j), ExpectedV(i, j), 1e-6);
            }
        }
    }

#endif // HAVE_EIGEN

BOOST_AUTO_TEST_SUITE_END()
