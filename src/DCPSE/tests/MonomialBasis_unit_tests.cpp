//
// Created by tommaso on 21/03/19.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <Space/Shape/Point.hpp>
#include <DCPSE/MonomialBasis.hpp>

BOOST_AUTO_TEST_SUITE(MonomialBasis_tests)

    BOOST_AUTO_TEST_CASE(Monomial_3D_test)
    {
        // Testing constructor from a Point
        Point<3, unsigned int> p({1, 2, 3});
        Monomial<3> mbe(p);
        BOOST_REQUIRE_EQUAL(mbe.order(), 6);
        BOOST_REQUIRE_EQUAL(mbe.getExponent(0), 1);
        BOOST_REQUIRE_EQUAL(mbe.getExponent(1), 2);
        BOOST_REQUIRE_EQUAL(mbe.getExponent(2), 3);

        // Testing copy constructor
        Monomial<3> mbe2(mbe);
        BOOST_REQUIRE_EQUAL(mbe2.order(), 6);
        BOOST_REQUIRE_EQUAL(mbe2.getExponent(0), 1);
        BOOST_REQUIRE_EQUAL(mbe2.getExponent(1), 2);
        BOOST_REQUIRE_EQUAL(mbe2.getExponent(2), 3);

        // Testing copy assignment operator
        Monomial<3> mbe3;
        mbe3 = mbe;
        BOOST_REQUIRE_EQUAL(mbe3.order(), 6);
        BOOST_REQUIRE_EQUAL(mbe3.getExponent(0), 1);
        BOOST_REQUIRE_EQUAL(mbe3.getExponent(1), 2);
        BOOST_REQUIRE_EQUAL(mbe3.getExponent(2), 3);

        // Testing default constructor
        Monomial<3> mbe0;
        BOOST_REQUIRE_EQUAL(mbe0.order(), 0);
        BOOST_REQUIRE_EQUAL(mbe0.getExponent(0), 0);
        BOOST_REQUIRE_EQUAL(mbe0.getExponent(1), 0);
        BOOST_REQUIRE_EQUAL(mbe0.getExponent(2), 0);

        // Testing one-by-one set values
        mbe0.setExponent(0, 5);
        mbe0.setExponent(1, 5);
        mbe0.setExponent(2, 5);
        BOOST_REQUIRE_EQUAL(mbe0.order(), 15);
        BOOST_REQUIRE_EQUAL(mbe0.getExponent(0), 5);
        BOOST_REQUIRE_EQUAL(mbe0.getExponent(1), 5);
        BOOST_REQUIRE_EQUAL(mbe0.getExponent(2), 5);
    }

    BOOST_AUTO_TEST_CASE(Monomial_evaluate_test)
    {
        // Test with uniform exponent and different coordinates
        {
            Monomial<3> monomial(Point<3, unsigned int>({1, 1, 1}));
            double res = monomial.evaluate(Point<3, double>({5, 3, 2}));
            double expected = 5 * 3 * 2;
            BOOST_REQUIRE_CLOSE(res, expected, 1e-6);
        }
        // Test with different exponents and uniform coordinates
        {
            Monomial<3> monomial(Point<3, unsigned int>({5, 3, 2}));
            double res = monomial.evaluate(Point<3, double>({2, 2, 2}));
            double expected = pow(2,5) * pow(2,3) * pow(2,2);
            BOOST_REQUIRE_CLOSE(res, expected, 1e-6);
        }
    }

    BOOST_AUTO_TEST_CASE(Monomial_derivative_test)
    {
        // Test differentiation along 3 dimensions independently
        {
            Monomial<3> monomial(Point<3, unsigned int>({6, 6, 6}));
            Monomial<3> derivative = monomial.getDerivative(Point<3, double>({3, 0, 0}));
            Monomial<3> expected(Point<3, unsigned int>({3, 6, 6}), 6*5*4);
            BOOST_REQUIRE_EQUAL(derivative, expected);
        }
        {
            Monomial<3> monomial(Point<3, unsigned int>({6, 6, 6}));
            Monomial<3> derivative = monomial.getDerivative(Point<3, double>({0, 3, 0}));
            Monomial<3> expected(Point<3, unsigned int>({6, 3, 6}), 6*5*4);
            BOOST_REQUIRE_EQUAL(derivative, expected);
        }
        {
            Monomial<3> monomial(Point<3, unsigned int>({6, 6, 6}));
            Monomial<3> derivative = monomial.getDerivative(Point<3, double>({0, 0, 3}));
            Monomial<3> expected(Point<3, unsigned int>({6, 6, 3}), 6*5*4);
            BOOST_REQUIRE_EQUAL(derivative, expected);
        }
        // Test for correctly differentiating beyond degree
        {
            Monomial<3> monomial(Point<3, unsigned int>({6, 6, 6}));
            Monomial<3> derivative = monomial.getDerivative(Point<3, double>({7, 0, 0}));
            Monomial<3> expected(Point<3, unsigned int>({0, 6, 6}), 0);
            BOOST_REQUIRE_EQUAL(derivative, expected);
        }
        {
            Monomial<3> monomial(Point<3, unsigned int>({6, 6, 6}));
            Monomial<3> derivative = monomial.getDerivative(Point<3, double>({0, 7, 0}));
            Monomial<3> expected(Point<3, unsigned int>({6, 0, 6}), 0);
            BOOST_REQUIRE_EQUAL(derivative, expected);
        }
        {
            Monomial<3> monomial(Point<3, unsigned int>({6, 6, 6}));
            Monomial<3> derivative = monomial.getDerivative(Point<3, double>({0, 0, 7}));
            Monomial<3> expected(Point<3, unsigned int>({6, 6, 0}), 0);
            BOOST_REQUIRE_EQUAL(derivative, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(MonomialBasis_2D_test)
    {
        {
            MonomialBasis<2> mb({1, 0}, 2);

            // Check that the basis contains all the elements it is supposed to contain
            std::vector<Monomial<2>> control;
            control.emplace_back(Point<2, unsigned int>({0, 0}));
            control.emplace_back(Point<2, unsigned int>({1, 0}));
            control.emplace_back(Point<2, unsigned int>({0, 1}));
            control.emplace_back(Point<2, unsigned int>({1, 1}));
            control.emplace_back(Point<2, unsigned int>({2, 0}));
            control.emplace_back(Point<2, unsigned int>({0, 2}));

            auto first = mb.getElements().begin();
            auto last = mb.getElements().end();
            int foundElements = 0;
            for (const auto &item : control)
            {
                auto pos = std::find(first, last, item);
                foundElements += (pos != last);
            }
            BOOST_REQUIRE_EQUAL(foundElements, control.size());
        }
        {
            MonomialBasis<2> mb({1, 1}, 2);
//            std::cout << mb << std::endl;
//            std::cout << mb.size() << std::endl;

            // Check that the basis contains all the elements it is supposed to contain
            std::vector<Monomial<2>> control;
            control.emplace_back(Point<2, unsigned int>({0, 0}));
            control.emplace_back(Point<2, unsigned int>({1, 0}));
            control.emplace_back(Point<2, unsigned int>({0, 1}));
            control.emplace_back(Point<2, unsigned int>({1, 1}));
            control.emplace_back(Point<2, unsigned int>({2, 0}));
            control.emplace_back(Point<2, unsigned int>({0, 2}));
            control.emplace_back(Point<2, unsigned int>({2, 1}));
            control.emplace_back(Point<2, unsigned int>({1, 2}));
            control.emplace_back(Point<2, unsigned int>({3, 0}));
            control.emplace_back(Point<2, unsigned int>({0, 3}));

            BOOST_REQUIRE_EQUAL(mb.size(), control.size());

            auto first = mb.getElements().begin();
            auto last = mb.getElements().end();
            int foundElements = 0;
            for (const auto &item : control)
            {
                auto pos = std::find(first, last, item);
                foundElements += (pos != last);
            }
            BOOST_REQUIRE_EQUAL(foundElements, control.size());
        }
    }

BOOST_AUTO_TEST_SUITE_END()