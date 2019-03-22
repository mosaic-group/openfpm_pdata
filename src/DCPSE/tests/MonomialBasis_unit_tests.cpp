//
// Created by tommaso on 21/03/19.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <Space/Shape/Point.hpp>
#include <DCPSE/MonomialBasisElement.hpp>
#include <DCPSE/MonomialBasis.hpp>

BOOST_AUTO_TEST_SUITE(MonomialBasis_tests)

    BOOST_AUTO_TEST_CASE(MonomialBasisElement_3D_test)
    {
        // Testing constructor from a Point
        Point<3, unsigned int> p({1, 2, 3});
        MonomialBasisElement<3> mbe(p);
        BOOST_REQUIRE_EQUAL(mbe.order(), 6);
        BOOST_REQUIRE_EQUAL(mbe.getExponent(0), 1);
        BOOST_REQUIRE_EQUAL(mbe.getExponent(1), 2);
        BOOST_REQUIRE_EQUAL(mbe.getExponent(2), 3);

        // Testing copy constructor
        MonomialBasisElement<3> mbe2(mbe);
        BOOST_REQUIRE_EQUAL(mbe2.order(), 6);
        BOOST_REQUIRE_EQUAL(mbe2.getExponent(0), 1);
        BOOST_REQUIRE_EQUAL(mbe2.getExponent(1), 2);
        BOOST_REQUIRE_EQUAL(mbe2.getExponent(2), 3);

        // Testing copy assignment operator
        MonomialBasisElement<3> mbe3;
        mbe3 = mbe;
        BOOST_REQUIRE_EQUAL(mbe3.order(), 6);
        BOOST_REQUIRE_EQUAL(mbe3.getExponent(0), 1);
        BOOST_REQUIRE_EQUAL(mbe3.getExponent(1), 2);
        BOOST_REQUIRE_EQUAL(mbe3.getExponent(2), 3);

        // Testing default constructor
        MonomialBasisElement<3> mbe0;
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

    BOOST_AUTO_TEST_CASE(MonomialBasis_2D_test)
    {
        {
            MonomialBasis<2> mb({1, 0}, 2);

            // Check that the basis contains all the elements it is supposed to contain
            std::vector<MonomialBasisElement<2>> control;
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

            // Check that the basis contains all the elements it is supposed to contain
            std::vector<MonomialBasisElement<2>> control;
            control.emplace_back(Point<2, unsigned int>({1, 0}));
            control.emplace_back(Point<2, unsigned int>({0, 1}));
            control.emplace_back(Point<2, unsigned int>({1, 1}));
            control.emplace_back(Point<2, unsigned int>({2, 0}));
            control.emplace_back(Point<2, unsigned int>({0, 2}));
            control.emplace_back(Point<2, unsigned int>({2, 1}));
            control.emplace_back(Point<2, unsigned int>({1, 2}));
            control.emplace_back(Point<2, unsigned int>({3, 0}));
            control.emplace_back(Point<2, unsigned int>({0, 3}));

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