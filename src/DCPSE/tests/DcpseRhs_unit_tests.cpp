//
// Created by tommaso on 28/03/19.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <DCPSE/DcpseRhs.hpp>
#include "../../openfpm_numerics/src/DMatrix/EMatrix.hpp"

BOOST_AUTO_TEST_SUITE(DcpseRhs_tests)

// If EIGEN is not present, EMatrix is not available and we don't need to build this test
#ifdef HAVE_EIGEN

    BOOST_AUTO_TEST_CASE(DcpseRhs_dx_test)
    {
        unsigned int r = 2;
        Point<2, unsigned int> m({1,0});
        MonomialBasis<2> mb(m.asArray(), r);
        std::cout << mb << std::endl;
        EMatrix<double, Eigen::Dynamic, 1> b(mb.size());

        DcpseRhs<2> rhs(mb, m);
        rhs.getVector<double>(b);
//        for (unsigned int i = 0; i < mb.size(); ++i)
//        {
//            std::cout << b(i) << std::endl;
//        }

        // Validation
        const MonomialBasis<2> &Dmb = mb.getDerivative(m);
        std::cout << Dmb << std::endl;
        auto p0 = Point<2, double>({0,0});
        BOOST_REQUIRE_CLOSE(b(0), -Dmb.getElement(0).evaluate(p0), 1e-16);
        BOOST_REQUIRE_CLOSE(b(1), -Dmb.getElement(1).evaluate(p0), 1e-16);
        BOOST_REQUIRE_CLOSE(b(2), -Dmb.getElement(2).evaluate(p0), 1e-16);
        BOOST_REQUIRE_CLOSE(b(3), -Dmb.getElement(3).evaluate(p0), 1e-16);
        BOOST_REQUIRE_CLOSE(b(4), -Dmb.getElement(4).evaluate(p0), 1e-16);
        BOOST_REQUIRE_CLOSE(b(5), -Dmb.getElement(5).evaluate(p0), 1e-16);
    }

    BOOST_AUTO_TEST_CASE(DcpseRhs_dxdy_test)
    {
        unsigned int r = 2;
        Point<2, unsigned int> m({1,1});
        MonomialBasis<2> mb(m.asArray(), r);
        std::cout << mb << std::endl;
        EMatrix<double, Eigen::Dynamic, Eigen::Dynamic> b(mb.size(), 1);

        DcpseRhs<2> rhs(mb, m);
        rhs.getVector<double>(b);
//        for (unsigned int i = 0; i < mb.size(); ++i)
//        {
//            std::cout << b(i) << std::endl;
//        }

        // Validation
        const MonomialBasis<2> &Dmb = mb.getDerivative(m);
        std::cout << Dmb << std::endl;
        auto p0 = Point<2, double>({0,0});
        BOOST_REQUIRE_CLOSE(b(0), Dmb.getElement(0).evaluate(p0), 1e-16);
        BOOST_REQUIRE_CLOSE(b(1), Dmb.getElement(1).evaluate(p0), 1e-16);
        BOOST_REQUIRE_CLOSE(b(2), Dmb.getElement(2).evaluate(p0), 1e-16);
        BOOST_REQUIRE_CLOSE(b(3), Dmb.getElement(3).evaluate(p0), 1e-16);
        BOOST_REQUIRE_CLOSE(b(4), Dmb.getElement(4).evaluate(p0), 1e-16);
        BOOST_REQUIRE_CLOSE(b(5), Dmb.getElement(5).evaluate(p0), 1e-16);
        BOOST_REQUIRE_CLOSE(b(6), Dmb.getElement(6).evaluate(p0), 1e-16);
        BOOST_REQUIRE_CLOSE(b(7), Dmb.getElement(7).evaluate(p0), 1e-16);
        BOOST_REQUIRE_CLOSE(b(8), Dmb.getElement(8).evaluate(p0), 1e-16);
        BOOST_REQUIRE_CLOSE(b(9), Dmb.getElement(9).evaluate(p0), 1e-16);
    }

    BOOST_AUTO_TEST_CASE(DcpseRhs_laplacian_test)
    {
        unsigned int r = 2;
        Point<2, unsigned int> m({2,2});
        MonomialBasis<2> mb(m.asArray(), r);
        std::cout << mb << std::endl;
        EMatrix<double, Eigen::Dynamic, Eigen::Dynamic> b(mb.size(), 1);

        DcpseRhs<2> rhs(mb, m);
        rhs.getVector<double>(b);

        // Validation
        const MonomialBasis<2> &Dmb = mb.getDerivative(m);
        std::cout << Dmb << std::endl;
        auto p0 = Point<2, double>({0,0});
        BOOST_REQUIRE_CLOSE(b(0), Dmb.getElement(0).evaluate(p0), 1e-16);
        BOOST_REQUIRE_CLOSE(b(1), Dmb.getElement(1).evaluate(p0), 1e-16);
        BOOST_REQUIRE_CLOSE(b(2), Dmb.getElement(2).evaluate(p0), 1e-16);
        BOOST_REQUIRE_CLOSE(b(3), Dmb.getElement(3).evaluate(p0), 1e-16);
        BOOST_REQUIRE_CLOSE(b(4), Dmb.getElement(4).evaluate(p0), 1e-16);
        BOOST_REQUIRE_CLOSE(b(5), Dmb.getElement(5).evaluate(p0), 1e-16);
        BOOST_REQUIRE_CLOSE(b(6), Dmb.getElement(6).evaluate(p0), 1e-16);
        BOOST_REQUIRE_CLOSE(b(7), Dmb.getElement(7).evaluate(p0), 1e-16);
        BOOST_REQUIRE_CLOSE(b(8), Dmb.getElement(8).evaluate(p0), 1e-16);
        BOOST_REQUIRE_CLOSE(b(9), Dmb.getElement(9).evaluate(p0), 1e-16);
    }

#endif // HAVE_EIGEN

BOOST_AUTO_TEST_SUITE_END()