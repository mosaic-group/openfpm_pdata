//
// Created by tommaso on 26/03/19.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include "Vector/vector_dist.hpp"
#include <Space/Shape/Point.hpp>
#include "DCPSE/SupportBuilder.hpp"

BOOST_AUTO_TEST_SUITE(Support_tests)

    BOOST_AUTO_TEST_CASE(SupportBuilder_2D_1_0_2spacing_test)
    {
        // Here build some easy domain and get some points around a given one
        size_t edgeSemiSize = 100;
        const size_t sz[2] = {2 * edgeSemiSize, 2 * edgeSemiSize};
        Box<2, double> box({-1.1, -1.1}, {1.1, 1.1});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing[2];
        spacing[0] = 1.0 / (sz[0] - 1);
        spacing[1] = 1.0 / (sz[1] - 1);
        Ghost<2, double> ghost(0.1);

        vector_dist<2, double, aggregate<double>> domain(0, box, bc, ghost);
        auto it = domain.getGridIterator(sz);
        size_t pointId = 0;
        size_t counter = 0;
        double minNormOne = 999;
        while (it.isNext())
        {
            domain.add();
            auto key = it.get();
            mem_id k0 = key.get(0) - edgeSemiSize;
            double x = k0 * spacing[0];
            domain.getLastPos()[0] = x;
            mem_id k1 = key.get(1) - edgeSemiSize;
            double y = k1 * spacing[1];
            domain.getLastPos()[1] = y;
            domain.template getLastProp<0>() = 0.0;

            ++counter;
            ++it;
        }
        // Now get iterator to point of interest
        auto itPoint = domain.getIterator();
        size_t foo = (sqrt(counter) + counter) / 2;
        for (int i = 0; i < foo; ++i)
        {
            ++itPoint;
        }
        // Get spatial position from point iterator
        vect_dist_key_dx p = itPoint.get();
        const auto pos = domain.getPos(p.getKey());
        std::cout << "p=(" << pos[0] << "," << pos[1] << ")" << std::endl;
//        BOOST_REQUIRE_CLOSE(pos[0], 0, 1e-16);
//        BOOST_REQUIRE_CLOSE(pos[1], 0, 1e-16);

        // Now that domain is built and populated, let's test SupportBuilder
        // We use (0,0) as initial point
        SupportBuilder<2, double, aggregate<double>> supportBuilder(domain, {1,0}, 2*spacing[0]);
        auto support = supportBuilder.getSupport(itPoint, 6);
//        for (const auto &off : support.getOffsets())
//        {
//            std::cout << off.toString() << std::endl;
//        }
        BOOST_REQUIRE_GE(support.size(), 6);
    }

    BOOST_AUTO_TEST_CASE(SupportBuilder_2D_2_2_2spacing_test)
    {
        // Here build some easy domain and get some points around a given one
        size_t edgeSemiSize = 25;
        const size_t sz[2] = {2 * edgeSemiSize, 2 * edgeSemiSize};
        Box<2, double> box({-1.2, -1.2}, {1.2, 1.2});
//        Box<2, double> innerDomain({-1.0, -1.0}, {1.0, 1.0});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing[2];
        spacing[0] = 1.0 / (sz[0] - 1);
        spacing[1] = 1.0 / (sz[1] - 1);
        Ghost<2, double> ghost(3.0 * spacing[0]);

        vector_dist<2, double, aggregate<double>> domain(0, box, bc, ghost);
        auto it = domain.getGridIterator(sz);
        size_t pointId = 0;
        size_t counter = 0;
        while (it.isNext())
        {
            domain.add();
            auto key = it.get();
            mem_id k0 = key.get(0) - edgeSemiSize;
            double x = k0 * spacing[0];
            domain.getLastPos()[0] = x;
            mem_id k1 = key.get(1) - edgeSemiSize;
            double y = k1 * spacing[1];
            domain.getLastPos()[1] = y;
            domain.template getLastProp<0>() = 0.0;
            if (abs(domain.getLastPos()[0]) + abs(domain.getLastPos()[1]) < 1e-16) // i.e. if we are at (0,0)
            {
                pointId = counter;
            }
            ++counter;
            ++it;
        }
        // Now get iterator to point of interest
        auto itPoint = domain.getDomainIterator();
        size_t foo = (sqrt(counter) + counter) / 2;
        for (int i = 0; i < foo; ++i)
        {
            ++itPoint;
        }
        // Get spatial position from point iterator
        vect_dist_key_dx p = itPoint.get();
        const auto pos = domain.getPos(p.getKey());
        std::cout << "p=(" << pos[0] << "," << pos[1] << ")" << std::endl;
//        BOOST_REQUIRE_CLOSE(pos[0], 0, 1e-16);
//        BOOST_REQUIRE_CLOSE(pos[1], 0, 1e-16);

        // Now that domain is built and populated, let's test SupportBuilder
        // We use (0,0) as initial point
        SupportBuilder<2, double, aggregate<double>> supportBuilder(domain, {2,2}, 2*spacing[0]);
        auto supportPoints = supportBuilder.getSupport(itPoint, 20);
//        for (const auto &k : supportPoints)
//        {
//            Point<2, double> pt = domain.getPos(k);
//            std::cout << pt.toString() << std::endl;
//        }
        BOOST_REQUIRE_GE(supportPoints.size(), 20);
    }

//    BOOST_AUTO_TEST_CASE(Support_CopyConstructor_test)
//    {
//
//    }

BOOST_AUTO_TEST_SUITE_END()