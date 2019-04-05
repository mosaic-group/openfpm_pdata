//
// Created by tommaso on 1/04/19.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <Vector/vector_dist.hpp>
#include <DCPSE/Dcpse.hpp>

template<typename T>
void check_small_or_close(T value, T expected, T tolerance)
{
    if (fabs(expected) < tolerance)
    {
        BOOST_CHECK_SMALL(value, tolerance);
    } else
    {
        BOOST_CHECK_CLOSE(value, expected, tolerance);
    }
}

template<typename T>
void check_small_or_close_abs(T value, T expected, T absTolerance)
{
    BOOST_CHECK_SMALL(value - expected, absTolerance);
}


BOOST_AUTO_TEST_SUITE(Dcpse_tests)

// If EIGEN is not present, EMatrix is not available and we don't need to build this test
#ifdef HAVE_EIGEN

    BOOST_AUTO_TEST_CASE(Dcpse_2D_test)
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        BOOST_TEST_MESSAGE("Init vars...");

        // Here build some easy domain and get some points around a given one
        size_t edgeSemiSize = 20;
        const size_t sz[2] = {2 * edgeSemiSize, 2 * edgeSemiSize};
        Box<2, double> box({0, 0}, {2 * M_PI, 2 * M_PI});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing[2];
        spacing[0] = 1.0 / (sz[0] - 1);
        spacing[1] = 1.0 / (sz[1] - 1);
        Ghost<2, double> ghost(0.1);

        double rCut = 2 * spacing[0];

        BOOST_TEST_MESSAGE("Init vector_dist...");
        vector_dist<2, double, aggregate<double, double, double>> domain(0, box, bc, ghost);

        if (rank == 0)
        {
            BOOST_TEST_MESSAGE("Init domain...");
            auto it = domain.getGridIterator(sz);
            size_t pointId = 0;
            size_t counter = 0;
            double minNormOne = 999;
            while (it.isNext())
            {
                domain.add();
                auto key = it.get();
                mem_id k0 = key.get(0);
                double x = k0 * spacing[0];
                domain.getLastPos()[0] = x;
                mem_id k1 = key.get(1);
                double y = k1 * spacing[1];
                domain.getLastPos()[1] = y;
                // Here fill the function value
                domain.template getLastProp<0>() = sin(x);
//            domain.template getLastProp<0>() = x * x;
//            domain.template getLastProp<0>() = x;
                // Here fill the validation value for Df/Dx
                domain.template getLastProp<2>() = cos(x);
//            domain.template getLastProp<2>() = 2 * x;
//            domain.template getLastProp<2>() = 1;

                ++counter;
                ++it;
            }
            BOOST_TEST_MESSAGE("Sync domain across processors...");
        }
        domain.map(); // Send particles to the right processors

        BOOST_TEST_MESSAGE("Getting ghost...");
        domain.ghost_get();

        // Setup finished, actual test below...
//        std::cout << "rCut = " << rCut << std::endl;
        BOOST_TEST_MESSAGE("DCPSE init & compute coefficients...");
        Dcpse<2, double, double, double, double> dcpse(domain, Point<2, unsigned int>({1, 0}), 2, rCut);
        BOOST_TEST_MESSAGE("DCPSE compute diff operator...");
        dcpse.template computeDifferentialOperator<0, 1>(domain);

        // Now check against the validation values
        BOOST_TEST_MESSAGE("Validating against ground truth...");
//        const double TOL = 1e-6;
        const double avgSpacing = spacing[0] + spacing[1];
        const double TOL = 2 * avgSpacing * avgSpacing;
        auto itVal = domain.getDomainIterator();
        bool check = true;
        double computedValue;
        double validationValue;
        while (itVal.isNext())
        {
            auto key = itVal.get();
            computedValue = domain.template getProp<1>(key);
            validationValue = domain.template getProp<2>(key);
            bool locCheck = (fabs(computedValue - validationValue) < TOL);
            check = check && locCheck;
            if (!locCheck)
            {
                std::cout << "FAILED CHECK :: pos=" << Point<2, double>(domain.getPos(key)).toString()
                          << ", computedValue=" << computedValue
                          << ", validationValue=" << validationValue
                          << ", difference=" << fabs(computedValue - validationValue)
                          << ", tolerance=" << TOL
                          << std::endl;
//                break;
            }
            ++itVal;
        }
        BOOST_REQUIRE(check);
    }

    BOOST_AUTO_TEST_CASE(Dcpse_2D_perturbed_test)
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        BOOST_TEST_MESSAGE("Init vars...");

        // Here build some easy domain and get some points around a given one
        size_t edgeSemiSize = 20;
        const size_t sz[2] = {2 * edgeSemiSize, 2 * edgeSemiSize};
        Box<2, double> box({0, 0}, {2 * M_PI, 2 * M_PI});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing[2];
        spacing[0] = 1.0 / (sz[0] - 1);
        spacing[1] = 1.0 / (sz[1] - 1);
        Ghost<2, double> ghost(0.1);

        double sigma2 = spacing[0] * spacing[1] / (2 * 4);
        std::cout
                << "sigma2 = " << sigma2
                << ", spacing[0] = " << spacing[0]
                << std::endl;
        double rCut = 2 * (spacing[0] + sqrt(sigma2));

        BOOST_TEST_MESSAGE("Init vector_dist...");
        vector_dist<2, double, aggregate<double, double, double>> domain(0, box, bc, ghost);

        if (rank == 0)
        {
            BOOST_TEST_MESSAGE("Init domain...");
//            std::random_device rd{};
//            std::mt19937 rng{rd()};
            std::mt19937 rng{6666666};

            std::normal_distribution<> gaussian{0, sigma2};

            auto it = domain.getGridIterator(sz);
            size_t pointId = 0;
            size_t counter = 0;
            double minNormOne = 999;
            while (it.isNext())
            {
                domain.add();
                auto key = it.get();
                mem_id k0 = key.get(0);
                double x = k0 * spacing[0];
                domain.getLastPos()[0] = x + gaussian(rng);
                mem_id k1 = key.get(1);
                double y = k1 * spacing[1];
                domain.getLastPos()[1] = y + gaussian(rng);
                // Here fill the function value
                domain.template getLastProp<0>() = sin(x);
//            domain.template getLastProp<0>() = x * x;
//            domain.template getLastProp<0>() = x;
                // Here fill the validation value for Df/Dx
                domain.template getLastProp<2>() = cos(x);
//            domain.template getLastProp<2>() = 2 * x;
//            domain.template getLastProp<2>() = 1;

                ++counter;
                ++it;
            }
            BOOST_TEST_MESSAGE("Sync domain across processors...");

        }
        domain.map(); // Send particles to the right processors

        // Setup finished, actual test below...
//        std::cout << "rCut = " << rCut << std::endl;
        BOOST_TEST_MESSAGE("DCPSE init & compute coefficients...");
        Dcpse<2, double, double, double, double> dcpse(domain, Point<2, unsigned int>({1, 0}), 2, rCut, 2);
        BOOST_TEST_MESSAGE("DCPSE compute diff operator...");
        dcpse.template computeDifferentialOperator<0, 1>(domain);

        // Now check against the validation values
        BOOST_TEST_MESSAGE("Validating against ground truth...");
//        const double TOL = 1e-6;
        const double avgSpacing = spacing[0] + spacing[1] + 2 * sqrt(sigma2);
        const double TOL = 2 * avgSpacing * avgSpacing;
        auto itVal = domain.getDomainIterator();
        bool check = true;
        double computedValue;
        double validationValue;
        while (itVal.isNext())
        {
            auto key = itVal.get();
            computedValue = domain.template getProp<1>(key);
            validationValue = domain.template getProp<2>(key);
            bool locCheck = (fabs(computedValue - validationValue) < TOL);
            check = check && locCheck;
            if (!locCheck)
            {
                std::cout << "FAILED CHECK :: pos=" << Point<2, double>(domain.getPos(key)).toString()
                          << ", computedValue=" << computedValue
                          << ", validationValue=" << validationValue
                          << ", difference=" << fabs(computedValue - validationValue)
                          << ", tolerance=" << TOL
                          << std::endl;
//                break;
            }
            ++itVal;
        }
        BOOST_REQUIRE(check);
    }

    BOOST_AUTO_TEST_CASE(Dcpse_3D_test)
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        BOOST_TEST_MESSAGE("Init vars...");

        // Here build some easy domain and get some points around a given one
        size_t edgeSemiSize = 20;
        const unsigned int DIM = 3;
        const size_t sz[DIM] = {2 * edgeSemiSize, 2 * edgeSemiSize, 2 * edgeSemiSize};
        Box<DIM, double> box({0, 0, 0}, {2 * M_PI, 2 * M_PI, 2 * M_PI});
        size_t bc[DIM] = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};
        double spacing[DIM];
        spacing[0] = 1.0 / (sz[0] - 1);
        spacing[1] = 1.0 / (sz[1] - 1);
        spacing[2] = 1.0 / (sz[2] - 1);
        Ghost<DIM, double> ghost(0.1);

        double rCut = 2 * spacing[0];

        BOOST_TEST_MESSAGE("Init vector_dist...");
        vector_dist<DIM, double, aggregate<double, double, double>> domain(0, box, bc, ghost);

        if (rank == 0)
        {
            BOOST_TEST_MESSAGE("Init domain...");
            auto it = domain.getGridIterator(sz);
            size_t pointId = 0;
            size_t counter = 0;
            double minNormOne = 999;
            while (it.isNext())
            {
                domain.add();
                auto key = it.get();
                mem_id k0 = key.get(0);
                double x = k0 * spacing[0];
                domain.getLastPos()[0] = x;
                mem_id k1 = key.get(1);
                double y = k1 * spacing[1];
                domain.getLastPos()[1] = y;
                mem_id k2 = key.get(2);
                double z = k2 * spacing[2];
                domain.getLastPos()[2] = z;
                // Here fill the function value
                domain.template getLastProp<0>() = sin(z);
//            domain.template getLastProp<0>() = x * x;
//            domain.template getLastProp<0>() = x;
                // Here fill the validation value for Df/Dx
                domain.template getLastProp<2>() = cos(z);
//            domain.template getLastProp<2>() = 2 * x;
//            domain.template getLastProp<2>() = 1;

                ++counter;
                ++it;
            }
            BOOST_TEST_MESSAGE("Sync domain across processors...");
        }
        domain.map(); // Send particles to the right processors

        // Setup finished, actual test below...
        BOOST_TEST_MESSAGE("DCPSE init & compute coefficients...");
        Dcpse<DIM, double, double, double, double> dcpse(domain, Point<DIM, unsigned int>({0, 0, 1}), 2, rCut);
        BOOST_TEST_MESSAGE("DCPSE compute diff operator...");
        dcpse.template computeDifferentialOperator<0, 1>(domain);

        // Now check against the validation values
        BOOST_TEST_MESSAGE("Validating against ground truth...");
//        const double TOL = 1e-6;
        const double avgSpacing = spacing[0] + spacing[1] + spacing[2];
        const double TOL = avgSpacing * avgSpacing;
        auto itVal = domain.getDomainIterator();
        bool check = true;
        double computedValue;
        double validationValue;
        while (itVal.isNext())
        {
            auto key = itVal.get();
            computedValue = domain.template getProp<1>(key);
            validationValue = domain.template getProp<2>(key);
            check = check && (fabs(computedValue - validationValue) < TOL);
            if (!check)
            {
                std::cout << "FAILED CHECK :: pos=" << Point<DIM, double>(domain.getPos(key)).toString()
                          << ", computedValue=" << computedValue
                          << ", validationValue=" << validationValue
                          << ", tolerance=" << TOL
                          << std::endl;
                break;
            }
            ++itVal;
        }
        BOOST_REQUIRE(check);
    }

    BOOST_AUTO_TEST_CASE(Dcpse_3D_2_test)
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        BOOST_TEST_MESSAGE("Init vars...");

        // Here build some easy domain and get some points around a given one
        size_t edgeSemiSize = 10;
        const unsigned int DIM = 3;
        const size_t sz[DIM] = {2 * edgeSemiSize, 2 * edgeSemiSize, 2 * edgeSemiSize};
        Box<DIM, double> box({0, 0, 0}, {2 * M_PI, 2 * M_PI, 2 * M_PI});
        size_t bc[DIM] = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};
        double spacing[DIM];
        spacing[0] = 1.0 / (sz[0] - 1);
        spacing[1] = 1.0 / (sz[1] - 1);
        spacing[2] = 1.0 / (sz[2] - 1);
        Ghost<DIM, double> ghost(0.1);

        double rCut = 2 * spacing[0];

        BOOST_TEST_MESSAGE("Init vector_dist...");
        vector_dist<DIM, double, aggregate<double, double, double>> domain(0, box, bc, ghost);

        if (rank == 0)
        {
            BOOST_TEST_MESSAGE("Init domain...");
            auto it = domain.getGridIterator(sz);
            size_t pointId = 0;
            size_t counter = 0;
            double minNormOne = 999;
            while (it.isNext())
            {
                domain.add();
                auto key = it.get();
                mem_id k0 = key.get(0);
                double x = k0 * spacing[0];
                domain.getLastPos()[0] = x;
                mem_id k1 = key.get(1);;
                double y = k1 * spacing[1];
                domain.getLastPos()[1] = y;
                mem_id k2 = key.get(2);
                double z = k2 * spacing[2];
                domain.getLastPos()[2] = z;
                // Here fill the function value
                domain.template getLastProp<0>() = x * x * sin(z);
                // Here fill the validation value for Df/Dx
                domain.template getLastProp<2>() = -2 * sin(z);

                ++counter;
                ++it;
            }
            BOOST_TEST_MESSAGE("Sync domain across processors...");
        }
        domain.map(); // Send particles to the right processors

        // Setup finished, actual test below...
        unsigned int r = 2;
        BOOST_TEST_MESSAGE("DCPSE init & compute coefficients...");
        Dcpse<DIM, double, double, double, double> dcpse(domain, Point<DIM, unsigned int>({2, 0, 2}), r, rCut);
        BOOST_TEST_MESSAGE("DCPSE compute diff operator...");
        dcpse.template computeDifferentialOperator<0, 1>(domain);

        // Now check against the validation values
        BOOST_TEST_MESSAGE("Validating against ground truth...");
//        const double TOL = 1e-6;
        const double avgSpacing = spacing[0] + spacing[1] + spacing[2];
        const double TOL = pow(avgSpacing, r);
        auto itVal = domain.getDomainIterator();
        bool check = true;
        double computedValue;
        double validationValue;
        while (itVal.isNext())
        {
            auto key = itVal.get();
            computedValue = domain.template getProp<1>(key);
            validationValue = domain.template getProp<2>(key);
            check = check && (fabs(computedValue - validationValue) < TOL);
            if (!check)
            {
                std::cout << "FAILED CHECK :: pos=" << Point<DIM, double>(domain.getPos(key)).toString()
                          << ", computedValue=" << computedValue
                          << ", validationValue=" << validationValue
                          << ", tolerance=" << TOL
                          << std::endl;
                break;
            }
            ++itVal;
        }
        BOOST_REQUIRE(check);
    }

#endif // HAVE_EIGEN

BOOST_AUTO_TEST_SUITE_END()
