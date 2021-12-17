//
// Created by landfried on 15.12.21.
//

#ifndef OPENFPM_PDATA_DOMAIN_HPP
#define OPENFPM_PDATA_DOMAIN_HPP

#include <Vector/vector_dist.hpp>

template <int dimension>
struct BoundaryCondition {
    size_t periodic[dimension]{};
    size_t non_periodic[dimension]{};

    BoundaryCondition() {
        for (int i = 0; i < dimension; i++) {
            periodic[i] = PERIODIC;
            non_periodic[i] = NON_PERIODIC;
        }
    }
};

template <int dimension, typename T>
Box<dimension, T> getDomain(T min, T max) {
    T domainMin[dimension];
    T domainMax[dimension];
    std::fill(std::begin(domainMin), std::end(domainMin), min);
    std::fill(std::begin(domainMax), std::end(domainMax), max);
    return Box<dimension, T>(domainMin, domainMax);
}


#endif //OPENFPM_PDATA_DOMAIN_HPP
