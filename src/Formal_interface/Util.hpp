//
// Created by landfried on 14.11.22.
//

#ifndef OPENFPM_PDATA_UTIL_HPP
#define OPENFPM_PDATA_UTIL_HPP

#include "Vector/vector_dist.hpp"


template<typename positionType, unsigned int dimension>
positionType abs2(Point<dimension, positionType> point) {
    return point.distance2(point.zero_p());
}

#endif //OPENFPM_PDATA_UTIL_HPP
