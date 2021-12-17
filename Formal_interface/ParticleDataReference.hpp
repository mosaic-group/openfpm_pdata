//
// Created by landfried on 09.12.21.
//

#ifndef OPENFPM_PDATA_PARTICLEDATAREFERENCE_HPP
#define OPENFPM_PDATA_PARTICLEDATAREFERENCE_HPP


#include <Vector/vector_dist.hpp>


template <int dimension, typename PositionType, typename ParticleType>
class ParticleDataReference {

public:
    vector_dist<dimension, PositionType, ParticleType>& vd_ref;

    ParticleDataReference(vector_dist<dimension, PositionType, ParticleType>& vd_in) : vd_ref(vd_in) {}

    template<unsigned int id> inline auto getProp(vect_dist_key_dx p) -> decltype(vd_ref.template getProp<id>(p)) {
        return vd_ref.template getProp<id>(p);
    }

    inline auto getPos(vect_dist_key_dx p) -> decltype(vd_ref.getPos(p)) {
        return vd_ref.getPos(p);
    }
};






#endif //OPENFPM_PDATA_PARTICLEDATAREFERENCE_HPP
