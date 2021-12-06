//
// Created by landfried on 06.12.21.
//

#ifndef OPENFPM_PDATA_PARTICLEDATA_HPP
#define OPENFPM_PDATA_PARTICLEDATA_HPP

#include <Vector/vector_dist.hpp>


template <typename ParticleType>
class ParticleData {

    float r_cut = 0.5;
    Box<1,float> box;
    size_t bc[1]={PERIODIC};
    Ghost<1,float> ghost;

public:
    vector_dist<1, float, ParticleType> vd;

    ParticleData() : box({0.0},{5.0}), ghost(r_cut), vd(0,box,bc,ghost) {}

    template<unsigned int id> inline auto getProp(vect_dist_key_dx p) -> decltype(vd.template getProp<id>(p)) {
        return vd.template getProp<id>(p);
    }

    inline auto getPos(vect_dist_key_dx p) -> decltype(vd.getPos(p)) {
        return vd.getPos(p);
    }
};




#endif //OPENFPM_PDATA_PARTICLEDATA_HPP
