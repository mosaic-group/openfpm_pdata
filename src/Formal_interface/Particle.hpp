//
// Created by landfried on 06.12.21.
//

#ifndef OPENFPM_PDATA_PARTICLE_HPP
#define OPENFPM_PDATA_PARTICLE_HPP


#include "Vector/vector_dist.hpp"
#include "memory_ly/Encap.hpp"
#include "memory_ly/memory_conf.hpp"
#include "ParticleData.hpp"
#include "DataContainer.hpp"
#include "OperationProxy.hpp"

template <typename ParticleSignatureType>
class Particle {

    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

    using DataContainerType = typename DataContainerFactory<ParticleSignatureType>::ContainerType;
    using DataKeyType = typename DataContainerFactory<ParticleSignatureType>::KeyType;

    PropertyType propertyAggregate;

protected:
    DataKeyType key;

    DataContainerType& dataContainer;

public:
    Particle(DataContainerType& dataContainer_in, DataKeyType key_in) : dataContainer(dataContainer_in), key(key_in) {}

    template<unsigned int id>
    inline auto property() -> decltype(dataContainer.template property<id>(key)) {
        return dataContainer.template property<id>(key);
    }

    template<unsigned int id>
    inline auto property_vec() {
        return dataContainer.template property_vec<id>(key);
    }

    inline auto position() -> decltype(dataContainer.position(key)) {
        return dataContainer.position(key);
    }


    inline auto position_vec() -> decltype(dataContainer.position_vec(key)) {
        return dataContainer.position_vec(key);
    }

    auto getID() -> decltype(key.getKey()) {
        return key.getKey();
    }

    bool operator== (Particle<ParticleSignatureType> rhs) {
        return getID() == rhs.getID();
    }

    bool operator!= (Particle<ParticleSignatureType> rhs) {
        return getID() != rhs.getID();
    }




};





#endif //OPENFPM_PDATA_PARTICLE_HPP
