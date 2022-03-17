//
// Created by landfried on 14.03.22.
//

#ifndef OPENFPM_PDATA_INSTANCE_HPP
#define OPENFPM_PDATA_INSTANCE_HPP

#include "iostream"
#include "ParticleData.hpp"

template <typename ParticleMethodType, typename SimulationParametersType>
class Instance {

protected:

    ParticleData<ParticleMethodType, SimulationParametersType> &particleData;

    void addParticle() {
        particleData.getOpenFPMContainer().add();
    }

    auto position() {
        return particleData.getOpenFPMContainer().getLastPos();
    }

    template<int dimension>
    auto property() {
        return particleData.getOpenFPMContainer().template getLastProp<dimension>();
    }

public:

    Instance(ParticleData<ParticleMethodType, SimulationParametersType> &particleData_in) : particleData(particleData_in) {

    }

    virtual void freePlacement() {}


};


#endif //OPENFPM_PDATA_INSTANCE_HPP
