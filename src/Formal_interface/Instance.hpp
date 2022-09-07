//
// Created by landfried on 14.03.22.
//

#ifndef OPENFPM_PDATA_INSTANCE_HPP
#define OPENFPM_PDATA_INSTANCE_HPP

#include "iostream"
#include "ParticleData.hpp"
#include "../../openfpm_numerics/src/Draw/DrawParticles.hpp"

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
    auto property() -> decltype(particleData.getOpenFPMContainer().template getLastProp<dimension>()) {
        return particleData.getOpenFPMContainer().template getLastProp<dimension>();
    }

    template <unsigned int dimension, typename PositionType>
    auto boxIterator(Point<dimension, PositionType> pointMin, Point<dimension, PositionType> pointMax, size_t (& sz)[dimension]) {
        Box<dimension, PositionType> drawBox(pointMin, pointMax);
        Box<dimension, PositionType> domain (particleData.simulationParameters.domainMin, particleData.simulationParameters.domainMax);

        auto boxIterator = DrawParticles::DrawBox(particleData.getOpenFPMContainer(),sz,domain,drawBox);
        return boxIterator;
    }

public:

    Instance(ParticleData<ParticleMethodType, SimulationParametersType> &particleData_in) : particleData(particleData_in) {

    }

    virtual void freePlacement() {}

    virtual void shapePlacement() {}

};


#endif //OPENFPM_PDATA_INSTANCE_HPP
