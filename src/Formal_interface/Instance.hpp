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

    template <unsigned int dimension, typename PositionType>
    auto skinIterator(Point<dimension, PositionType> pointMin, Point<dimension, PositionType> pointMax, PositionType width,  size_t (& sz)[dimension]) {
        Box<dimension, PositionType> drawOuterBox(pointMin, pointMax);
        Box<dimension, PositionType> drawInnerBox(pointMin + width, pointMax - width);
        Box<dimension, PositionType> domain (particleData.simulationParameters.domainMin, particleData.simulationParameters.domainMax);

        auto boxIterator = DrawParticles::DrawSkin(particleData.getOpenFPMContainer(),sz,domain,drawInnerBox, drawOuterBox);
        return boxIterator;
    }

    template <unsigned int dimension, typename PositionType>
    auto skinOpenIterator(Point<dimension, PositionType> pointMin, Point<dimension, PositionType> pointMax, PositionType width,  size_t (& sz)[dimension]) {
        Point<dimension, PositionType> upperCutoutPoint = pointMax - width;
        upperCutoutPoint[dimension - 1] = pointMax[dimension - 1];

        Box<dimension, PositionType> drawOuterBox(pointMin, pointMax);
        Box<dimension, PositionType> drawInnerBox(pointMin + width, upperCutoutPoint);
        Box<dimension, PositionType> domain (particleData.simulationParameters.domainMin, particleData.simulationParameters.domainMax);

        auto boxIterator = DrawParticles::DrawSkin(particleData.getOpenFPMContainer(),sz,domain,drawInnerBox, drawOuterBox);
        return boxIterator;
    }

public:

    Instance(ParticleData<ParticleMethodType, SimulationParametersType> &particleData_in) : particleData(particleData_in) {

    }

    virtual void freePlacement() {}

    virtual void shapePlacement() {}

    virtual void setGlobalVariable() {}

};


#endif //OPENFPM_PDATA_INSTANCE_HPP
