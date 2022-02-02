//
// Created by landfried on 06.12.21.
//

#ifndef OPENFPM_PDATA_PARTICLEDATA_HPP
#define OPENFPM_PDATA_PARTICLEDATA_HPP

#include <Vector/vector_dist.hpp>
#include "GlobalVar.hpp"
#include "DataContainer.hpp"
template <typename ParticleMethodType, typename SimulationParametersType>
class ParticleData {

    using ParticleSignatureType = typename ParticleMethodType::ParticleSignature;
    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

    SimulationParametersType simulationParameters;

    using ContainerFactoryType = DataContainerFactory<ParticleSignatureType>;
    ContainerFactoryType dataContainerFactory;
    using DataContainerType = typename ContainerFactoryType::ContainerType;
    using DataStructureType = typename ContainerFactoryType::ContainerType::DataStructureType;

public:

    DataContainerType dataContainer;

    ParticleData() : dataContainer(dataContainerFactory.template createContainer<SimulationParametersType>(simulationParameters))
    {
        dataContainer.printType();
    }

    // return reference to underlying OpenFPM data structure
    // vector_dist or grid_dist_id
    DataStructureType& getOpenFPMContainer() {
        return dataContainer.getContainer();
    }

    // return reference to abstract data container
    // can be free or mesh-based particle container
    DataContainerType& getDataContainer() {
        return dataContainer;
    }

    auto getParticleIterator() -> decltype(getOpenFPMContainer().getDomainIterator()){
        return getOpenFPMContainer().getDomainIterator();
    }


};


#endif //OPENFPM_PDATA_PARTICLEDATA_HPP
