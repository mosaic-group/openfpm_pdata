//
// Created by landfried on 06.12.21.
//

#ifndef OPENFPM_PDATA_PARTICLEDATA_HPP
#define OPENFPM_PDATA_PARTICLEDATA_HPP

#include "Vector/vector_dist.hpp"
#include "GlobalVar.hpp"
#include "DataContainer.hpp"
#include "DataContainerFactory.hpp"
#include "Constants.hpp"
#include "DomainIterator.hpp"

template <typename ParticleMethodType, typename SimulationParametersType>
class ParticleData {

    using ParticleSignatureType = typename ParticleMethodType::ParticleSignature;
    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;


    using ContainerFactoryType = DataContainerFactory<ParticleSignatureType>;
    ContainerFactoryType dataContainerFactory;
    using DataContainerType = typename ContainerFactoryType::ContainerType;
    using DataStructureType = typename ContainerFactoryType::ContainerType::DataStructureType;

    DomainIterator<ParticleMethodType, SimulationParametersType> domainIterator;



    /**
    * helper function to call ghost_get for all properties up to N
    * @tparam prp pass std::make_integer_sequence<int, N>{}
    */
    template <int ... prp>
    void ghost_get_N (std::integer_sequence<int, prp...>)
    {
        getOpenFPMContainer().template ghost_get<prp...>();
    }


public:

    SimulationParametersType simulationParameters;

    DataContainerType dataContainer;

    ParticleData() : dataContainer(dataContainerFactory.template createContainer<SimulationParametersType>(simulationParameters))
    {
//        dataContainer.printType();
//        getOpenFPMContainer().map();
//        ghost_get_all();
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

    auto getParticleIterator() {
        return domainIterator.getDomainIterator(dataContainer);
    }



    /**
     * calls ghost_get for all properties
     */
    void ghost_get_all() {
        // get number of properties
        const int max_prop = ParticleSignatureType::properties::max_prop;
        // call ghost_get for all properties
        ghost_get_N(std::make_integer_sequence<int, max_prop>{});
    }

    /**
 * Executes dynamic load balancing
 */
    void dynamicLoadBalancing() {
        if (simulationParameters.dynamicLoadBalancing) {
            dataContainer.dynamicLoadBalancing();
        }
    }

};


#endif //OPENFPM_PDATA_PARTICLEDATA_HPP
