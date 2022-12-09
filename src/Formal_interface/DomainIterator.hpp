//
// Created by landfried on 07.02.22.
//

#ifndef OPENFPM_PDATA_DOMAINITERATOR_HPP
#define OPENFPM_PDATA_DOMAINITERATOR_HPP

#include <iostream>
#include "Constants.hpp"
#include "DataContainerFactory.hpp"


/**
 * Primary template
 * Wraps access to domain iterator, which depends on the data structure type and interaction type
 * @tparam ParticleMethodType
 * @tparam SimulationParametersType
 * @tparam DataStructureType can be FREE (vector_Dist) or MESH (grid_dist_id), which have different method names
 * @tparam InteractionType determines if iterator has to iterate through ghost particles or not
 */
template <typename ParticleMethodType, typename SimulationParametersType,
        typename DataStructureType = typename ParticleMethodType::ParticleSignature::dataStructure,
        int InteractionType = SimulationParametersType::interactionType>
class DomainIterator {
    using DataContainerType = typename DataContainerFactory<typename ParticleMethodType::ParticleSignature>::ContainerType;
public:
    auto getDomainIterator(DataContainerType &dataContainer) = 0;

    bool executeGhostPut() {
        return false;
    }
};




// Free particles
// = vector_dist

template <typename ParticleMethodType, typename SimulationParametersType>
class DomainIterator <ParticleMethodType, SimulationParametersType, FREE_PARTICLES, INTERACTION_PULL> {
    using DataContainerType = typename DataContainerFactory<typename ParticleMethodType::ParticleSignature>::ContainerType;
public:
    DomainIterator() {std::cout << "DomainIterator: FREE, PULL" << std::endl;}
    auto getDomainIterator(DataContainerType &dataContainer) {
        return dataContainer.getContainer().getDomainIterator();
    }

    bool executeGhostPut() {
        return false;
    }
};

template <typename ParticleMethodType, typename SimulationParametersType>
class DomainIterator <ParticleMethodType, SimulationParametersType, FREE_PARTICLES, INTERACTION_PUSH> {
    using DataContainerType = typename DataContainerFactory<typename ParticleMethodType::ParticleSignature>::ContainerType;
public:
    DomainIterator() {std::cout << "DomainIterator: FREE, PUSH" << std::endl;}

    auto getDomainIterator(DataContainerType &dataContainer) {
        return dataContainer.getContainer().getDomainAndGhostIterator();
    }

    bool executeGhostPut() {
        return false;
    }
};

template <typename ParticleMethodType, typename SimulationParametersType>
class DomainIterator <ParticleMethodType, SimulationParametersType, FREE_PARTICLES, INTERACTION_PULL_PUSH> {
    using DataContainerType = typename DataContainerFactory<typename ParticleMethodType::ParticleSignature>::ContainerType;
public:
    DomainIterator() {std::cout << "DomainIterator: FREE, PULL_PUSH" << std::endl;}
    auto getDomainIterator(DataContainerType &dataContainer) {
        return dataContainer.getContainer().getDomainAndGhostIterator();
    }

    bool executeGhostPut() {
        return false;
    }
};

template <typename ParticleMethodType, typename SimulationParametersType>
class DomainIterator <ParticleMethodType, SimulationParametersType, FREE_PARTICLES, INTERACTION_SYMMETRIC> {
    using DataContainerType = typename DataContainerFactory<typename ParticleMethodType::ParticleSignature>::ContainerType;
public:
    DomainIterator() {std::cout << "DomainIterator: FREE, SYMMETRIC" << std::endl;}
    auto getDomainIterator(DataContainerType &dataContainer) {
        return dataContainer.getContainer().getDomainAndGhostIterator();
    }

    bool executeGhostPut() {
        return false;
    }
};




// Mesh particles
// = grid_dist_id

template <typename ParticleMethodType, typename SimulationParametersType>
class DomainIterator <ParticleMethodType, SimulationParametersType, MESH_PARTICLES, INTERACTION_PULL> {
    using DataContainerType = typename DataContainerFactory<typename ParticleMethodType::ParticleSignature>::ContainerType;
public:
    DomainIterator() {std::cout << "DomainIterator: MESH, PULL" << std::endl;}

    auto getDomainIterator(DataContainerType &dataContainer) {
        return dataContainer.getContainer().getDomainIterator();
    }

    bool executeGhostPut() {
        return false;
    }
};

template <typename ParticleMethodType, typename SimulationParametersType>
class DomainIterator <ParticleMethodType, SimulationParametersType, MESH_PARTICLES, INTERACTION_PUSH> {
    using DataContainerType = typename DataContainerFactory<typename ParticleMethodType::ParticleSignature>::ContainerType;
public:
    DomainIterator() {std::cout << "DomainIterator: MESH, PUSH" << std::endl;}

    auto getDomainIterator(DataContainerType &dataContainer) {
        return dataContainer.getContainer().getDomainGhostIterator();
    }

    bool executeGhostPut() {
        return false;
    }
};

template <typename ParticleMethodType, typename SimulationParametersType>
class DomainIterator <ParticleMethodType, SimulationParametersType, MESH_PARTICLES, INTERACTION_PULL_PUSH> {
    using DataContainerType = typename DataContainerFactory<typename ParticleMethodType::ParticleSignature>::ContainerType;
public:
    DomainIterator() {std::cout << "DomainIterator: MESH, PULL_PUSH" << std::endl;}

    auto getDomainIterator(DataContainerType &dataContainer) {
        return dataContainer.getContainer().getDomainGhostIterator();
    }

    bool executeGhostPut() {
        return false;
    }
};

template <typename ParticleMethodType, typename SimulationParametersType>
class DomainIterator <ParticleMethodType, SimulationParametersType, MESH_PARTICLES, INTERACTION_SYMMETRIC> {
    using DataContainerType = typename DataContainerFactory<typename ParticleMethodType::ParticleSignature>::ContainerType;
public:
    DomainIterator() {std::cout << "DomainIterator: MESH, SYMMETRIC" << std::endl;}

    auto getDomainIterator(DataContainerType &dataContainer) {
        return dataContainer.getContainer().getDomainIterator();
    }

    bool executeGhostPut() {
        return true;
    }
};

#endif //OPENFPM_PDATA_DOMAINITERATOR_HPP
