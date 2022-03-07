//
// Created by landfried on 13.01.22.
//

#ifndef OPENFPM_PDATA_INTERACTION_HPP
#define OPENFPM_PDATA_INTERACTION_HPP

#include <iostream>
#include "ParticleData.hpp"
#include "Particle.hpp"
#include "Constants.hpp"
#include "CellList.hpp"

/**
 *
 * @tparam NeighborhoodDetermination
 * @tparam ParticleMethodType
 * @tparam SimulationParametersType
 */
template <typename NeighborhoodDetermination, typename ParticleMethodType, typename SimulationParametersType>
class Interaction_Impl {
public:

    explicit Interaction_Impl(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {}

    void executeInteraction(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {}
};



template <typename ParticleMethodType, typename SimulationParametersType>
class Interaction_Impl<NEIGHBORHOOD_ALLPARTICLES, ParticleMethodType, SimulationParametersType> {

    using ParticleSignatureType = typename ParticleMethodType::ParticleSignature;
    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

    ParticleMethodType particleMethod;
    SimulationParametersType simulationParameters;


public:

    explicit Interaction_Impl(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {
    }


    void executeInteraction(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {


        // iterate through all particles
        auto iteratorAll = particleData.getParticleIterator();

        while (iteratorAll.isNext())
        {
            auto p = iteratorAll.get();
            Particle<ParticleSignatureType> particle(particleData.dataContainer, p);

            // iterate through all particles as neighbors
            auto iteratorNeighbors = particleData.getParticleIterator();

            while (iteratorNeighbors.isNext()) {
                auto n = iteratorNeighbors.get();
                Particle<ParticleSignatureType> neighbor(particleData.dataContainer, n);

                if (particle != neighbor) {
                    particleMethod.interact(particle, neighbor);
                }
                ++iteratorNeighbors;
            }
            ++iteratorAll;
        }
    }
};

/*

template <typename ParticleMethodType, typename SimulationParametersType>
class Interaction_Impl<NEIGHBORHHOD_CELLLIST, ParticleMethodType, SimulationParametersType> {

    using ParticleSignatureType = typename ParticleMethodType::ParticleSignature;
    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

    ParticleMethodType particleMethod;
    SimulationParametersType simulationParameters;

    using CellListType = typename vector_dist<dimension, PositionType, PropertyType>::CellList_type;

    CellListType cellList;

public:

    explicit Interaction_Impl(ParticleData<ParticleMethodType, SimulationParametersType> &particleData): cellList(
            createCellList(particleData)) {



        auto & vcl = create_vcluster();
        std::cout << "interaction impl cell list constructed " << vcl.getProcessUnitID() << std::endl;

    }

    CellListType createCellList(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {
        auto & vcl = create_vcluster();

        std::cout << "create cell list " << vcl.getProcessUnitID() << std::endl;

        // symmetric cell list
        if (simulationParameters.interactionType == INTERACTION_SYMMETRIC) {
            std::cout << "symmetric " << vcl.getProcessUnitID() << std::endl;
            return particleData.getOpenFPMContainer().getCellListSym(simulationParameters.cellWidth);

        }

        // unsymmetric cell list
        std::cout << "unsymmetric" << std::endl;
        return particleData.getOpenFPMContainer().getCellList(simulationParameters.cellWidth);
    }


    void executeInteraction(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {

        particleData.getOpenFPMContainer().map();
        particleData.ghost_get_all();

        // update symmetric cell list
        if (simulationParameters.interactionType == INTERACTION_SYMMETRIC) {
//            particleData.getOpenFPMContainer().template updateCellListSym(cellList);
            particleData.getOpenFPMContainer().updateCellListSym(cellList);
        }
        else {
            particleData.getOpenFPMContainer().updateCellList(cellList);
        }




        // iterate through all particles
        auto iteratorAll = particleData.getParticleIterator();
        while (iteratorAll.isNext())
        {
            auto p = iteratorAll.get();
            Particle<ParticleSignatureType> particle(particleData.dataContainer, p);

            // iterate through all particles in neighbor cells
            auto iteratorNeighbors = cellList.template getNNIteratorSym<NO_CHECK>(cellList.getCell(
                    particleData.getOpenFPMContainer().getPos(p)), p.getKey(), particleData.getOpenFPMContainer().getPosVector());
//            auto iteratorNeighbors = cellList.template getNNIterator<NO_CHECK>(cellList.getCell(
//                    particleData.getOpenFPMContainer().getPos(p)));

            while (iteratorNeighbors.isNext()) {
                auto n = iteratorNeighbors.get();
                Particle<ParticleSignatureType> neighbor(particleData.dataContainer, n);

                if (particle != neighbor) {
                    particleMethod.interact(particle, neighbor);
                }
                ++iteratorNeighbors;
            }
            ++iteratorAll;
        }
    }

};

*/


template <typename ParticleMethodType, typename SimulationParametersType>
class Interaction_Impl<NEIGHBORHHOD_CELLLIST, ParticleMethodType, SimulationParametersType> {

    using ParticleSignatureType = typename ParticleMethodType::ParticleSignature;
    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

    ParticleMethodType particleMethod;
    SimulationParametersType simulationParameters;

    CellList_Impl<ParticleMethodType, SimulationParametersType> cellListImpl;

public:

    explicit Interaction_Impl(ParticleData<ParticleMethodType, SimulationParametersType> &particleData): cellListImpl(particleData) {}

    void executeInteraction(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {
      cellListImpl.executeInteraction(particleData);
    }

};



template <typename ParticleMethodType, typename SimulationParametersType>
class Interaction_Impl<NEIGHBORHOOD_MESH, ParticleMethodType, SimulationParametersType> {

    using ParticleSignatureType = typename ParticleMethodType::ParticleSignature;
    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

    ParticleMethodType particleMethod;
    SimulationParametersType simulationParameters;

    std::vector<std::array<int, dimension>> stencil;

public:

    explicit Interaction_Impl(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {

        createStencil();

        if (simulationParameters.interactionType == INTERACTION_SYMMETRIC)
            makeSymmetricStencil();
    }


    void createStencil() {

        // calculate farthest neighbor distance (in mesh nodes)
        int neighborDistance = round(simulationParameters.cutoff_radius / simulationParameters.meshSpacing);
        std::cout << "neighborDistance " << neighborDistance << std::endl;

        // width of nD cube
        int diameter = (2 * neighborDistance) + 1;

        // offset to shift nD cube to the center node
        int offset = floor(diameter / 2);

        // create nD cube
        for (int i = 0; i < pow(diameter, dimension); ++i) {
            std::array<int, dimension> neighbor;

            // compute stencil node
            for (int component = 0; component < dimension; ++component) {
                neighbor[component] = ((int)round(i / pow(diameter, component)) % diameter) - offset;
            }

            // add to stencil
            stencil.push_back(neighbor);
        }

        // remove all neighbors outside of cutoff radius
        PositionType cutoffRadius2 = simulationParameters.cutoff_radius * simulationParameters.cutoff_radius;
        PositionType spacing = simulationParameters.meshSpacing;
        auto it_rcut = std::remove_if(stencil.begin(), stencil.end(), [&cutoffRadius2, &spacing](std::array<int, dimension> node)
        {
            // calculate squared distance
            // from center to neighbor node
            float distance2 = 0;
            for (int i = 0; i < dimension; ++i) {
                distance2 += node[i] * node[i] * spacing * spacing;
            }
            return (distance2 > cutoffRadius2);
        });
        stencil.erase(it_rcut, stencil.end());

        // remove center particle
        auto it_center = std::remove_if(stencil.begin(), stencil.end(), [](std::array<int, dimension> node)
            { return std::all_of(node.begin(), node.end(), [](int comp){return comp == 0;}); });
        stencil.erase(it_center, stencil.end());

    }


    void makeSymmetricStencil() {

        // remove half of the particles to create symmetric stencil
        auto it_symm = std::remove_if(stencil.begin(), stencil.end(), [](std::array<int, dimension> node)
        {
            for (int i = 0; i < dimension; ++i) {
                // remove node if component > 0
                if (node[i] > 0)
                    return true;
                // keep node if component < 0
                if (node[i] < 0)
                    return false;
                // if component == 0, check next dimension
            }

            // remove center node [0, 0, 0]
            return true;
        });
        stencil.erase(it_symm, stencil.end());

    }


    void executeInteraction(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {


/*
        // flag all domain particles

        auto & vcl = create_vcluster();
        size_t pid = vcl.getProcessUnitID();

        auto iteratorFlag = particleData.getOpenFPMContainer().getDomainIterator();
        while (iteratorFlag.isNext()) {

            auto p = iteratorFlag.get();
            particleData.getOpenFPMContainer().template getProp<0>(p) = pid;
            ++iteratorFlag;
        }
*/

        // iterate through all particles
        auto iteratorAll = particleData.getParticleIterator();
        while (iteratorAll.isNext())
        {
            auto p = iteratorAll.get();
            Particle<ParticleSignatureType> particle(particleData.dataContainer, p);

//             auto & vcl = create_vcluster();
//            std::cout << vcl.getProcessUnitID() << " ~ interact " << p.to_string() << std::endl;
//            std::cout << "find neighbors for " << p.to_string() << std::endl;

            // iterate through all neighbor mesh nodes
            for (std::array<int, dimension> node : stencil) {
                auto n = p;

                // move to stencil position
                for (int component = 0; component < dimension; ++component) {
                    n = n.move(component, node[component]);
                }



                if (simulationParameters.interactionType == INTERACTION_SYMMETRIC) {

                    // check if neighbor is in on this processor domain
                    int check_count = 0;
                    auto &gdb_ext = particleData.dataContainer.getContainer().getLocalGridsInfo();
                    for (size_t s = 0; s < gdb_ext.size(); s++) {
                        Box<2, long int> bx = gdb_ext.get(s).Dbox;
                        bx += gdb_ext.get(s).origin;
                        if (bx.isInside(n.getKey().toPoint())) {
                            check_count++;
                            break;
                        }
                    }

                    // skip if not in processor domain
                    if (check_count == 0) {
                        continue;
                    }
                }


                // execute
                Particle<ParticleSignatureType> neighbor(particleData.dataContainer, n);
                particleMethod.interact(particle, neighbor);
            }

//            std::cout << "  ### interaction executed " << std::endl;


            ++iteratorAll;
        }
    }
};


#endif //OPENFPM_PDATA_INTERACTION_HPP
