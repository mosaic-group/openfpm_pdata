//
// Created by landfried on 01.02.22.
//

#ifndef OPENFPM_PDATA_CONSTANTS_HPP
#define OPENFPM_PDATA_CONSTANTS_HPP


// SimulationParameters

// Type of interaction
constexpr int INTERACTION_PULL = 0;
constexpr int INTERACTION_PUSH = 1;
constexpr int INTERACTION_PULL_PUSH = 2;
constexpr int INTERACTION_SYMMETRIC = 3;

// Neighborhood determination
struct NEIGHBORHOOD_ALLPARTICLES {};
struct NEIGHBORHHOD_CELLLIST {};
struct NEIGHBORHOOD_MESH {};

// Initial placement of free particles
struct INITIALCONDITION_NONE {};
struct INITIALCONDITION_MESH {};
struct INITIALCONDITION_RANDOM {};

#endif //OPENFPM_PDATA_CONSTANTS_HPP
