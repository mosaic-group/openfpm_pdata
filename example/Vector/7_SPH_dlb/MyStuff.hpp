#ifndef OPENFPM_PDATA_MYSTUFF_HPP
#define OPENFPM_PDATA_MYSTUFF_HPP

#include "./AbstractStrategyModels.hpp"
#include "./AbstractDecompositionStrategy.hpp"
#include "./AbstractDistributionStrategy.hpp"

using val_t = double;
using Box3D = Box<3, val_t>;
constexpr unsigned int REBALANCING_TIME_STEPS_THRESHOLD = 200;

// Properties

// A constant to indicate boundary particles
#define BOUNDARY 0

// A constant to indicate fluid particles
#define FLUID 1

// FLUID or BOUNDARY
const size_t type = 0;

// Type of the vector containing particles
constexpr unsigned int SPACE_N_DIM = 3;
using SpaceType = double;
typedef vector_dist<SPACE_N_DIM,
                    SpaceType,
                    aggregate<size_t,
                              double,
                              double,
                              double,
                              double,
                              double[3],
                              double[3],
                              double[3]>>
    particles;

struct MyComputationalCostsModel : ModelComputationalCosts {
  template <typename DistributionStrategy, typename vector>
  void addToComputation(DistributionStrategy& dist,
                        vector& vd,
                        size_t v,
                        size_t p) {
    if (vd.template getProp<type>(p) == FLUID) {
      dist.addComputationCost(v, 4);
    } else {
      dist.addComputationCost(v, 3);
    }
  }

  // todo where to use it
  template <typename DistributionStrategy>
  void init(DistributionStrategy& dist) {
    for (size_t i = 0; i < dist.getNOwnerSubSubDomains(); i++) {
      dist.setComputationCost(dist.getOwnerSubSubDomain(i), 1);
    }
  }

  template <typename particles,
            typename DecompositionStrategy,
            typename DistributionStrategy>
  void calculate(particles& vd,
                 DecompositionStrategy& dec,
                 DistributionStrategy& dist) {
    CellDecomposer_sm<SPACE_N_DIM, SpaceType, shift<SPACE_N_DIM, SpaceType>>
        cdsm;
    cdsm.setDimensions(dec.getDomain(), dist.getGrid().getSize(), 0);
    for (auto it = vd.getDomainIterator(); !it.hasEnded(); ++it) {
      Point<SPACE_N_DIM, SpaceType> p = vd.getPos(it.get());
      const size_t v = cdsm.getCell(p);
      addToComputation(dist, vd, v, it.get().getKey());
    }
  }

  template <typename DecompositionStrategy, typename DistributionStrategy>
  void computeCommunicationAndMigrationCosts(DecompositionStrategy& dec,
                                             DistributionStrategy& dist,
                                             const size_t ts = 1) {
    float migration;
    size_t norm;
    std::tie(migration, norm) = dec.computeCommunicationCosts(dist.getGhost());
    dist.setMigrationCosts(migration, norm, ts);
  }
};

struct MyDecompositionModel : ModelDecompose {};

struct MyDistributionModel : ModelDistribute {
  val_t toll() { return 1.01; }
  template <typename DistributionStrategy>
  void applyModel(DistributionStrategy& dist, size_t v) {
    const size_t id = v;
    const size_t weight = dist.getSubSubDomainComputationCost(v) *
                          dist.getSubSubDomainComputationCost(v);
    dist.setComputationCost(id, weight);
  }

  template <typename DistributionStrategy, typename Graph>
  void finalize(DistributionStrategy& dist, Graph& graph) {
    for (auto i = 0; i < dist.getNOwnerSubSubDomains(); i++) {
      // apply model to all the sub-sub-domains
      applyModel(dist, dist.getOwnerSubSubDomain(i));
    }

    dist.setDistTol(graph, toll());
  }
};

using MyDecompositionStrategy =
    AbstractDecompositionStrategy<SPACE_N_DIM, SpaceType>;

using MyDistributionStrategy =
    AbstractDistributionStrategy<SPACE_N_DIM, SpaceType>;

using ParmetisGraph = Parmetis<Graph_CSR<nm_v<SPACE_N_DIM>, nm_e>>;

#endif  // OPENFPM_PDATA_MYSTUFF_HPP
