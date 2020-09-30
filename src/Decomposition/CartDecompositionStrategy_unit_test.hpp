#ifndef SRC_DECOMPOSITION_CARTDECOMPOSITIONSTRATEGY_UNIT_TEST_HPP
#define SRC_DECOMPOSITION_CARTDECOMPOSITIONSTRATEGY_UNIT_TEST_HPP

#include "Decomposition/AbstractStrategyModels.hpp"
#include "Decomposition/CartDecompositionStrategy.hpp"
#include "Decomposition/CartDistributionStrategy.hpp"
#include "util/generic.hpp"

#define SUB_UNIT_FACTOR 1024

// Properties

// A constant to indicate boundary particles
#define BOUNDARY 0

// A constant to indicate fluid particles
#define FLUID 1

// FLUID or BOUNDARY
const size_t type = 0;

// Type of the vector containing particles
constexpr unsigned int SPACE_N_DIM = 3;
using domain_type = double;

using AbstractDecStrategy =
    AbstractDecompositionStrategy<SPACE_N_DIM, domain_type>;

using MyDecompositionStrategy =
    CartDecompositionStrategy<SPACE_N_DIM, domain_type>;

using MyDistributionStrategy =
    CartDistributionStrategy<SPACE_N_DIM, domain_type>;

using ParmetisGraph = Parmetis<Graph_CSR<nm_v<SPACE_N_DIM>, nm_e>>;

struct MyComputationalCostsModel : ModelComputationalCosts {
  template <typename DecompositionStrategy, typename DistributionStrategy>
  void computeCommunicationAndMigrationCosts(DecompositionStrategy &dec,
                                             DistributionStrategy &dist,
                                             const size_t ts = 1) {
    domain_type migration;
    size_t norm;
    std::tie(migration, norm) = dec.computeCommunicationCosts(dist.dist.getGhost());
    dist.dist.setMigrationCosts(migration, norm, ts);
  }
};

struct MyDecompositionModel : ModelDecompose {};

struct MyDistributionModel : ModelDistribute {
  val_t toll() { return 1.01; }
  
  template <typename DistributionStrategy>
  void applyModel(DistributionStrategy &dist, size_t v) {
    const size_t id = v;
    const size_t weight = dist.getSubSubDomainComputationCost(v) *
                          dist.getSubSubDomainComputationCost(v);
    dist.setComputationCost(id, weight);
  }

  template <typename DistributionStrategy, typename Graph>
  void finalize(DistributionStrategy &dist, Graph &graph) {
    for (auto i = 0; i < dist.getNOwnerSubSubDomains(); i++) {
      // apply model to all the sub-sub-domains
      applyModel(dist, dist.getOwnerSubSubDomain(i));
    }

    dist.setDistTol(graph, toll());
  }
};

void CartDecomposition_non_periodic_test(const unsigned int nProcs) {
  Vcluster<> &vcl = create_vcluster();

  // specify
  // - how we want to add the computational cost ...
  MyComputationalCostsModel mcc;

  // - how we want to decompose ...
  MyDecompositionStrategy dec(vcl);
  MyDecompositionModel mde;

  // - how we want to distribute ...
  MyDistributionStrategy dist(vcl);
  MyDistributionModel mdi;

  // ... and our shared information

  //! Convert the graph to parmetis format
  ParmetisGraph parmetis_graph(vcl, vcl.getProcessingUnits());

  // Physical domain
  Box<SPACE_N_DIM, domain_type> box({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
  size_t div[SPACE_N_DIM];
  size_t div_sub[SPACE_N_DIM];

  // Get the number of processor and calculate the number of sub-domain
  // for each processor (SUB_UNIT_FACTOR=64)
  size_t n_proc = vcl.getProcessingUnits();
  size_t n_sub = n_proc * SUB_UNIT_FACTOR * 4 * 4 * 4;

  // Set the number of sub-domains on each dimension (in a scalable way)
  for (int i = 0; i < SPACE_N_DIM; i++) {
    div[i] = openfpm::math::round_big_2(pow(n_sub, 1.0 / 3));
  }

  // create a sub_distribution grid
  for (int i = 0; i < SPACE_N_DIM; i++) {
    div_sub[i] = div[i] / 4;
  }

  grid_sm<SPACE_N_DIM, void> gsub(div_sub);

  // Define ghost
  Ghost<SPACE_N_DIM, domain_type> g(0.01);

  // Boundary conditions
  size_t bc[] = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};

  // init
  dec.setParameters(div, box, bc, gsub);
  dist.setParameters(dec.getGrid(), g, gsub);
  dist.createCartGraph(bc, box);

  //////////////////////////////////////////////////////////////////// decompose
  dec.dec.reset();
  dist.dist.reset(parmetis_graph);
  if (dec.dec.shouldSetCosts()) {
    mcc.computeCommunicationAndMigrationCosts(dec, dist);
  }
  dec.decompose(mde, parmetis_graph, dist.dist.getVtxdist());

  /////////////////////////////////////////////////////////////////// distribute
  dist.dist.distribute(parmetis_graph);

  ///////////////////////////////////////////////////////////////////////  merge
  dec.merge(dist.dist.getGraph(), dist.dist.getGhost(), dist.getGrid());

  ///////////////////////////////////////////////////////////////////// finalize
  dist.dist.onEnd();
  dec.dec.onEnd(dist.dist.getGhost());

  // For each calculated ghost box
  for (size_t i = 0; i < dec.dec.getNIGhostBox(); ++i) {
    SpaceBox<SPACE_N_DIM, domain_type> b = dec.dec.getIGhostBox(i);
    size_t proc = dec.dec.getIGhostBoxProcessor(i);

    // sample one point inside the box
    Point<SPACE_N_DIM, domain_type> p = b.rnd();

    // Check that ghost_processorsID return that processor number
    const openfpm::vector<size_t> &pr =
        dec.dec.ghost_processorID<AbstractDecStrategy::processor_id>(p, UNIQUE);
    bool found = isIn(pr, proc);

    printMe(vcl);
    std::cout << "assert " << found << " == true" << std::endl;
  }

  // Check the consistency
  bool val = dec.dec.check_consistency();

  printMe(vcl);
  std::cout << "assert " << val << " == true" << std::endl;
}

#endif // SRC_DECOMPOSITION_CARTDECOMPOSITIONSTRATEGY_UNIT_TEST_HPP
