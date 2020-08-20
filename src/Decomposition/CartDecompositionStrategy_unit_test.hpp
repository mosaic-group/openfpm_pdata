#ifndef SRC_DECOMPOSITION_MYDECOMPOSITIONSTRATEGY_UNIT_TEST_HPP
#define SRC_DECOMPOSITION_MYDECOMPOSITIONSTRATEGY_UNIT_TEST_HPP

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
using SpaceType = double;
typedef vector_dist<
    SPACE_N_DIM, SpaceType,
    aggregate<size_t, double, double, double, double, double[SPACE_N_DIM],
              double[SPACE_N_DIM], double[SPACE_N_DIM]>>
    particles;

using MyDecompositionStrategy =
    AbstractDecompositionStrategy<SPACE_N_DIM, SpaceType>;

using MyDistributionStrategy =
    AbstractDistributionStrategy<SPACE_N_DIM, SpaceType>;

using ParmetisGraph = Parmetis<Graph_CSR<nm_v<SPACE_N_DIM>, nm_e>>;

struct MyComputationalCostsModel : ModelComputationalCosts {
  template <typename DistributionStrategy, typename vector>
  void addToComputation(DistributionStrategy &dist, vector &vd, size_t v,
                        size_t p) {
    if (vd.template getProp<type>(p) == FLUID) {
      dist.addComputationCost(v, 4);
    } else {
      dist.addComputationCost(v, 3);
    }
  }

  // todo where to use it
  template <typename DistributionStrategy>
  void init(DistributionStrategy &dist) {
    for (size_t i = 0; i < dist.getNOwnerSubSubDomains(); i++) {
      dist.setComputationCost(dist.getOwnerSubSubDomain(i), 1);
    }
  }

  template <typename particles, typename DecompositionStrategy,
            typename DistributionStrategy>
  void calculate(particles &vd, DecompositionStrategy &dec,
                 DistributionStrategy &dist) {
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
  void computeCommunicationAndMigrationCosts(DecompositionStrategy &dec,
                                             DistributionStrategy &dist,
                                             const size_t ts = 1) {
    SpaceType migration;
    size_t norm;
    std::tie(migration, norm) = dec.computeCommunicationCosts(dist.getGhost());
    dist.setMigrationCosts(migration, norm, ts);
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
  Box<SPACE_N_DIM, SpaceType> box({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
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
  Ghost<SPACE_N_DIM, SpaceType> g(0.01);

  // Boundary conditions
  size_t bc[] = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};

  // init
  dec.setParameters(div, box, bc, gsub);
  dist.setParameters(dec.getGrid(), g, gsub);
  dist.createCartGraph(bc, box);

  //////////////////////////////////////////////////////////////////// decompose
  dec.reset();
  dist.reset(parmetis_graph);
  if (dec.shouldSetCosts()) {
    mcc.computeCommunicationAndMigrationCosts(dec, dist);
  }
  dec.decompose(mde, parmetis_graph, dist.getVtxdist());

  /////////////////////////////////////////////////////////////////// distribute
  dist.distribute(parmetis_graph);

  ///////////////////////////////////////////////////////////////////////  merge
  dec.merge(dist.getGraph(), dist.getGhost(), dist.getGrid());

  ///////////////////////////////////////////////////////////////////// finalize
  dist.onEnd();
  dec.onEnd(dist.getGhost());

  // For each calculated ghost box
  printVar(dec.getNIGhostBox());
  for (size_t i = 0; i < dec.getNIGhostBox(); ++i) {
    SpaceBox<SPACE_N_DIM, SpaceType> b = dec.getIGhostBox(i);
    size_t proc = dec.getIGhostBoxProcessor(i);

    // sample one point inside the box
    Point<SPACE_N_DIM, SpaceType> p = b.rnd();

    // Check that ghost_processorsID return that processor number
    const openfpm::vector<size_t> &pr =
        dec.ghost_processorID<MyDecompositionStrategy::processor_id>(p, UNIQUE);

    printMe(vcl);
    std::cout << "ghost_processorsID " << proc << std::endl;
    for (auto x : pr) {
      std::cout << x << ", ";
    }
    std::cout << std::endl;

    bool found = isIn(pr, proc);
    /* question why is it useful? if (!found) {
      const openfpm::vector<size_t> pr2 =
          dec.ghost_processorID<MyDecompositionStrategy::processor_id>(p);
    } */

    printMe(vcl);
    std::cout << "assert " << found << " == true" << std::endl;
  }

  // Check the consistency
  bool val = dec.check_consistency();

  printMe(vcl);
  std::cout << "assert " << val << " == true" << std::endl;
}

#endif // SRC_DECOMPOSITION_MYDECOMPOSITIONSTRATEGY_UNIT_TEST_HPP
