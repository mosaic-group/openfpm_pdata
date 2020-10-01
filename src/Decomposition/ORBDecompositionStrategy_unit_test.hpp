#ifndef SRC_DECOMPOSITION_ORBDECOMPOSITIONSTRATEGY_UNIT_TEST_HPP
#define SRC_DECOMPOSITION_ORBDECOMPOSITIONSTRATEGY_UNIT_TEST_HPP

#include "Decomposition/AbstractStrategyModels.hpp"
#include "Decomposition/OrbDecompositionStrategy.hpp"
#include "Decomposition/OrbDistributionStrategy.hpp"
#include "util/generic.hpp"

#define SUB_UNIT_FACTOR 1024

// Properties

// A constant to indicate boundary particles
#define BOUNDARY 0

// A constant to indicate fluid particles
#define FLUID 1

// FLUID or BOUNDARY
const size_t type = 0;

// n of points sampled
#define N_POINTS 1024

// Type of the vector containing particles
constexpr unsigned int SPACE_N_DIM = 3;
using domain_type = double;

using AbstractDecStrategy =
    AbstractDecompositionStrategy<SPACE_N_DIM, domain_type>;

using MyDecompositionStrategy =
    OrbDecompositionStrategy<SPACE_N_DIM, domain_type>;

using MyDistributionStrategy =
    OrbDistributionStrategy<SPACE_N_DIM, domain_type>;

using ParmetisGraph = Parmetis<Graph_CSR<nm_v<SPACE_N_DIM>, nm_e>>;

using MyPoint = Point<SPACE_N_DIM, domain_type>;
using MyPoints = openfpm::vector<MyPoint>;

typedef vector_dist<
    SPACE_N_DIM, SpaceType,
    aggregate<size_t, double, double, double, double, double[SPACE_N_DIM],
              double[SPACE_N_DIM], double[SPACE_N_DIM]>>
    particles;

struct MyComputationalCostsModel : ModelComputationalCosts {
  template <typename DistributionStrategy, typename vector>
  void addToComputation(DistributionStrategy &dist, size_t v) {
    dist.addComputationCost(v, 2);
  }

  template <typename DistributionStrategy>
  void init(DistributionStrategy &dist) {
    for (size_t i = 0; i < dist.getNOwnerSubSubDomains(); i++) {
      dist.setComputationCost(dist.getOwnerSubSubDomain(i), 1);
    }
  }

  template <typename particles, typename DecompositionStrategy,
            typename DistributionStrategy>
  void computeCommunicationAndMigrationCosts(MyPoints &vd, DecompositionStrategy &dec, DistributionStrategy &dist) {
    CellDecomposer_sm<SPACE_N_DIM, domain_type, shift<SPACE_N_DIM, domain_type>>
        cdsm;  // question can really use this ?
    cdsm.setDimensions(dec.getDomain(), dist.getGrid().getSize(), 0);
    for (auto it = vd.getIterator(); !it.hasEnded(); ++it) {
      auto p = it.get();
      const size_t v = cdsm.getCell(p);
      addToComputation(dist, v);
    }
  }
};

struct MyDecompositionModel : ModelDecompose {};

void OrbDecomposition_non_periodic_test(const unsigned int nProcs) {
  Vcluster<> &vcl = create_vcluster();

  // specify
  // - how we want to add the computational cost ...
  MyComputationalCostsModel mcc;

  // - how we want to decompose ...
  MyDecompositionStrategy dec(vcl);
  MyDecompositionModel mde;

  // - how we want to distribute ...
  MyDistributionStrategy dist(vcl);

  // ... and our shared information

  //! Convert the graph to parmetis format
  ParmetisGraph parmetis_graph(vcl, vcl.getProcessingUnits());

  // Physical domain
  Box<SPACE_N_DIM, domain_type> box({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
  
  // fill particles
  MyPoints vp(N_POINTS);

  // create the random generator engine
	std::srand(create_vcluster().getProcessUnitID());
  std::default_random_engine eg;
  std::uniform_real_distribution<float> ud(0.0, 1.0);

	auto vp_it = vp.getIterator();
	while (vp_it.isNext()) {
		auto key = vp_it.get();

		vp.get<p::x>(key)[0] = ud(eg);
		vp.get<p::x>(key)[1] = ud(eg);
		vp.get<p::x>(key)[2] = ud(eg);

		++vp_it;
	}

  // Define ghost
  Ghost<SPACE_N_DIM, domain_type> g(0.01);

  // Boundary conditions
  size_t bc[] = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};

  // init
  mcc.init(dist);
  dec.setParameters(box, bc);
  dist.setParameters(g);

  //////////////////////////////////////////////////////////////////// decompose
  dec.dec.reset();
  dist.dist.reset(parmetis_graph);
  if (dec.dec.shouldSetCosts()) {
    mcc.computeCommunicationAndMigrationCosts(vp, dec, dist);
  }
  dec.decompose(mde, parmetis_graph, dist.getVtxdist());

  /////////////////////////////////////////////////////////////////// distribute
  dist.dist.distribute(parmetis_graph);

  ///////////////////////////////////////////////////////////////////////  merge
  dec.merge(dist.dist.getGraph(), dist.dist.getGhost());

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

#endif // SRC_DECOMPOSITION_ORBDECOMPOSITIONSTRATEGY_UNIT_TEST_HPP
