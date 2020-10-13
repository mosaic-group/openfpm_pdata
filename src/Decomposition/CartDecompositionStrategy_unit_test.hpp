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

void CartDecomposition_non_periodic_test(const unsigned int nProcs) {
  Vcluster<> &vcl = create_vcluster();
  CartDecompositionStrategy<SPACE_N_DIM, domain_type> dec(vcl);
  CartDistributionStrategy<SPACE_N_DIM, domain_type> dist(vcl);

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

  ///////////////////////////////////////////////////////////////////////// init
  dec.setParameters(div, box, bc, g, gsub);
  dist.setParameters(dec.inner().getGraph(), box, dec.getGrid(), bc, gsub);

  //////////////////////////////////////////////////////////////////// decompose
  dec.decompose();

  /////////////////////////////////////////////////////////////////// distribute
  dist.distribute(dec.inner().getGraph());

  //////////////////////////////////////////////////////////////////////// merge
  dec.merge(dist.getGrid());

  ///////////////////////////////////////////////////////////////////// finalize
  dist.onEnd();
  dec.onEnd();

  // For each calculated ghost box
  for (size_t i = 0; i < dec.getNIGhostBox(); ++i) {
    SpaceBox<SPACE_N_DIM, domain_type> b = dec.getIGhostBox(i);
    size_t proc = dec.getIGhostBoxProcessor(i);

    // sample one point inside the box
    Point<SPACE_N_DIM, domain_type> p = b.rnd();

    // Check that ghost_processorsID return that processor number
    const openfpm::vector<size_t> &pr =
        dec.ghost_processorID<CartDecompositionStrategy<SPACE_N_DIM, domain_type>::processor_id>(p, UNIQUE);
    bool found = isIn(pr, proc);

    printMe(vcl);
    std::cout << "assert " << found << " == true" << std::endl;
  }

  // Check the consistency
  bool val = dec.check_consistency();

  printMe(vcl);
  std::cout << "assert " << val << " == true" << std::endl;
}

#endif // SRC_DECOMPOSITION_CARTDECOMPOSITIONSTRATEGY_UNIT_TEST_HPP
