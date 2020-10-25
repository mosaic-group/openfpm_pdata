#ifndef SRC_DECOMPOSITION_ORB_DECOMPOSITIONSTRATEGY_UNIT_TEST_HPP
#define SRC_DECOMPOSITION_ORB_DECOMPOSITIONSTRATEGY_UNIT_TEST_HPP

#include "Decomposition/AbstractStrategyModels.hpp"
#include "Decomposition/ORBDecompositionStrategy.hpp"
#include "Decomposition/Distribution/SequentialDistributionStrategy.hpp"
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
using domain_type = float;

// n of points sampled
#define N_POINTS 1024

typedef vector_dist<
    SPACE_N_DIM, domain_type,
    aggregate<size_t, domain_type, domain_type, domain_type, domain_type, domain_type[SPACE_N_DIM],
              domain_type[SPACE_N_DIM], domain_type[SPACE_N_DIM]>> particles;

void ORBDecomposition_non_periodic_test(const unsigned int nProcs) {
  Vcluster<> &vcl = create_vcluster();
  OrbDecompositionStrategy<SPACE_N_DIM, domain_type> dec(vcl);
  SequentialDistributionStrategy<SPACE_N_DIM, domain_type> dist(vcl);

  // Physical domain
  Box<SPACE_N_DIM, domain_type> box({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});

  std::srand(create_vcluster().getProcessUnitID());
  std::default_random_engine eg;
  std::uniform_real_distribution<domain_type> ud(0.0, 1.0);

  typedef Point<SPACE_N_DIM, domain_type> p;

  // create a random local vector of particles
  openfpm::vector<Point<SPACE_N_DIM, domain_type>> vp(N_POINTS);

  // fill the particles
  auto vp_it = vp.getIterator();
  while (vp_it.isNext()) {
    auto key = vp_it.get();

    vp.get<p::x>(key)[0] = ud(eg);
    vp.get<p::x>(key)[1] = ud(eg);
    vp.get<p::x>(key)[2] = ud(eg);

    ++vp_it;
  }

  ///////////////////////////////////////////////////////////////////////// init
  dec.setParameters(box);

  //////////////////////////////////////////////////////////////////// decompose
  dec.decompose(vp);

  /////////////////////////////////////////////////////////////////// distribute
  dist.distribute(dec.inner().getGraph());

  //////////////////////////////////////////////////////////////////////// merge
  // no need (we're NOT using sub(sub)domains)

  ASSERT: i-th processing unit has i-th vertex
  auto myRank = vcl.rank();
  auto myLeaf = dec.inner().getGraph().template vertex_p<nm_v_proc_id>(myRank);
  printMe(vcl);
  std::cout << "assert " << myRank << " == " << myLeaf << std::endl;

  // todo more asserts
}

#endif // SRC_DECOMPOSITION_ORB_DECOMPOSITIONSTRATEGY_UNIT_TEST_HPP
