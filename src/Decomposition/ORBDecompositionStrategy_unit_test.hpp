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
using MyPoint = Point<SPACE_N_DIM, domain_type>;
using MyPoints = openfpm::vector<MyPoint>;

typedef vector_dist<
    SPACE_N_DIM, SpaceType,
    aggregate<size_t, double, double, double, double, double[SPACE_N_DIM],
              double[SPACE_N_DIM], double[SPACE_N_DIM]>>
    particles;

void OrbDecomposition_non_periodic_test(const unsigned int nProcs) {
  Vcluster<> &vcl = create_vcluster();
  ORBDecompositionStrategy<SPACE_N_DIM, domain_type> dec(vcl);
  ORBDistributionStrategy<SPACE_N_DIM, domain_type> dist(vcl);

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

  // Boundary conditions
  size_t bc[] = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};

  ///////////////////////////////////////////////////////////////////////// init
  dec.setParameters(box, bc);

  /* todo
  1 creo un vector di particles
  2 chiamo bisect (fino a che # foglie >= # processori)
  3 trasformo il Tree (alias di Graph_CSR_s) dell'ORB a Parmetis e lo passo a ParmetisDistribution
  4 distribute()
  5 NON c'è proprio il passo del merge (in quanto non ci sono sub-sub-domains) ...
  6 faccio qualche assert
  */

  //////////////////////////////////////////////////////////////////// decompose
  dec.decompose(vp);

  /////////////////////////////////////////////////////////////////// distribute
  // trasformo il Tree (alias di Graph_CSR_s) dell'ORB a Parmetis e lo passo a ParmetisDistribution
  dist.distribute(dec.getTree());

  ///////////////////////////////////////////////////////////////////// finalize
  dist.onEnd();
  dec.onEnd(d);

  // todo faccio qualche assert
}

#endif // SRC_DECOMPOSITION_ORBDECOMPOSITIONSTRATEGY_UNIT_TEST_HPP
