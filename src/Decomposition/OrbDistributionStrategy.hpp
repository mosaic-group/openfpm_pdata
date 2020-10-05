#ifndef SRC_DECOMPOSITION_ORB_DISTRIBUTION_STRATEGY_HPP
#define SRC_DECOMPOSITION_ORB_DISTRIBUTION_STRATEGY_HPP

#include "Decomposition/AbstractDistributionStrategy.hpp"

template <unsigned int dim, typename domain_type, typename Memory = HeapMemory,
          template <typename> class layout_base = memory_traits_lin,
          typename AbstractDistStrategy =
              AbstractDistributionStrategy<dim, domain_type>>
class OrbDistributionStrategy {

  using Box = SpaceBox<dim, domain_type>;

public:
  OrbDistributionStrategy(Vcluster<> &v_cl) : dist(v_cl) {}

  ~OrbDistributionStrategy() {}

  void setParameters(const Ghost<dim, domain_type> &ghost) {
    dist.setGhost(ghost);
  }

  void distribute(DecompositionGraph& gp) {
    // todo: "TrivialDistribuzion" che assegna ad ogni nodo un processore in maniera requenziale

    dist.distribute();
  }

  void reset(DecompositionGraph& gp) {
    // todo
  }

  void onEnd() {
    // todo
  }

private:
  AbstractDistStrategy inner;
};
#endif // SRC_DECOMPOSITION_ORB_DISTRIBUTION_STRATEGY_HPP
