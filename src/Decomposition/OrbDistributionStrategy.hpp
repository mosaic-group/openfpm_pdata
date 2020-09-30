#ifndef SRC_DECOMPOSITION_ORB_DISTRIBUTION_STRATEGY_HPP
#define SRC_DECOMPOSITION_ORB_DISTRIBUTION_STRATEGY_HPP

#include "Decomposition/AbstractDistributionStrategy.hpp"

template <unsigned int dim, typename domain_type,
          typename AbstractDistStrategy =
              AbstractDistributionStrategy<dim, domain_type>, typename Memory = HeapMemory,
          template <typename> class layout_base = memory_traits_lin>
class OrbDistributionStrategy {

  using Box = SpaceBox<dim, domain_type>;
  using DGrid = grid_sm<dim, void>;

public:
  OrbDistributionStrategy(Vcluster<> &v_cl) : dist(v_cl) {}

  ~OrbDistributionStrategy() {}

  void setParameters(const Ghost<dim, domain_type> &ghost) {
    dist.ghost = ghost;
  }

// todo private:

  AbstractDistStrategy dist;
};
#endif // SRC_DECOMPOSITION_ORB_DISTRIBUTION_STRATEGY_HPP
