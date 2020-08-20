#ifndef SRC_DECOMPOSITION_CART_DISTRIBUTION_STRATEGY_HPP
#define SRC_DECOMPOSITION_CART_DISTRIBUTION_STRATEGY_HPP

#include "Decomposition/AbstractDistributionStrategy.hpp"

template <unsigned int dim, typename SpaceType,
          typename AbstractDistributionStrategy =
              AbstractDistributionStrategy<dim, SpaceType>>
class CartDistributionStrategy {

public:
  CartDistributionStrategy(Vcluster<> &v_cl) : dec(v_cl) {}

  ~CartDistributionStrategy() {}

private:
  AbstractDistributionStrategy dec;
};
#endif // SRC_DECOMPOSITION_CART_DISTRIBUTION_STRATEGY_HPP
