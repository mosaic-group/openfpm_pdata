#ifndef SRC_DECOMPOSITION_ABSTRACT_DECOMPOSITION_STRATEGY_HPP
#define SRC_DECOMPOSITION_ABSTRACT_DECOMPOSITION_STRATEGY_HPP

#include "Decomposition/AbstractDecompositionStrategy.hpp"

template <unsigned int dim,
          typename SpaceType,
          typename AbstractDecompositionStrategy = AbstractDecompositionStrategy<dim, SpaceType>>
class CartDecompositionStrategy {

public:
  CartDecompositionStrategy(Vcluster<>& v_cl)
    : dec(v_cl) {}

  ~CartDecompositionStrategy() {}
private:
  AbstractDecompositionStrategy dec;
};
#endif  // SRC_DECOMPOSITION_ABSTRACT_DECOMPOSITION_STRATEGY_HPP
