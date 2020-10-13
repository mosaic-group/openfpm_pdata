#ifndef SRC_DECOMPOSITION_SEQUENTIAL_DISTRIBUTION_STRATEGY_HPP
#define SRC_DECOMPOSITION_SEQUENTIAL_DISTRIBUTION_STRATEGY_HPP

#include "Decomposition/AbstractDistributionStrategy.hpp"

/**
 * assegna ad ogni nodo un processore in maniera sequenziale
 */
template <unsigned int dim, typename domain_type, typename Memory = HeapMemory,
          template <typename> class layout_base = memory_traits_lin, typename DecompositionGraph = Graph_CSR<nm_v<dim>, nm_e>,
          typename AbstractDistStrategy =
              AbstractDistributionStrategy<dim, domain_type, Parmetis<DecompositionGraph>>>
class SequentialDistributionStrategy {

  using Box = SpaceBox<dim, domain_type>;

public:
  SequentialDistributionStrategy(Vcluster<> &v_cl) : _inner(v_cl) {}

  ~SequentialDistributionStrategy() {}

  void distribute(DecompositionGraph& gp) {
    auto i = inner().getVcluster().rank();
    gp.vertex(i).template get<nm_v_proc_id>() = i;
    inner().distribute();
  }

  void reset(DecompositionGraph& gp) {
    // todo
  }

  AbstractDistStrategy &inner() {
    return _inner;
  }

private:
  AbstractDistStrategy _inner;
};
#endif // SRC_DECOMPOSITION_SEQUENTIAL_DISTRIBUTION_STRATEGY_HPP
