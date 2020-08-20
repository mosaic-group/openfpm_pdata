#ifndef SRC_DECOMPOSITION_CART_DECOMPOSITION_STRATEGY_HPP
#define SRC_DECOMPOSITION_CART_DECOMPOSITION_STRATEGY_HPP

#include "Decomposition/AbstractDecompositionStrategy.hpp"

template <unsigned int dim, typename SpaceType,
          typename AbstractDecompositionStrategy =
              AbstractDecompositionStrategy<dim, SpaceType>>
class CartDecompositionStrategy {

public:
  CartDecompositionStrategy(Vcluster<> &v_cl) : dec(v_cl) {}

  ~CartDecompositionStrategy() {}

  /*! \brief Copy constructor
   *
   * \param cart object to copy
   *
   */
  CartDecompositionStrategy(
      const CartDecompositionStrategy<dim, SpaceType, AbstractDecompositionStrategy>& cart) {
    this->operator=(cart);
  }

  /*! \brief Copy the element
   *
   * \param cart element to copy
   *
   * \return itself
   *
   */
  CartDecompositionStrategy<dim, SpaceType> &
  operator=(const AbstractDecompositionStrategy &cart) {
    // todo
  }

  template <typename Model, typename Graph>
  void decompose(Model m, Graph &graph, openfpm::vector<rid> &vtxdist) {
    // todo decompose
  }

  template <typename Graph>
  void merge(Graph &graph, Ghost<dim, SpaceType> &ghost) {
    // todo decompose
  }

private:
  AbstractDecompositionStrategy dec;
};
#endif // SRC_DECOMPOSITION_CART_DECOMPOSITION_STRATEGY_HPP
