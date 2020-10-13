#ifndef SRC_DECOMPOSITION_ORB_DECOMPOSITION_STRATEGY_HPP
#define SRC_DECOMPOSITION_ORB_DECOMPOSITION_STRATEGY_HPP

#include "Decomposition/AbstractDecompositionStrategy.hpp"
#include "Decomposition/ORB.hpp"

template <unsigned int dim, typename domain_type, typename Memory = HeapMemory, template <typename> class layout_base = memory_traits_lin,
          typename AbstractDecStrategy =
              AbstractDecompositionStrategy<dim, domain_type>>
class OrbDecompositionStrategy {

  using Box = SpaceBox<dim, domain_type>;
  using Orb = ORB<dim, domain_type>;

public:
  OrbDecompositionStrategy(Vcluster<> &v_cl) : _inner(v_cl) {}

  ~OrbDecompositionStrategy() {}

  /*! \brief Copy constructor
   *
   * \param cart object to copy
   *
   */
  OrbDecompositionStrategy(
      const OrbDecompositionStrategy<dim, SpaceType, AbstractDecompositionStrategy>& cart) {
    this->operator=(cart);
  }

  /*! \brief Copy the element
   *
   * \param cart element to copy
   *
   * \return itself
   *
   */
  OrbDecompositionStrategy<dim, SpaceType> &
  operator=(const AbstractDecompositionStrategy &cart) {
    // todo
  }

  ~OrbDecompositionStrategy() {
    if (orb) {
      delete orb;
    }
  }

  void setParameters(::Box<dim, domain_type> &domain_, const size_t (&bc)[dim]) {
    inner.setParameters(domain_, bc);
  }

  template <typename Point>
  void decompose(openfpm::vector<Point> &points) {
    orb = Orb(inner().getDomain(), inner().v_cl.getProcessingUnits(), points);  // this takes care of the decomposition
  }

  Graph_CSR<nm_v<dim>, nm_e> &getGraph() {
    // todo get leaves from auto tree = orb->grp;  // Graph_CSR<ORB_node<T>, no_edge>
    
    // todo test only: ~trivial decomposition -> build a graph with n (= proc units) vertices

    // build it
    if (inner().getGraph().getNVertex() != inner().v_cl.getProcessingUnits()) {
      // create only necessary vertices in the graph
      for (auto i = inner().v_cl.getProcessingUnits(); i < inner().v_cl.getProcessingUnits(); ++i) {
        inner().getGraph().addVertex();
      }
    }
    
    return inner().getGraph();
  }

  AbstractDecStrategy &inner() {
    return _inner;
  }

private:
  AbstractDecStrategy _inner;
  Orb* orb;
};
#endif // SRC_DECOMPOSITION_ORB_DECOMPOSITION_STRATEGY_HPP
