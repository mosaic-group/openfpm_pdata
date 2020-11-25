#ifndef SRC_DECOMPOSITION_ORB_DECOMPOSITION_STRATEGY_HPP
#define SRC_DECOMPOSITION_ORB_DECOMPOSITION_STRATEGY_HPP

#include "Decomposition/dec_optimizer.hpp"
#include "DLB/DLB.hpp"
#include "Decomposition/Decomposition.hpp"
#include "Decomposition/Domain_icells_cart.hpp"
#include "Decomposition/common.hpp"
#include "Decomposition/ie_ghost.hpp"
#include "Decomposition/ie_loc_ghost.hpp"
#include "Decomposition/nn_processor.hpp"
#include "Graph/CartesianGraphFactory.hpp"
#include "GraphMLWriter/GraphMLWriter.hpp"
#include "NN/CellList/CellDecomposer.hpp"
#include "NN/CellList/CellList.hpp"
#include "Space/Ghost.hpp"
#include "Space/Shape/Box.hpp"
#include "Space/Shape/Point.hpp"
#include "SubdomainGraphNodes.hpp"
#include "VCluster/VCluster.hpp"
#include "Vector/map_vector.hpp"
#include "config.h"
#include "data_type/aggregate.hpp"
#include "util/generic.hpp"
#include "util/mathutil.hpp"
#include "util/se_util.hpp"

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

  /*! \brief Copy constructor
   *
   * \param cart object to copy
   *
   */
  OrbDecompositionStrategy(
      const OrbDecompositionStrategy<dim, domain_type, AbstractDecStrategy>& cart) {
    this->operator=(cart);
  }

  /*! \brief Copy the element
   *
   * \param cart element to copy
   *
   * \return itself
   *
   */
  OrbDecompositionStrategy<dim, domain_type> &
  operator=(const AbstractDecStrategy &cart) {
    // todo
  }

  ~OrbDecompositionStrategy() {
    if (orb) {
      delete orb;
    }
  }


  void setParameters(const size_t (&div_)[dim],
    ::Box<dim, domain_type> &domain_,
    const size_t (&bc)[dim],
    const Ghost<dim, domain_type> &ghost,
    const grid_sm<dim, void> &sec_dist = grid_sm<dim, void>()) {
    inner().setDomain(domain_);
  }

  /*! \brief Stub method to homogenize the interface do not use
   *
   */
  void createCartGraph() {
  }

  void decompose(openfpm::vector<Point<dim, domain_type>> points) {
    orb = new Orb(inner().getDomain(), inner().getVcluster().getProcessingUnits(), points);  // this takes care of the decomposition
  }

  Graph_CSR<nm_v<dim>, nm_e> &getGraph() {
    // todo get leaves from auto tree = orb->grp;  // Graph_CSR<ORB_node<T>, no_edge>
    // todo test only: ~trivial decomposition

    // build it
    if (inner().getGraph().getNVertex() < inner().v_cl.getProcessingUnits()) {
      // create only necessary vertices in the graph
      for (auto i = inner().getGraph().getNVertex(); i < inner().v_cl.getProcessingUnits(); ++i) {
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
