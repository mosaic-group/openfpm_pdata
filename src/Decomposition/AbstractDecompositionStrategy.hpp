#ifndef SRC_DECOMPOSITION_ABSTRACT_DECOMPOSITION_STRATEGY_HPP
#define SRC_DECOMPOSITION_ABSTRACT_DECOMPOSITION_STRATEGY_HPP

#include <cmath>
#include <initializer_list>
#include <unordered_map>
#include <utility>
#include <vector>

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

template <unsigned int dim, typename domain_type, typename Memory = HeapMemory,
          template <typename> class layout_base = memory_traits_lin,
          typename DGrid = grid_sm<dim, void>>
class AbstractDecompositionStrategy {
  using Box = SpaceBox<dim, domain_type>;
  using SubDomains =
      openfpm::vector<Box, Memory, typename layout_base<Box>::type,
                      layout_base>;

public:
  /*! \brief Abstract decomposition constructor
   *
   * \param v_cl Virtual cluster, used internally to handle or pipeline
   * comm0unication
   *
   */
  AbstractDecompositionStrategy(Vcluster<> &v_cl) : ref_cnt(0) {}

  /*! \brief Copy constructor
   *
   * \param cart object to copy
   *
   */
  AbstractDecompositionStrategy(
      const AbstractDecompositionStrategy<dim, domain_type, Memory, layout_base, DGrid>
          &cart) : ref_cnt(0) {
    this->operator=(cart);
  }

  /*! \brief Copy the element
   *
   * \param cart element to copy
   *
   * \return itself
   *
   */
  AbstractDecompositionStrategy<dim, domain_type, Memory, layout_base, DGrid> &
  operator=(const AbstractDecompositionStrategy &cart) {
    // todo
  }

  //! Destructor
  ~AbstractDecompositionStrategy() {}

  //! Increment the reference counter
  void incRef() { ref_cnt++; }

  //! Decrement the reference counter
  void decRef() { ref_cnt--; }

  //! Return the reference counter
  long int ref() { return ref_cnt; }

  bool shouldSetCosts() { return !costBeenSet; }

  /*! \brief Return the box of the physical domain
   *
   * \return The physical domain box
   *
   */
  const ::Box<dim, domain_type> &getDomain() const { return domain; }

  /*! \brief Delete the decomposition and reset the data-structure */
  void reset() {}

  void decompose() {}

  void onEnd() {}

// todo private:

  //! (rectangular) domain to decompose
  ::Box<dim, domain_type> domain;

  //! reference counter of the object in case is shared between object
  long int ref_cnt;

  bool costBeenSet = false;

  //! Runtime virtual cluster machine
  Vcluster<> &v_cl;
};
#endif // SRC_DECOMPOSITION_ABSTRACT_DECOMPOSITION_STRATEGY_HPP