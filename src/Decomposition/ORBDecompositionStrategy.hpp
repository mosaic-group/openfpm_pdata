#ifndef SRC_DECOMPOSITION_ORB_DECOMPOSITION_STRATEGY_HPP
#define SRC_DECOMPOSITION_ORB_DECOMPOSITION_STRATEGY_HPP

#include "Decomposition/AbstractDecompositionStrategy.hpp"
#include "Decomposition/ORB.hpp"

template <unsigned int dim, typename domain_type,
          typename AbstractDecStrategy =
              AbstractDecompositionStrategy<dim, domain_type>, typename Memory = HeapMemory, template <typename> class layout_base = memory_traits_lin>
class OrbDecompositionStrategy {

  using Box = SpaceBox<dim, domain_type>;
  using Orb = ORB<dim, domain_type>;

public:
  OrbDecompositionStrategy(Vcluster<> &v_cl) : dec(v_cl) {}

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
      delete orb;  // todo needed ??
    }
  }

  void setParameters(::Box<dim, domain_type> &domain, const size_t (&bc)[dim]) {
    std::copy(bc, bc + dim, dec.bc);  // boundary conditions

    dec.domain = domain;  // domain
  }

  template <typename Point>
  void decompose(openfpm::vector<Point> &points) {
    orb = Orb(dec.domain, dec.v_cl.getProcessingUnits(), points);
  }

  Orb::Tree& getTree() {
    return orb->grp;
  }

// todo private:
  AbstractDecStrategy dec;
  Orb* orb;
};
#endif // SRC_DECOMPOSITION_ORB_DECOMPOSITION_STRATEGY_HPP
