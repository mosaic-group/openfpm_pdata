#ifndef OPENFPM_PDATA_ABSTRACTDISTRIBUTIONSTRATEGY_HPP
#define OPENFPM_PDATA_ABSTRACTDISTRIBUTIONSTRATEGY_HPP

#include "Decomposition/Domain_NN_calculator_cart.hpp"
#include "Graph/CartesianGraphFactory.hpp"
#include "Graph/ids.hpp"
#include "SubdomainGraphNodes.hpp"

/*! \brief Class that distribute sub-sub-domains across processors
 */
template <unsigned int dim, typename T>
class AbstractDistributionStrategy : public domain_nn_calculator_cart<dim> {
  //! It simplify to access the SpaceBox element
  using Box = SpaceBox<dim, T>;

public:
  //! Vcluster
  Vcluster<>& v_cl;

  /*! Constructor
   *
   * \param v_cl Vcluster to use as communication object in this class
   */
  AbstractDistributionStrategy(Vcluster<>& v_cl) : v_cl(v_cl) {}

  /*! \brief Return the global id of the owned sub-sub-domain
   *
   * \param id in the list of owned sub-sub-domains
   *
   * \return the global id
   *
   */
  size_t getOwnerSubSubDomain(size_t id) const { return 0; }

  /*! \brief Return the total number of sub-sub-domains this processor own
   *
   * \return the total number of sub-sub-domains owned by this processor
   *
   */
  size_t getNOwnerSubSubDomains() const { return 0; }

  /*! \brief Returns total number of sub-sub-domains in the distribution graph
   *
   * \return the total number of sub-sub-domains
   *
   */
  size_t getNSubSubDomains() const { return 0; }

  /*! \brief Set migration cost of the vertex id
   *
   * \param id of the vertex to update
   * \param migration cost of the migration
   */
  void setMigrationCost(size_t id, size_t migration) {}

  /*! \brief Set communication cost of the edge id
   *
   * \param v_id Id of the source vertex of the edge
   * \param e i child of the vertex
   * \param communication Communication value
   */
  void setCommunicationCost(size_t v_id, size_t e, size_t communication) {}

  /*! \brief Returns total number of neighbors of the sub-sub-domain id
   *
   * \param id id of the sub-sub-domain
   *
   * \return the number of neighborhood sub-sub-domains for each sub-domain
   *
   */
  size_t getNSubSubDomainNeighbors(size_t id) { return 0; }

  /*! \brief Set the tolerance for each partition
   *
   * \param tol tolerance
   *
   */
  void setDistTol(double tol) {}

  void setMigrationCosts(const float migration,
                         const size_t norm,
                         const size_t ts) {
    for (auto i = 0; i < getNSubSubDomains(); i++) {
      setMigrationCost(i, norm * migration);

      for (auto s = 0; s < getNSubSubDomainNeighbors(i); s++) {
        // We have to remove getSubSubDomainComputationCost(i) otherwise the
        // graph is not directed
        setCommunicationCost(i, s, 1 * ts);
      }
    }
  }

  template <typename DecompositionStrategy, typename Model>
  void distribute(DecompositionStrategy& dec, Model m) {}

  void onEnd() {
    domain_nn_calculator_cart<dim>::reset();
    domain_nn_calculator_cart<dim>::setParameters(proc_box);
  }

private:
  //! Processor domain bounding box
  ::Box<dim, size_t> proc_box;
};

#endif  // OPENFPM_PDATA_ABSTRACTDISTRIBUTIONSTRATEGY_HPP
