#ifndef OPENFPM_PDATA_ABSTRACTDISTRIBUTIONSTRATEGY_HPP
#define OPENFPM_PDATA_ABSTRACTDISTRIBUTIONSTRATEGY_HPP

#include "SubdomainGraphNodes.hpp"
#include "Graph/ids.hpp"
#include "Graph/CartesianGraphFactory.hpp"

/*! \brief Class that distribute sub-sub-domains across processors
 */
template<unsigned int dim, typename T>
class AbstractDistributionStrategy {
public:

    //! Vcluster
    Vcluster<> & v_cl;

    /*! Constructor
     *
     * \param v_cl Vcluster to use as communication object in this class
     */
    AbstractDistributionStrategy(Vcluster<> & v_cl) : v_cl(v_cl) {}

    /*! \brief Return the global id of the owned sub-sub-domain
	 *
	 * \param id in the list of owned sub-sub-domains
	 *
	 * \return the global id
	 *
	 */
    size_t getOwnerSubSubDomain(size_t id) const
    {
      return 0;
    }

    /*! \brief Return the total number of sub-sub-domains this processor own
	 *
	 * \return the total number of sub-sub-domains owned by this processor
	 *
	 */
    size_t getNOwnerSubSubDomains() const
    {
      return 0;
    }

    /*! \brief Set the tolerance for each partition
	 *
	 * \param tol tolerance
	 *
	 */
    void setDistTol(double tol) {}

    template<typename DecompositionStrategy, typename Model>
    void distribute(DecompositionStrategy & dec, Model m) {}
};

#endif //OPENFPM_PDATA_ABSTRACTDISTRIBUTIONSTRATEGY_HPP
